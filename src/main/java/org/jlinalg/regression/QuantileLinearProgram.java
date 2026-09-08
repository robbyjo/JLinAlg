/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.PivotedQrFactor;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/**
 * Deterministic Mehrotra predictor-corrector solver for the actual quantile LP:
 * min tau*sum(u)+(1-tau)*sum(v), X*beta+u-v=y, u,v >= 0.
 * Its dual maximizes y'z subject to X'z=0 and {@code tau-1 <= z <= tau}.
 * The barrier is an LP algorithm, not a smoothed pinball objective. Numerical
 * convergence requires primal/dual feasibility, complementarity and duality gap.
 * Dense storage is O(n*p+p*p), with weighted pivoted QR for Newton systems.
 * No coefficient uniqueness, inferential covariance or standard errors are assumed.
 * Original-unit reconstruction is also certified. Large fitted offsets can lose
 * residual precision in double arithmetic and return converged=false even after
 * the internally scaled LP converges; callers may explicitly choose a tolerance
 * appropriate to their data's representable precision.
 */
public final class QuantileLinearProgram {
    private QuantileLinearProgram() { }

    public record Options(int maximumIterations, double tolerance) {
        public Options {
            if (maximumIterations < 1 || !Double.isFinite(tolerance) || !(tolerance > 0 && tolerance < 1))
                throw new IllegalArgumentException("positive iteration budget and tolerance in (0,1) required");
        }
        public static Options defaults() { return new Options(200, 1e-9); }
    }

    /** Residuals and gap are dimensionless after design/effect scaling; objectives use original units. */
    public record Certificate(double primalObjective, double dualObjective,
            double primalResidual, double dualResidual, double boundViolation,
            double complementarity, double relativeGap) { }

    /** The dual weights independently certify an optimum; coefficients need not be unique. */
    public record Result(QuantileRegressionResult fit, double[] dualWeights, Certificate certificate) {
        public Result { dualWeights = dualWeights.clone(); }
        @Override public double[] dualWeights() { return dualWeights.clone(); }
    }

    public static Result solve(double[] response, double[][] predictors, double quantile) {
        return solve(response, predictors, quantile, Options.defaults());
    }

    public static Result solve(double[] response, double[][] predictors, double quantile, Options options) {
        if (response == null || predictors == null || options == null || response.length == 0
                || predictors.length != response.length || !(quantile > 0 && quantile < 1))
            throw new IllegalArgumentException("finite effects, design and a quantile in (0,1) required");
        int n = response.length, p = predictors[0] == null ? 0 : predictors[0].length;
        if (p < 1 || p > n) throw new IllegalArgumentException("the exact LP requires full column rank and n >= p >= 1");
        double[] x = new double[Math.multiplyExact(n, p)], columnScale = new double[p];
        double responseScale = 0;
        for (int i = 0; i < n; i++) {
            if (!Double.isFinite(response[i]) || predictors[i] == null || predictors[i].length != p)
                throw new IllegalArgumentException("nonfinite response or inconsistent predictor rows");
            responseScale = Math.max(responseScale, Math.abs(response[i]));
            for (int j = 0; j < p; j++) {
                if (!Double.isFinite(predictors[i][j])) throw new IllegalArgumentException("predictors must be finite");
                columnScale[j] = Math.max(columnScale[j], Math.abs(predictors[i][j]));
            }
        }
        for (int j = 0; j < p; j++) if (!(columnScale[j] > 0))
            throw new IllegalArgumentException("rank-deficient design: zero column");
        for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) x[i*p+j] = predictors[i][j] / columnScale[j];
        if (responseScale == 0) responseScale = 1;
        double[] scaledY = new double[n];
        for (int i = 0; i < n; i++) scaledY[i] = response[i] / responseScale;
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            ComputeBackend backend = context.backend();
            PivotedQrFactor initialQr = backend.dgeqp3(x.clone(), n, p);
            requireRank(initialQr, p);
            double[] origin = initialQr.solveLeastSquares(scaledY), y = new double[n];
            double residualScale = 0;
            for (int i = 0; i < n; i++) {
                y[i] = scaledY[i] - rowDot(x, i, origin, p);
                residualScale = Math.max(residualScale, Math.abs(y[i]));
            }
            // Subtract an OLS fitted vector without changing the LP estimand. This
            // prevents a large intercept from hiding the scale of the residual loss.
            if (residualScale < 32 * Math.ulp(1.0)) residualScale = 1;
            for (int i = 0; i < n; i++) y[i] /= residualScale;
            double unit = responseScale * residualScale;
            if (!(unit > 0) || !Double.isFinite(unit)) throw new ArithmeticException("response scale cannot be represented");
            double[] beta = new double[p], u = new double[n], v = new double[n], z = new double[n];
            for (int i = 0; i < n; i++) { u[i] = Math.max(y[i], 0) + 1; v[i] = Math.max(-y[i], 0) + 1; }
            State state = new State(y, x, beta, u, v, z, quantile, n, p);
            int iterations = 0;
            for (; iterations < options.maximumIterations(); iterations++) {
                if (certified(state, options.tolerance())) break;
                double mu = state.mu();
                NewtonSystem system = system(state, backend);
                Direction affine = direction(state, 0, null, system);
                double alphaP = primalStep(state, affine), alphaD = dualStep(state, affine);
                double muAffine = 0;
                for (int i = 0; i < n; i++) {
                    muAffine += (u[i] + alphaP*affine.u[i]) * (quantile-z[i]-alphaD*affine.z[i]);
                    muAffine += (v[i] + alphaP*affine.v[i]) * (1-quantile+z[i]+alphaD*affine.z[i]);
                }
                muAffine /= 2*n;
                double sigma = Math.pow(Math.max(0, Math.min(1, muAffine/mu)), 3);
                Direction corrected = direction(state, sigma*mu, affine, system);
                alphaP = Math.min(1, .995*primalStep(state, corrected));
                alphaD = Math.min(1, .995*dualStep(state, corrected));
                if (!(alphaP > 0 && alphaD > 0)) throw new ArithmeticException("LP step lost strict feasibility");
                for (int j = 0; j < p; j++) beta[j] += alphaP*corrected.beta[j];
                for (int i = 0; i < n; i++) {
                    u[i] += alphaP*corrected.u[i]; v[i] += alphaP*corrected.v[i]; z[i] += alphaD*corrected.z[i];
                    if (!(u[i] > 0 && v[i] > 0 && quantile-z[i] > 0 && 1-quantile+z[i] > 0)
                            || !Double.isFinite(u[i]+v[i])) throw new ArithmeticException("LP numerical feasibility failure");
                }
            }
            double[] coefficients = new double[p], fitted = new double[n], residual = new double[n];
            for (int j = 0; j < p; j++) {
                coefficients[j] = (origin[j]*responseScale + beta[j]*unit) / columnScale[j];
                if (!Double.isFinite(coefficients[j])) throw new ArithmeticException("coefficient overflow");
            }
            double primal = 0, primalError = 0, comp = 0;
            for (int i = 0; i < n; i++) {
                fitted[i] = dot(predictors[i], coefficients);
                residual[i] = response[i] - fitted[i];
                if (!Double.isFinite(residual[i])) throw new ArithmeticException("nonfinite fitted residual");
                primal += pinball(residual[i], quantile);
                primalError = Math.max(primalError, Math.abs(residual[i]/unit - u[i] + v[i]));
                comp = Math.max(comp, Math.abs(residual[i]/unit) * (residual[i] >= 0 ? quantile-z[i] : 1-quantile+z[i]));
            }
            double dual = dot(response, z), denominator = Math.max(unit, Math.max(Math.abs(primal), Math.abs(dual)));
            Certificate certificate = new Certificate(primal, dual, primalError, dualResidual(state),
                boundViolation(state), comp, Math.abs(primal-dual)/denominator);
            if (!finite(certificate)) throw new ArithmeticException("nonfinite exact quantile optimality certificate");
            boolean converged = certified(state, options.tolerance())
                && primalError <= options.tolerance()*10 && comp <= options.tolerance()*10
                && certificate.relativeGap() <= options.tolerance()*10;
            return new Result(new QuantileRegressionResult(coefficients, fitted, residual, primal,
                quantile, iterations, converged), z, certificate);
        }
    }

    /** Both predictor and corrector share exactly the same weighted QR factorization. */
    private static NewtonSystem system(State s, ComputeBackend backend) {
        double[] w = new double[s.n]; double largest = 0;
        for (int i = 0; i < s.n; i++) {
            w[i] = 1/(s.u[i]/(s.tau-s.z[i])+s.v[i]/(1-s.tau+s.z[i]));
            if (!(w[i] > 0) || !Double.isFinite(w[i])) throw new ArithmeticException("LP Newton weights are not representable");
            largest = Math.max(largest,w[i]);
        }
        double[] weighted = new double[s.n*s.p];
        for (int i = 0; i < s.n; i++) for (int j = 0; j < s.p; j++)
            weighted[i*s.p+j] = Math.sqrt(w[i]/largest)*s.x[i*s.p+j];
        PivotedQrFactor qr = backend.dgeqp3(weighted,s.n,s.p);
        requireRank(qr,s.p);
        return new NewtonSystem(w,largest,qr);
    }

    private static Direction direction(State s, double target, Direction affine, NewtonSystem system) {
        int n = s.n, p = s.p;
        double[] w = system.weights, h = new double[n], cu = new double[n], cv = new double[n];
        double largestWeight = system.largestWeight;
        for (int i = 0; i < n; i++) {
            double a = s.tau-s.z[i], b = 1-s.tau+s.z[i];
            cu[i] = target-s.u[i]*a; cv[i] = target-s.v[i]*b;
            if (affine != null) { cu[i] += affine.u[i]*affine.z[i]; cv[i] -= affine.v[i]*affine.z[i]; }
            h[i] = s.y[i]-rowDot(s.x,i,s.beta,p)-s.u[i]+s.v[i]-cu[i]/a+cv[i]/b;
            if (!Double.isFinite(h[i]))
                throw new ArithmeticException("LP Newton scaling is not representable");
        }
        double[] rhs = new double[p];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                rhs[j] += s.x[i*p+j] * ((w[i]/largestWeight)*h[i] + s.z[i]/largestWeight);
            }
        }
        double[] dbeta = solveNormal(system.qr,rhs), du = new double[n], dv = new double[n], dz = new double[n];
        for (int i = 0; i < n; i++) {
            dz[i] = w[i]*(h[i]-rowDot(s.x,i,dbeta,p));
            du[i] = (cu[i]+s.u[i]*dz[i])/(s.tau-s.z[i]);
            dv[i] = (cv[i]-s.v[i]*dz[i])/(1-s.tau+s.z[i]);
            if (!Double.isFinite(du[i]) || !Double.isFinite(dv[i]) || !Double.isFinite(dz[i]))
                throw new ArithmeticException("nonfinite LP Newton direction");
        }
        return new Direction(dbeta,du,dv,dz);
    }

    /** Solve R'R in pivoted coordinates without forming the squared-condition normal matrix. */
    private static double[] solveNormal(PivotedQrFactor qr, double[] rhs) {
        int p = rhs.length; double[] packed = qr.packed(), first = new double[p], second = new double[p], out = new double[p];
        int[] pivot = qr.pivot();
        for (int i = 0; i < p; i++) {
            double value = rhs[pivot[i]];
            for (int j = 0; j < i; j++) value -= packed[j*p+i]*first[j];
            first[i] = value/packed[i*p+i];
        }
        for (int i = p-1; i >= 0; i--) {
            double value = first[i];
            for (int j = i+1; j < p; j++) value -= packed[i*p+j]*second[j];
            second[i] = value/packed[i*p+i]; out[pivot[i]] = second[i];
        }
        return out;
    }

    private static void requireRank(PivotedQrFactor qr, int p) {
        if (qr.rank() != p) throw new IllegalArgumentException("rank-deficient or numerically singular exact quantile design/system");
        double[] r = qr.packed(); double largest = 0, smallest = Double.POSITIVE_INFINITY;
        for (int j = 0; j < p; j++) { largest = Math.max(largest,Math.abs(r[j*p+j])); smallest = Math.min(smallest,Math.abs(r[j*p+j])); }
        if (!(smallest > 1e-12*largest)) throw new IllegalArgumentException("exact quantile system is too ill-conditioned for a reliable solution");
    }
    private static double primalStep(State s, Direction d) {
        double step = 1;
        for (int i = 0; i < s.n; i++) {
            if (d.u[i] < 0) step = Math.min(step,-s.u[i]/d.u[i]);
            if (d.v[i] < 0) step = Math.min(step,-s.v[i]/d.v[i]);
        }
        return step;
    }
    private static double dualStep(State s, Direction d) {
        double step = 1;
        for (int i = 0; i < s.n; i++) {
            if (d.z[i] > 0) step = Math.min(step,(s.tau-s.z[i])/d.z[i]);
            if (d.z[i] < 0) step = Math.min(step,-(1-s.tau+s.z[i])/d.z[i]);
        }
        return step;
    }
    private static boolean certified(State s, double tolerance) {
        double primal = 0, error = 0, complementarity = 0;
        for (int i = 0; i < s.n; i++) {
            double residual = s.y[i]-rowDot(s.x,i,s.beta,s.p);
            primal += pinball(residual,s.tau);
            error = Math.max(error,Math.abs(residual-s.u[i]+s.v[i]));
            complementarity = Math.max(complementarity,Math.max(s.u[i]*(s.tau-s.z[i]),s.v[i]*(1-s.tau+s.z[i])));
        }
        double dual = dot(s.y,s.z), gap = Math.abs(primal-dual)/Math.max(1,Math.max(Math.abs(primal),Math.abs(dual)));
        return Double.isFinite(gap) && gap <= tolerance && error <= tolerance
            && dualResidual(s) <= tolerance && boundViolation(s) <= tolerance && complementarity <= tolerance;
    }
    private static double dualResidual(State s) {
        double error = 0;
        for (int j = 0; j < s.p; j++) {
            double value = 0, norm = 0;
            for (int i = 0; i < s.n; i++) { value += s.x[i*s.p+j]*s.z[i]; norm += Math.abs(s.x[i*s.p+j]); }
            error = Math.max(error,Math.abs(value)/Math.max(1,norm));
        }
        return error;
    }
    private static double boundViolation(State s) {
        double error = 0;
        for (double z : s.z) error = Math.max(error,Math.max(z-s.tau,s.tau-1-z));
        return error;
    }
    private static boolean finite(Certificate c) {
        return Double.isFinite(c.primalObjective()) && Double.isFinite(c.dualObjective()) && Double.isFinite(c.primalResidual())
            && Double.isFinite(c.dualResidual()) && Double.isFinite(c.complementarity()) && Double.isFinite(c.relativeGap());
    }
    private static double pinball(double residual,double tau) { return residual >= 0 ? tau*residual : (tau-1)*residual; }
    private static double rowDot(double[] x,int row,double[] beta,int p) {
        double value = 0; for (int j = 0; j < p; j++) value += x[row*p+j]*beta[j]; return value;
    }
    private static double dot(double[] a,double[] b) {
        double value = 0, correction = 0;
        for (int i = 0; i < a.length; i++) { double term = a[i]*b[i]-correction, next = value+term; correction = (next-value)-term; value = next; }
        return value;
    }
    private record Direction(double[] beta,double[] u,double[] v,double[] z) { }
    private record NewtonSystem(double[] weights,double largestWeight,PivotedQrFactor qr) { }
    private record State(double[] y,double[] x,double[] beta,double[] u,double[] v,double[] z,double tau,int n,int p) {
        double mu() {
            double value = 0; for (int i = 0; i < n; i++) value += u[i]*(tau-z[i])+v[i]*(1-tau+z[i]); return value/(2*n);
        }
    }
}
