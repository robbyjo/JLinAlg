/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;
import java.util.function.DoubleUnaryOperator;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;

/** Shared dense GLS; all matrices are row major. Sampling covariance is known. */
final class MetaGls {
    private MetaGls() { }

    static double[] matrix(double[][] a, int n) {
        if (a == null || a.length != n) throw new IllegalArgumentException("covariance must be n by n");
        for (double[] row : a) if (row == null || row.length != n)
            throw new IllegalArgumentException("covariance must be square");
        double[] out = new double[n * n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            if (!Double.isFinite(a[i][j]) || Math.abs(a[i][j] - a[j][i]) >
                    1e-12 * Math.max(1.0, Math.max(Math.abs(a[i][j]), Math.abs(a[j][i]))))
                throw new IllegalArgumentException("covariance must be finite and symmetric");
            out[i * n + j] = a[i][j];
        }
        return out;
    }

    static double[] design(double[][] mods, List<String> names, int n, boolean intercept) {
        if (mods == null || mods.length != n || names == null)
            throw new IllegalArgumentException("moderator rows and names are required");
        int m = names.size(), p = m + (intercept ? 1 : 0);
        if (p < 1 || n <= p) throw new IllegalArgumentException("more studies than coefficients required");
        for (String name : names) if (name == null || name.isBlank())
            throw new IllegalArgumentException("moderator names must not be blank");
        double[] x = new double[n * p];
        for (int i = 0; i < n; i++) {
            if (mods[i] == null || mods[i].length != m) throw new IllegalArgumentException("moderator widths must match names");
            if (intercept) x[i * p] = 1;
            for (int j = 0; j < m; j++) {
                if (!Double.isFinite(mods[i][j])) throw new IllegalArgumentException("moderators must be finite");
                x[i * p + j + (intercept ? 1 : 0)] = mods[i][j];
            }
        }
        return x;
    }

    static Fit fit(double[] y, double[] x, int p, double[] v, ComputeBackend backend) {
        int n = y.length;
        CholeskyFactor vf = backend.dpotrf(v.clone(), n);
        double[] wx = vf.solve(x, p), wy = vf.solve(y);
        double[] info = new double[p * p], rhs = new double[p], xtx = new double[p * p];
        for (int a = 0; a < p; a++) for (int i = 0; i < n; i++) {
            rhs[a] += x[i * p + a] * wy[i];
            for (int b = 0; b < p; b++) {
                info[a * p + b] += x[i * p + a] * wx[i * p + b];
                xtx[a * p + b] += x[i * p + a] * x[i * p + b];
            }
        }
        CholeskyFactor f = backend.dpotrf(info, p);
        double[] beta = f.solve(rhs), bread = f.solve(identity(p), p), e = y.clone();
        for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) e[i] -= x[i * p + j] * beta[j];
        double[] we = vf.solve(e);
        double q = 0;
        for (int i = 0; i < n; i++) q += e[i] * we[i];
        double ml = -0.5 * (n * Math.log(2 * Math.PI) + vf.logDeterminant() + q);
        // Orthonormal error contrasts: metafor's REML includes the |X'X| constant.
        double reml = -0.5 * ((n - p) * Math.log(2 * Math.PI) + vf.logDeterminant()
            + f.logDeterminant() - backend.dpotrf(xtx, p).logDeterminant() + q);
        return new Fit(beta, bread, e, wx, vf, q, ml, reml);
    }

    static double[] identity(int n) {
        double[] a = new double[n * n];
        for (int i = 0; i < n; i++) a[i * n + i] = 1;
        return a;
    }

    static double[] addDiagonal(double[] v, int n, double tau) {
        double[] out = v.clone();
        for (int i = 0; i < n; i++) out[i * n + i] += tau;
        return out;
    }

    static double tau(double[] y, double[] x, int p, double[] v,
                      MetaAnalysisOptions options, ComputeBackend backend) {
        if (options.method() == MetaAnalysisMethod.FIXED_EFFECT) return 0;
        int n = y.length;
        Fit zero = fit(y, x, p, v, backend);
        if (options.tauSquaredEstimator() == TauSquaredEstimator.DERSIMONIAN_LAIRD) {
            double[] w = zero.factor.solve(identity(n), n);
            double trace = 0;
            for (int i = 0; i < n; i++) {
                trace += w[i * n + i];
                for (int a = 0; a < p; a++) for (int b = 0; b < p; b++)
                    trace -= zero.wx[i * p + a] * zero.bread[a * p + b] * zero.wx[i * p + b];
            }
            return Math.max(0, (zero.q - (n - p)) / trace);
        }
        double mean = java.util.Arrays.stream(y).average().orElseThrow(), upper = 0;
        for (double value : y) upper += (value - mean) * (value - mean);
        upper = Math.max(1e-8, upper / (n - 1));
        if (options.tauSquaredEstimator() == TauSquaredEstimator.PAULE_MANDEL) {
            if (zero.q <= n - p) return 0;
            while (fit(y, x, p, addDiagonal(v, n, upper), backend).q > n - p) {
                upper *= 4;
                if (!Double.isFinite(upper)) throw new ArithmeticException("cannot bracket heterogeneity");
            }
            double lower = 0;
            for (int i = 0; i < options.maximumIterations(); i++) {
                double mid = (lower + upper) / 2;
                if (fit(y, x, p, addDiagonal(v, n, mid), backend).q > n - p) lower = mid;
                else upper = mid;
                if (upper - lower <= options.tolerance() * Math.max(1, upper)) return (lower + upper) / 2;
            }
            throw new ArithmeticException("Paule-Mandel did not converge");
        }
        DoubleUnaryOperator objective = t -> -fit(y, x, p, addDiagonal(v, n, t), backend).reml;
        double prev = objective.applyAsDouble(0), current = objective.applyAsDouble(upper);
        while (current < prev) {
            prev = current; upper *= 4;
            if (!Double.isFinite(upper)) throw new ArithmeticException("cannot bracket REML");
            current = objective.applyAsDouble(upper);
        }
        double lower = 0, ratio = (Math.sqrt(5) - 1) / 2;
        double a = upper - ratio * upper, b = ratio * upper;
        double fa = objective.applyAsDouble(a), fb = objective.applyAsDouble(b);
        for (int i = 0; i < options.maximumIterations(); i++) {
            if (fa < fb) { upper = b; b = a; fb = fa; a = upper - ratio * (upper - lower); fa = objective.applyAsDouble(a); }
            else { lower = a; a = b; fa = fb; b = lower + ratio * (upper - lower); fb = objective.applyAsDouble(b); }
            if (upper - lower <= options.tolerance() * Math.max(1, upper)) {
                double candidate = (lower + upper) / 2;
                return objective.applyAsDouble(0) <= objective.applyAsDouble(candidate) ? 0 : candidate;
            }
        }
        throw new ArithmeticException("REML did not converge");
    }

    record Fit(double[] beta, double[] bread, double[] residual, double[] wx,
               CholeskyFactor factor, double q, double ml, double reml) { }
}
