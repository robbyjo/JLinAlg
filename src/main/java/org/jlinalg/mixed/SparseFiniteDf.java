/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.function.Function;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.inference.DegreesOfFreedomMethod;

/**
 * Finite-denominator-DF inference driven by sparse REML likelihood evaluations.
 * The fitting layer supplies sparse solves, not caller-written derivatives.
 * Numerical derivatives and information inverses have parameter/fixed-effect
 * dimensions and never materialize observation-scale covariance matrices.
 */
public final class SparseFiniteDf {
    private SparseFiniteDf() { }

    /** One sparse likelihood evaluation in absolute, linear covariance coordinates. */
    record Point(double logLikelihood, double restrictedLogDeterminant,
                 double[] covariance, double residualVariance) { }
    record JointState(double[] unadjustedCovariance,
            double[][] covarianceGradient,
            double[] varianceParameterCovariance) { }
    record Inference(double[] covariance, double[] degreesOfFreedom,
            JointState jointState) { }

    /** Box-boundary estimates do not have the unconstrained information law. */
    static void requireInteriorVariances(double[] variances, double minimum, double maximum) {
        double residual = variances[variances.length - 1];
        for (double variance : variances) {
            // Resolve boundaries at the optimizer's physical variance scale,
            // not by shrinking finite differences until floating-point noise fits.
            double tolerance = 1e-6 * Math.max(variance, residual);
            if (variance - minimum <= tolerance || maximum - variance <= tolerance)
                throw new IllegalArgumentException("finite-DF inference requires variances away from their physical bounds");
        }
    }

    /**
     * Differentiates the sparse likelihood and fixed-effect covariance. Only
     * parameter-by-parameter and fixed-effect-by-fixed-effect matrices are held.
     * Satterthwaite uses observed REML information. KR uses expected information
     * -1/2 d2 log|V|_REML and C_A = C - sum_ij W_ij C_ij. The latter identity
     * is the usual 2 C sum W_ij(Q_ij-P_i C P_j) C adjustment, valid because
     * these are linear (absolute covariance-entry) coordinates, not Cholesky ones.
     */
    static Inference compute(double[] parameters, Function<double[], Point> evaluate,
            int columns, DegreesOfFreedomMethod method, ComputeBackend backend) {
        return compute(parameters, null, evaluate, columns, method, backend);
    }

    static Inference compute(double[] parameters, double[] parameterScales,
            Function<double[], Point> evaluate, int columns,
            DegreesOfFreedomMethod method, ComputeBackend backend) {
        int count = parameters.length;
        Point center = evaluate.apply(parameters.clone());
        double[] step = new double[count];
        Point[] plus = new Point[count], minus = new Point[count];
        for (int i = 0; i < count; i++) {
            step[i] = 1e-3 * (parameterScales == null
                ? Math.max(Math.abs(parameters[i]), Math.abs(parameters[count - 1]) * .1)
                : parameterScales[i]);
            boolean valid = false;
            for (int attempt = 0; attempt < 12; attempt++) {
                try {
                    plus[i] = shifted(evaluate, parameters, i, step[i], i, 0);
                    minus[i] = shifted(evaluate, parameters, i, -step[i], i, 0);
                    valid = true;
                    break;
                } catch (IllegalArgumentException | IllegalStateException failure) {
                    step[i] *= .5;
                }
            }
            if (!valid) throw new IllegalArgumentException(
                "finite-DF derivatives require an interior identifiable covariance estimate");
        }
        boolean kr = method == DegreesOfFreedomMethod.KENWARD_ROGER;
        double[] information = new double[count * count];
        double[][] gradient = new double[count][columns * columns];
        double[][] curvature = kr ? new double[count * count][] : null;
        for (int i = 0; i < count; i++) {
            for (int k = 0; k < columns * columns; k++)
                gradient[i][k] = (plus[i].covariance()[k] - minus[i].covariance()[k]) / (2 * step[i]);
            for (int j = 0; j <= i; j++) {
                Point pp, pm, mp, mm;
                double denominator;
                if (i == j) {
                    pp = plus[i]; pm = center; mp = center; mm = minus[i];
                    denominator = step[i] * step[i];
                } else {
                    pp = shifted(evaluate, parameters, i, step[i], j, step[j]);
                    pm = shifted(evaluate, parameters, i, step[i], j, -step[j]);
                    mp = shifted(evaluate, parameters, i, -step[i], j, step[j]);
                    mm = shifted(evaluate, parameters, i, -step[i], j, -step[j]);
                    denominator = 4 * step[i] * step[j];
                }
                double hessian = kr
                    ? .5 * (pp.restrictedLogDeterminant() - pm.restrictedLogDeterminant()
                        - mp.restrictedLogDeterminant() + mm.restrictedLogDeterminant()) / denominator
                    : (pp.logLikelihood() - pm.logLikelihood() - mp.logLikelihood()
                        + mm.logLikelihood()) / denominator;
                information[i * count + j] = information[j * count + i] = -hessian;
                if (kr) {
                    double[] derivative = new double[columns * columns];
                    for (int k = 0; k < derivative.length; k++)
                        derivative[k] = (pp.covariance()[k] - pm.covariance()[k]
                            - mp.covariance()[k] + mm.covariance()[k]) / denominator;
                    curvature[i * count + j] = curvature[j * count + i] = derivative;
                }
            }
        }
        double[] parameterCovariance;
        try {
            parameterCovariance = backend.dpotrf(information, count)
                .solve(MatrixOps.identity(count), count);
        } catch (IllegalArgumentException | IllegalStateException failure) {
            throw new IllegalArgumentException("finite-DF information is not positive definite; "
                + "boundary or unidentifiable covariance parameters", failure);
        }
        double[] covariance = center.covariance().clone();
        if (kr) for (int i = 0; i < count * count; i++)
            for (int k = 0; k < covariance.length; k++)
                covariance[k] -= parameterCovariance[i] * curvature[i][k];
        double[] degrees = new double[columns];
        for (int c = 0; c < columns; c++) {
            double[] selected = new double[count];
            for (int i = 0; i < count; i++) selected[i] = gradient[i][c * columns + c];
            // For a one-dimensional KR contrast A1=A2; its moment-matched DF
            // reduces to 2/A2, using unadjusted C and expected-information W.
            degrees[c] = satterthwaite(center.covariance()[c * columns + c],
                selected, parameterCovariance);
            if (!(degrees[c] > 0) || !Double.isFinite(degrees[c])
                    || !(covariance[c * columns + c] > 0))
                throw new IllegalArgumentException("finite-DF inference is not estimable");
        }
        return new Inference(covariance, degrees,
            kr ? new JointState(center.covariance().clone(), gradient,
                parameterCovariance.clone()) : null);
    }

    private static Point shifted(Function<double[], Point> evaluate, double[] parameters,
            int i, double left, int j, double right) {
        double[] shifted = parameters.clone();
        shifted[i] += left;
        shifted[j] += right;
        Point result = evaluate.apply(shifted);
        if (!Double.isFinite(result.logLikelihood()))
            throw new IllegalArgumentException("nonfinite derivative likelihood");
        return result;
    }

    /** Satterthwaite DF from a scalar contrast variance and its derivatives. */
    public static double satterthwaite(double variance,
            double[] varianceDerivatives, double[] varianceParameterCovariance) {
        if (!(variance > 0.0) || !Double.isFinite(variance)
                || varianceDerivatives == null || varianceParameterCovariance == null
                || varianceParameterCovariance.length != varianceDerivatives.length
                    * varianceDerivatives.length)
            throw new IllegalArgumentException("finite-DF inputs are invalid");
        double uncertainty = 0.0;
        int n = varianceDerivatives.length;
        for (int row = 0; row < n; row++) for (int column = 0; column < n; column++)
            uncertainty += varianceDerivatives[row] * varianceParameterCovariance[row * n + column]
                * varianceDerivatives[column];
        if (!(uncertainty > 0.0) || !Double.isFinite(uncertainty))
            return Double.POSITIVE_INFINITY;
        return 2.0 * variance * variance / uncertainty;
    }

    /** Simple variance addition, not a Kenward-Roger inference procedure. */
    public static double adjustedVariance(double modelVariance,
            double adjustment) {
        if (!(modelVariance > 0.0) || !Double.isFinite(modelVariance)
                || adjustment < 0.0 || !Double.isFinite(adjustment))
            throw new IllegalArgumentException("variance adjustment is invalid");
        return modelVariance + adjustment;
    }
}
