/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;

/** GLS meta-analysis for correlated or multilevel study effects. */
public final class MetaCorrelatedAnalysis {
    private MetaCorrelatedAnalysis() { }

    public static MetaCorrelatedResult fit(List<MetaStudy> studies,
            double[][] samplingCovariance, MetaAnalysisOptions options,
            BackendPolicy backendPolicy) {
        MetaMath.Data data = MetaMath.data(studies);
        if (options == null || backendPolicy == null)
            throw new IllegalArgumentException("options and backend policy are required");
        int n = studies.size();
        double[] covariance = validateCovariance(samplingCovariance, n);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double tau = options.method() == MetaAnalysisMethod.FIXED_EFFECT ? 0.0
                : estimateTau(data.effects(), covariance, options, backend);
            Fit fit = fit(data.effects(), covariance, tau, backend);
            double df = n - 1.0;
            AssociationStatistics statistics = MetaAnalysis.statistics(options,
                new double[] {fit.beta}, new double[] {Math.sqrt(fit.variance)}, df);
            double[] weights = new double[n];
            double sum = 0.0;
            for (int i = 0; i < n; i++) { weights[i] = fit.inverse[i]; sum += weights[i]; }
            if (sum != 0.0) for (int i = 0; i < n; i++) weights[i] /= sum;
            return new MetaCorrelatedResult(studies.stream().map(MetaStudy::name).toList(),
                statistics, new double[] {fit.variance}, tau, fit.q, df, weights,
                context.provenance());
        }
    }

    private static double estimateTau(double[] effects, double[] covariance,
            MetaAnalysisOptions options, ComputeBackend backend) {
        double upper = 0.0;
        for (double value : effects) upper = Math.max(upper, value * value);
        upper = Math.max(1e-8, upper);
        double left = 0.0, right = upper;
        double best = fit(effects, covariance, 0.0, backend).q;
        if (best <= effects.length - 1.0) return 0.0;
        for (int i = 0; i < options.maximumIterations(); i++) {
            double middle = 0.5 * (left + right);
            if (fit(effects, covariance, middle, backend).q > effects.length - 1.0)
                left = middle;
            else right = middle;
            if (right - left <= options.tolerance() * Math.max(1.0, right)) break;
        }
        return 0.5 * (left + right);
    }

    private static Fit fit(double[] effects, double[] samplingCovariance,
            double tau, ComputeBackend backend) {
        int n = effects.length;
        double[] covariance = samplingCovariance.clone();
        for (int i = 0; i < n; i++) covariance[i * n + i] += tau;
        CholeskyFactor factor = backend.dpotrf(covariance, n);
        double[] ones = new double[n]; Arrays.fill(ones, 1.0);
        double[] inverseOnes = factor.solve(ones);
        double[] inverseEffects = factor.solve(effects);
        double denominator = backend.ddot(n, ones, 0, 1, inverseOnes, 0, 1);
        double beta = backend.ddot(n, ones, 0, 1, inverseEffects, 0, 1) / denominator;
        double[] residual = effects.clone();
        for (int i = 0; i < n; i++) residual[i] -= beta;
        double[] inverseResidual = factor.solve(residual);
        return new Fit(beta, 1.0 / denominator,
            backend.ddot(n, residual, 0, 1, inverseResidual, 0, 1), inverseOnes);
    }

    private static double[] validateCovariance(double[][] matrix, int n) {
        if (matrix == null || matrix.length != n)
            throw new IllegalArgumentException("sampling covariance must be n by n");
        double[] result = new double[n * n];
        for (int i = 0; i < n; i++) {
            if (matrix[i] == null || matrix[i].length != n)
                throw new IllegalArgumentException("sampling covariance must be square");
            for (int j = 0; j < n; j++) {
                double value = matrix[i][j];
                if (!Double.isFinite(value)
                        || Math.abs(value - matrix[j][i]) > 1e-10)
                    throw new IllegalArgumentException("sampling covariance must be finite and symmetric");
                result[i * n + j] = value;
            }
        }
        return result;
    }

    private record Fit(double beta, double variance, double q, double[] inverse) { }
}
