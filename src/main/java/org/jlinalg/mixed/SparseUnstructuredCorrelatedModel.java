/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;

/** Sparse grouped mixed model with an estimated unstructured random covariance. */
public final class SparseUnstructuredCorrelatedModel {
    private SparseUnstructuredCorrelatedModel() { }

    /**
     * Estimates a lower-Cholesky covariance shape by deterministic coordinate
     * refits. Sparse equation factorization remains inside each likelihood fit.
     */
    public static Result fit(double[] response, double[] fixedEffects, int rows, int columns,
                             List<String> groups, List<String> effectNames,
                             double[][] observationEffectDesign, RemlOptions options,
                             BackendPolicy backendPolicy) {
        if (response == null || fixedEffects == null || groups == null || observationEffectDesign == null
                || options == null || backendPolicy == null || effectNames == null || effectNames.isEmpty())
            throw new IllegalArgumentException("unstructured sparse inputs are invalid");
        int effects = effectNames.size(), parameters = effects * (effects + 1) / 2;
        double[] point = new double[parameters]; int position = 0; for (int row = 0; row < effects; row++) for (int column = 0; column <= row; column++) point[position++] = row == column ? 0.0 : 0.0;
        double[] lower = new double[parameters], upper = new double[parameters]; Arrays.fill(lower, -4.0); Arrays.fill(upper, 4.0);
        double bestValue = Double.NEGATIVE_INFINITY; SparseLinearMixedModelResult best = null; int evaluations = 0; double step = 0.5;
        int evaluationLimit = Math.max(24, options.maximumIterations() * 4);
        while (evaluations < evaluationLimit && step > options.relativeTolerance() * 10.0) {
            boolean improved = false;
            for (int parameter = 0; parameter < parameters; parameter++) {
                for (int direction : new int[] {-1, 1}) {
                    double[] candidate = point.clone(); candidate[parameter] = Math.max(lower[parameter], Math.min(upper[parameter], candidate[parameter] + direction * step));
                    SparseLinearMixedModelResult fit = refit(response, fixedEffects, rows, columns, groups, effectNames, observationEffectDesign, candidate, options, backendPolicy); evaluations++;
                    if (fit.logLikelihood() > bestValue) { bestValue = fit.logLikelihood(); best = fit; point = candidate; improved = true; }
                    if (evaluations >= evaluationLimit) break;
                }
            }
            if (!improved) step *= 0.5;
        }
        if (best == null) { best = refit(response, fixedEffects, rows, columns, groups, effectNames, observationEffectDesign, point, options, backendPolicy); bestValue = best.logLikelihood(); evaluations++; }
        return new Result(best, covariance(point, effects), evaluations, step <= options.relativeTolerance() * 10.0);
    }

    private static SparseLinearMixedModelResult refit(double[] response, double[] fixedEffects, int rows, int columns, List<String> groups, List<String> effectNames, double[][] design, double[] point, RemlOptions options, BackendPolicy backendPolicy) {
        return SparseCorrelatedLinearMixedModel.fit(response, fixedEffects, rows, columns, List.of(SparseCorrelatedRandomEffectBlock.of("unstructured", groups, effectNames, design, covariance(point, effectNames.size()))), options, backendPolicy);
    }
    private static double[][] covariance(double[] point, int effects) { double[][] lower = new double[effects][effects]; int position = 0; for (int row = 0; row < effects; row++) for (int column = 0; column <= row; column++) lower[row][column] = row == column ? Math.exp(point[position++]) : point[position++]; double[][] result = new double[effects][effects]; for (int row = 0; row < effects; row++) for (int column = 0; column < effects; column++) for (int k = 0; k <= Math.min(row, column); k++) result[row][column] += lower[row][k] * lower[column][k]; return result; }

    public record Result(SparseLinearMixedModelResult fit, double[][] covarianceShape,
                         int evaluations, boolean converged) {
        public Result { covarianceShape = copy(covarianceShape); }
        public double[][] covarianceShape() { return copy(covarianceShape); }
        private static double[][] copy(double[][] value) { double[][] result = new double[value.length][]; for (int row = 0; row < value.length; row++) result[row] = value[row].clone(); return result; }
    }
}
