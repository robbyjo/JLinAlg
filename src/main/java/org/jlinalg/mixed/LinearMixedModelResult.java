/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.reml.RemlResult;

/** Exact Gaussian REML fit with conditional random-effect estimates. */
public final class LinearMixedModelResult {
    private final RemlResult reml;
    private final List<RandomEffectEstimates> randomEffects;
    private final Map<String, RandomEffectEstimates> randomEffectsByName;
    private final double[] conditionalFittedValues;
    private final double[] conditionalResiduals;

    LinearMixedModelResult(
            RemlResult reml,
            List<RandomEffectEstimates> randomEffects,
            double[] conditionalFittedValues,
            double[] conditionalResiduals) {
        this.reml = Objects.requireNonNull(reml, "reml");
        this.randomEffects = List.copyOf(randomEffects);
        Map<String, RandomEffectEstimates> byName = new LinkedHashMap<>();
        for (RandomEffectEstimates value : randomEffects) {
            if (byName.put(value.termName(), value) != null) {
                throw new IllegalArgumentException(
                    "duplicate random-effect term: " + value.termName());
            }
        }
        this.randomEffectsByName = Map.copyOf(byName);
        this.conditionalFittedValues = conditionalFittedValues.clone();
        this.conditionalResiduals = conditionalResiduals.clone();
    }

    /**
     * Adapts an exact sparse coefficient-space fit to the dense LMM result
     * contract without materializing an observation covariance matrix.
     */
    public static LinearMixedModelResult fromSparse(
            SparseLinearMixedModelResult sparse,
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns) {
        Objects.requireNonNull(sparse, "sparse");
        if (response == null || response.length != rows
                || fixedEffects == null
                || fixedEffects.length != rows * columns) {
            throw new IllegalArgumentException(
                "response and fixed-effect dimensions are invalid");
        }
        double[] beta = sparse.beta();
        if (beta.length != columns) {
            throw new IllegalArgumentException(
                "sparse fixed-effect count does not match the design");
        }
        double[] marginalFitted = new double[rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                marginalFitted[row] += fixedEffects[row * columns + column]
                    * beta[column];
            }
        }
        double[] marginalResiduals = new double[rows];
        for (int row = 0; row < rows; row++) {
            marginalResiduals[row] = response[row] - marginalFitted[row];
        }
        RemlResult reml = RemlResult.fromCoefficientSpace(
            sparse.componentNames(), sparse.varianceComponents(),
            sparse.associationStatistics(), sparse.fixedEffectCovariance(),
            marginalFitted, marginalResiduals, sparse.logLikelihood(),
            sparse.varianceEstimation(), rows, columns,
            sparse.functionEvaluations(), sparse.converged(), sparse.backend());
        return new LinearMixedModelResult(reml, sparse.randomEffects(),
            sparse.conditionalFittedValues(), sparse.conditionalResiduals());
    }

    public RemlResult reml() { return reml; }
    public AssociationStatistics associationStatistics() {
        return reml.associationStatistics();
    }
    public double[] beta() { return reml.beta(); }
    /** lme4-compatible alias for fixed effects. */
    public double[] fixef() { return reml.beta(); }
    public double[] standardErrors() { return reml.standardErrors(); }
    public double[] tStatistics() { return reml.tStatistics(); }
    public double[] pValues() { return reml.pValues(); }
    public List<RandomEffectEstimates> randomEffects() { return randomEffects; }
    /** lme4-compatible named random-effect table. */
    public Map<String, RandomEffectEstimates> ranef() {
        return randomEffectsByName;
    }

    /** Scalar variance-component summaries, including residual variance. */
    public List<VarianceComponentSummary> varCorr() {
        double[] variances = reml.varianceComponents();
        List<String> names = reml.componentNames();
        java.util.ArrayList<VarianceComponentSummary> result =
            new java.util.ArrayList<>(variances.length);
        for (int index = 0; index < variances.length; index++)
            result.add(VarianceComponentSummary.of(
                names.get(index), variances[index]));
        return List.copyOf(result);
    }

    public RandomEffectEstimates randomEffects(String termName) {
        RandomEffectEstimates value = randomEffectsByName.get(termName);
        if (value == null) {
            throw new IllegalArgumentException(
                "unknown random-effect term: " + termName);
        }
        return value;
    }

    public double[] conditionalFittedValues() {
        return conditionalFittedValues.clone();
    }

    public double[] fittedValues() { return conditionalFittedValues(); }

    public double[] conditionalResiduals() {
        return conditionalResiduals.clone();
    }
    public double[] residuals() { return conditionalResiduals(); }
}
