/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;
import java.util.Objects;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmPqlResult;

/** PQL generalized additive fit with REML working-model smoothing updates. */
public final class GeneralizedGamResult {
    private final GlmFamily family;
    private final GlmmPqlResult workingModel;
    private final int parametricColumns;
    private final List<SmoothTermEstimate> smoothTerms;
    private final double totalEffectiveDegreesOfFreedom;

    GeneralizedGamResult(
            GlmFamily family,
            GlmmPqlResult workingModel,
            int parametricColumns,
            List<SmoothTermEstimate> smoothTerms,
            double totalEffectiveDegreesOfFreedom) {
        this.family = Objects.requireNonNull(family, "family");
        this.workingModel = Objects.requireNonNull(workingModel, "workingModel");
        this.parametricColumns = parametricColumns;
        this.smoothTerms = List.copyOf(smoothTerms);
        this.totalEffectiveDegreesOfFreedom = totalEffectiveDegreesOfFreedom;
    }

    public String family() { return family.name(); }
    public boolean converged() { return workingModel.converged(); }
    public String convergenceMessage() {
        return workingModel.convergenceMessage();
    }
    public int iterations() { return workingModel.iterations(); }
    public double conditionalDeviance() {
        return workingModel.conditionalDeviance();
    }
    public double[] linearPredictor() {
        return workingModel.linearPredictor();
    }
    public double[] fittedMeans() { return workingModel.fittedMeans(); }
    public double[] parametricCoefficients() {
        return prefix(workingModel.beta(), parametricColumns);
    }
    public double[] parametricStandardErrors() {
        return prefix(workingModel.standardErrors(), parametricColumns);
    }
    public double[] parametricTStatistics() {
        return prefix(workingModel.tStatistics(), parametricColumns);
    }
    public double[] parametricPValues() {
        return prefix(workingModel.pValues(), parametricColumns);
    }
    public List<SmoothTermEstimate> smoothTerms() { return smoothTerms; }
    public double totalEffectiveDegreesOfFreedom() {
        return totalEffectiveDegreesOfFreedom;
    }
    public double residualDegreesOfFreedom() {
        return fittedMeans().length - totalEffectiveDegreesOfFreedom;
    }
    public GlmmPqlResult workingModel() { return workingModel; }

    /** Predicts the additive predictor, including an optional link-scale offset. */
    public double[] predictLinearPredictor(
            double[] parametricDesign,
            int rows,
            List<double[]> smoothCovariates,
            double[] offset) {
        if (rows < 1 || parametricDesign == null
                || parametricDesign.length != rows * parametricColumns
                || smoothCovariates == null
                || smoothCovariates.size() != smoothTerms.size()
                || (offset != null && offset.length != rows)) {
            throw new IllegalArgumentException("prediction dimensions are invalid");
        }
        double[] beta = parametricCoefficients();
        double[] result = offset == null ? new double[rows] : offset.clone();
        for (int row = 0; row < rows; row++) {
            if (!Double.isFinite(result[row])) {
                throw new IllegalArgumentException("offset must be finite");
            }
            for (int column = 0; column < parametricColumns; column++) {
                result[row] += parametricDesign[row * parametricColumns + column]
                    * beta[column];
            }
        }
        for (int term = 0; term < smoothTerms.size(); term++) {
            double[] contribution = smoothTerms.get(term).predict(
                smoothCovariates.get(term));
            if (contribution.length != rows) {
                throw new IllegalArgumentException(
                    "each smooth covariate must have one value per prediction row");
            }
            for (int row = 0; row < rows; row++) {
                result[row] += contribution[row];
            }
        }
        return result;
    }

    /** Predicts conditional response means. */
    public double[] predictMean(
            double[] parametricDesign,
            int rows,
            List<double[]> smoothCovariates,
            double[] offset) {
        double[] result = predictLinearPredictor(
            parametricDesign, rows, smoothCovariates, offset);
        for (int row = 0; row < rows; row++) {
            result[row] = family.inverseLink(result[row]);
        }
        return result;
    }

    private static double[] prefix(double[] values, int length) {
        return java.util.Arrays.copyOf(values, length);
    }
}
