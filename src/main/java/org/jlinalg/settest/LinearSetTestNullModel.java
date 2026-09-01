/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.Arrays;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.OlsResult;

/**
 * Covariate-only Gaussian null model reused by Burden, SKAT, and SKAT-O.
 *
 * <p>The response and each variant are residualized against the same fixed
 * design. The retained information inverse avoids a model refit per set.</p>
 */
public final class LinearSetTestNullModel implements GaussianSetTestNullModel {
    private final double[] responseResidual;
    private final double[] fixedDesign;
    private final double[] informationInverse;
    private final int observations;
    private final int fixedColumns;
    private final int residualDegreesOfFreedom;
    private final double residualSumSquares;
    private final double residualVariance;
    private final BackendPolicy backendPolicy;

    private LinearSetTestNullModel(
            double[] responseResidual, double[] fixedDesign,
            double[] informationInverse, int observations, int fixedColumns,
            int residualDegreesOfFreedom, double residualSumSquares,
            double residualVariance, BackendPolicy backendPolicy) {
        this.responseResidual = responseResidual;
        this.fixedDesign = fixedDesign;
        this.informationInverse = informationInverse;
        this.observations = observations;
        this.fixedColumns = fixedColumns;
        this.residualDegreesOfFreedom = residualDegreesOfFreedom;
        this.residualSumSquares = residualSumSquares;
        this.residualVariance = residualVariance;
        this.backendPolicy = backendPolicy;
    }

    public static LinearSetTestNullModel prepare(
            double[] response, double[][] fixedDesign,
            OlsOptions options, BackendPolicy backendPolicy) {
        if (response == null || fixedDesign == null || options == null
                || backendPolicy == null)
            throw new IllegalArgumentException("null-model inputs are required");
        OlsResult fit = Ols.fit(
            response, fixedDesign, options, backendPolicy);
        if (fit.rankDeficient())
            throw new IllegalArgumentException(
                "set-test fixed design must have full column rank");
        double variance = fit.residualVariance();
        if (!(variance > 0) || !Double.isFinite(variance))
            throw new IllegalArgumentException(
                "set-test null model requires positive residual variance");
        double[] inverse = fit.covariance();
        for (int index = 0; index < inverse.length; index++)
            inverse[index] /= variance;
        return new LinearSetTestNullModel(fit.residuals(),
            MatrixOps.rowMajor(fixedDesign, response.length), inverse,
            response.length, fixedDesign[0].length,
            fit.residualDegreesOfFreedom(), fit.residualSumOfSquares(),
            variance, backendPolicy);
    }

    /** Rows are variants and columns are aligned samples. */
    public double[][] residualize(double[][] variants) {
        if (variants == null || variants.length == 0)
            throw new IllegalArgumentException("variant rows are required");
        int variantCount = variants.length;
        double[] sampleByVariant = new double[observations * variantCount];
        for (int variant = 0; variant < variantCount; variant++) {
            if (variants[variant] == null
                    || variants[variant].length != observations)
                throw new IllegalArgumentException(
                    "variant row length must equal null-model observations");
            MatrixOps.requireFinite(variants[variant], "variant row");
            for (int sample = 0; sample < observations; sample++)
                sampleByVariant[sample * variantCount + variant] =
                    variants[variant][sample];
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] cross = MatrixOps.transposeMultiply(backend,
                fixedDesign, observations, fixedColumns,
                sampleByVariant, variantCount);
            double[] coefficients = MatrixOps.multiply(backend,
                informationInverse, fixedColumns, fixedColumns,
                cross, variantCount);
            double[] fitted = MatrixOps.multiply(backend,
                fixedDesign, observations, fixedColumns,
                coefficients, variantCount);
            double[][] result = new double[variantCount][observations];
            for (int sample = 0; sample < observations; sample++)
                for (int variant = 0; variant < variantCount; variant++)
                    result[variant][sample] =
                        sampleByVariant[sample * variantCount + variant]
                        - fitted[sample * variantCount + variant];
            return result;
        }
    }

    /**
     * Returns standardized Gaussian scores and their covariance under the
     * fixed-effect null model.
     */
    @Override
    public SetTestScoreState score(double[][] variantRows) {
        double[][] residual = residualize(variantRows);
        double scale = Math.sqrt(residualVariance);
        double[] scores = new double[residual.length];
        double[] information = new double[residual.length * residual.length];
        for (int left = 0; left < residual.length; left++) {
            scores[left] = dot(residual[left], responseResidual) / scale;
            for (int right = 0; right <= left; right++) {
                double value = dot(residual[left], residual[right]);
                information[left * residual.length + right] = value;
                information[right * residual.length + left] = value;
            }
        }
        return new SetTestScoreState(scores, information, residual.length);
    }

    public double[] responseResiduals() { return responseResidual.clone(); }
    @Override public double degreesOfFreedom() {
        return residualDegreesOfFreedom;
    }
    public int observations() { return observations; }
    public int fixedColumns() { return fixedColumns; }
    public int residualDegreesOfFreedom() { return residualDegreesOfFreedom; }
    public double residualSumSquares() { return residualSumSquares; }
    public double residualVariance() { return residualVariance; }
    public BackendPolicy backendPolicy() { return backendPolicy; }

    double[] responseResidualView() { return responseResidual; }

    private static double dot(double[] left, double[] right) {
        double value = 0;
        for (int index = 0; index < left.length; index++)
            value += left[index] * right[index];
        return value;
    }
}
