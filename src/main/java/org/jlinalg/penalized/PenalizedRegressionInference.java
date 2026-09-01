/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import java.util.ArrayList;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.OlsResult;

/** Explicit inference procedures for penalized Gaussian regression. */
public final class PenalizedRegressionInference {
    private PenalizedRegressionInference() { }

    /**
     * Fits ridge and computes its model-based covariance. Tests use Student t
     * with residual DF based on the ridge hat-matrix trace.
     */
    public static RidgeRegressionResult ridge(
            double[] response,
            double[][] predictors,
            double lambda,
            ElasticNetOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] design = MatrixOps.rowMajor(predictors, response.length);
        int columns = predictors[0].length;
        MatrixOps.validateModelData(
            response, design, response.length, columns);
        if (options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "options and backendPolicy are required");
        }
        ElasticNetOptions ridgeOptions = PenalizedRegression.copyWithAlpha(
            options, 0.0);
        PenalizedRegressionResult fit = PenalizedRegression.fit(
            response, predictors, lambda, ridgeOptions);
        PenalizedRegression.PreparedData data = PenalizedRegression.prepare(
            response, design, response.length, columns, ridgeOptions);

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] gram = weightedGram(data);
            double[] penalized = gram.clone();
            for (int column = 0; column < columns; column++) {
                penalized[column * columns + column] +=
                    lambda * data.penaltyFactors()[column];
            }
            double[] inverse = backend.dpotrf(penalized, columns).solve(
                MatrixOps.identity(columns), columns);
            double[] influence = MatrixOps.multiply(
                backend, inverse, columns, columns, gram, columns);
            double effective = ridgeOptions.fitIntercept() ? 1.0 : 0.0;
            for (int column = 0; column < columns; column++) {
                effective += influence[column * columns + column];
            }
            double residualDf = response.length - effective;
            if (!(residualDf > 0.0)) {
                throw new IllegalArgumentException(
                    "ridge effective model dimension leaves no residual DF");
            }
            double residualVariance = fit.weightedResidualSumOfSquares()
                / residualDf;
            double[] covariance = MatrixOps.multiply(
                backend, influence, columns, columns, inverse, columns);
            for (int row = 0; row < columns; row++) {
                for (int column = 0; column < columns; column++) {
                    covariance[row * columns + column] *= residualVariance
                        / response.length
                        / data.scales()[row] / data.scales()[column];
                }
            }
            symmetrize(covariance, columns);
            double[] standardErrors = new double[columns];
            for (int column = 0; column < columns; column++) {
                standardErrors[column] = Math.sqrt(Math.max(
                    0.0, covariance[column * columns + column]));
            }
            AssociationStatistics association = AssociationStatistics.studentT(
                fit.coefficients(), standardErrors, residualDf,
                DegreesOfFreedomMethod.EFFECTIVE_RESIDUAL);
            return new RidgeRegressionResult(
                fit, association, covariance, effective,
                residualDf, residualVariance);
        }
    }

    /**
     * Refits the nonzero penalized predictors by OLS. Inference is conditional
     * on the selected active set and is not selection-adjusted.
     */
    public static PostSelectionOlsResult refitActiveSet(
            double[] response,
            double[][] predictors,
            PenalizedRegressionResult penalizedFit,
            ElasticNetOptions options,
            double activeTolerance,
            BackendPolicy backendPolicy) {
        if (response == null || penalizedFit == null
                || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "response, penalizedFit, options, and backendPolicy are required");
        }
        if (!Double.isFinite(activeTolerance) || activeTolerance < 0.0) {
            throw new IllegalArgumentException(
                "activeTolerance must be finite and nonnegative");
        }
        double[] source = MatrixOps.rowMajor(predictors, response.length);
        int sourceColumns = predictors[0].length;
        MatrixOps.validateModelData(
            response, source, response.length, sourceColumns);
        double[] penalizedBeta = penalizedFit.coefficients();
        if (penalizedBeta.length != sourceColumns) {
            throw new IllegalArgumentException(
                "penalized fit and predictors have different column counts");
        }
        List<Integer> active = new ArrayList<>();
        for (int column = 0; column < sourceColumns; column++) {
            if (Math.abs(penalizedBeta[column]) > activeTolerance) {
                active.add(column);
            }
        }
        int refitColumns = active.size() + (options.fitIntercept() ? 1 : 0);
        if (refitColumns == 0) {
            throw new IllegalArgumentException(
                "active set is empty and the model has no intercept");
        }
        double[][] refitDesign = new double[response.length][refitColumns];
        for (int row = 0; row < response.length; row++) {
            int target = 0;
            if (options.fitIntercept()) {
                refitDesign[row][target++] = 1.0;
            }
            for (int column : active) {
                refitDesign[row][target++] = source[row * sourceColumns + column];
            }
        }
        OlsResult ols = Ols.fit(response, refitDesign,
            OlsOptions.defaults(), backendPolicy);
        int[] indices = new int[active.size()];
        for (int index = 0; index < active.size(); index++) {
            indices[index] = active.get(index);
        }
        return new PostSelectionOlsResult(
            indices, options.fitIntercept(), ols);
    }

    private static double[] weightedGram(
            PenalizedRegression.PreparedData data) {
        int columns = data.columns();
        double[] result = new double[columns * columns];
        for (int row = 0; row < columns; row++) {
            for (int column = 0; column <= row; column++) {
                double value = 0.0;
                for (int observation = 0;
                        observation < data.rows(); observation++) {
                    value += data.weights()[observation]
                        * data.workingPredictors()[observation * columns + row]
                        * data.workingPredictors()[observation * columns + column]
                        / data.rows();
                }
                result[row * columns + column] = value;
                result[column * columns + row] = value;
            }
        }
        return result;
    }

    private static void symmetrize(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++) {
            for (int column = row + 1; column < dimension; column++) {
                double value = 0.5 * (matrix[row * dimension + column]
                    + matrix[column * dimension + row]);
                matrix[row * dimension + column] = value;
                matrix[column * dimension + row] = value;
            }
        }
    }
}
