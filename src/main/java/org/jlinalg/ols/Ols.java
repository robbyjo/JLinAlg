/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.ols;

import java.util.Arrays;
import jdistlib.T;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.LeastSquaresSolver.Solution;
import org.jlinalg.internal.CompleteCases;
import org.jlinalg.internal.MatrixOps;

/** Ordinary least-squares fitting with pivoted QR and an optional SVD fallback. */
public final class Ols {
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);

    private Ols() {
    }

    /** Fits OLS from a conventional rectangular Java matrix. */
    public static OlsResult fit(double[] response, double[][] design) {
        return fit(response, design, OlsOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits OLS from a conventional rectangular Java matrix. */
    public static OlsResult fit(
            double[] response,
            double[][] design,
            OlsOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, design, null, null, options, backendPolicy);
    }

    /** Fits weighted OLS with an optional additive offset. */
    public static OlsResult fit(
            double[] response,
            double[][] design,
            double[] weights,
            double[] offset,
            OlsOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajorUnchecked(design, response.length);
        int columns = design[0].length;
        return fit(response, rowMajor, response.length, columns,
            weights, offset, options, backendPolicy);
    }

    /** Fits OLS from a contiguous row-major design matrix. */
    public static OlsResult fit(
            double[] response,
            double[] design,
            int rows,
            int columns,
            OlsOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, design, rows, columns,
            null, null, options, backendPolicy);
    }

    /** Contiguous row-major weighted OLS with an optional additive offset. */
    public static OlsResult fit(
            double[] response,
            double[] design,
            int rows,
            int columns,
            double[] weights,
            double[] offset,
            OlsOptions options,
            BackendPolicy backendPolicy) {
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        CompleteCases.Data complete = CompleteCases.prepare(
            response, design, rows, columns, weights, offset,
            options.missingDataPolicy());
        int effectiveRows = complete.response().length;
        double[] effectiveWeights = prepareWeights(
            complete.weights(), effectiveRows);
        double[] effectiveOffset = complete.offset() == null
            ? new double[effectiveRows] : complete.offset();
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return fit(complete.response(), complete.design(),
                effectiveRows, columns, effectiveWeights, effectiveOffset,
                complete.retainedRows(), complete.originalRows(), options,
                context.backend(), context.provenance());
        }
    }

    private static OlsResult fit(
            double[] response,
            double[] design,
            int rows,
            int columns,
            double[] weights,
            double[] offset,
            int[] retainedRows,
            int originalRows,
            OlsOptions options,
            ComputeBackend backend,
            BackendProvenance provenance) {
        MatrixOps.validateModelData(response, design, rows, columns);

        double[] target = new double[rows];
        double[] weightedDesign = new double[design.length];
        double[] weightedTarget = new double[rows];
        for (int row = 0; row < rows; row++) {
            target[row] = response[row] - offset[row];
            double scale = Math.sqrt(weights[row]);
            weightedTarget[row] = scale * target[row];
            for (int column = 0; column < columns; column++) {
                weightedDesign[row * columns + column] =
                    scale * design[row * columns + column];
            }
        }

        Solution solution = LeastSquaresSolver.solve(
            weightedDesign, weightedTarget, rows, columns,
            options.rankDeficiencyStrategy() == RankDeficiencyStrategy.MINIMUM_NORM,
            backend);
        int degreesOfFreedom = rows - solution.rank();
        if (degreesOfFreedom < 1) {
            throw new IllegalArgumentException(
                "OLS requires at least one residual degree of freedom");
        }

        double[] fitted = MatrixOps.multiply(
            backend, design, rows, columns, solution.coefficients());
        for (int row = 0; row < rows; row++) {
            fitted[row] += offset[row];
        }
        double[] residuals = MatrixOps.subtract(response, fitted);
        double rss = 0.0;
        double logWeightSum = 0.0;
        for (int row = 0; row < rows; row++) {
            rss += weights[row] * residuals[row] * residuals[row];
            logWeightSum += Math.log(weights[row]);
        }
        rss = Math.max(0.0, rss);
        double residualVariance = rss / degreesOfFreedom;
        double[] covariance = solution.unscaledCovariance().clone();
        for (int index = 0; index < covariance.length; index++) {
            covariance[index] *= residualVariance;
        }

        double[] standardErrors = new double[columns];
        double[] tStatistics = new double[columns];
        double[] pValues = new double[columns];
        double[] confidenceLower = new double[columns];
        double[] confidenceUpper = new double[columns];
        double alpha = 1.0 - options.confidenceLevel();
        double critical = T.quantile(1.0 - alpha / 2.0,
            degreesOfFreedom, true, false);

        for (int column = 0; column < columns; column++) {
            standardErrors[column] = Math.sqrt(Math.max(
                0.0, covariance[column * columns + column]));
            double estimate = solution.coefficients()[column];
            double standardError = standardErrors[column];
            if (standardError == 0.0) {
                tStatistics[column] = estimate == 0.0
                    ? Double.NaN : Math.copySign(Double.POSITIVE_INFINITY, estimate);
            } else {
                tStatistics[column] = estimate / standardError;
            }
            pValues[column] = Double.isNaN(tStatistics[column])
                ? Double.NaN
                : Math.min(1.0, 2.0 * T.cumulative(
                    Math.abs(tStatistics[column]), degreesOfFreedom, false, false));
            double margin = critical * standardError;
            confidenceLower[column] = estimate - margin;
            confidenceUpper[column] = estimate + margin;
        }

        double maximumLikelihoodVariance = rss / rows;
        double logLikelihood = maximumLikelihoodVariance == 0.0
            ? Double.POSITIVE_INFINITY
            : -0.5 * (rows * (LOG_TWO_PI
                + Math.log(maximumLikelihoodVariance) + 1.0) - logWeightSum);

        return new OlsResult(
            solution.coefficients(), fitted, residuals, covariance,
            standardErrors, tStatistics, pValues,
            confidenceLower, confidenceUpper,
            rows, columns, solution.rank(), degreesOfFreedom,
            rss, residualVariance, logLikelihood,
            solution.minimumNorm(), solution.tolerance(),
            retainedRows, originalRows, provenance);
    }

    private static double[] prepareWeights(double[] weights, int rows) {
        if (weights == null) {
            double[] result = new double[rows];
            Arrays.fill(result, 1.0);
            return result;
        }
        for (double value : weights) {
            if (!(value > 0.0) || !Double.isFinite(value)) {
                throw new IllegalArgumentException(
                    "weights must be finite and strictly positive");
            }
        }
        return weights;
    }

}
