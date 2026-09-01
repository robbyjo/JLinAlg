/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Cox proportional-hazards regression by Newton partial likelihood. */
public final class CoxRegression {
    private CoxRegression() { }

    public static CoxResult fit(
            CoxSurvivalData survival, double[][] covariates) {
        return fit(survival, covariates, null,
            CoxOptions.defaults(), BackendPolicy.PREFERRED);
    }

    public static CoxResult fit(
            CoxSurvivalData survival,
            double[][] covariates,
            double[] offset,
            CoxOptions options,
            BackendPolicy backendPolicy) {
        if (survival == null || covariates == null || covariates.length == 0
                || options == null || backendPolicy == null)
            throw new IllegalArgumentException("Cox model inputs are required");
        int rows = survival.observations();
        if (covariates.length != rows || covariates[0] == null
                || covariates[0].length == 0)
            throw new IllegalArgumentException(
                "covariate rows must match survival observations");
        int columns = covariates[0].length;
        double[] design = MatrixOps.rowMajor(covariates, rows);
        validateColumns(design, rows, columns);
        double[] modelOffset = offset == null ? new double[rows]
            : MatrixOps.finiteCopy(offset, "offset");
        if (modelOffset.length != rows)
            throw new IllegalArgumentException(
                "one Cox offset is required per observation");

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] beta = new double[columns];
            CoxPartialLikelihood.Evaluation evaluation =
                CoxPartialLikelihood.evaluate(survival, design, columns,
                    beta, modelOffset, options.ties());
            boolean converged = false;
            int iterations = 0;
            String message = "maximum iterations reached";
            for (int iteration = 1;
                    iteration <= options.maximumIterations(); iteration++) {
                iterations = iteration;
                if (CoxMath.maximumAbsolute(evaluation.score())
                        <= options.scoreTolerance()) {
                    converged = true;
                    message = "partial-likelihood score tolerance reached";
                    break;
                }
                double[] step = CoxMath.solvePositive(backend,
                    evaluation.information(), columns, evaluation.score(),
                    options.informationRidge());
                double scale = 1;
                CoxPartialLikelihood.Evaluation candidate = null;
                double[] next = null;
                for (int halving = 0;
                        halving <= options.maximumStepHalvings(); halving++) {
                    next = beta.clone();
                    for (int column = 0; column < columns; column++)
                        next[column] += scale * step[column];
                    candidate = CoxPartialLikelihood.evaluate(
                        survival, design, columns, next, modelOffset,
                        options.ties());
                    if (candidate.logLikelihood()
                            >= evaluation.logLikelihood() - 1e-12) break;
                    scale *= 0.5;
                }
                if (candidate.logLikelihood()
                        < evaluation.logLikelihood() - 1e-10) {
                    message = "step halving failed to improve partial likelihood";
                    break;
                }
                double relative = 0;
                for (int column = 0; column < columns; column++)
                    relative = Math.max(relative,
                        Math.abs(next[column] - beta[column])
                            / (1 + Math.abs(beta[column])));
                beta = next;
                evaluation = candidate;
                if (relative <= options.relativeTolerance()
                        && CoxMath.maximumAbsolute(evaluation.score())
                            <= Math.sqrt(options.scoreTolerance())) {
                    converged = true;
                    message = "coefficient and score tolerances reached";
                    break;
                }
            }
            double[] covariance = CoxMath.inversePositive(backend,
                evaluation.information(), columns,
                options.informationRidge());
            return new CoxResult(beta, covariance,
                CoxPartialLikelihood.baseline(survival, design, columns,
                    beta, modelOffset, options.ties()),
                evaluation.logLikelihood(), options, iterations, converged,
                message, context.provenance());
        }
    }

    private static void validateColumns(
            double[] design, int rows, int columns) {
        for (int column = 0; column < columns; column++) {
            double minimum = Double.POSITIVE_INFINITY;
            double maximum = Double.NEGATIVE_INFINITY;
            for (int row = 0; row < rows; row++) {
                double value = design[row * columns + column];
                minimum = Math.min(minimum, value);
                maximum = Math.max(maximum, value);
            }
            if (!(maximum > minimum))
                throw new IllegalArgumentException(
                    "Cox covariates must not include an intercept or constant column");
        }
    }
}
