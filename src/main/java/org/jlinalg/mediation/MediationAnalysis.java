/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mediation;

import java.util.Objects;
import jdistlib.Normal;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.OlsResult;

/**
 * Gaussian linear mediation using ordinary least squares and analytic
 * product-of-coefficients inference.
 *
 * <p>This implementation intentionally does not simulate draws, bootstrap
 * rows, or sample from an approximate posterior. The indirect-effect standard
 * error is the first-order Sobel/delta-method estimate from the two path
 * coefficient variances. Its confidence interval and p-value use the
 * asymptotic standard-normal approximation.</p>
 */
public final class MediationAnalysis {
    private MediationAnalysis() {
    }

    /** Fits a mediation model without additional covariates. */
    public static MediationResult fit(
            double[] outcome, double[] treatment, double[] mediator) {
        return fit(outcome, treatment, mediator, null,
            OlsOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a mediation model with additional covariates. */
    public static MediationResult fit(
            double[] outcome, double[] treatment, double[] mediator,
            double[][] covariates) {
        return fit(outcome, treatment, mediator, covariates,
            OlsOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /**
     * Fits a mediation model with caller-selected OLS and compute-backend
     * policies. Covariates should not contain an intercept column; one is
     * added to each component model.
     */
    public static MediationResult fit(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            OlsOptions options,
            BackendPolicy backendPolicy) {
        Objects.requireNonNull(options, "options");
        Objects.requireNonNull(backendPolicy, "backendPolicy");
        PreparedData data = prepare(
            outcome, treatment, mediator, covariates,
            options.missingDataPolicy());

        double[][] mediatorDesign = design(
            data.treatment(), data.mediator(), data.covariates(), false);
        double[][] outcomeDesign = design(
            data.treatment(), data.mediator(), data.covariates(), true);
        double[][] totalDesign = design(
            data.treatment(), null, data.covariates(), false);

        OlsResult mediatorModel = Ols.fit(
            data.mediator(), mediatorDesign, options, backendPolicy);
        OlsResult outcomeModel = Ols.fit(
            data.outcome(), outcomeDesign, options, backendPolicy);
        OlsResult totalModel = Ols.fit(
            data.outcome(), totalDesign, options, backendPolicy);
        requireFullRank(mediatorModel, "mediator");
        requireFullRank(outcomeModel, "outcome");
        requireFullRank(totalModel, "total-effect");

        MediationEffect aPath = coefficientEffect(
            "a", mediatorModel, 1);
        MediationEffect bPath = coefficientEffect(
            "b", outcomeModel, 2);
        MediationEffect directEffect = coefficientEffect(
            "direct", outcomeModel, 1);
        MediationEffect totalEffect = coefficientEffect(
            "total", totalModel, 1);
        MediationEffect indirectEffect = indirectEffect(
            aPath.estimate(), bPath.estimate(),
            mediatorModel.covariance()[1 * mediatorModel.parameters() + 1],
            outcomeModel.covariance()[2 * outcomeModel.parameters() + 2],
            options.confidenceLevel());

        return new MediationResult(
            mediatorModel, outcomeModel, totalModel,
            aPath, bPath, indirectEffect, directEffect, totalEffect,
            data.outcome().length, data.originalObservations(),
            data.retainedRows());
    }

    private static MediationEffect coefficientEffect(
            String name, OlsResult fit, int column) {
        return new MediationEffect(
            name,
            fit.coefficients()[column],
            fit.standardErrors()[column],
            fit.tStatistics()[column],
            fit.pValues()[column],
            fit.confidenceLower()[column],
            fit.confidenceUpper()[column],
            fit.residualDegreesOfFreedom());
    }

    private static MediationEffect indirectEffect(
            double a, double b, double varianceA, double varianceB,
            double confidenceLevel) {
        double estimate = a * b;
        double variance = b * b * varianceA + a * a * varianceB;
        double standardError = Math.sqrt(Math.max(0.0, variance));
        double statistic;
        double pValue;
        if (standardError == 0.0) {
            statistic = estimate == 0.0
                ? Double.NaN : Math.copySign(Double.POSITIVE_INFINITY, estimate);
            pValue = estimate == 0.0 ? Double.NaN : 0.0;
        } else {
            statistic = estimate / standardError;
            pValue = Math.min(1.0, 2.0 * Normal.cumulative(
                Math.abs(statistic), 0.0, 1.0, false, false));
        }
        double critical = Normal.quantile(
            0.5 + confidenceLevel / 2.0, 0.0, 1.0, true, false);
        return new MediationEffect(
            "indirect", estimate, standardError, statistic, pValue,
            estimate - critical * standardError,
            estimate + critical * standardError,
            Double.POSITIVE_INFINITY);
    }

    private static double[][] design(
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            boolean includeMediator) {
        int rows = treatment.length;
        int covariateCount = covariates == null ? 0 : covariates[0].length;
        int columns = 2 + covariateCount + (includeMediator ? 1 : 0);
        double[][] design = new double[rows][columns];
        for (int row = 0; row < rows; row++) {
            design[row][0] = 1.0;
            design[row][1] = treatment[row];
            int column = 2;
            if (includeMediator) design[row][column++] = mediator[row];
            if (covariates != null) {
                System.arraycopy(covariates[row], 0,
                    design[row], column, covariateCount);
            }
        }
        return design;
    }

    private static PreparedData prepare(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            MissingDataPolicy missingDataPolicy) {
        if (outcome == null || treatment == null || mediator == null
                || outcome.length == 0
                || outcome.length != treatment.length
                || outcome.length != mediator.length) {
            throw new IllegalArgumentException(
                "outcome, treatment, and mediator must have equal positive lengths");
        }
        int rows = outcome.length;
        int covariateCount = validateCovariates(covariates, rows);
        boolean[] retained = new boolean[rows];
        int retainedCount = 0;
        for (int row = 0; row < rows; row++) {
            boolean finite = Double.isFinite(outcome[row])
                && Double.isFinite(treatment[row])
                && Double.isFinite(mediator[row]);
            for (int column = 0; finite && column < covariateCount; column++) {
                finite = Double.isFinite(covariates[row][column]);
            }
            if (!finite && missingDataPolicy == MissingDataPolicy.ERROR) {
                throw new IllegalArgumentException(
                    "mediation data contain a non-finite value in row " + row);
            }
            retained[row] = finite;
            if (finite) retainedCount++;
        }
        if (retainedCount == 0) {
            throw new IllegalArgumentException("no complete mediation observations remain");
        }

        double[] compactOutcome = new double[retainedCount];
        double[] compactTreatment = new double[retainedCount];
        double[] compactMediator = new double[retainedCount];
        double[][] compactCovariates = covariates == null
            ? null : new double[retainedCount][covariateCount];
        int[] retainedRows = new int[retainedCount];
        int target = 0;
        for (int row = 0; row < rows; row++) {
            if (!retained[row]) continue;
            compactOutcome[target] = outcome[row];
            compactTreatment[target] = treatment[row];
            compactMediator[target] = mediator[row];
            if (compactCovariates != null) {
                System.arraycopy(covariates[row], 0,
                    compactCovariates[target], 0, covariateCount);
            }
            retainedRows[target++] = row;
        }
        return new PreparedData(
            compactOutcome, compactTreatment, compactMediator,
            compactCovariates, rows, retainedRows);
    }

    private static int validateCovariates(double[][] covariates, int rows) {
        if (covariates == null) return 0;
        if (covariates.length != rows || covariates[0] == null
                || covariates[0].length == 0) {
            throw new IllegalArgumentException(
                "covariates must have one non-empty row per observation");
        }
        int columns = covariates[0].length;
        for (double[] row : covariates) {
            if (row == null || row.length != columns) {
                throw new IllegalArgumentException("covariates must be rectangular");
            }
        }
        return columns;
    }

    private static void requireFullRank(OlsResult fit, String model) {
        if (fit.rank() != fit.parameters()) {
            throw new IllegalArgumentException(
                "mediation requires a full-rank " + model + " model");
        }
    }

    private record PreparedData(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            int originalObservations,
            int[] retainedRows) {
    }
}
