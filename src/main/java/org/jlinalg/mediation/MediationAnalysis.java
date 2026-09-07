/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mediation;

import java.util.Objects;
import java.util.ArrayList;
import java.util.List;
import jdistlib.Normal;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.OlsResult;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.pedigree.SparsePedigreeMixedModel;
import org.jlinalg.reml.RemlOptions;

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
    private static final double DEFAULT_CONFIDENCE_LEVEL = 0.95;

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

    /**
     * Fits mediation with one or more independent sparse Gaussian random
     * effects shared by the mediator, outcome, and total-effect models.
     * Sparse REML structure is prepared once and reused across all three
     * fits.
     */
    public static MediationMixedResult fitMixed(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            List<RandomEffectTerm> randomEffects) {
        return fitMixed(outcome, treatment, mediator, null, randomEffects,
            RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits sparse REML mediation with optional fixed-effect covariates. */
    public static MediationMixedResult fitMixed(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        Objects.requireNonNull(randomEffects, "randomEffects");
        if (randomEffects.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one ordinary random-effect term is required");
        }
        Objects.requireNonNull(options, "options");
        Objects.requireNonNull(backendPolicy, "backendPolicy");
        return fitMixed(outcome, treatment, mediator, covariates,
            randomEffects, options, DEFAULT_CONFIDENCE_LEVEL, backendPolicy);
    }

    /** Fits sparse REML mediation with an explicit confidence level. */
    public static MediationMixedResult fitMixed(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            double confidenceLevel,
            BackendPolicy backendPolicy) {
        Objects.requireNonNull(randomEffects, "randomEffects");
        if (randomEffects.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one ordinary random-effect term is required");
        }
        Objects.requireNonNull(options, "options");
        Objects.requireNonNull(backendPolicy, "backendPolicy");
        validateConfidenceLevel(confidenceLevel);
        PreparedData data = prepare(
            outcome, treatment, mediator, covariates,
            MissingDataPolicy.ERROR);
        List<RandomEffectTerm> compactTerms = compactTerms(
            randomEffects, data, "random-effect");
        try (SparseLinearMixedModel.Prepared prepared =
                SparseLinearMixedModel.prepare(
                    data.outcome().length, compactTerms,
                    options, backendPolicy)) {
            return fitPreparedMixed(
                data, prepared, confidenceLevel);
        }
    }

    /**
     * Fits sparse REML mediation with pedigree random effects and optional
     * independent random effects. The same terms are shared by all three
     * component models; pedigree precision remains sparse in coefficient
     * space.
     */
    public static MediationMixedResult fitPedigree(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            List<PedigreeRandomEffectTerm> pedigreeEffects) {
        return fitPedigree(outcome, treatment, mediator, null,
            pedigreeEffects, List.of(), RemlOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    /** Fits pedigree-aware sparse REML mediation with fixed covariates. */
    public static MediationMixedResult fitPedigree(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            List<PedigreeRandomEffectTerm> pedigreeEffects,
            List<RandomEffectTerm> ordinaryEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        Objects.requireNonNull(pedigreeEffects, "pedigreeEffects");
        if (pedigreeEffects.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one pedigree random-effect term is required");
        }
        Objects.requireNonNull(options, "options");
        Objects.requireNonNull(backendPolicy, "backendPolicy");
        return fitPedigree(outcome, treatment, mediator, covariates,
            pedigreeEffects, ordinaryEffects, options,
            DEFAULT_CONFIDENCE_LEVEL, backendPolicy);
    }

    /** Fits pedigree-aware sparse REML mediation with an explicit confidence level. */
    public static MediationMixedResult fitPedigree(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            List<PedigreeRandomEffectTerm> pedigreeEffects,
            List<RandomEffectTerm> ordinaryEffects,
            RemlOptions options,
            double confidenceLevel,
            BackendPolicy backendPolicy) {
        Objects.requireNonNull(pedigreeEffects, "pedigreeEffects");
        if (pedigreeEffects.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one pedigree random-effect term is required");
        }
        Objects.requireNonNull(options, "options");
        Objects.requireNonNull(backendPolicy, "backendPolicy");
        validateConfidenceLevel(confidenceLevel);
        PreparedData data = prepare(
            outcome, treatment, mediator, covariates,
            MissingDataPolicy.ERROR);
        List<PedigreeRandomEffectTerm> compactPedigree = compactPedigreeTerms(
            pedigreeEffects, data);
        List<RandomEffectTerm> compactOrdinary = ordinaryEffects == null
            ? List.of() : compactTerms(ordinaryEffects, data, "random-effect");
        try (SparseLinearMixedModel.Prepared prepared =
                SparsePedigreeMixedModel.prepare(
                    data.outcome().length, compactPedigree, compactOrdinary,
                    options, backendPolicy)) {
            return fitPreparedMixed(
                data, prepared, confidenceLevel);
        }
    }

    private static MediationMixedResult fitPreparedMixed(
            PreparedData data,
            SparseLinearMixedModel.Prepared prepared,
            double confidenceLevel) {
        double[] mediatorDesign = designRowMajor(
            data.treatment(), data.mediator(), data.covariates(), false);
        double[] outcomeDesign = designRowMajor(
            data.treatment(), data.mediator(), data.covariates(), true);
        double[] totalDesign = designRowMajor(
            data.treatment(), null, data.covariates(), false);
        int covariates = data.covariates() == null
            ? 0 : data.covariates()[0].length;
        int mediatorColumns = 2 + covariates;
        int outcomeColumns = mediatorColumns + 1;

        SparseLinearMixedModelResult mediatorModel = prepared.fit(
            data.mediator(), mediatorDesign, mediatorColumns);
        prepared.warmStart(mediatorModel.varianceComponents());
        SparseLinearMixedModelResult outcomeModel = prepared.fit(
            data.outcome(), outcomeDesign, outcomeColumns);
        prepared.warmStart(outcomeModel.varianceComponents());
        SparseLinearMixedModelResult totalModel = prepared.fit(
            data.outcome(), totalDesign, mediatorColumns);

        MediationEffect aPath = mixedCoefficientEffect(
            "a", mediatorModel, 1, confidenceLevel);
        MediationEffect bPath = mixedCoefficientEffect(
            "b", outcomeModel, 2, confidenceLevel);
        MediationEffect directEffect = mixedCoefficientEffect(
            "direct", outcomeModel, 1, confidenceLevel);
        MediationEffect totalEffect = mixedCoefficientEffect(
            "total", totalModel, 1, confidenceLevel);
        MediationEffect indirectEffect = indirectEffect(
            aPath.estimate(), bPath.estimate(),
            diagonal(mediatorModel.fixedEffectCovariance(), mediatorColumns, 1),
            diagonal(outcomeModel.fixedEffectCovariance(), outcomeColumns, 2),
            confidenceLevel);
        return new MediationMixedResult(
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

    static MediationEffect indirectEffect(
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

    private static MediationEffect mixedCoefficientEffect(
            String name, SparseLinearMixedModelResult fit, int column,
            double confidenceLevel) {
        double[] degrees = fit.associationStatistics().degreesOfFreedom();
        double df = degrees[column];
        double estimate = fit.beta()[column];
        double standardError = fit.standardErrors()[column];
        double statistic = fit.tStatistics()[column];
        double critical = jdistlib.T.quantile(
            0.5 + confidenceLevel / 2.0, df, true, false);
        return new MediationEffect(
            name, estimate, standardError, statistic,
            fit.pValues()[column],
            estimate - critical * standardError,
            estimate + critical * standardError, df);
    }

    private static double diagonal(double[] covariance, int columns, int index) {
        return covariance[index * columns + index];
    }

    private static void validateConfidenceLevel(double confidenceLevel) {
        if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0)
                || !Double.isFinite(confidenceLevel)) {
            throw new IllegalArgumentException(
                "confidenceLevel must be strictly between zero and one");
        }
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

    private static double[] designRowMajor(
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            boolean includeMediator) {
        double[][] matrix = design(
            treatment, mediator, covariates, includeMediator);
        double[] result = new double[matrix.length * matrix[0].length];
        for (int row = 0; row < matrix.length; row++) {
            System.arraycopy(matrix[row], 0, result,
                row * matrix[row].length, matrix[row].length);
        }
        return result;
    }

    private static List<RandomEffectTerm> compactTerms(
            List<RandomEffectTerm> terms,
            PreparedData data,
            String kind) {
        List<RandomEffectTerm> result = new ArrayList<>(terms.size());
        for (RandomEffectTerm term : terms) {
            if (term == null || term.observations() != data.originalObservations()) {
                throw new IllegalArgumentException(
                    kind + " rows must equal mediation observations");
            }
            int rows = data.retainedRows().length;
            int[] retained = data.retainedRows();
            if (rows == data.originalObservations()) {
                result.add(term);
                continue;
            }
            if (term.sparse()) {
                int[] originalStarts = term.rowPointers();
                int[] originalColumns = term.columnIndices();
                double[] originalValues = term.sparseValues();
                int[] starts = new int[rows + 1];
                List<Integer> columns = new ArrayList<>();
                List<Double> values = new ArrayList<>();
                for (int row = 0; row < rows; row++) {
                    starts[row] = columns.size();
                    int source = retained[row];
                    for (int index = originalStarts[source];
                            index < originalStarts[source + 1]; index++) {
                        columns.add(originalColumns[index]);
                        values.add(originalValues[index]);
                    }
                }
                starts[rows] = columns.size();
                int[] compactColumns = columns.stream().mapToInt(Integer::intValue).toArray();
                double[] compactValues = values.stream().mapToDouble(Double::doubleValue).toArray();
                result.add(RandomEffectTerm.ofSparseCsr(
                    term.name(), rows, term.coefficients(), starts,
                    compactColumns, compactValues, term.coefficientNames()));
            } else {
                double[] original = term.design();
                double[] compact = new double[rows * term.coefficients()];
                for (int row = 0; row < rows; row++) {
                    System.arraycopy(original, retained[row] * term.coefficients(),
                        compact, row * term.coefficients(), term.coefficients());
                }
                result.add(RandomEffectTerm.of(
                    term.name(), compact, rows, term.coefficients(),
                    term.coefficientNames()));
            }
        }
        return List.copyOf(result);
    }

    private static List<PedigreeRandomEffectTerm> compactPedigreeTerms(
            List<PedigreeRandomEffectTerm> terms, PreparedData data) {
        List<PedigreeRandomEffectTerm> result = new ArrayList<>(terms.size());
        for (PedigreeRandomEffectTerm term : terms) {
            if (term == null) {
                throw new IllegalArgumentException(
                    "pedigree random-effect terms must not be null");
            }
            RandomEffectTerm compact = compactTerms(
                List.of(term.randomEffect()), data, "pedigree random-effect")
                .get(0);
            result.add(new PedigreeRandomEffectTerm(compact, term.precision()));
        }
        return List.copyOf(result);
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
