/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Cox model with Gaussian frailties and Laplace-profiled variances. */
public final class CoxMixedModel {
    private CoxMixedModel() { }

    public static CoxMixedResult fit(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<CoxRandomEffectTerm> randomEffects) {
        return fit(survival, fixedEffects, randomEffects, null,
            CoxMixedOptions.defaults(), BackendPolicy.PREFERRED);
    }

    public static CoxMixedResult fit(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<CoxRandomEffectTerm> randomEffects,
            double[] offset,
            CoxMixedOptions options,
            BackendPolicy backendPolicy) {
        if (survival == null || fixedEffects == null
                || fixedEffects.length == 0 || randomEffects == null
                || randomEffects.isEmpty() || options == null
                || backendPolicy == null)
            throw new IllegalArgumentException("mixed Cox inputs are required");
        int rows = survival.observations();
        if (fixedEffects.length != rows || fixedEffects[0] == null
                || fixedEffects[0].length == 0)
            throw new IllegalArgumentException(
                "fixed-effect rows must match survival observations");
        if (options.initialVariances().length != randomEffects.size())
            throw new IllegalArgumentException(
                "one initial variance is required per Cox frailty term");
        validateTerms(randomEffects, rows);
        int fixedColumns = fixedEffects[0].length;
        double[] fixed = MatrixOps.rowMajor(fixedEffects, rows);
        validateFixedColumns(fixed, rows, fixedColumns);
        double[] modelOffset = offset == null ? new double[rows]
            : MatrixOps.finiteCopy(offset, "offset");
        if (modelOffset.length != rows)
            throw new IllegalArgumentException(
                "one Cox offset is required per observation");

        Combined combined = combine(fixed, rows, fixedColumns, randomEffects);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] logVariances = java.util.Arrays.stream(
                options.initialVariances()).map(Math::log).toArray();
            ModeFit current = mode(survival, combined, randomEffects,
                variances(logVariances), modelOffset, options,
                new double[combined.columns()], backend);
            double logStep = 1;
            boolean varianceConverged = false;
            int varianceIterations = 0;
            for (int iteration = 1;
                    iteration <= options.maximumVarianceIterations(); iteration++) {
                varianceIterations = iteration;
                boolean improved = false;
                for (int term = 0; term < randomEffects.size(); term++) {
                    ModeFit best = current;
                    double bestLog = logVariances[term];
                    for (double direction : new double[] {-1, 1}) {
                        double candidateLog = clamp(logVariances[term]
                            + direction * logStep,
                            Math.log(options.minimumVariance()),
                            Math.log(options.maximumVariance()));
                        if (candidateLog == logVariances[term]) continue;
                        double[] candidateLogs = logVariances.clone();
                        candidateLogs[term] = candidateLog;
                        ModeFit candidate = mode(survival, combined,
                            randomEffects, variances(candidateLogs), modelOffset,
                            options, current.coefficients(), backend);
                        if (candidate.laplaceLogLikelihood()
                                > best.laplaceLogLikelihood() + 1e-9) {
                            best = candidate;
                            bestLog = candidateLog;
                        }
                    }
                    if (best != current) {
                        current = best;
                        logVariances[term] = bestLog;
                        improved = true;
                    }
                }
                if (!improved) logStep *= 0.5;
                if (logStep <= options.logVarianceTolerance()) {
                    varianceConverged = true;
                    break;
                }
            }
            double[] jointCovariance = CoxMath.inversePositive(backend,
                current.penalizedInformation(), combined.columns(),
                options.coxOptions().informationRidge());
            double[] fixedCovariance = new double[fixedColumns * fixedColumns];
            for (int row = 0; row < fixedColumns; row++)
                for (int column = 0; column < fixedColumns; column++)
                    fixedCovariance[row * fixedColumns + column] =
                        jointCovariance[row * combined.columns() + column];
            double[] beta = java.util.Arrays.copyOf(
                current.coefficients(), fixedColumns);
            List<CoxRandomEffectEstimates> estimates = new ArrayList<>();
            double[] fittedVariances = variances(logVariances);
            int start = fixedColumns;
            for (int term = 0; term < randomEffects.size(); term++) {
                CoxRandomEffectTerm value = randomEffects.get(term);
                estimates.add(new CoxRandomEffectEstimates(value.name(),
                    value.coefficientNames(), fittedVariances[term],
                    java.util.Arrays.copyOfRange(current.coefficients(),
                        start, start + value.coefficients())));
                start += value.coefficients();
            }
            boolean converged = varianceConverged && current.converged();
            String message = converged
                ? "Laplace variance and conditional-mode tolerances reached"
                : (!current.converged() ? current.message()
                    : "maximum frailty-variance iterations reached");
            return new CoxMixedResult(beta, fixedCovariance, estimates,
                CoxPartialLikelihood.baseline(survival, combined.design(),
                    combined.columns(), current.coefficients(), modelOffset,
                    options.coxOptions().ties()),
                current.partialLogLikelihood(),
                current.penalizedLogLikelihood(),
                current.laplaceLogLikelihood(), options,
                varianceIterations, converged, message,
                context.provenance());
        }
    }

    private static ModeFit mode(
            CoxSurvivalData survival,
            Combined combined,
            List<CoxRandomEffectTerm> terms,
            double[] variances,
            double[] offset,
            CoxMixedOptions options,
            double[] initial,
            ComputeBackend backend) {
        double[] coefficients = initial.clone();
        PenalizedEvaluation evaluation = penalized(survival, combined,
            terms, variances, coefficients, offset, options, backend);
        boolean converged = false;
        String message = "maximum conditional-mode iterations reached";
        for (int iteration = 0;
                iteration < options.coxOptions().maximumIterations(); iteration++) {
            if (CoxMath.maximumAbsolute(evaluation.score())
                    <= options.coxOptions().scoreTolerance()) {
                converged = true;
                message = "conditional-mode score tolerance reached";
                break;
            }
            double[] step = CoxMath.solvePositive(backend,
                evaluation.information(), combined.columns(),
                evaluation.score(), options.coxOptions().informationRidge());
            double scale = 1;
            PenalizedEvaluation candidate = null;
            double[] next = null;
            for (int halving = 0;
                    halving <= options.coxOptions().maximumStepHalvings(); halving++) {
                next = coefficients.clone();
                for (int column = 0; column < next.length; column++)
                    next[column] += scale * step[column];
                candidate = penalized(survival, combined, terms, variances,
                    next, offset, options, backend);
                if (candidate.penalizedLogLikelihood()
                        >= evaluation.penalizedLogLikelihood() - 1e-12) break;
                scale *= 0.5;
            }
            if (candidate.penalizedLogLikelihood()
                    < evaluation.penalizedLogLikelihood() - 1e-10) {
                message = "conditional-mode step halving failed";
                break;
            }
            double relative = 0;
            for (int column = 0; column < next.length; column++)
                relative = Math.max(relative,
                    Math.abs(next[column] - coefficients[column])
                        / (1 + Math.abs(coefficients[column])));
            coefficients = next;
            evaluation = candidate;
            if (relative <= options.coxOptions().relativeTolerance()
                    && CoxMath.maximumAbsolute(evaluation.score())
                        <= Math.sqrt(options.coxOptions().scoreTolerance())) {
                converged = true;
                message = "conditional-mode coefficient/score tolerances reached";
                break;
            }
        }
        double laplace = laplace(evaluation, combined.fixedColumns(), terms,
            variances, options, backend);
        return new ModeFit(coefficients,
            evaluation.partialLogLikelihood(),
            evaluation.penalizedLogLikelihood(), laplace,
            evaluation.information(), converged, message);
    }

    private static PenalizedEvaluation penalized(
            CoxSurvivalData survival, Combined combined,
            List<CoxRandomEffectTerm> terms, double[] variances,
            double[] coefficients, double[] offset,
            CoxMixedOptions options, ComputeBackend backend) {
        CoxPartialLikelihood.Evaluation partial =
            CoxPartialLikelihood.evaluate(survival, combined.design(),
                combined.columns(), coefficients, offset,
                options.coxOptions().ties());
        double[] score = partial.score().clone();
        double[] information = partial.information().clone();
        double penalty = 0;
        int start = combined.fixedColumns();
        for (int term = 0; term < terms.size(); term++) {
            CoxRandomEffectTerm value = terms.get(term);
            int dimension = value.coefficients();
            double inverseVariance = 1 / variances[term];
            for (int row = 0; row < dimension; row++) {
                double product = 0;
                for (int column = 0; column < dimension; column++) {
                    double precision = inverseVariance
                        * value.precisionView()[row * dimension + column];
                    product += precision * coefficients[start + column];
                    information[(start + row) * combined.columns()
                        + start + column] += precision;
                }
                score[start + row] -= product;
                penalty += 0.5 * coefficients[start + row] * product;
            }
            start += dimension;
        }
        return new PenalizedEvaluation(partial.logLikelihood(),
            partial.logLikelihood() - penalty, score, information);
    }

    private static double laplace(
            PenalizedEvaluation evaluation, int fixedColumns,
            List<CoxRandomEffectTerm> terms, double[] variances,
            CoxMixedOptions options, ComputeBackend backend) {
        int randomColumns = terms.stream()
            .mapToInt(CoxRandomEffectTerm::coefficients).sum();
        double[] randomInformation = new double[randomColumns * randomColumns];
        int totalColumns = fixedColumns + randomColumns;
        for (int row = 0; row < randomColumns; row++)
            for (int column = 0; column < randomColumns; column++)
                randomInformation[row * randomColumns + column] =
                    evaluation.information()[(fixedColumns + row)
                        * totalColumns + fixedColumns + column];
        double logDeterminantInformation = CoxMath.factor(backend,
            randomInformation, randomColumns,
            options.coxOptions().informationRidge()).logDeterminant();
        double logDeterminantPrecision = 0;
        for (int term = 0; term < terms.size(); term++) {
            CoxRandomEffectTerm value = terms.get(term);
            logDeterminantPrecision += CoxMath.factor(backend,
                value.precisionView(), value.coefficients(),
                options.coxOptions().informationRidge()).logDeterminant()
                - value.coefficients() * Math.log(variances[term]);
        }
        return evaluation.penalizedLogLikelihood()
            + 0.5 * logDeterminantPrecision
            - 0.5 * logDeterminantInformation;
    }

    private static Combined combine(
            double[] fixed, int rows, int fixedColumns,
            List<CoxRandomEffectTerm> terms) {
        int columns = fixedColumns + terms.stream()
            .mapToInt(CoxRandomEffectTerm::coefficients).sum();
        double[] design = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(fixed, row * fixedColumns,
                design, row * columns, fixedColumns);
            int destination = fixedColumns;
            for (CoxRandomEffectTerm term : terms) {
                System.arraycopy(term.designView(), row * term.coefficients(),
                    design, row * columns + destination, term.coefficients());
                destination += term.coefficients();
            }
        }
        return new Combined(design, columns, fixedColumns);
    }

    private static void validateTerms(
            List<CoxRandomEffectTerm> terms, int rows) {
        Set<String> names = new HashSet<>();
        for (CoxRandomEffectTerm term : terms)
            if (term == null || term.observations() != rows
                    || !names.add(term.name()))
                throw new IllegalArgumentException(
                    "frailty terms must be nonnull, named uniquely, and row-aligned");
    }

    private static void validateFixedColumns(
            double[] design, int rows, int columns) {
        for (int column = 0; column < columns; column++) {
            double minimum = Double.POSITIVE_INFINITY;
            double maximum = Double.NEGATIVE_INFINITY;
            for (int row = 0; row < rows; row++) {
                minimum = Math.min(minimum, design[row * columns + column]);
                maximum = Math.max(maximum, design[row * columns + column]);
            }
            if (!(maximum > minimum))
                throw new IllegalArgumentException(
                    "Cox fixed effects must not include an intercept/constant");
        }
    }

    private static double[] variances(double[] logVariances) {
        return java.util.Arrays.stream(logVariances).map(Math::exp).toArray();
    }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    private record Combined(double[] design, int columns, int fixedColumns) { }
    private record PenalizedEvaluation(
        double partialLogLikelihood, double penalizedLogLikelihood,
        double[] score, double[] information) { }
    private record ModeFit(
        double[] coefficients,
        double partialLogLikelihood,
        double penalizedLogLikelihood,
        double laplaceLogLikelihood,
        double[] penalizedInformation,
        boolean converged,
        String message) { }
}
