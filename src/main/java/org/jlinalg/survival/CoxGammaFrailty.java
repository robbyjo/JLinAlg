/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.Gamma;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Cox model with multiplicative shared gamma frailty and Laplace integration. */
public final class CoxGammaFrailty {
    private CoxGammaFrailty() { }

    public static CoxGammaFrailtyResult fit(CoxSurvivalData survival,
            double[][] fixedEffects, List<String> groupIds) {
        return fit(survival, fixedEffects, groupIds, null,
            CoxGammaFrailtyOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Profiles the gamma variance on the log scale. */
    public static CoxGammaFrailtyResult fit(CoxSurvivalData survival,
            double[][] fixedEffects, List<String> groupIds, double[] offset,
            CoxGammaFrailtyOptions options, BackendPolicy backendPolicy) {
        Prepared prepared = prepare(survival, fixedEffects, groupIds, offset,
            options);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            double logTheta = Math.log(options.initialVariance());
            Mode current = mode(prepared, Math.exp(logTheta), null,
                context.backend());
            double step = 1.0;
            boolean varianceConverged = false;
            int iterations = 0;
            for (int iteration = 1;
                    iteration <= options.maximumVarianceIterations(); iteration++) {
                iterations = iteration;
                double lower = clamp(logTheta - step,
                    Math.log(options.minimumVariance()),
                    Math.log(options.maximumVariance()));
                double upper = clamp(logTheta + step,
                    Math.log(options.minimumVariance()),
                    Math.log(options.maximumVariance()));
                Mode best = current;
                double selected = logTheta;
                if (lower != logTheta) {
                    Mode candidate = mode(prepared, Math.exp(lower), current,
                        context.backend());
                    if (candidate.laplace() > best.laplace()) {
                        best = candidate;
                        selected = lower;
                    }
                }
                if (upper != logTheta) {
                    Mode candidate = mode(prepared, Math.exp(upper), current,
                        context.backend());
                    if (candidate.laplace() > best.laplace()) {
                        best = candidate;
                        selected = upper;
                    }
                }
                if (selected == logTheta) step *= 0.5;
                else {
                    logTheta = selected;
                    current = best;
                }
                if (step <= options.logVarianceTolerance()) {
                    varianceConverged = true;
                    break;
                }
            }
            return result(prepared, Math.exp(logTheta), current, iterations,
                varianceConverged, context);
        }
    }

    /** Fits conditional modes at one caller-fixed gamma variance. */
    public static CoxGammaFrailtyResult fitAtVariance(
            CoxSurvivalData survival, double[][] fixedEffects,
            List<String> groupIds, double[] offset, double variance,
            CoxGammaFrailtyOptions options, BackendPolicy backendPolicy) {
        if (options == null || variance < options.minimumVariance()
                || variance > options.maximumVariance())
            throw new IllegalArgumentException(
                "gamma variance must lie within configured bounds");
        Prepared prepared = prepare(survival, fixedEffects, groupIds, offset,
            options);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            Mode mode = mode(prepared, variance, null, context.backend());
            return result(prepared, variance, mode, 0, true, context);
        }
    }

    private static CoxGammaFrailtyResult result(Prepared prepared,
            double theta, Mode mode, int varianceIterations,
            boolean varianceConverged, BackendContext context) {
        int fixed = prepared.fixedColumns();
        int total = fixed + prepared.groups().size();
        double[] joint = CoxMath.inversePositive(context.backend(),
            mode.information(), total,
            prepared.options().coxOptions().informationRidge());
        double[] covariance = new double[fixed * fixed];
        for (int row = 0; row < fixed; row++)
            for (int column = 0; column < fixed; column++)
                covariance[row * fixed + column] =
                    joint[row * total + column];
        double[] beta = java.util.Arrays.copyOf(mode.coefficients(), fixed);
        double[] logFrailties = java.util.Arrays.copyOfRange(
            mode.coefficients(), fixed, total);
        boolean converged = varianceConverged && mode.converged();
        return new CoxGammaFrailtyResult(beta, covariance, theta,
            prepared.groups(), logFrailties,
            CoxPartialLikelihood.baseline(prepared.survival(),
                prepared.design(), total, mode.coefficients(),
                prepared.offset(), prepared.options().coxOptions().ties(),
                prepared.rightPlan(), prepared.countingPlan()),
            mode.partial(), mode.laplace(), varianceIterations, converged,
            converged
                ? "gamma variance and conditional-mode tolerances reached"
                : (!mode.converged() ? mode.message()
                    : "maximum gamma-variance iterations reached"),
            context.provenance());
    }

    private static Mode mode(Prepared prepared, double theta, Mode start,
            ComputeBackend backend) {
        int fixed = prepared.fixedColumns();
        int random = prepared.groups().size();
        int total = fixed + random;
        double[] coefficients = start == null
            ? new double[total] : start.coefficients().clone();
        Penalized evaluation = evaluate(prepared, coefficients, theta);
        boolean converged = false;
        String message = "maximum gamma conditional-mode iterations reached";
        for (int iteration = 0; iteration
                < prepared.options().coxOptions().maximumIterations(); iteration++) {
            if (CoxMath.maximumAbsolute(evaluation.score())
                    <= prepared.options().coxOptions().scoreTolerance()) {
                converged = true;
                message = "gamma conditional-mode score tolerance reached";
                break;
            }
            double[] step = CoxMath.solvePositive(backend,
                evaluation.information(), total, evaluation.score(),
                prepared.options().coxOptions().informationRidge());
            double scale = 1.0;
            Penalized candidate = null;
            double[] next = null;
            for (int halving = 0; halving
                    <= prepared.options().coxOptions().maximumStepHalvings();
                    halving++) {
                next = coefficients.clone();
                for (int column = 0; column < total; column++)
                    next[column] += scale * step[column];
                candidate = evaluate(prepared, next, theta);
                if (candidate.logPosterior()
                        >= evaluation.logPosterior() - 1e-12) break;
                scale *= 0.5;
            }
            if (candidate.logPosterior() < evaluation.logPosterior() - 1e-10) {
                message = "gamma conditional-mode step halving failed";
                break;
            }
            double relative = 0.0;
            for (int column = 0; column < total; column++)
                relative = Math.max(relative, Math.abs(next[column]
                    - coefficients[column]) / (1.0 + Math.abs(coefficients[column])));
            coefficients = next;
            evaluation = candidate;
            if (relative <= prepared.options().coxOptions().relativeTolerance()
                    && CoxMath.maximumAbsolute(evaluation.score()) <= Math.sqrt(
                        prepared.options().coxOptions().scoreTolerance())) {
                converged = true;
                message = "gamma conditional-mode coefficient/score tolerances reached";
                break;
            }
        }
        double[] randomInformation = new double[random * random];
        for (int row = 0; row < random; row++)
            for (int column = 0; column < random; column++)
                randomInformation[row * random + column] =
                    evaluation.information()[(fixed + row) * total
                        + fixed + column];
        double logDeterminant = CoxMath.factor(backend, randomInformation,
            random, prepared.options().coxOptions().informationRidge())
            .logDeterminant();
        double laplace = evaluation.logPosterior()
            + 0.5 * random * Math.log(2.0 * Math.PI)
            - 0.5 * logDeterminant;
        return new Mode(coefficients, evaluation.information(),
            evaluation.partial(), laplace, converged, message);
    }

    private static Penalized evaluate(
            Prepared prepared, double[] coefficients, double theta) {
        int fixed = prepared.fixedColumns();
        int total = fixed + prepared.groups().size();
        CoxPartialLikelihood.Evaluation partial =
            CoxPartialLikelihood.evaluate(prepared.survival(),
                prepared.design(), total, coefficients, prepared.offset(),
                prepared.options().coxOptions().ties(), prepared.rightPlan(),
                prepared.countingPlan());
        double[] score = partial.score().clone();
        double[] information = partial.information().clone();
        double logPrior = 0.0;
        double shape = 1.0 / theta;
        for (int group = 0; group < prepared.groups().size(); group++) {
            int index = fixed + group;
            double logFrailty = coefficients[index];
            double frailty = Math.exp(Math.min(700.0, logFrailty));
            logPrior += Gamma.density(frailty, shape, theta, true)
                + logFrailty;
            score[index] += shape - frailty / theta;
            information[index * total + index] += frailty / theta;
        }
        return new Penalized(partial.logLikelihood(),
            partial.logLikelihood() + logPrior, score, information);
    }

    private static Prepared prepare(CoxSurvivalData survival,
            double[][] fixedEffects, List<String> groupIds, double[] offset,
            CoxGammaFrailtyOptions options) {
        if (survival == null || fixedEffects == null || groupIds == null
                || options == null || fixedEffects.length != survival.observations()
                || groupIds.size() != survival.observations())
            throw new IllegalArgumentException(
                "one gamma-frailty group is required per survival row");
        int rows = survival.observations();
        int fixed = fixedEffects[0].length;
        double[] fixedValues = MatrixOps.rowMajor(fixedEffects, rows);
        Map<String, Integer> groupIndex = new LinkedHashMap<>();
        for (String group : groupIds) {
            if (group == null || group.isBlank())
                throw new IllegalArgumentException(
                    "gamma-frailty group IDs must be nonblank");
            groupIndex.computeIfAbsent(group, ignored -> groupIndex.size());
        }
        List<String> groups = new ArrayList<>(groupIndex.keySet());
        int total = fixed + groups.size();
        double[] design = new double[rows * total];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(fixedValues, row * fixed, design,
                row * total, fixed);
            design[row * total + fixed + groupIndex.get(groupIds.get(row))] = 1.0;
        }
        double[] offsets = offset == null ? new double[rows] : offset.clone();
        if (offsets.length != rows)
            throw new IllegalArgumentException("one Cox offset is required per row");
        CoxRiskSetPlan right = CoxRiskSetPlan.prepare(survival);
        return new Prepared(survival, design, fixed, groups, offsets, options,
            right, right == null ? CoxCountingProcessPlan.prepare(survival) : null);
    }

    private static double clamp(double value, double lower, double upper) {
        return Math.max(lower, Math.min(upper, value));
    }

    private record Prepared(CoxSurvivalData survival, double[] design,
        int fixedColumns, List<String> groups, double[] offset,
        CoxGammaFrailtyOptions options, CoxRiskSetPlan rightPlan,
        CoxCountingProcessPlan countingPlan) { }
    private record Penalized(double partial, double logPosterior,
        double[] score, double[] information) { }
    private record Mode(double[] coefficients, double[] information,
        double partial, double laplace, boolean converged, String message) { }
}
