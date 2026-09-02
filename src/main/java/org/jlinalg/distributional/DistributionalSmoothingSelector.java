/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PenalizedPredictor;
import org.jlinalg.gam.QuadraticPenalizedPredictor;
import org.jlinalg.gam.QuadraticSmoothTerm;

/** AIC outer selection of every parameter-specific quadratic smooth penalty. */
public final class DistributionalSmoothingSelector {
    private DistributionalSmoothingSelector() { }

    /** Selects independent smoothing parameters for every term and margin. */
    public static DistributionalSmoothingResult fit(
            double[] response,
            List<double[][]> parametricDesigns,
            List<List<QuadraticSmoothTerm>> smoothTerms,
            List<List<double[]>> initialSmoothing,
            DistributionalFamily family,
            DistributionalOptions fitOptions,
            DistributionalSmoothingOptions smoothingOptions,
            BackendPolicy backendPolicy) {
        validate(response, parametricDesigns, smoothTerms, family,
            fitOptions, smoothingOptions, backendPolicy);
        List<List<double[]>> logs = logs(initialSmoothing, smoothTerms);
        List<Coordinate> coordinates = coordinates(logs);
        Candidate best = evaluate(response, parametricDesigns, smoothTerms,
            exponentiate(logs), family, fitOptions, backendPolicy, null);
        if (best.fit() == null) {
            throw new IllegalArgumentException("initial distributional smoothing fit failed");
        }
        int evaluations = 1;
        double step = smoothingOptions.initialLogStep();
        for (int sweep = 0; sweep < smoothingOptions.maximumSweeps(); sweep++) {
            boolean improved = false;
            for (Coordinate coordinate : coordinates) {
                double[] vector = logs.get(coordinate.parameter())
                    .get(coordinate.term());
                double original = vector[coordinate.penalty()];
                double selected = original;
                Candidate coordinateBest = best;
                for (double direction : new double[] {-1.0, 1.0}) {
                    double trial = Math.max(smoothingOptions.minimumLogSmoothing(),
                        Math.min(smoothingOptions.maximumLogSmoothing(),
                            original + direction * step));
                    if (trial == original) continue;
                    vector[coordinate.penalty()] = trial;
                    Candidate candidate = evaluate(response, parametricDesigns,
                        smoothTerms, exponentiate(logs), family, fitOptions,
                        backendPolicy, best.fit());
                    evaluations++;
                    if (candidate.aic() < coordinateBest.aic()) {
                        coordinateBest = candidate;
                        selected = trial;
                    }
                }
                vector[coordinate.penalty()] = selected;
                if (coordinateBest.aic() < best.aic()
                        - smoothingOptions.tolerance()
                            * (1.0 + Math.abs(best.aic()))) {
                    best = coordinateBest;
                    improved = true;
                }
            }
            if (!improved) step *= 0.5;
            if (step <= smoothingOptions.tolerance()) break;
        }
        List<List<double[]>> smoothing = exponentiate(logs);
        best = evaluate(response, parametricDesigns, smoothTerms, smoothing,
            family, fitOptions, backendPolicy, best.fit());
        evaluations++;
        return new DistributionalSmoothingResult(best.fit(), best.predictors(),
            smoothing, best.aic(), evaluations);
    }

    private static Candidate evaluate(
            double[] response,
            List<double[][]> parametricDesigns,
            List<List<QuadraticSmoothTerm>> smoothTerms,
            List<List<double[]>> smoothing,
            DistributionalFamily family,
            DistributionalOptions options,
            BackendPolicy backendPolicy,
            DistributionalResult warm) {
        try {
            List<PenalizedPredictor> predictors = new ArrayList<>(family.parameterCount());
            for (int parameter = 0; parameter < family.parameterCount(); parameter++) {
                if (smoothTerms.get(parameter).isEmpty()) {
                    predictors.add(PenalizedPredictor.linear(
                        parametricDesigns.get(parameter)));
                } else {
                    predictors.add(QuadraticPenalizedPredictor.compile(
                        parametricDesigns.get(parameter), smoothTerms.get(parameter),
                        smoothing.get(parameter), backendPolicy));
                }
            }
            List<double[]> starting = null;
            if (warm != null) {
                starting = new ArrayList<>(warm.parameters().size());
                for (DistributionalParameterResult parameter : warm.parameters()) {
                    starting.add(parameter.coefficients());
                }
                for (int parameter = 0; parameter < predictors.size(); parameter++) {
                    if (starting.get(parameter).length != predictors.get(parameter).columns()) {
                        starting = null;
                        break;
                    }
                }
            }
            DistributionalResult fit = DistributionalModel.fit(response,
                predictors, family, options, backendPolicy, starting);
            double edf = 0.0;
            for (DistributionalParameterResult parameter : fit.parameters()) {
                edf += parameter.effectiveDegreesOfFreedom();
            }
            double aic = fit.converged()
                ? -2.0 * fit.logLikelihood() + 2.0 * edf
                : Double.POSITIVE_INFINITY;
            return new Candidate(fit, List.copyOf(predictors), aic);
        } catch (IllegalArgumentException | IllegalStateException exception) {
            return new Candidate(null, List.of(), Double.POSITIVE_INFINITY);
        }
    }

    private static List<List<double[]>> logs(
            List<List<double[]>> supplied,
            List<List<QuadraticSmoothTerm>> terms) {
        List<List<double[]>> result = new ArrayList<>(terms.size());
        for (int parameter = 0; parameter < terms.size(); parameter++) {
            List<double[]> parameterResult = new ArrayList<>(terms.get(parameter).size());
            if (supplied != null && supplied.size() != terms.size()) {
                throw new IllegalArgumentException("one smoothing collection is required per parameter");
            }
            for (int term = 0; term < terms.get(parameter).size(); term++) {
                int count = terms.get(parameter).get(term).penaltyCount();
                double[] values = new double[count];
                if (supplied != null) {
                    if (supplied.get(parameter).size() != terms.get(parameter).size()
                            || supplied.get(parameter).get(term).length != count) {
                        throw new IllegalArgumentException("initial smoothing dimensions are invalid");
                    }
                    for (int penalty = 0; penalty < count; penalty++) {
                        double value = supplied.get(parameter).get(term)[penalty];
                        if (!(value > 0.0) || !Double.isFinite(value)) {
                            throw new IllegalArgumentException("smoothing parameters must be positive");
                        }
                        values[penalty] = Math.log(value);
                    }
                }
                parameterResult.add(values);
            }
            result.add(parameterResult);
        }
        return result;
    }

    private static List<Coordinate> coordinates(List<List<double[]>> values) {
        List<Coordinate> result = new ArrayList<>();
        for (int parameter = 0; parameter < values.size(); parameter++) {
            for (int term = 0; term < values.get(parameter).size(); term++) {
                for (int penalty = 0; penalty < values.get(parameter).get(term).length;
                        penalty++) {
                    result.add(new Coordinate(parameter, term, penalty));
                }
            }
        }
        return List.copyOf(result);
    }

    private static List<List<double[]>> exponentiate(List<List<double[]>> logs) {
        List<List<double[]>> result = new ArrayList<>(logs.size());
        for (List<double[]> parameter : logs) {
            List<double[]> terms = new ArrayList<>(parameter.size());
            for (double[] term : parameter) {
                double[] values = new double[term.length];
                for (int penalty = 0; penalty < term.length; penalty++) {
                    values[penalty] = Math.exp(term[penalty]);
                }
                terms.add(values);
            }
            result.add(List.copyOf(terms));
        }
        return List.copyOf(result);
    }

    private static void validate(
            double[] response,
            List<double[][]> parametricDesigns,
            List<List<QuadraticSmoothTerm>> smoothTerms,
            DistributionalFamily family,
            DistributionalOptions fitOptions,
            DistributionalSmoothingOptions smoothingOptions,
            BackendPolicy backendPolicy) {
        if (response == null || parametricDesigns == null || smoothTerms == null
                || family == null || fitOptions == null || smoothingOptions == null
                || backendPolicy == null
                || parametricDesigns.size() != family.parameterCount()
                || smoothTerms.size() != family.parameterCount()) {
            throw new IllegalArgumentException("invalid distributional smoothing inputs");
        }
        for (double[][] design : parametricDesigns) {
            if (design == null || design.length != response.length) {
                throw new IllegalArgumentException("parametric designs must match response");
            }
        }
    }

    private record Coordinate(int parameter, int term, int penalty) { }
    private record Candidate(
            DistributionalResult fit,
            List<PenalizedPredictor> predictors,
            double aic) { }
}
