/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import java.util.ArrayList;
import java.util.TreeSet;
import jdistlib.math.MultivariableFunction;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.OptimizationResult;

/** Small bounded-optimization adapter around JDistlib BOBYQA. */
final class BoundedOptimizer {
    private static final double GOLDEN_RATIO = 0.6180339887498949;

    private BoundedOptimizer() { }

    static Result minimize(
            double[] initial,
            double[] lower,
            double[] upper,
            MultivariableFunction objective,
            int maximumEvaluations,
            double tolerance) {
        if (initial.length == 0) {
            double value = objective.eval(initial);
            return new Result(initial.clone(), value, 1, Double.isFinite(value) && value < Double.MAX_VALUE / 16);
        }
        if (initial.length == 1) {
            return scanAndRefine(initial[0], lower[0], upper[0], objective, maximumEvaluations, tolerance);
        }
        try {
            int interpolationPoints = 2 * initial.length + 1;
            OptimizationResult result = Bobyqa.bobyqa(
                initial, lower, upper, objective, interpolationPoints,
                0.5, tolerance, maximumEvaluations, true);
            boolean valid = result.mX != null
                && result.mX.length == initial.length
                && Double.isFinite(result.mF) && result.mF < Double.MAX_VALUE / 16;
            // OptimizationResult.isMinimum selects the optimization direction,
            // not successful termination. A finite result alone is not convergence.
            boolean check = valid && result.numFunctionCalls + 2 * initial.length < maximumEvaluations;
            boolean converged = check && stationary(result.mX, lower, upper, objective, result.mF, tolerance);
            return new Result(valid ? result.mX : initial.clone(),
                valid ? result.mF : objective.eval(initial),
                result.numFunctionCalls + (check ? 2 * initial.length : 0), converged);
        } catch (RuntimeException exception) {
            return coordinateSearch(initial, lower, upper, objective,
                maximumEvaluations, tolerance);
        }
    }

    /** Explore the bounded interval before refining every sampled basin. ARMA
     * likelihoods need not be unimodal. Always retain the caller's seed and
     * boundary values; a worse local solution cannot replace them. */
    private static Result scanAndRefine(double initial, double lower, double upper,
            MultivariableFunction objective, int maximumEvaluations, double tolerance) {
        int divisions = Math.min(64, Math.max(2, maximumEvaluations / 4));
        TreeSet<Double> grid = new TreeSet<>();
        grid.add(Math.max(lower, Math.min(upper, initial)));
        for (int i = 0; i <= divisions; i++) grid.add(lower + (upper - lower) * i / divisions);
        double[] points = grid.stream().mapToDouble(Double::doubleValue).toArray();
        double[] values = new double[points.length];
        double bestPoint = points[0], bestValue = Double.POSITIVE_INFINITY;
        int evaluations = 0;
        for (int i = 0; i < points.length; i++) {
            values[i] = objective.eval(new double[] {points[i]}); evaluations++;
            if (values[i] < bestValue) { bestPoint = points[i]; bestValue = values[i]; }
        }
        ArrayList<Integer> basins = new ArrayList<>();
        for (int i = 0; i < points.length; i++) {
            if ((i == 0 || values[i] <= values[i - 1])
                    && (i == points.length - 1 || values[i] <= values[i + 1])) basins.add(i);
        }
        boolean resolved = true;
        for (int basin = 0; basin < basins.size(); basin++) {
            int budget = (maximumEvaluations - evaluations) / (basins.size() - basin);
            if (budget < 3) { resolved = false; break; }
            int i = basins.get(basin);
            Result candidate = goldenSection(points[Math.max(0, i - 1)],
                points[Math.min(points.length - 1, i + 1)], objective, budget, tolerance);
            evaluations += candidate.evaluations();
            resolved &= candidate.converged();
            if (candidate.objective() < bestValue) {
                bestValue = candidate.objective(); bestPoint = candidate.parameters()[0];
            }
        }
        boolean check = resolved && Double.isFinite(bestValue) && bestValue < Double.MAX_VALUE / 16
            && evaluations + 2 <= maximumEvaluations;
        boolean converged = check && stationary(new double[] {bestPoint}, new double[] {lower},
            new double[] {upper}, objective, bestValue, tolerance);
        return new Result(new double[] {bestPoint}, bestValue, evaluations + (check ? 2 : 0), converged);
    }

    private static boolean stationary(double[] point, double[] lower, double[] upper,
            MultivariableFunction objective, double value, double tolerance) {
        double threshold = Math.max(1e-5, Math.sqrt(tolerance)) * (1 + Math.sqrt(Math.abs(value)));
        boolean result = true;
        for (int i = 0; i < point.length; i++) {
            double h = 1e-5 * (1 + Math.abs(point[i]));
            double[] plus = point.clone(), minus = point.clone();
            plus[i] = Math.min(upper[i], point[i] + h);
            minus[i] = Math.max(lower[i], point[i] - h);
            double fp = objective.eval(plus), fm = objective.eval(minus);
            double score = (fp - fm) / (plus[i] - minus[i]);
            if (point[i] <= lower[i] + h) score = Math.min(0, score);
            if (point[i] >= upper[i] - h) score = Math.max(0, score);
            result &= Double.isFinite(score) && Math.abs(score) <= threshold;
        }
        return result;
    }

    private static Result goldenSection(
            double lower,
            double upper,
            MultivariableFunction objective,
            int maximumEvaluations,
            double tolerance) {
        double left = lower;
        double right = upper;
        double x1 = right - GOLDEN_RATIO * (right - left);
        double x2 = left + GOLDEN_RATIO * (right - left);
        double f1 = objective.eval(x1);
        double f2 = objective.eval(x2);
        int evaluations = 2;
        while (evaluations < maximumEvaluations
                && right - left > tolerance * (1.0 + Math.abs(left + right))) {
            if (f1 <= f2) {
                right = x2;
                x2 = x1;
                f2 = f1;
                x1 = right - GOLDEN_RATIO * (right - left);
                f1 = objective.eval(x1);
            } else {
                left = x1;
                x1 = x2;
                f1 = f2;
                x2 = left + GOLDEN_RATIO * (right - left);
                f2 = objective.eval(x2);
            }
            evaluations++;
        }
        boolean first = f1 <= f2;
        return new Result(new double[] {first ? x1 : x2}, first ? f1 : f2,
            evaluations, Double.isFinite(first ? f1 : f2)
                && (first ? f1 : f2) < Double.MAX_VALUE / 16
                && right - left <= tolerance * (1.0 + Math.abs(left + right)));
    }

    private static Result coordinateSearch(
            double[] initial,
            double[] lower,
            double[] upper,
            MultivariableFunction objective,
            int maximumEvaluations,
            double tolerance) {
        double[] point = initial.clone();
        double value = objective.eval(point);
        int evaluations = 1;
        double step = 0.5;
        while (evaluations < maximumEvaluations && step > tolerance) {
            boolean improved = false;
            for (int dimension = 0;
                    dimension < point.length
                        && evaluations < maximumEvaluations; dimension++) {
                double original = point[dimension];
                for (int direction : new int[] {-1, 1}) {
                    if (evaluations >= maximumEvaluations) break;
                    point[dimension] = Math.max(lower[dimension],
                        Math.min(upper[dimension], original + direction * step));
                    double candidate = objective.eval(point);
                    evaluations++;
                    if (candidate < value) {
                        value = candidate;
                        original = point[dimension];
                        improved = true;
                    } else {
                        point[dimension] = original;
                    }
                }
            }
            if (!improved) {
                step *= 0.5;
            }
        }
        boolean check = step <= tolerance && Double.isFinite(value) && value < Double.MAX_VALUE / 16
            && evaluations + 2 * point.length <= maximumEvaluations;
        boolean converged = check && stationary(point, lower, upper, objective, value, tolerance);
        return new Result(point, value, evaluations + (check ? 2 * point.length : 0), converged);
    }

    record Result(double[] parameters, double objective,
                  int evaluations, boolean converged) { }
}
