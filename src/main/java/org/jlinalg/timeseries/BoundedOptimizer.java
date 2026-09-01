/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

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
            return new Result(initial.clone(), objective.eval(initial), 1, true);
        }
        if (initial.length == 1) {
            return goldenSection(
                lower[0], upper[0], objective, maximumEvaluations, tolerance);
        }
        try {
            int interpolationPoints = 2 * initial.length + 1;
            OptimizationResult result = Bobyqa.bobyqa(
                initial, lower, upper, objective, interpolationPoints,
                0.5, tolerance, maximumEvaluations, true);
            boolean converged = result.mX != null
                && result.mX.length == initial.length
                && Double.isFinite(result.mF);
            return new Result(converged ? result.mX : initial.clone(),
                converged ? result.mF : objective.eval(initial),
                result.numFunctionCalls, converged);
        } catch (RuntimeException exception) {
            return coordinateSearch(initial, lower, upper, objective,
                maximumEvaluations, tolerance);
        }
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
            evaluations, right - left <= Math.sqrt(tolerance));
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
        return new Result(point, value, evaluations, step <= tolerance);
    }

    record Result(double[] parameters, double objective,
                  int evaluations, boolean converged) { }
}
