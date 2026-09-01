/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import java.util.Arrays;
import jdistlib.math.MultivariableFunction;
import org.jlinalg.internal.MatrixOps;

/** Conditional Gaussian AR, MA, ARMA, ARIMA, and seasonal ARIMA fitting. */
public final class Arima {
    private static final double TRANSFORM_BOUND = 3.8;
    private static final double MINIMUM_VARIANCE = 1e-14;

    private Arima() { }

    public static ArimaResult fit(double[] series, ArimaOrder order) {
        return fit(series, order, ArimaOptions.defaults());
    }

    public static ArimaResult fit(
            double[] series, ArimaOrder order, ArimaOptions options) {
        double[] values = MatrixOps.finiteCopy(series, "series");
        if (order == null || options == null) {
            throw new IllegalArgumentException("order and options are required");
        }
        SeasonalArimaOrder seasonal = options.seasonalOrder();
        int totalDifferences = order.differences() + seasonal.differences();
        if (options.includeDrift() && totalDifferences != 1) {
            throw new IllegalArgumentException(
                "drift is supported only with exactly one differencing operator");
        }
        boolean locationIncluded = totalDifferences == 0
            ? options.includeMean() : options.includeDrift();
        double[] differenced = ArimaMath.difference(values, order, seasonal);
        int maximumArLag = Math.max(order.autoregressive(),
            seasonal.autoregressive() * seasonal.period()
                + order.autoregressive());
        int parameterCount = order.autoregressive() + order.movingAverage()
            + seasonal.autoregressive() + seasonal.movingAverage()
            + (locationIncluded ? 1 : 0);
        if (differenced.length - maximumArLag
                <= Math.max(2, parameterCount)) {
            throw new IllegalArgumentException(
                "series is too short for the requested ARIMA order");
        }

        double[] initial = new double[parameterCount];
        double[] lower = new double[parameterCount];
        double[] upper = new double[parameterCount];
        Arrays.fill(lower, -TRANSFORM_BOUND);
        Arrays.fill(upper, TRANSFORM_BOUND);
        if (locationIncluded) {
            int index = parameterCount - 1;
            double mean = mean(differenced);
            double scale = Math.max(1.0, standardDeviation(differenced, mean));
            initial[index] = mean;
            lower[index] = mean - 10.0 * scale;
            upper[index] = mean + 10.0 * scale;
        }

        MultivariableFunction objective = parameters -> {
            ArimaMath.Coefficients coefficients = ArimaMath.decode(
                parameters, order, seasonal, locationIncluded,
                options.includeDrift());
            ArimaMath.Innovations innovations = ArimaMath.innovations(
                differenced, coefficients);
            double variance = innovations.rss() / innovations.used();
            return variance > MINIMUM_VARIANCE && Double.isFinite(variance)
                ? 0.5 * innovations.used() * Math.log(variance)
                : Double.MAX_VALUE / 4.0;
        };
        BoundedOptimizer.Result optimized = multiStart(
            initial, lower, upper, objective, parameterCount
                - (locationIncluded ? 1 : 0),
            options.maximumFunctionEvaluations(),
            options.optimizationTolerance());
        ArimaMath.Coefficients coefficients = ArimaMath.decode(
            optimized.parameters(), order, seasonal, locationIncluded,
            options.includeDrift());
        ArimaMath.Innovations innovations = ArimaMath.innovations(
            differenced, coefficients);
        double innovationVariance = Math.max(MINIMUM_VARIANCE,
            innovations.rss() / innovations.used());
        double logLikelihood = -0.5 * innovations.used()
            * (Math.log(2.0 * Math.PI) + 1.0 + Math.log(innovationVariance));
        int estimatedParameters = parameterCount + 1;
        double aic = -2.0 * logLikelihood + 2.0 * estimatedParameters;
        double denominator = innovations.used() - estimatedParameters - 1.0;
        double aicc = denominator > 0.0
            ? aic + 2.0 * estimatedParameters * (estimatedParameters + 1.0)
                / denominator
            : Double.POSITIVE_INFINITY;
        double bic = -2.0 * logLikelihood
            + Math.log(innovations.used()) * estimatedParameters;
        return new ArimaResult(order, seasonal,
            coefficients.ar(), coefficients.ma(), coefficients.seasonalAr(),
            coefficients.seasonalMa(), coefficients.location(),
            coefficients.drift(), innovationVariance, innovations.values(),
            differenced, values, logLikelihood, aic, aicc, bic,
            innovations.used(), optimized.evaluations(), optimized.converged(),
            optimized.converged() ? "conditional likelihood optimized"
                : "optimizer stopped before its tolerance was reached");
    }

    /** Returns a stationary ARMA correlation matrix with unit diagonal. */
    public static double[] correlationMatrix(
            int observations, double[] autoregressive, double[] movingAverage) {
        double[] ar = MatrixOps.finiteCopy(autoregressive, "autoregressive");
        double[] ma = MatrixOps.finiteCopy(movingAverage, "movingAverage");
        return ArimaMath.correlationMatrix(observations, ar, ma);
    }

    private static double mean(double[] values) {
        double result = 0.0;
        for (double value : values) {
            result += value;
        }
        return result / values.length;
    }

    private static double standardDeviation(double[] values, double mean) {
        double sum = 0.0;
        for (double value : values) {
            double centered = value - mean;
            sum += centered * centered;
        }
        return Math.sqrt(sum / values.length);
    }

    static BoundedOptimizer.Result multiStart(
            double[] initial,
            double[] lower,
            double[] upper,
            MultivariableFunction objective,
            int dynamicParameters,
            int maximumEvaluations,
            double tolerance) {
        if (dynamicParameters == 0 || initial.length == 1) {
            return BoundedOptimizer.minimize(initial, lower, upper, objective,
                maximumEvaluations, tolerance);
        }
        double[][] starts = new double[5][];
        starts[0] = initial.clone();
        double[] levels = {0.55, -0.55, 0.3, -0.3};
        for (int start = 1; start < starts.length; start++) {
            starts[start] = initial.clone();
            for (int parameter = 0;
                    parameter < dynamicParameters; parameter++) {
                starts[start][parameter] = levels[
                    (start - 1 + parameter) % levels.length];
            }
        }
        int perStart = Math.max(20, maximumEvaluations / starts.length);
        BoundedOptimizer.Result best = null;
        int evaluations = 0;
        boolean converged = false;
        for (double[] start : starts) {
            BoundedOptimizer.Result candidate = BoundedOptimizer.minimize(
                start, lower, upper, objective, perStart, tolerance);
            evaluations += candidate.evaluations();
            converged |= candidate.converged();
            if (best == null || candidate.objective() < best.objective()) {
                best = candidate;
            }
        }
        return new BoundedOptimizer.Result(best.parameters(), best.objective(),
            evaluations, converged);
    }
}
