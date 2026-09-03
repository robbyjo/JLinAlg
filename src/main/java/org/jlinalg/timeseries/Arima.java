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
        boolean closedFormAutoregression = !seasonal.present()
            && order.autoregressive() > 0 && order.movingAverage() == 0
            && initializeAutoregression(
                differenced, order.autoregressive(), locationIncluded,
                initial, lower, upper);

        double[] innovationWorkspace = new double[differenced.length];
        MultivariableFunction objective = parameters -> {
            ArimaMath.Coefficients coefficients = ArimaMath.decode(
                parameters, order, seasonal, locationIncluded,
                options.includeDrift());
            int used = differenced.length - coefficients.effectiveAr().length;
            double variance = ArimaMath.innovationRss(
                differenced, coefficients, innovationWorkspace) / used;
            return variance > MINIMUM_VARIANCE && Double.isFinite(variance)
                ? 0.5 * used * Math.log(variance)
                : Double.MAX_VALUE / 4.0;
        };
        BoundedOptimizer.Result optimized = closedFormAutoregression
            ? new BoundedOptimizer.Result(
                initial.clone(), objective.eval(initial), 1, true)
            : multiStart(initial, lower, upper, objective, parameterCount
                - (locationIncluded ? 1 : 0),
                options.maximumFunctionEvaluations(),
                options.optimizationTolerance(), options.optimizationStarts());
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

    private static boolean initializeAutoregression(
            double[] series,
            int autoregressive,
            boolean locationIncluded,
            double[] initial,
            double[] lower,
            double[] upper) {
        int columns = autoregressive + (locationIncluded ? 1 : 0);
        double[] gram = new double[columns * columns];
        double[] right = new double[columns];
        double[] row = new double[columns];
        for (int time = autoregressive; time < series.length; time++) {
            for (int lag = 0; lag < autoregressive; lag++) {
                row[lag] = series[time - lag - 1];
            }
            if (locationIncluded) row[columns - 1] = 1.0;
            for (int first = 0; first < columns; first++) {
                right[first] += row[first] * series[time];
                for (int second = 0; second <= first; second++) {
                    gram[first * columns + second] += row[first] * row[second];
                }
            }
        }
        for (int first = 0; first < columns; first++) {
            for (int second = 0; second < first; second++) {
                gram[second * columns + first] = gram[first * columns + second];
            }
        }
        double[] solution = solve(gram, right, columns);
        if (solution == null) return false;
        double[] ar = Arrays.copyOf(solution, autoregressive);
        double[] encoded = ArimaMath.encodeAutoregressive(ar);
        if (encoded == null) return false;
        for (int index = 0; index < encoded.length; index++) {
            if (encoded[index] < lower[index] || encoded[index] > upper[index]) {
                return false;
            }
        }
        double location = 0.0;
        if (locationIncluded) {
            double denominator = 1.0;
            for (double coefficient : ar) denominator -= coefficient;
            if (Math.abs(denominator) <= 1e-10) return false;
            location = solution[columns - 1] / denominator;
            int index = initial.length - 1;
            if (!Double.isFinite(location)
                    || location < lower[index] || location > upper[index]) {
                return false;
            }
        }
        System.arraycopy(encoded, 0, initial, 0, autoregressive);
        if (locationIncluded) initial[initial.length - 1] = location;
        return true;
    }

    private static double[] solve(double[] matrix, double[] right, int size) {
        double[] values = matrix.clone();
        double[] result = right.clone();
        for (int pivot = 0; pivot < size; pivot++) {
            int selected = pivot;
            for (int row = pivot + 1; row < size; row++) {
                if (Math.abs(values[row * size + pivot])
                        > Math.abs(values[selected * size + pivot])) {
                    selected = row;
                }
            }
            if (!(Math.abs(values[selected * size + pivot]) > 1e-12)) return null;
            if (selected != pivot) {
                for (int column = pivot; column < size; column++) {
                    double temporary = values[pivot * size + column];
                    values[pivot * size + column] = values[selected * size + column];
                    values[selected * size + column] = temporary;
                }
                double temporary = result[pivot];
                result[pivot] = result[selected];
                result[selected] = temporary;
            }
            double diagonal = values[pivot * size + pivot];
            for (int row = pivot + 1; row < size; row++) {
                double factor = values[row * size + pivot] / diagonal;
                for (int column = pivot + 1; column < size; column++) {
                    values[row * size + column] -=
                        factor * values[pivot * size + column];
                }
                result[row] -= factor * result[pivot];
            }
        }
        for (int row = size - 1; row >= 0; row--) {
            for (int column = row + 1; column < size; column++) {
                result[row] -= values[row * size + column] * result[column];
            }
            result[row] /= values[row * size + row];
        }
        return result;
    }
    static BoundedOptimizer.Result multiStart(
            double[] initial,
            double[] lower,
            double[] upper,
            MultivariableFunction objective,
            int dynamicParameters,
            int maximumEvaluations,
            double tolerance) {
        return multiStart(initial, lower, upper, objective, dynamicParameters,
            maximumEvaluations, tolerance, 5);
    }

    private static BoundedOptimizer.Result multiStart(
            double[] initial,
            double[] lower,
            double[] upper,
            MultivariableFunction objective,
            int dynamicParameters,
            int maximumEvaluations,
            double tolerance,
            int startCount) {
        if (dynamicParameters == 0 || initial.length == 1) {
            return BoundedOptimizer.minimize(initial, lower, upper, objective,
                maximumEvaluations, tolerance);
        }
        double[][] starts = new double[startCount][];
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
