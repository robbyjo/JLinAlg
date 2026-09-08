/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import jdistlib.math.MultivariableFunction;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Exact stationary Gaussian likelihood for complete or incomplete ARMA series. */
public final class ExactArma {
    private static final double BOUND = 3.8;
    private ExactArma() { }

    public static ExactArmaResult fit(double[] series, ArimaOrder order) {
        return fitPanel(List.of(series), order, true, BackendPolicy.PREFERRED);
    }

    public static ExactArmaResult fit(
            double[] series, ArimaOrder order, boolean includeMean,
            BackendPolicy backendPolicy) {
        return fitPanel(List.of(series), order, includeMean, backendPolicy);
    }

    /** Fits shared ARMA parameters to independent series (block-diagonal covariance). */
    public static ExactArmaResult fitPanel(
            List<double[]> series,
            ArimaOrder order,
            boolean includeMean,
            BackendPolicy backendPolicy) {
        if (series == null || series.isEmpty() || order == null
                || backendPolicy == null || order.differences() != 0) {
            throw new IllegalArgumentException("stationary series, order, and backend are required");
        }
        List<ObservedSeries> observed = series.stream().map(ExactArma::observed).toList();
        int count = observed.stream().mapToInt(value -> value.values().length).sum();
        int dynamic = order.autoregressive() + order.movingAverage();
        int parameters = dynamic + (includeMean ? 1 : 0);
        if (count <= parameters + series.size()) {
            throw new IllegalArgumentException("too few observed values for exact ARMA");
        }
        double[] initial = new double[parameters];
        double[] lower = new double[parameters];
        double[] upper = new double[parameters];
        Arrays.fill(lower, -BOUND);
        Arrays.fill(upper, BOUND);
        if (includeMean) {
            double mean = observed.stream().flatMapToDouble(value -> Arrays.stream(value.values()))
                .average().orElseThrow();
            double scale = Math.max(1.0, Math.sqrt(observed.stream()
                .flatMapToDouble(value -> Arrays.stream(value.values()))
                .map(value -> (value - mean) * (value - mean)).average().orElse(1)));
            initial[parameters - 1] = mean;
            lower[parameters - 1] = mean - 10.0 * scale;
            upper[parameters - 1] = mean + 10.0 * scale;
        }
        boolean conditionallyInitialized = initializeFromConditional(
            series, observed, order, includeMean, initial, lower, upper);
        SeasonalArimaOrder seasonal = SeasonalArimaOrder.none();
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            MultivariableFunction objective = point -> evaluate(
                point, order, seasonal, includeMean, observed, backend).negativeLogLikelihood();
            BoundedOptimizer.Result optimized = conditionallyInitialized
                ? BoundedOptimizer.minimize(
                    initial, lower, upper, objective, 5_000, 1e-8)
                : Arima.multiStart(
                    initial, lower, upper, objective, dynamic, 5_000, 1e-8);
            Likelihood likelihood = evaluate(optimized.parameters(), order,
                seasonal, includeMean, observed, backend);
            if (!Double.isFinite(likelihood.variance()) || !(likelihood.variance() > 0)
                    || likelihood.negativeLogLikelihood() >= Double.MAX_VALUE / 16)
                throw new IllegalArgumentException("no finite exact likelihood: degenerate or numerically invalid innovation variance");
            double[] rawCovariance = inverseHessian(
                optimized.parameters(), objective, backend);
            double[] coefficientCovariance = deltaCovariance(
                optimized.parameters(), rawCovariance, order,
                seasonal, includeMean, backend);
            if (!optimized.converged()) Arrays.fill(coefficientCovariance, Double.NaN);
            double[] standardErrors = new double[parameters];
            for (int index = 0; index < parameters; index++) {
                standardErrors[index] = Math.sqrt(Math.max(0.0,
                    coefficientCovariance[index * parameters + index]));
            }
            ArimaMath.Coefficients coefficients = ArimaMath.decode(
                optimized.parameters(), order, seasonal, includeMean, false);
            double innovationVariance = likelihood.variance();
            int likelihoodParameters = parameters + 1;
            double logLikelihood = -likelihood.negativeLogLikelihood();
            return new ExactArmaResult(order, coefficients.ar(), coefficients.ma(),
                coefficients.location(), innovationVariance, coefficientCovariance,
                standardErrors, logLikelihood,
                -2.0 * logLikelihood + 2.0 * likelihoodParameters,
                -2.0 * logLikelihood + Math.log(count) * likelihoodParameters,
                count, series.size(), optimized.evaluations(), optimized.converged(),
                context.provenance());
        }
    }

    private static boolean initializeFromConditional(
            List<double[]> series,
            List<ObservedSeries> observed,
            ArimaOrder order,
            boolean includeMean,
            double[] initial,
            double[] lower,
            double[] upper) {
        if (series.size() != 1 || order.autoregressive() == 0
                || order.movingAverage() != 0
                || observed.get(0).values().length != series.get(0).length) {
            return false;
        }
        ArimaResult conditional;
        try {
            conditional = Arima.fit(series.get(0), order,
                ArimaOptions.builder().includeMean(includeMean).build());
        } catch (IllegalArgumentException exception) {
            return false; // A conditional seed is optional; its sample-size rules differ.
        }
        double[] encoded = ArimaMath.encodeAutoregressive(
            conditional.autoregressive());
        if (encoded == null) return false;
        for (int index = 0; index < encoded.length; index++) {
            if (encoded[index] < lower[index] || encoded[index] > upper[index]) {
                return false;
            }
        }
        System.arraycopy(encoded, 0, initial, 0, encoded.length);
        if (includeMean) {
            int index = initial.length - 1;
            initial[index] = Math.max(lower[index],
                Math.min(upper[index], conditional.location()));
        }
        return true;
    }
    private static Likelihood evaluate(
            double[] parameters, ArimaOrder order, SeasonalArimaOrder seasonal,
            boolean includeMean, List<ObservedSeries> series, ComputeBackend backend) {
        ArimaMath.Coefficients coefficients = ArimaMath.decode(
            parameters, order, seasonal, includeMean, false);
        double quadratic = 0.0;
        double logDeterminant = 0.0;
        int observations = 0;
        ArimaStateSpace model;
        try {
            model = new ArimaStateSpace(coefficients.effectiveAr(), coefficients.effectiveMa(), new double[] {1});
        } catch (IllegalArgumentException exception) {
            return new Likelihood(Double.MAX_VALUE / 8.0, Double.NaN);
        }
        for (ObservedSeries value : series) {
            try {
                ArimaStateSpace.Evaluation contribution = model.filter(value.original(), coefficients.location(), 0, false);
                quadratic += contribution.quadratic();
                logDeterminant += contribution.logDeterminant();
                observations += contribution.observations();
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return new Likelihood(Double.MAX_VALUE / 8.0, Double.NaN);
            }
        }
        double variance = quadratic / observations;
        if (!(variance > 0) || !Double.isFinite(variance)) {
            return new Likelihood(Double.MAX_VALUE / 8.0, Double.NaN);
        }
        double nll = 0.5 * (observations
            * (Math.log(2.0 * Math.PI) + 1.0 + Math.log(variance)) + logDeterminant);
        return new Likelihood(nll, variance);
    }

    static double[] toeplitzLikelihood(
            double[] values, double[] correlation) {
        int size = values.length;
        double[] previous = new double[size];
        double[] current = new double[size];
        double variance = correlation[0];
        double quadratic = values[0] * values[0] / variance;
        double logDeterminant = Math.log(variance);
        for (int order = 1; order < size; order++) {
            double numerator = correlation[order];
            for (int lag = 1; lag < order; lag++) {
                numerator -= previous[lag] * correlation[order - lag];
            }
            double reflection = numerator / variance;
            if (!Double.isFinite(reflection) || Math.abs(reflection) >= 1.0) {
                return null;
            }
            for (int lag = 1; lag < order; lag++) {
                current[lag] = previous[lag]
                    - reflection * previous[order - lag];
            }
            current[order] = reflection;
            variance *= 1.0 - reflection * reflection;
            if (!(variance > 1e-14) || !Double.isFinite(variance)) return null;
            double innovation = values[order];
            for (int lag = 1; lag <= order; lag++) {
                innovation -= current[lag] * values[order - lag];
            }
            quadratic += innovation * innovation / variance;
            logDeterminant += Math.log(variance);
            double[] temporary = previous;
            previous = current;
            current = temporary;
        }
        return new double[] {quadratic, logDeterminant};
    }
    private static ObservedSeries observed(double[] series) {
        if (series == null || series.length < 3) {
            throw new IllegalArgumentException("each series must have at least three positions");
        }
        int count = 0;
        for (double value : series) {
            if (Double.isInfinite(value)) throw new IllegalArgumentException("infinite observation; use NaN for missing");
            if (!Double.isNaN(value)) count++;
        }
        if (count < 3) throw new IllegalArgumentException("each series needs three observed values");
        double[] values = new double[count];
        int target = 0;
        for (int index = 0; index < series.length; index++) {
            if (Double.isFinite(series[index])) {
                values[target++] = series[index];
            }
        }
        return new ObservedSeries(values, series.clone());
    }

    private static double[] inverseHessian(
            double[] point, MultivariableFunction objective, ComputeBackend backend) {
        int size = point.length;
        if (size == 0) return new double[0];
        double[] hessian = hessian(point, objective, 1e-4);
        double[] alternate = hessian(point, objective, 2e-4);
        double[] normalized = new double[size * size];
        double error = 0;
        boolean valid = true;
        for (int i = 0; i < size; i++) valid &= hessian[i * size + i] > 0;
        if (valid) for (int i = 0; i < size; i++) {
            double rowError = 0;
            for (int j = 0; j < size; j++) {
                double scale = Math.sqrt(hessian[i * size + i]) * Math.sqrt(hessian[j * size + j]);
                normalized[i * size + j] = hessian[i * size + j] / scale;
                rowError += Math.abs(hessian[i * size + j] - alternate[i * size + j]) / scale;
            }
            error = Math.max(error, rowError);
        }
        if (valid && Double.isFinite(error)) {
            // Require positive information separated from differentiation error.
            // Subtracting this guard tests identification; never add a ridge to
            // manufacture finite SEs along an AR/MA cancellation direction.
            double guard = Math.max(1e-6, 4 * error);
            for (int i = 0; i < size; i++) normalized[i * size + i] -= guard;
            try {
                backend.dpotrf(normalized, size);
                return backend.dpotrf(hessian, size).solve(MatrixOps.identity(size), size);
            } catch (IllegalArgumentException | IllegalStateException ignored) { }
        }
        double[] result = new double[size * size];
        Arrays.fill(result, Double.NaN);
        return result;
    }

    private static double[] hessian(double[] point, MultivariableFunction objective, double step) {
        int size = point.length;
        double[] hessian = new double[size * size];
        double center = objective.eval(point);
        for (int first = 0; first < size; first++) {
            double h1 = step * (1.0 + Math.abs(point[first]));
            double[] plus = point.clone();
            double[] minus = point.clone();
            plus[first] += h1;
            minus[first] -= h1;
            hessian[first * size + first] =
                (objective.eval(plus) - 2.0 * center + objective.eval(minus)) / (h1 * h1);
            for (int second = 0; second < first; second++) {
                double h2 = step * (1.0 + Math.abs(point[second]));
                double[] pp = point.clone(); pp[first] += h1; pp[second] += h2;
                double[] pm = point.clone(); pm[first] += h1; pm[second] -= h2;
                double[] mp = point.clone(); mp[first] -= h1; mp[second] += h2;
                double[] mm = point.clone(); mm[first] -= h1; mm[second] -= h2;
                double value = (objective.eval(pp) - objective.eval(pm)
                    - objective.eval(mp) + objective.eval(mm)) / (4.0 * h1 * h2);
                hessian[first * size + second] = value;
                hessian[second * size + first] = value;
            }
        }
        return hessian;
    }

    private static double[] deltaCovariance(
            double[] point, double[] rawCovariance, ArimaOrder order,
            SeasonalArimaOrder seasonal, boolean includeMean, ComputeBackend backend) {
        int size = point.length;
        if (size == 0) return new double[0];
        double[] jacobian = new double[size * size];
        for (int parameter = 0; parameter < size; parameter++) {
            double step = 1e-5 * (1.0 + Math.abs(point[parameter]));
            double[] plus = point.clone(); plus[parameter] += step;
            double[] minus = point.clone(); minus[parameter] -= step;
            double[] first = reported(plus, order, seasonal, includeMean);
            double[] second = reported(minus, order, seasonal, includeMean);
            for (int row = 0; row < size; row++) {
                jacobian[row * size + parameter] = (first[row] - second[row]) / (2.0 * step);
            }
        }
        double[] temporary = MatrixOps.multiply(
            backend, jacobian, size, size, rawCovariance, size);
        double[] result = new double[size * size];
        backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
            size, size, size, 1.0, temporary, jacobian, 0.0, result);
        return result;
    }

    private static double[] reported(
            double[] point, ArimaOrder order, SeasonalArimaOrder seasonal,
            boolean includeMean) {
        ArimaMath.Coefficients value = ArimaMath.decode(
            point, order, seasonal, includeMean, false);
        double[] result = new double[point.length];
        int position = 0;
        System.arraycopy(value.ar(), 0, result, position, value.ar().length);
        position += value.ar().length;
        System.arraycopy(value.ma(), 0, result, position, value.ma().length);
        position += value.ma().length;
        if (includeMean) result[position] = value.location();
        return result;
    }

    private record ObservedSeries(double[] values, double[] original) { }
    private record Likelihood(double negativeLogLikelihood, double variance) { }
}
