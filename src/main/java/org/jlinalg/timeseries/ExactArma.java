/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
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
            double scale = Math.max(1.0, Math.abs(mean));
            initial[parameters - 1] = mean;
            lower[parameters - 1] = mean - 10.0 * scale;
            upper[parameters - 1] = mean + 10.0 * scale;
        }
        SeasonalArimaOrder seasonal = SeasonalArimaOrder.none();
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            MultivariableFunction objective = point -> evaluate(
                point, order, seasonal, includeMean, observed, backend).negativeLogLikelihood();
            BoundedOptimizer.Result optimized = Arima.multiStart(
                initial, lower, upper, objective, dynamic, 5_000, 1e-8);
            Likelihood likelihood = evaluate(optimized.parameters(), order,
                seasonal, includeMean, observed, backend);
            double[] rawCovariance = inverseHessian(
                optimized.parameters(), objective, backend);
            double[] coefficientCovariance = deltaCovariance(
                optimized.parameters(), rawCovariance, order,
                seasonal, includeMean, backend);
            double[] standardErrors = new double[parameters];
            for (int index = 0; index < parameters; index++) {
                standardErrors[index] = Math.sqrt(Math.max(0.0,
                    coefficientCovariance[index * parameters + index]));
            }
            ArimaMath.Coefficients coefficients = ArimaMath.decode(
                optimized.parameters(), order, seasonal, includeMean, false);
            double innovationVariance = likelihood.variance()
                / ArimaMath.marginalVariancePerInnovation(
                    coefficients.effectiveAr(), coefficients.effectiveMa());
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

    private static Likelihood evaluate(
            double[] parameters, ArimaOrder order, SeasonalArimaOrder seasonal,
            boolean includeMean, List<ObservedSeries> series, ComputeBackend backend) {
        ArimaMath.Coefficients coefficients = ArimaMath.decode(
            parameters, order, seasonal, includeMean, false);
        double quadratic = 0.0;
        double logDeterminant = 0.0;
        int observations = 0;
        for (ObservedSeries value : series) {
            double[] full = ArimaMath.correlationMatrix(
                value.originalLength(), coefficients.effectiveAr(), coefficients.effectiveMa());
            int size = value.values().length;
            double[] covariance = new double[size * size];
            double[] centered = new double[size];
            for (int row = 0; row < size; row++) {
                centered[row] = value.values()[row] - coefficients.location();
                for (int column = 0; column < size; column++) {
                    covariance[row * size + column] = full[
                        value.indices()[row] * value.originalLength() + value.indices()[column]];
                }
            }
            try {
                CholeskyFactor factor = backend.dpotrf(covariance, size);
                double[] solved = factor.solve(centered);
                quadratic += backend.ddot(size, centered, 0, 1, solved, 0, 1);
                logDeterminant += factor.logDeterminant();
                observations += size;
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return new Likelihood(Double.MAX_VALUE / 8.0, Double.NaN);
            }
        }
        double variance = quadratic / observations;
        if (!(variance > 1e-14) || !Double.isFinite(variance)) {
            return new Likelihood(Double.MAX_VALUE / 8.0, Double.NaN);
        }
        double nll = 0.5 * (observations
            * (Math.log(2.0 * Math.PI) + 1.0 + Math.log(variance)) + logDeterminant);
        return new Likelihood(nll, variance);
    }

    private static ObservedSeries observed(double[] series) {
        if (series == null || series.length < 3) {
            throw new IllegalArgumentException("each series must have at least three positions");
        }
        int count = 0;
        for (double value : series) if (Double.isFinite(value)) count++;
        if (count < 3) throw new IllegalArgumentException("each series needs three observed values");
        double[] values = new double[count];
        int[] indices = new int[count];
        int target = 0;
        for (int index = 0; index < series.length; index++) {
            if (Double.isFinite(series[index])) {
                values[target] = series[index];
                indices[target++] = index;
            }
        }
        return new ObservedSeries(values, indices, series.length);
    }

    private static double[] inverseHessian(
            double[] point, MultivariableFunction objective, ComputeBackend backend) {
        int size = point.length;
        if (size == 0) return new double[0];
        double[] hessian = new double[size * size];
        double center = objective.eval(point);
        for (int first = 0; first < size; first++) {
            double h1 = 1e-4 * (1.0 + Math.abs(point[first]));
            double[] plus = point.clone();
            double[] minus = point.clone();
            plus[first] += h1;
            minus[first] -= h1;
            hessian[first * size + first] =
                (objective.eval(plus) - 2.0 * center + objective.eval(minus)) / (h1 * h1);
            for (int second = 0; second < first; second++) {
                double h2 = 1e-4 * (1.0 + Math.abs(point[second]));
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
        double ridge = 1e-10;
        for (int attempt = 0; attempt < 12; attempt++) {
            double[] regularized = hessian.clone();
            for (int index = 0; index < size; index++)
                regularized[index * size + index] += ridge;
            try {
                return backend.dpotrf(regularized, size)
                    .solve(MatrixOps.identity(size), size);
            } catch (IllegalArgumentException | IllegalStateException ignored) {
                ridge *= 10.0;
            }
        }
        double[] result = new double[size * size];
        Arrays.fill(result, Double.NaN);
        return result;
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

    private record ObservedSeries(double[] values, int[] indices, int originalLength) { }
    private record Likelihood(double negativeLogLikelihood, double variance) { }
}
