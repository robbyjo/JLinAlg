/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import jdistlib.math.MultivariableFunction;

/** Exact diffuse Gaussian likelihood for integrated (including seasonal) ARIMA. */
public final class DiffuseArima {
    private DiffuseArima() { }

    /**
     * Fits the original observation sequence, marginalizing NaN observations.
     * The stationary ARMA block has its unconditional covariance; integrated
     * states have symbolic infinite variance. As in R's arima ML convention,
     * diffuse observations contribute neither a quadratic nor a determinant.
     */
    public static Result fit(double[] series, ArimaOrder order, ArimaOptions options) {
        ArimaOptions controls = options == null ? ArimaOptions.defaults() : options;
        SeasonalArimaOrder seasonal = controls.seasonalOrder();
        if (series == null || order == null
                || order.differences() + seasonal.differences() < 1)
            throw new IllegalArgumentException("an integrated ARIMA order is required");
        double[] values = series.clone();
        int observed = 0;
        for (double value : values) {
            if (Double.isInfinite(value)) throw new IllegalArgumentException("infinite observations are invalid; use NaN for missing");
            if (!Double.isNaN(value)) observed++;
        }
        if (controls.includeDrift() && order.differences() + seasonal.differences() != 1)
            throw new IllegalArgumentException("drift requires exactly one differencing operator");
        double[] difference = ArimaMath.differencingPolynomial(order, seasonal);
        int diffuse = difference.length - 1;
        int dynamic = order.autoregressive() + order.movingAverage()
            + seasonal.autoregressive() + seasonal.movingAverage();
        int parameters = dynamic + (controls.includeDrift() ? 1 : 0);
        if (values.length <= diffuse || observed - diffuse <= parameters + 1)
            throw new IllegalArgumentException("too few observed values for the requested diffuse ARIMA");
        double[] initial = new double[dynamic], lower = new double[dynamic], upper = new double[dynamic];
        Arrays.fill(lower, -3.8); Arrays.fill(upper, 3.8);
        // location is the mean of the differenced process. A seasonal difference
        // of period s maps a linear trend with slope location/s to that mean.
        double driftDivisor = seasonal.differences() > 0 ? seasonal.period() : 1;
        MultivariableFunction objective = point -> {
            try {
                ArimaMath.Coefficients c = ArimaMath.decode(point, order, seasonal, false, controls.includeDrift());
                ArimaStateSpace model = new ArimaStateSpace(c.effectiveAr(), c.effectiveMa(), difference);
                double nll = (controls.includeDrift() ? model.profileSlope(values, false)
                    : model.filter(values, 0, 0, false)).negativeLogLikelihood();
                return Double.isFinite(nll) ? nll : Double.MAX_VALUE / 8;
            } catch (IllegalArgumentException exception) { return Double.MAX_VALUE / 8; }
        };
        BoundedOptimizer.Result optimized = Arima.multiStart(initial, lower, upper, objective, dynamic,
            controls.maximumFunctionEvaluations(), controls.optimizationTolerance(), controls.optimizationStarts());
        ArimaMath.Coefficients c = ArimaMath.decode(optimized.parameters(), order, seasonal,
            false, controls.includeDrift());
        ArimaStateSpace model = new ArimaStateSpace(c.effectiveAr(), c.effectiveMa(), difference);
        ArimaStateSpace.Evaluation likelihood = controls.includeDrift() ? model.profileSlope(values, true)
            : model.filter(values, 0, 0, true);
        if (!Double.isFinite(likelihood.negativeLogLikelihood()))
            throw new IllegalArgumentException("diffuse states are unidentified or the innovation variance is degenerate");
        int count = likelihood.observations(), k = parameters + 1;
        double ll = -likelihood.negativeLogLikelihood(), aic = -2 * ll + 2 * k;
        double aicc = count > k + 1 ? aic + 2.0 * k * (k + 1) / (count - k - 1) : Double.POSITIVE_INFINITY;
        ArimaResult result = new ArimaResult(order, seasonal, c.ar(), c.ma(), c.seasonalAr(), c.seasonalMa(),
            likelihood.fittedSlope() * driftDivisor, c.drift(), likelihood.variance(), likelihood.innovations(),
            ArimaMath.difference(values, order, seasonal), values, ll, aic, aicc, -2 * ll + Math.log(count) * k,
            count, optimized.evaluations(), optimized.converged(), optimized.converged()
                ? "exact diffuse likelihood optimized" : "optimizer stopped before convergence",
            likelihood.state());
        return new Result(result, diffuse, count, true);
    }

    public record Result(ArimaResult fit, int diffuseStateCount, int likelihoodObservations,
                         boolean diffuseLikelihood) { }
}
