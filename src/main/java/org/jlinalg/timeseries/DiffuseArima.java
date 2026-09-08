/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

/** Integrated ARIMA fit with explicit diffuse initial-state bookkeeping. */
public final class DiffuseArima {
    private DiffuseArima() { }

    /**
     * Fits the differenced likelihood and reports the number of diffuse
     * observations removed before stationary innovation evaluation.
     */
    public static Result fit(double[] series, ArimaOrder order, ArimaOptions options) {
        if (series == null || order == null || order.differences() < 1)
            throw new IllegalArgumentException("an integrated ARIMA order is required");
        ArimaResult result = Arima.fit(series, order, options == null ? ArimaOptions.defaults() : options);
        int diffuse = order.differences() + (options == null ? 0 : options.seasonalOrder().differences());
        return new Result(result, diffuse, result.effectiveObservations(), true);
    }

    public record Result(ArimaResult fit, int diffuseStateCount, int likelihoodObservations,
                         boolean diffuseLikelihood) { }
}
