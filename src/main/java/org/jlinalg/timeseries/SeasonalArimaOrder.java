/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

/** Multiplicative seasonal ARIMA order {@code (P,D,Q)[period]}. */
public record SeasonalArimaOrder(
        int autoregressive,
        int differences,
        int movingAverage,
        int period) {
    private static final SeasonalArimaOrder NONE =
        new SeasonalArimaOrder(0, 0, 0, 1);

    public SeasonalArimaOrder {
        if (autoregressive < 0 || differences < 0 || movingAverage < 0) {
            throw new IllegalArgumentException(
                "seasonal ARIMA orders must be nonnegative");
        }
        if (period < 1) {
            throw new IllegalArgumentException("seasonal period must be positive");
        }
        if (period == 1
                && (autoregressive > 0 || differences > 0 || movingAverage > 0)) {
            throw new IllegalArgumentException(
                "a nonzero seasonal order requires period greater than one");
        }
    }

    public static SeasonalArimaOrder none() { return NONE; }

    public static SeasonalArimaOrder of(
            int ar, int differences, int ma, int period) {
        return new SeasonalArimaOrder(ar, differences, ma, period);
    }

    public boolean present() {
        return autoregressive > 0 || differences > 0 || movingAverage > 0;
    }
}
