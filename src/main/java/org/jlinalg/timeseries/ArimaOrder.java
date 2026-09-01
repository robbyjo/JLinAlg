/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

/** Nonseasonal ARIMA order {@code (p,d,q)}. */
public record ArimaOrder(int autoregressive, int differences, int movingAverage) {
    public ArimaOrder {
        if (autoregressive < 0 || differences < 0 || movingAverage < 0) {
            throw new IllegalArgumentException("ARIMA orders must be nonnegative");
        }
    }

    public static ArimaOrder ar(int order) { return new ArimaOrder(order, 0, 0); }
    public static ArimaOrder ma(int order) { return new ArimaOrder(0, 0, order); }
    public static ArimaOrder arma(int ar, int ma) { return new ArimaOrder(ar, 0, ma); }
    public static ArimaOrder arima(int ar, int differences, int ma) {
        return new ArimaOrder(ar, differences, ma);
    }
}
