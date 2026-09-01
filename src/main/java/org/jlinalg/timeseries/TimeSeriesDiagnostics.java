/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import jdistlib.ChiSquare;
import org.jlinalg.internal.MatrixOps;

/** Autocorrelation, partial autocorrelation, and Ljung-Box diagnostics. */
public final class TimeSeriesDiagnostics {
    private TimeSeriesDiagnostics() { }

    /** Returns autocorrelations at lags zero through {@code maximumLag}. */
    public static double[] autocorrelation(double[] values, int maximumLag) {
        double[] data = MatrixOps.finiteCopy(values, "values");
        validateLag(data.length, maximumLag);
        double mean = 0.0;
        for (double value : data) {
            mean += value;
        }
        mean /= data.length;
        double denominator = 0.0;
        for (double value : data) {
            double centered = value - mean;
            denominator += centered * centered;
        }
        if (!(denominator > 0.0)) {
            throw new IllegalArgumentException(
                "autocorrelation is undefined for a constant series");
        }
        double[] result = new double[maximumLag + 1];
        result[0] = 1.0;
        for (int lag = 1; lag <= maximumLag; lag++) {
            double numerator = 0.0;
            for (int index = lag; index < data.length; index++) {
                numerator += (data[index] - mean) * (data[index - lag] - mean);
            }
            result[lag] = numerator / denominator;
        }
        return result;
    }

    /** Durbin-Levinson partial autocorrelations at lags zero through max lag. */
    public static double[] partialAutocorrelation(
            double[] values, int maximumLag) {
        double[] acf = autocorrelation(values, maximumLag);
        double[] result = new double[maximumLag + 1];
        result[0] = 1.0;
        double[] previous = new double[maximumLag + 1];
        for (int order = 1; order <= maximumLag; order++) {
            double numerator = acf[order];
            double denominator = 1.0;
            for (int lag = 1; lag < order; lag++) {
                numerator -= previous[lag] * acf[order - lag];
                denominator -= previous[lag] * acf[lag];
            }
            double reflection = denominator > 1e-14
                ? numerator / denominator : Double.NaN;
            result[order] = reflection;
            double[] updated = previous.clone();
            updated[order] = reflection;
            for (int lag = 1; lag < order; lag++) {
                updated[lag] = previous[lag]
                    - reflection * previous[order - lag];
            }
            previous = updated;
        }
        return result;
    }

    /** Computes the Ljung-Box statistic with fitted-parameter DF reduction. */
    public static LjungBoxResult ljungBox(
            double[] residuals, int lags, int fittedArmaParameters) {
        if (fittedArmaParameters < 0 || lags <= fittedArmaParameters) {
            throw new IllegalArgumentException(
                "lags must exceed the nonnegative fitted parameter count");
        }
        double[] acf = autocorrelation(residuals, lags);
        int observations = residuals.length;
        double statistic = 0.0;
        for (int lag = 1; lag <= lags; lag++) {
            statistic += acf[lag] * acf[lag] / (observations - lag);
        }
        statistic *= observations * (observations + 2.0);
        int degrees = lags - fittedArmaParameters;
        double pValue = ChiSquare.cumulative(
            statistic, degrees, false, false);
        return new LjungBoxResult(statistic, lags, degrees, pValue);
    }

    private static void validateLag(int observations, int maximumLag) {
        if (maximumLag < 1 || maximumLag >= observations) {
            throw new IllegalArgumentException(
                "maximum lag must lie between one and observations minus one");
        }
    }
}
