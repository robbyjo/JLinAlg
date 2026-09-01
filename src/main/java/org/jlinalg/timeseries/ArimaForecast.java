/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

/** Forecast means, standard errors, and pointwise normal intervals. */
public final class ArimaForecast {
    private final double[] means;
    private final double[] standardErrors;
    private final double[] lowerBounds;
    private final double[] upperBounds;
    private final double confidenceLevel;

    ArimaForecast(
            double[] means,
            double[] standardErrors,
            double[] lowerBounds,
            double[] upperBounds,
            double confidenceLevel) {
        this.means = means.clone();
        this.standardErrors = standardErrors.clone();
        this.lowerBounds = lowerBounds.clone();
        this.upperBounds = upperBounds.clone();
        this.confidenceLevel = confidenceLevel;
    }

    public double[] means() { return means.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] lowerBounds() { return lowerBounds.clone(); }
    public double[] upperBounds() { return upperBounds.clone(); }
    public double confidenceLevel() { return confidenceLevel; }
}
