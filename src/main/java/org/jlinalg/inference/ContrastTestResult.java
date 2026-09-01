/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

/** Immutable joint linear-contrast Wald test. */
public final class ContrastTestResult {
    private final double[] estimates;
    private final double[] covariance;
    private final int numeratorDegreesOfFreedom;
    private final double denominatorDegreesOfFreedom;
    private final double statistic;
    private final double pValue;
    private final StatisticDistribution distribution;

    ContrastTestResult(
            double[] estimates, double[] covariance,
            int numeratorDegreesOfFreedom,
            double denominatorDegreesOfFreedom,
            double statistic, double pValue,
            StatisticDistribution distribution) {
        this.estimates = estimates.clone();
        this.covariance = covariance.clone();
        this.numeratorDegreesOfFreedom = numeratorDegreesOfFreedom;
        this.denominatorDegreesOfFreedom = denominatorDegreesOfFreedom;
        this.statistic = statistic;
        this.pValue = pValue;
        this.distribution = distribution;
    }

    public double[] estimates() { return estimates.clone(); }
    public double[] covariance() { return covariance.clone(); }
    public int numeratorDegreesOfFreedom() { return numeratorDegreesOfFreedom; }
    public double denominatorDegreesOfFreedom() {
        return denominatorDegreesOfFreedom;
    }
    /** F statistic for finite denominator DF; chi-square otherwise. */
    public double statistic() { return statistic; }
    public double pValue() { return pValue; }
    public StatisticDistribution distribution() { return distribution; }
}
