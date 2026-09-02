/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

/** Quasi-likelihood Tweedie variance model with log link and fixed power. */
public final class TweedieLogFamily implements GlmFamily {
    private static final double MINIMUM = 1e-12;
    private final double power;

    /** Power is commonly in (1,2) for semicontinuous mass-at-zero outcomes. */
    public TweedieLogFamily(double power) {
        if (!Double.isFinite(power) || power == 1.0 || power == 2.0) {
            throw new IllegalArgumentException(
                "Tweedie power must be finite and different from one and two");
        }
        this.power = power;
    }
    public double power() { return power; }
    @Override public String name() { return "Tweedie(log,p=" + power + ")"; }
    @Override public void validateResponse(double response, double priorWeight) {
        if (!Double.isFinite(response) || response < 0.0) {
            throw new IllegalArgumentException("Tweedie responses must be finite and nonnegative");
        }
    }
    @Override public double initialMean(double response) {
        return Math.max(0.1, response + 0.1);
    }
    @Override public double link(double mean) {
        return Math.log(Math.max(MINIMUM, mean));
    }
    @Override public double inverseLink(double predictor) {
        return Math.max(MINIMUM, Math.exp(Math.min(700.0, predictor)));
    }
    @Override public double meanDerivative(double predictor) {
        return inverseLink(predictor);
    }
    @Override public double variance(double mean) {
        return Math.max(MINIMUM, Math.pow(mean, power));
    }
    @Override public double unitDeviance(double response, double mean) {
        double first = response == 0.0 ? 0.0
            : Math.pow(response, 2.0 - power)
                / ((1.0 - power) * (2.0 - power));
        double second = response * Math.pow(mean, 1.0 - power)
            / (1.0 - power);
        double third = Math.pow(mean, 2.0 - power) / (2.0 - power);
        return 2.0 * (first - second + third);
    }
    @Override public double logLikelihood(
            double response, double mean, double priorWeight, double dispersion) {
        return Double.NaN;
    }
    @Override public boolean fixedDispersion() { return false; }
}
