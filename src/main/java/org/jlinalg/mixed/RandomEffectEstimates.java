/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import java.util.List;
import java.util.Objects;

/** Conditional modes and prediction-error variances for one random term. */
public final class RandomEffectEstimates {
    private final String termName;
    private final List<String> coefficientNames;
    private final double variance;
    private final double[] estimates;
    private final double[] predictionErrorVariances;

    RandomEffectEstimates(
            String termName,
            List<String> coefficientNames,
            double variance,
            double[] estimates,
            double[] predictionErrorVariances) {
        this.termName = Objects.requireNonNull(termName, "termName");
        this.coefficientNames = List.copyOf(coefficientNames);
        this.variance = variance;
        this.estimates = estimates.clone();
        this.predictionErrorVariances = predictionErrorVariances.clone();
    }

    public String termName() { return termName; }
    public List<String> coefficientNames() { return coefficientNames; }
    public double variance() { return variance; }
    public double[] estimates() { return estimates.clone(); }
    public double[] predictionErrorVariances() {
        return predictionErrorVariances.clone();
    }

    public double estimate(String coefficientName) {
        int index = coefficientNames.indexOf(coefficientName);
        if (index < 0) {
            throw new IllegalArgumentException(
                "unknown random-effect coefficient: " + coefficientName);
        }
        return estimates[index];
    }
}
