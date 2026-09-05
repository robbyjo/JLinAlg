/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

import java.util.Map;
import java.util.Objects;
import java.util.Optional;

/**
 * Uniform immutable result returned by the basic statistical-test facade.
 * Named maps preserve test-specific degrees of freedom, estimates, and null
 * values without changing the common result shape.
 */
public record StatisticalTestResult(
        String method,
        String statisticName,
        double statistic,
        Map<String, Double> parameters,
        double pValue,
        Map<String, Double> estimates,
        Map<String, Double> nullValues,
        Optional<ConfidenceInterval> confidenceInterval,
        Alternative alternative,
        String pValueMethod) {

    public StatisticalTestResult {
        method = requireText(method, "method");
        statisticName = requireText(statisticName, "statistic name");
        parameters = Map.copyOf(Objects.requireNonNull(parameters, "parameters"));
        estimates = Map.copyOf(Objects.requireNonNull(estimates, "estimates"));
        nullValues = Map.copyOf(Objects.requireNonNull(nullValues, "null values"));
        confidenceInterval = Objects.requireNonNull(
            confidenceInterval, "confidence interval");
        alternative = Objects.requireNonNull(alternative, "alternative");
        pValueMethod = requireText(pValueMethod, "p-value method");
        if ((!Double.isNaN(pValue) && (pValue < 0.0 || pValue > 1.0))
                || Double.isInfinite(pValue)) {
            throw new IllegalArgumentException("p-value must be in [0, 1] or NaN");
        }
    }

    private static String requireText(String value, String name) {
        if (value == null || value.isBlank()) {
            throw new IllegalArgumentException(name + " is required");
        }
        return value;
    }
}
