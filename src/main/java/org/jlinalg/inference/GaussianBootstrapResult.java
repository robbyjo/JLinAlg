/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Parametric-bootstrap draws and percentile summaries for a Gaussian model. */
public final class GaussianBootstrapResult {
    private final int requestedSimulations;
    private final long randomSeed;
    private final double confidenceLevel;
    private final double[] observedFixedEffects;
    private final List<String> varianceComponentNames;
    private final double[] observedVarianceComponents;
    private final double[][] fixedEffectReplicates;
    private final double[][] varianceComponentReplicates;
    private final List<BootstrapFailure> failures;

    public GaussianBootstrapResult(
            int requestedSimulations,
            long randomSeed,
            double confidenceLevel,
            double[] observedFixedEffects,
            List<String> varianceComponentNames,
            double[] observedVarianceComponents,
            double[][] fixedEffectReplicates,
            double[][] varianceComponentReplicates,
            List<BootstrapFailure> failures) {
        this.requestedSimulations = requestedSimulations;
        this.randomSeed = randomSeed;
        this.confidenceLevel = confidenceLevel;
        this.observedFixedEffects = observedFixedEffects.clone();
        this.varianceComponentNames = List.copyOf(varianceComponentNames);
        this.observedVarianceComponents = observedVarianceComponents.clone();
        this.fixedEffectReplicates = deepCopy(fixedEffectReplicates);
        this.varianceComponentReplicates = deepCopy(
            varianceComponentReplicates);
        this.failures = List.copyOf(failures);
        if (this.varianceComponentNames.size()
                != this.observedVarianceComponents.length
                || this.fixedEffectReplicates.length
                    != this.varianceComponentReplicates.length)
            throw new IllegalArgumentException(
                "bootstrap result dimensions are invalid");
    }

    public int requestedSimulations() { return requestedSimulations; }
    public int successfulSimulations() { return fixedEffectReplicates.length; }
    public int failedSimulations() { return failures.size(); }
    public long randomSeed() { return randomSeed; }
    public double confidenceLevel() { return confidenceLevel; }
    public double[][] fixedEffectReplicates() {
        return deepCopy(fixedEffectReplicates);
    }
    public double[][] varianceComponentReplicates() {
        return deepCopy(varianceComponentReplicates);
    }
    public List<BootstrapFailure> failures() { return failures; }

    public List<BootstrapParameterSummary> fixedEffectSummaries() {
        List<BootstrapParameterSummary> result =
            new ArrayList<>(observedFixedEffects.length);
        for (int parameter = 0; parameter < observedFixedEffects.length;
                parameter++)
            result.add(summary("beta[" + parameter + "]",
                observedFixedEffects[parameter], column(
                    fixedEffectReplicates, parameter)));
        return List.copyOf(result);
    }

    public List<BootstrapParameterSummary> varianceComponentSummaries() {
        List<BootstrapParameterSummary> result =
            new ArrayList<>(observedVarianceComponents.length);
        for (int parameter = 0;
                parameter < observedVarianceComponents.length; parameter++)
            result.add(summary(varianceComponentNames.get(parameter),
                observedVarianceComponents[parameter], column(
                    varianceComponentReplicates, parameter)));
        return List.copyOf(result);
    }

    private BootstrapParameterSummary summary(
            String name, double estimate, double[] values) {
        if (values.length < 2)
            return new BootstrapParameterSummary(name, estimate,
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                values.length);
        double mean = Arrays.stream(values).average().orElse(Double.NaN);
        double sumSquares = 0;
        for (double value : values) {
            double difference = value - mean;
            sumSquares += difference * difference;
        }
        double standardError = Math.sqrt(sumSquares / (values.length - 1));
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        double alpha = (1 - confidenceLevel) / 2;
        return new BootstrapParameterSummary(name, estimate, mean,
            mean - estimate, standardError, quantile(sorted, alpha),
            quantile(sorted, 1 - alpha), values.length);
    }

    private static double quantile(double[] sorted, double probability) {
        double position = probability * (sorted.length - 1);
        int lower = (int) Math.floor(position);
        int upper = (int) Math.ceil(position);
        if (lower == upper) return sorted[lower];
        double fraction = position - lower;
        return sorted[lower] * (1 - fraction) + sorted[upper] * fraction;
    }

    private static double[] column(double[][] matrix, int column) {
        double[] result = new double[matrix.length];
        for (int row = 0; row < matrix.length; row++)
            result[row] = matrix[row][column];
        return result;
    }

    private static double[][] deepCopy(double[][] values) {
        double[][] result = new double[values.length][];
        for (int index = 0; index < values.length; index++)
            result[index] = values[index].clone();
        return result;
    }
}
