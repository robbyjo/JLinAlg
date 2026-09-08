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
        if (requestedSimulations < 2 || !(confidenceLevel > 0 && confidenceLevel < 1)
                || observedFixedEffects == null || observedVarianceComponents == null
                || varianceComponentNames == null || fixedEffectReplicates == null
                || varianceComponentReplicates == null || failures == null
                || (long) fixedEffectReplicates.length + failures.size() != requestedSimulations) {
            throw new IllegalArgumentException("bootstrap accounting or inputs are invalid");
        }
        validateReplicates(fixedEffectReplicates, observedFixedEffects.length);
        validateReplicates(varianceComponentReplicates, observedVarianceComponents.length);
        java.util.HashSet<Integer> failureIds = new java.util.HashSet<>();
        for (BootstrapFailure failure : failures) {
            if (failure == null || failure.simulation() >= requestedSimulations
                    || !failureIds.add(failure.simulation())) {
                throw new IllegalArgumentException("bootstrap failure indices must be distinct and in range");
            }
        }
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
        double scale = 0;
        for (double value : values) scale = Math.max(scale, Math.abs(value));
        if (scale == 0) scale = 1;
        double normalizedSum = 0, correction = 0;
        for (double value : values) {
            double adjusted = value / scale - correction;
            double next = normalizedSum + adjusted;
            correction = (next - normalizedSum) - adjusted;
            normalizedSum = next;
        }
        double normalizedMean = normalizedSum / values.length;
        double mean = normalizedMean * scale;
        double sumSquares = 0;
        for (double value : values) {
            double difference = value / scale - normalizedMean;
            sumSquares += difference * difference;
        }
        double standardError = Math.sqrt(sumSquares / (values.length - 1)) * scale;
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

    private static void validateReplicates(double[][] values, int columns) {
        for (double[] row : values) {
            if (row == null || row.length != columns) {
                throw new IllegalArgumentException("bootstrap replicate dimensions are invalid");
            }
            for (double value : row) if (!Double.isFinite(value)) {
                throw new IllegalArgumentException("successful bootstrap replicates must be finite");
            }
        }
    }
}
