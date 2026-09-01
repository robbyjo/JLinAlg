/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.IntStream;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.reml.RemlResult;

/** Immutable estimates and animal predictions from pedigree REML. */
public final class PedigreeRemlResult {
    private final RemlResult reml;
    private final List<String> individualIds;
    private final Map<String, Integer> indexById;
    private final double[] breedingValues;
    private final double[] predictionErrorVariances;
    private final double[] reliabilities;

    PedigreeRemlResult(
            RemlResult reml,
            List<String> individualIds,
            double[] breedingValues,
            double[] predictionErrorVariances,
            double[] reliabilities) {
        this.reml = Objects.requireNonNull(reml, "reml");
        this.individualIds = List.copyOf(individualIds);
        this.indexById = IntStream.range(0, individualIds.size()).boxed()
            .collect(java.util.stream.Collectors.toUnmodifiableMap(
                individualIds::get, index -> index));
        this.breedingValues = breedingValues.clone();
        this.predictionErrorVariances = predictionErrorVariances.clone();
        this.reliabilities = reliabilities.clone();
    }

    /** Returns the underlying two-component Gaussian REML result. */
    public RemlResult reml() { return reml; }

    /** Fixed-effect beta, SE, t statistic, denominator DF, and p-value table. */
    public AssociationStatistics associationStatistics() {
        return reml.associationStatistics();
    }
    public double[] beta() { return reml.beta(); }
    public double[] standardErrors() { return reml.standardErrors(); }
    public double[] tStatistics() { return reml.tStatistics(); }
    public double[] pValues() { return reml.pValues(); }

    /** Returns the additive genetic variance estimate. */
    public double additiveGeneticVariance() {
        return reml.varianceComponents()[0];
    }

    /** Returns the independent residual variance estimate. */
    public double residualVariance() {
        return reml.varianceComponents()[1];
    }

    /** Returns narrow-sense heritability on the modeled observation scale. */
    public double heritability() {
        double genetic = additiveGeneticVariance();
        return genetic / (genetic + residualVariance());
    }

    /** Returns pedigree identifiers in prediction-array order. */
    public List<String> individualIds() { return individualIds; }

    /** Returns predicted additive breeding values for all pedigree individuals. */
    public double[] breedingValues() { return breedingValues.clone(); }

    /** Returns breeding-value prediction error variances. */
    public double[] predictionErrorVariances() {
        return predictionErrorVariances.clone();
    }

    /** Returns reliabilities, {@code 1 - PEV / genetic variance}, by individual. */
    public double[] reliabilities() { return reliabilities.clone(); }

    /** Returns one individual's predicted additive breeding value. */
    public double breedingValue(String individualId) {
        return breedingValues[indexOf(individualId)];
    }

    /** Returns one individual's prediction reliability. */
    public double reliability(String individualId) {
        return reliabilities[indexOf(individualId)];
    }

    private int indexOf(String individualId) {
        Integer index = indexById.get(individualId);
        if (index == null) {
            throw new IllegalArgumentException(
                "individual is absent from the pedigree result: " + individualId);
        }
        return index;
    }
}
