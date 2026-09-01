/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.jlinalg.glmm.GlmmPqlResult;
import org.jlinalg.inference.AssociationStatistics;

/** Dense pedigree GLMM estimates from the first-order PQL implementation. */
public final class PedigreeGlmmPqlResult {
    private final GlmmPqlResult glmm;
    private final List<String> observationIndividualIds;
    private final Map<String, Double> observedBreedingValues;

    PedigreeGlmmPqlResult(
            GlmmPqlResult glmm, List<String> observationIndividualIds) {
        this.glmm = Objects.requireNonNull(glmm, "glmm");
        this.observationIndividualIds = List.copyOf(observationIndividualIds);
        double[] predictions = glmm.randomLinearPredictor();
        Map<String, double[]> accumulators = new LinkedHashMap<>();
        for (int index = 0; index < predictions.length; index++) {
            double[] accumulator = accumulators.computeIfAbsent(
                observationIndividualIds.get(index), ignored -> new double[2]);
            accumulator[0] += predictions[index];
            accumulator[1] += 1.0;
        }
        Map<String, Double> values = new LinkedHashMap<>();
        for (Map.Entry<String, double[]> entry : accumulators.entrySet()) {
            values.put(entry.getKey(), entry.getValue()[0] / entry.getValue()[1]);
        }
        this.observedBreedingValues = Map.copyOf(values);
    }

    /** Underlying generic PQL fit; its sole variance component is additive genetic. */
    public GlmmPqlResult glmm() { return glmm; }

    public AssociationStatistics associationStatistics() {
        return glmm.associationStatistics();
    }
    public double[] beta() { return glmm.beta(); }
    public double[] standardErrors() { return glmm.standardErrors(); }
    public double[] tStatistics() { return glmm.tStatistics(); }
    public double[] pValues() { return glmm.pValues(); }

    public double additiveGeneticVariance() {
        return glmm.varianceComponents()[0];
    }

    public List<String> observationIndividualIds() {
        return observationIndividualIds;
    }

    /**
     * Observation-supported additive prediction on the link scale. Ancestors
     * without observations are not included in this dense PQL result.
     */
    public Map<String, Double> observedBreedingValues() {
        return observedBreedingValues;
    }

    public double observedBreedingValue(String individualId) {
        Double value = observedBreedingValues.get(individualId);
        if (value == null) {
            throw new IllegalArgumentException(
                "individual has no observation-supported PQL prediction: "
                    + individualId);
        }
        return value;
    }
}
