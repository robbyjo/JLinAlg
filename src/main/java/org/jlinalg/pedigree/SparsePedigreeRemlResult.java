/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;

import java.util.List;
import java.util.LinkedHashMap;
import java.util.Map;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.mixed.SparseLinearMixedModelResult;

/** Animal-model REML fitted through sparse A-inverse equations. */
public final class SparsePedigreeRemlResult {
    private final SparseLinearMixedModelResult mixedModel;
    private final List<String> individualIds;
    private final double[] breedingValues;
    private final double[] predictionErrorVariances;
    private final double[] reliabilities;

    SparsePedigreeRemlResult(
            SparseLinearMixedModelResult mixedModel,
            List<String> individualIds,
            double[] breedingValues,
            double[] inbreedingCoefficients) {
        this.mixedModel = mixedModel;
        this.individualIds = List.copyOf(individualIds);
        this.breedingValues = breedingValues.clone();
        this.predictionErrorVariances = mixedModel
            .randomEffects("additive genetic").predictionErrorVariances();
        this.reliabilities = new double[predictionErrorVariances.length];
        double geneticVariance = additiveGeneticVariance();
        for (int index = 0; index < reliabilities.length; index++) {
            double priorVariance = geneticVariance
                * (1.0 + inbreedingCoefficients[index]);
            predictionErrorVariances[index] = Math.max(0.0,
                Math.min(priorVariance, predictionErrorVariances[index]));
            double value = priorVariance > 0.0
                ? 1.0 - predictionErrorVariances[index] / priorVariance : 0.0;
            reliabilities[index] = Math.max(0.0, Math.min(1.0, value));
        }
    }

    public SparseLinearMixedModelResult mixedModel() { return mixedModel; }
    public AssociationStatistics associationStatistics() {
        return mixedModel.associationStatistics();
    }
    public double[] beta() { return mixedModel.beta(); }
    public double[] standardErrors() { return mixedModel.standardErrors(); }
    public double[] tStatistics() { return mixedModel.tStatistics(); }
    public double[] pValues() { return mixedModel.pValues(); }
    public double additiveGeneticVariance() {
        return mixedModel.varianceComponents()[0];
    }
    public double residualVariance() {
        return mixedModel.varianceComponents()[1];
    }
    public double heritability() {
        double genetic = additiveGeneticVariance();
        return genetic / (genetic + residualVariance());
    }
    public List<String> individualIds() { return individualIds; }
    public double[] breedingValues() { return breedingValues.clone(); }
    /** Diagonal breeding-value prediction-error variances in pedigree order. */
    public double[] predictionErrorVariances() {
        return predictionErrorVariances.clone();
    }
    /** Selects diagonal PEVs without constructing a dense covariance block. */
    public Map<String, Double> predictionErrorVariances(
            List<String> selectedIndividuals) {
        Map<String, Double> result = new LinkedHashMap<>();
        for (String id : selectedIndividuals) {
            int index = indexOf(id);
            result.put(id, predictionErrorVariances[index]);
        }
        return java.util.Collections.unmodifiableMap(result);
    }
    /** Reliability, {@code 1 - PEV / additive variance}, in pedigree order. */
    public double[] reliabilities() { return reliabilities.clone(); }
    public double breedingValue(String individualId) {
        return breedingValues[indexOf(individualId)];
    }
    public double reliability(String individualId) {
        return reliabilities[indexOf(individualId)];
    }

    private int indexOf(String individualId) {
        int index = individualIds.indexOf(individualId);
        if (index < 0)
            throw new IllegalArgumentException(
                "individual is absent from pedigree result: " + individualId);
        return index;
    }
}
