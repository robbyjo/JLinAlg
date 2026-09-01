/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;

import java.util.List;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.mixed.SparseLinearMixedModelResult;

/** Animal-model REML fitted through sparse A-inverse equations. */
public final class SparsePedigreeRemlResult {
    private final SparseLinearMixedModelResult mixedModel;
    private final List<String> individualIds;
    private final double[] breedingValues;

    SparsePedigreeRemlResult(
            SparseLinearMixedModelResult mixedModel,
            List<String> individualIds,
            double[] breedingValues) {
        this.mixedModel = mixedModel;
        this.individualIds = List.copyOf(individualIds);
        this.breedingValues = breedingValues.clone();
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
    public double breedingValue(String individualId) {
        int index = individualIds.indexOf(individualId);
        if (index < 0)
            throw new IllegalArgumentException(
                "individual is absent from pedigree result: " + individualId);
        return breedingValues[index];
    }
}
