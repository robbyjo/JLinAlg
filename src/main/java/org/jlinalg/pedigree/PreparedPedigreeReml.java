/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;

/** Retained pedigree animal-model structure for fast response refits. */
public final class PreparedPedigreeReml {
    private final double[] fixedEffects;
    private final int observations;
    private final int fixedColumns;
    private final List<String> observationIndividualIds;
    private final int[] observationIndices;
    private final Pedigree pedigree;
    private final RemlOptions options;
    private final BackendPolicy backendPolicy;

    public PreparedPedigreeReml(
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (fixedEffects == null || fixedEffects.length == 0
                || fixedEffects[0] == null || fixedEffects[0].length == 0)
            throw new IllegalArgumentException("fixed-effect design is required");
        observations = fixedEffects.length;
        fixedColumns = fixedEffects[0].length;
        this.fixedEffects = MatrixOps.rowMajor(fixedEffects, observations);
        this.pedigree = Objects.requireNonNull(pedigree, "pedigree");
        if (observationIndividualIds == null
                || observationIndividualIds.size() != observations)
            throw new IllegalArgumentException(
                "one pedigree individual identifier is required per observation");
        this.observationIndividualIds = List.copyOf(observationIndividualIds);
        observationIndices = new int[observations];
        for (int row = 0; row < observations; row++) {
            String id = this.observationIndividualIds.get(row);
            if (id == null)
                throw new IllegalArgumentException(
                    "observation individual identifiers must not be null");
            observationIndices[row] = pedigree.indexOf(id);
        }
        this.options = Objects.requireNonNull(options, "options");
        this.backendPolicy = Objects.requireNonNull(backendPolicy, "backendPolicy");
    }

    public static PreparedPedigreeReml of(
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree) {
        return new PreparedPedigreeReml(fixedEffects,
            observationIndividualIds, pedigree, RemlOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    public PedigreeRemlResult fit(double[] response) {
        validateResponse(response);
        return PedigreeReml.fit(response, fixedEffects, observations,
            fixedColumns, observationIndividualIds, pedigree, options,
            backendPolicy);
    }

    public PedigreeRemlResult refit(
            PedigreeRemlResult previous, double[] response) {
        Objects.requireNonNull(previous, "previous");
        validateResponse(response);
        RemlOptions warm = options.toBuilder().initialVariances(
            previous.reml().varianceComponents()).build();
        return PedigreeReml.fit(response, fixedEffects, observations,
            fixedColumns, observationIndividualIds, pedigree, warm,
            backendPolicy);
    }

    public int observations() { return observations; }
    public int fixedEffectColumns() { return fixedColumns; }
    public List<String> observationIndividualIds() {
        return observationIndividualIds;
    }
    public Pedigree pedigree() { return pedigree; }
    double[] fixedEffectsView() { return fixedEffects; }
    int[] observationIndicesView() { return observationIndices; }
    BackendPolicy backendPolicy() { return backendPolicy; }

    private void validateResponse(double[] response) {
        if (response == null || response.length != observations)
            throw new IllegalArgumentException(
                "response length must equal prepared observations");
        MatrixOps.requireFinite(response, "response");
    }
}
