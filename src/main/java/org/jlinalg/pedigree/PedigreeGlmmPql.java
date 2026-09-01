/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmPql;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.glmm.GlmmPqlResult;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.VarianceComponent;

/** First-order PQL GLMM with a pedigree numerator-relationship random effect. */
public final class PedigreeGlmmPql {
    private static final String GENETIC_COMPONENT = "additive genetic";

    private PedigreeGlmmPql() { }

    /** Fits with unit prior weights, zero offset, and preferred acceleration. */
    public static PedigreeGlmmPqlResult fit(
            double[] response,
            double[][] fixedEffects,
            GlmFamily family,
            List<String> observationIndividualIds,
            Pedigree pedigree) {
        return fit(response, fixedEffects, family, observationIndividualIds,
            pedigree, null, null, GlmmPqlOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    /** Fits a dense pedigree PQL model from a conventional fixed-effect matrix. */
    public static PedigreeGlmmPqlResult fit(
            double[] response,
            double[][] fixedEffects,
            GlmFamily family,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, rowMajor, response.length, fixedEffects[0].length,
            family, observationIndividualIds, pedigree, priorWeights, offset,
            options, backendPolicy);
    }

    /** Fits from a contiguous row-major fixed-effect matrix. */
    public static PedigreeGlmmPqlResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            GlmFamily family,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        if (pedigree == null) {
            throw new IllegalArgumentException("pedigree is required");
        }
        int[] observationIndices = observationIndices(
            observationIndividualIds, rows, pedigree);
        double[] relationship = observationRelationship(
            observationIndices, pedigree);
        GlmmPqlResult result = GlmmPql.fit(
            response, fixedEffects, rows, columns, family,
            List.of(new VarianceComponent(
                GENETIC_COMPONENT, rows, relationship)),
            priorWeights, offset, options, backendPolicy);
        return new PedigreeGlmmPqlResult(result, observationIndividualIds);
    }

    private static int[] observationIndices(
            List<String> individualIds, int rows, Pedigree pedigree) {
        if (individualIds == null || individualIds.size() != rows) {
            throw new IllegalArgumentException(
                "one pedigree individual identifier is required per observation");
        }
        int[] result = new int[rows];
        for (int index = 0; index < rows; index++) {
            String id = individualIds.get(index);
            if (id == null) {
                throw new IllegalArgumentException(
                    "observation individual identifiers must not be null");
            }
            result[index] = pedigree.indexOf(id);
        }
        return result;
    }

    private static double[] observationRelationship(
            int[] observationIndices, Pedigree pedigree) {
        int rows = observationIndices.length;
        double[] result = new double[rows * rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column <= row; column++) {
                double value = pedigree.relationship(
                    observationIndices[row], observationIndices[column]);
                result[row * rows + column] = value;
                result[column * rows + row] = value;
            }
        }
        return result;
    }
}
