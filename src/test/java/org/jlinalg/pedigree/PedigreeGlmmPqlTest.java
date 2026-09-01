/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class PedigreeGlmmPqlTest {
    @Test
    void fitsBinomialPedigreePqlAndReturnsAssociationInference() {
        double[] response = {
            0, 0, 0, 0, 1,
            0, 0, 0, 1, 0,
            0, 0, 1, 1, 0,
            1, 1, 1, 0, 0,
            1, 1, 1, 1, 0,
            1, 1, 1, 0, 1
        };
        double[][] fixed = new double[response.length][1];
        for (double[] row : fixed) {
            row[0] = 1.0;
        }
        List<PedigreeIndividual> individuals = new ArrayList<>();
        List<String> observed = new ArrayList<>();
        for (int individual = 0; individual < 6; individual++) {
            String id = "animal-" + individual;
            individuals.add(PedigreeIndividual.founder(id));
            for (int record = 0; record < 5; record++) {
                observed.add(id);
            }
        }
        Pedigree pedigree = Pedigree.of(individuals);
        GlmmPqlOptions options = GlmmPqlOptions.builder()
            .maximumIterations(50)
            .relativeTolerance(1e-5)
            .remlOptions(RemlOptions.builder()
                .initialVariances(0.5)
                .scoreTolerance(1e-6)
                .build())
            .build();

        PedigreeGlmmPqlResult result = PedigreeGlmmPql.fit(
            response, fixed, GlmFamilies.binomial(), observed, pedigree,
            null, null, options, BackendPolicy.CPU);

        assertTrue(result.glmm().converged(),
            result.glmm().convergenceMessage());
        assertTrue(result.additiveGeneticVariance() > 0.01);
        assertEquals(6, result.observedBreedingValues().size());
        assertEquals(1, result.associationStatistics().beta().length);
        assertTrue(result.associationStatistics().pValues()[0] >= 0.0);
    }
}
