/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class SparsePedigreeRemlTest {
    @Test
    void sparseAInverseVarianceFitMatchesDenseAnimalModel() {
        double[] response = {0, 1, 2, 4, 5, 6, 8, 9, 10};
        double[][] fixed = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("A"),
            PedigreeIndividual.founder("B"),
            PedigreeIndividual.founder("C")));
        List<String> observed = List.of(
            "A", "A", "A", "B", "B", "B", "C", "C", "C");
        RemlOptions options = RemlOptions.builder()
            .initialVariances(10.0, 2.0)
            .maximumIterations(200).build();

        PedigreeRemlResult dense = PedigreeReml.fit(response, fixed,
            observed, pedigree, options, BackendPolicy.CPU);
        SparsePedigreeRemlResult sparse = SparsePedigreeReml.fit(
            response, fixed, observed, pedigree, options, BackendPolicy.CPU);

        assertArrayEquals(dense.reml().varianceComponents(),
            sparse.mixedModel().varianceComponents(), 2e-5);
        assertArrayEquals(dense.beta(), sparse.beta(), 1e-8);
        assertArrayEquals(dense.breedingValues(),
            sparse.breedingValues(), 2e-5);
        assertEquals(dense.heritability(), sparse.heritability(), 2e-6);
        assertTrue(sparse.mixedModel().equationNonzeroCount()
            <= pedigree.sparseRelationshipMatrixInverse().nonzeroCount());
    }

    @Test
    void composesPedigreeAndOrdinaryRandomTermsInOneSparseFit() {
        double[] response = {0, 1, 2, 4, 5, 6, 8, 9, 10};
        double[][] fixed = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("A"),
            PedigreeIndividual.founder("B"),
            PedigreeIndividual.founder("C")));
        List<String> animals = List.of(
            "A", "A", "A", "B", "B", "B", "C", "C", "C");
        PedigreeRandomEffectTerm additive = PedigreeRandomEffectTerm.of(
            "animal", animals, pedigree);
        org.jlinalg.mixed.RandomEffectTerm batch =
            org.jlinalg.mixed.RandomEffectTerm.randomIntercept("batch",
                List.of("x", "y", "x", "x", "y", "x", "x", "y", "x"));
        org.jlinalg.mixed.SparseLinearMixedModelResult result =
            SparsePedigreeMixedModel.fit(response, fixed, List.of(additive),
                List.of(batch), RemlOptions.builder()
                    .initialVariances(10, 1, 2).maximumIterations(200).build(),
                BackendPolicy.CPU);

        assertEquals(3, result.varianceComponents().length);
        assertEquals(2, result.randomEffects().size());
        assertEquals(pedigree.size(), result.randomEffects("animal").estimates().length);
        assertTrue(Double.isFinite(result.logLikelihood()));
    }
}
