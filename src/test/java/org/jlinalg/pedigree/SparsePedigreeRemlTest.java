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
            new PedigreeIndividual("C", "A", "B")));
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
        assertArrayEquals(dense.predictionErrorVariances(),
            sparse.predictionErrorVariances(), 2e-5);
        assertArrayEquals(dense.reliabilities(), sparse.reliabilities(), 2e-5);
        assertEquals(sparse.predictionErrorVariances()[1],
            sparse.predictionErrorVariances(List.of("B")).get("B"), 0.0);
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

    @Test
    void diagonalPevScalesBeyondFormerCoefficientCutoff() {
        int individuals = 300;
        double[] response = new double[individuals * 2];
        double[][] fixed = new double[response.length][1];
        java.util.ArrayList<PedigreeIndividual> members = new java.util.ArrayList<>();
        java.util.ArrayList<String> observed = new java.util.ArrayList<>();
        for (int individual = 0; individual < individuals; individual++) {
            String id = "member-" + individual;
            members.add(PedigreeIndividual.founder(id));
            for (int repeat = 0; repeat < 2; repeat++) {
                int row = 2 * individual + repeat;
                fixed[row][0] = 1.0;
                response[row] = individual % 11 + (repeat == 0 ? -0.25 : 0.25);
                observed.add(id);
            }
        }
        SparsePedigreeRemlResult result = SparsePedigreeReml.fit(response,
            fixed, observed, Pedigree.of(members), RemlOptions.builder()
                .initialVariances(5.0, 1.0).maximumIterations(40).build(),
            BackendPolicy.CPU);

        assertEquals(individuals, result.predictionErrorVariances().length);
        assertTrue(java.util.Arrays.stream(result.predictionErrorVariances())
            .allMatch(Double::isFinite));
        assertTrue(java.util.Arrays.stream(result.reliabilities())
            .allMatch(value -> value >= 0.0 && value <= 1.0));
    }
}
