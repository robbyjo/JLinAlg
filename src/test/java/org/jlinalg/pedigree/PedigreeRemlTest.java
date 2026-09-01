/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class PedigreeRemlTest {
    @Test
    void repeatedFounderRecordsMatchBalancedAnimalModel() {
        double[] response = {0.0, 1.0, 2.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0};
        double[][] fixed = new double[response.length][1];
        for (double[] row : fixed) {
            row[0] = 1.0;
        }
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("A"),
            PedigreeIndividual.founder("B"),
            PedigreeIndividual.founder("C")));
        List<String> observed = List.of(
            "A", "A", "A", "B", "B", "B", "C", "C", "C");

        PedigreeRemlResult result = PedigreeReml.fit(
            response, fixed, observed, pedigree,
            RemlOptions.builder().initialVariances(10.0, 2.0).build(),
            BackendPolicy.CPU);

        assertTrue(result.reml().converged(),
            result.reml().convergenceMessage());
        assertEquals(47.0 / 3.0, result.additiveGeneticVariance(), 1e-5);
        assertEquals(1.0, result.residualVariance(), 1e-6);
        assertEquals(47.0 / 50.0, result.heritability(), 1e-6);
        assertEquals(5.0, result.reml().fixedEffects()[0], 1e-10);
        assertEquals(List.of("additive genetic", "residual"),
            result.reml().componentNames());
        assertEquals(-result.breedingValue("C"),
            result.breedingValue("A"), 1e-9);
        assertEquals(0.0, result.breedingValue("B"), 1e-9);
        for (double reliability : result.reliabilities()) {
            assertTrue(reliability > 0.0 && reliability < 1.0);
        }
        org.junit.jupiter.api.Assertions.assertArrayEquals(
            result.beta(), result.fixef(), 0);
        assertEquals(result.breedingValue("A"), result.ranef().get("A"));
        assertEquals(2, result.varCorr().size());
        org.junit.jupiter.api.Assertions.assertArrayEquals(
            response, add(result.fittedValues(), result.residuals()), 1e-12);
    }

    @Test
    void requiresOneKnownPedigreeIdPerObservation() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("A")));
        double[] response = {1.0, 2.0, 3.0};
        double[][] fixed = {{1.0}, {1.0}, {1.0}};

        assertThrows(IllegalArgumentException.class, () -> PedigreeReml.fit(
            response, fixed, List.of("A", "A"), pedigree));
        assertThrows(IllegalArgumentException.class, () -> PedigreeReml.fit(
            response, fixed, List.of("A", "A", "missing"), pedigree));
    }

    private static double[] add(double[] left, double[] right) {
        double[] result = new double[left.length];
        for (int index = 0; index < result.length; index++)
            result[index] = left[index] + right[index];
        return result;
    }
}
