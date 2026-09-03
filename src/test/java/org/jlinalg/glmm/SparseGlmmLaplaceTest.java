/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.GammCovariances;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

/** Regression tests for coefficient-space sparse Laplace GLMM fitting. */
final class SparseGlmmLaplaceTest {
    @Test
    void groupedBinomialFitMatchesLme4LaplaceReference() {
        int groups = 20;
        int perGroup = 15;
        int observations = groups * perGroup;
        double[] response = new double[observations];
        double[] fixed = new double[observations * 2];
        List<String> group = new ArrayList<>(observations);
        long state = 1_234_567L;
        for (int row = 0; row < observations; row++) {
            int cluster = row / perGroup;
            double x = -1.0 + 2.0 * (row % perGroup) / (perGroup - 1.0);
            double eta = -0.2 + 0.8 * x + 1.1 * Math.sin(1.3 * cluster);
            double probability = 1.0 / (1.0 + Math.exp(-eta));
            state = 48_271L * state % 2_147_483_647L;
            response[row] = state / 2_147_483_647.0 < probability ? 1.0 : 0.0;
            fixed[row * 2] = 1.0;
            fixed[row * 2 + 1] = x;
            group.add(Integer.toString(cluster));
        }

        GlmmLaplaceResult fit = SparseGlmmLaplace.fit(
            response, fixed, observations, 2, GlmFamilies.binomial(),
            List.of(RandomEffectTerm.randomIntercept("group", group)),
            null, null,
            new GlmmLaplaceOptions(30, 100, 1e-7, 1.0,
                1e-8, 100.0, null), BackendPolicy.CPU);

        assertTrue(fit.converged());
        assertEquals(-0.08778265, fit.beta()[0], 0.01);
        assertEquals(0.7488305, fit.beta()[1], 0.04);
        assertEquals(0.8167182, fit.varianceComponents()[0], 0.02);
        assertEquals(-191.0959, fit.marginalLogLikelihood(), 0.05);
        assertEquals(0.38780058844855192, fit.fittedMeans()[0], 0.02);
        assertEquals(0.37318139192372202, fit.fittedMeans()[149], 0.02);
        assertEquals(0.54865646216866559, fit.fittedMeans()[299], 0.02);
    }

    @Test
    void preparedScanReusesStructureAndAcceptsWarmStart() {
        double[] response = {0, 1, 0, 1, 1, 0, 1, 0};
        double[] fixed = {
            1, -1.0, 1, -0.7, 1, -0.3, 1, 0.0,
            1, 0.2, 1, 0.5, 1, 0.8, 1, 1.1};
        RandomEffectTerm term = RandomEffectTerm.randomIntercept(
            "family", List.of("a", "a", "b", "b", "c", "c", "d", "d"));
        GlmmLaplaceOptions controls = new GlmmLaplaceOptions(
            8, 40, 1e-5, 1.0, 1e-6, 100.0, null);

        try (SparseGlmmLaplace.Prepared prepared = SparseGlmmLaplace.prepare(
                response.length, GlmFamilies.binomial(), List.of(term),
                controls, BackendPolicy.CPU)) {
            prepared.warmStart(0.5);
            GlmmLaplaceResult first = prepared.fit(response, fixed, 2);
            GlmmLaplaceResult second = prepared.fit(response, fixed, 2);
            assertEquals(first.beta()[1], second.beta()[1], 1e-5);
            assertTrue(Double.isFinite(second.marginalLogLikelihood()));
        }
    }

    @Test
    void pedigreePrecisionMatchesDenseRelationshipReference() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("sire"),
            PedigreeIndividual.founder("dam"),
            new PedigreeIndividual("a", "sire", "dam"),
            new PedigreeIndividual("b", "sire", "dam"),
            new PedigreeIndividual("c", "a", "b")));
        List<String> animals = new ArrayList<>();
        int observations = 75;
        double[] response = new double[observations];
        double[] fixed = new double[observations * 2];
        String[] observed = {"a", "b", "c"};
        long state = 8_675_309L;
        for (int row = 0; row < observations; row++) {
            String animal = observed[row % observed.length];
            double x = -1.0 + 2.0 * (row % 25) / 24.0;
            double animalEffect = switch (animal) {
                case "a" -> -0.7;
                case "b" -> 0.2;
                default -> 0.8;
            };
            double probability = 1.0 / (1.0
                + Math.exp(-(-0.1 + 0.6 * x + animalEffect)));
            state = 48_271L * state % 2_147_483_647L;
            response[row] = state / 2_147_483_647.0 < probability ? 1.0 : 0.0;
            fixed[row * 2] = 1.0;
            fixed[row * 2 + 1] = x;
            animals.add(animal);
        }
        GlmmLaplaceOptions controls = new GlmmLaplaceOptions(
            20, 80, 1e-6, 1.0, 1e-7, 100.0, null);
        VarianceComponent denseComponent = GammCovariances.pedigree(
            "animal", pedigree, animals);
        GlmmLaplaceResult dense = GlmmLaplace.fit(response, fixed,
            observations, 2, GlmFamilies.binomial(),
            List.of(denseComponent), null, null, controls, BackendPolicy.CPU);
        PedigreeRandomEffectTerm sparseTerm = PedigreeRandomEffectTerm.of(
            "animal", animals, pedigree);
        GlmmLaplaceResult sparse = SparseGlmmLaplace.fitWithPrecision(
            response, fixed, observations, 2, GlmFamilies.binomial(),
            List.of(sparseTerm.randomEffect()), List.of(sparseTerm.precision()),
            null, null, controls, BackendPolicy.CPU);

        assertEquals(dense.beta()[0], sparse.beta()[0], 2e-3);
        assertEquals(dense.beta()[1], sparse.beta()[1], 2e-3);
        assertEquals(dense.varianceComponents()[0],
            sparse.varianceComponents()[0], 2e-2);
        assertEquals(dense.marginalLogLikelihood(),
            sparse.marginalLogLikelihood(), 2e-2);
    }
}
