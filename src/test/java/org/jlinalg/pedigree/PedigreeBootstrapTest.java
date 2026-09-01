/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.BootstrapOptions;
import org.jlinalg.inference.GaussianBootstrapResult;
import org.jlinalg.mixed.MixedModelSimulationMode;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class PedigreeBootstrapTest {
    @Test
    void pedigreeSimulationUsesConditionalModesAndIsReproducible() {
        Fixture fixture = fixture();
        double[][] first = PedigreeSimulation.simulate(
            fixture.model(), fixture.fit(), 3, 71,
            MixedModelSimulationMode.MARGINAL);
        double[][] second = PedigreeSimulation.simulate(
            fixture.model(), fixture.fit(), 3, 71,
            MixedModelSimulationMode.MARGINAL);
        for (int simulation = 0; simulation < first.length; simulation++)
            assertArrayEquals(first[simulation], second[simulation], 0);

        double[][] conditional = PedigreeSimulation.simulate(
            fixture.model(), fixture.fit(), 2_000, 72,
            MixedModelSimulationMode.CONDITIONAL);
        double[] means = means(conditional);
        double fixed = fixture.fit().beta()[0];
        List<String> observed = fixture.model().observationIndividualIds();
        for (int row = 0; row < means.length; row++)
            assertEquals(fixed + fixture.fit().breedingValue(observed.get(row)),
                means[row], 0.12);
    }

    @Test
    void pedigreeParametricBootstrapUsesWarmRefits() {
        Fixture fixture = fixture();
        GaussianBootstrapResult result = PedigreeBootstrap.bootstrap(
            fixture.model(), fixture.fit(),
            new BootstrapOptions(12, 0.90, 808, 2));

        assertEquals(12, result.successfulSimulations()
            + result.failedSimulations());
        assertTrue(result.successfulSimulations() >= 6,
            result.failures().toString());
        assertEquals(List.of("additive genetic", "residual"),
            result.varianceComponentSummaries().stream()
                .map(value -> value.name()).toList());
        assertTrue(Double.isFinite(
            result.fixedEffectSummaries().get(0).standardError()));
    }

    private static Fixture fixture() {
        int animals = 10;
        int repeats = 4;
        List<PedigreeIndividual> entries = new ArrayList<>();
        List<String> observed = new ArrayList<>();
        double[][] fixed = new double[animals * repeats][1];
        double[] response = new double[animals * repeats];
        double[] residualPattern = {-0.35, -0.10, 0.15, 0.30};
        for (int animal = 0; animal < animals; animal++) {
            String id = "a" + animal;
            entries.add(PedigreeIndividual.founder(id));
            double breeding = (animal - 4.5) * 0.50;
            for (int repeat = 0; repeat < repeats; repeat++) {
                int row = animal * repeats + repeat;
                fixed[row][0] = 1;
                response[row] = 2 + breeding + residualPattern[repeat];
                observed.add(id);
            }
        }
        PreparedPedigreeReml model = new PreparedPedigreeReml(
            fixed, observed, Pedigree.of(entries),
            RemlOptions.builder().initialVariances(2, 0.2).build(),
            BackendPolicy.CPU);
        return new Fixture(model, model.fit(response));
    }

    private static double[] means(double[][] simulations) {
        double[] result = new double[simulations[0].length];
        for (double[] simulation : simulations)
            for (int row = 0; row < result.length; row++)
                result[row] += simulation[row] / simulations.length;
        return result;
    }

    private record Fixture(
        PreparedPedigreeReml model, PedigreeRemlResult fit) { }
}
