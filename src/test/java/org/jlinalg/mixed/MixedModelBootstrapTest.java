/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.BootstrapOptions;
import org.jlinalg.inference.GaussianBootstrapResult;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class MixedModelBootstrapTest {
    @Test
    void marginalAndConditionalSimulationAreReproducible() {
        Fixture fixture = fixture();
        double[][] first = MixedModelSimulation.simulate(
            fixture.model(), fixture.fit(), 3, 17,
            MixedModelSimulationMode.MARGINAL);
        double[][] second = MixedModelSimulation.simulate(
            fixture.model(), fixture.fit(), 3, 17,
            MixedModelSimulationMode.MARGINAL);
        for (int simulation = 0; simulation < first.length; simulation++)
            assertArrayEquals(first[simulation], second[simulation], 0);

        double[][] conditional = MixedModelSimulation.simulate(
            fixture.model(), fixture.fit(), 2_000, 91,
            MixedModelSimulationMode.CONDITIONAL);
        double[] means = means(conditional);
        assertArrayEquals(fixture.fit().conditionalFittedValues(), means, 0.12);

        double[][] marginal = MixedModelSimulation.simulate(
            fixture.model(), fixture.fit(), 2_000, 92,
            MixedModelSimulationMode.MARGINAL);
        double[] marginalMeans = means(marginal);
        for (double mean : marginalMeans)
            assertEquals(fixture.fit().beta()[0], mean, 0.18);
    }

    @Test
    void parametricBootstrapReturnsOrderedIntervalsAndFailures() {
        Fixture fixture = fixture();
        BootstrapOptions options = new BootstrapOptions(16, 0.90, 1234, 2);

        GaussianBootstrapResult result = MixedModelBootstrap.bootstrap(
            fixture.model(), fixture.fit(), options);
        GaussianBootstrapResult serial = MixedModelBootstrap.bootstrap(
            fixture.model(), fixture.fit(),
            new BootstrapOptions(16, 0.90, 1234, 1));

        assertEquals(16, result.requestedSimulations());
        assertEquals(16, result.successfulSimulations()
            + result.failedSimulations());
        assertTrue(result.successfulSimulations() >= 8,
            result.failures().toString());
        assertEquals(1, result.fixedEffectSummaries().size());
        assertEquals(2, result.varianceComponentSummaries().size());
        assertTrue(Double.isFinite(
            result.fixedEffectSummaries().get(0).standardError()));
        assertTrue(result.fixedEffectSummaries().get(0).lower()
            <= result.fixedEffectSummaries().get(0).upper());
        assertEquals(List.of("subject", "residual"),
            result.varianceComponentSummaries().stream()
                .map(value -> value.name()).toList());
        assertEquals(result.failures(), serial.failures());
        double[][] parallelDraws = result.fixedEffectReplicates();
        double[][] serialDraws = serial.fixedEffectReplicates();
        assertEquals(parallelDraws.length, serialDraws.length);
        for (int simulation = 0; simulation < parallelDraws.length; simulation++)
            assertArrayEquals(parallelDraws[simulation],
                serialDraws[simulation], 0);
    }

    private static Fixture fixture() {
        int groups = 10;
        int repeats = 4;
        double[][] fixed = new double[groups * repeats][1];
        double[] response = new double[groups * repeats];
        List<String> labels = new ArrayList<>();
        double[] residualPattern = {-0.45, -0.05, 0.10, 0.40};
        for (int group = 0; group < groups; group++) {
            double randomEffect = (group - 4.5) * 0.55;
            for (int repeat = 0; repeat < repeats; repeat++) {
                int row = group * repeats + repeat;
                fixed[row][0] = 1;
                response[row] = 3 + randomEffect + residualPattern[repeat];
                labels.add("g" + group);
            }
        }
        PreparedLinearMixedModel model = new PreparedLinearMixedModel(
            fixed, List.of(RandomEffectTerm.randomIntercept(
                "subject", labels)),
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
        PreparedLinearMixedModel model, LinearMixedModelResult fit) { }
}
