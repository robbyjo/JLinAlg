/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class NonlinearModelTest {
    private static NonlinearMeanFunction exponentialMean() {
        return (parameters, row) -> {
            double x = row * 0.25;
            double value = parameters[0] * Math.exp(parameters[1] * x);
            return new NonlinearMeanFunction.Evaluation(value,
                new double[] {Math.exp(parameters[1] * x), x * value});
        };
    }

    @Test
    void fixedModelRecoversAnalyticExponentialMean() {
        double[] response = new double[12];
        for (int row = 0; row < response.length; row++)
            response[row] = 2.0 * Math.exp(0.4 * row * 0.25);
        NonlinearFitResult result = NonlinearFixedModel.fit(response,
            new double[] {1.0, 0.1}, exponentialMean(),
            NonlinearModelOptions.defaults(), BackendPolicy.CPU);
        assertEquals(2.0, result.beta()[0], 1e-7);
        assertEquals(0.4, result.beta()[1], 1e-7);
        assertTrue(result.converged());
    }

    @Test
    void ordinarySparseMixedModelUsesNonlinearLinearization() {
        double[] response = new double[12];
        List<String> groups = List.of("a", "a", "a", "a", "b", "b", "b", "b", "c", "c", "c", "c");
        double[] random = {0.25, -0.15, 0.05};
        for (int row = 0; row < response.length; row++)
            response[row] = 1.5 * Math.exp(0.3 * row * 0.2)
                + random[row / 4];
        NonlinearMixedFitResult result = NonlinearMixedModel.fit(response,
            new double[] {1.0, 0.1}, exponentialMean(),
            List.of(RandomEffectTerm.randomIntercept("1|group", groups)),
            RemlOptions.defaults(), NonlinearModelOptions.defaults(), BackendPolicy.CPU);
        assertTrue(Double.isFinite(result.objective()));
        assertTrue(result.linearizedModel().randomEffects("1|group").estimates().length == 3);
    }

    @Test
    void pedigreeSparseMixedModelAcceptsNonlinearMean() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("sire"),
            PedigreeIndividual.founder("dam"),
            new PedigreeIndividual("child", "sire", "dam")));
        List<String> observed = List.of("sire", "sire", "dam", "dam", "child", "child", "child", "sire", "dam");
        double[] response = new double[observed.size()];
        for (int row = 0; row < response.length; row++)
            response[row] = 1.2 * Math.exp(0.25 * row * 0.15)
                + (observed.get(row).equals("child") ? 0.2 : 0.0);
        NonlinearMixedFitResult result = NonlinearMixedModel.fitPedigree(response,
            new double[] {1.0, 0.1}, exponentialMean(),
            List.of(PedigreeRandomEffectTerm.of("additive", observed, pedigree)),
            List.of(), RemlOptions.defaults(), NonlinearModelOptions.defaults(), BackendPolicy.CPU);
        assertTrue(Double.isFinite(result.objective()));
        assertEquals(3, result.linearizedModel().randomEffects("additive").estimates().length);
    }
}
