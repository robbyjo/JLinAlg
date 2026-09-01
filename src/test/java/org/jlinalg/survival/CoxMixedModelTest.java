/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.junit.jupiter.api.Test;

class CoxMixedModelTest {
    @Test
    void gaussianFrailtyProfilesPositiveVarianceAndHazardRatio() {
        Fixture fixture = clusteredFixture();
        CoxRandomEffectTerm frailty = CoxRandomEffectTerm.independent(
            RandomEffectTerm.randomIntercept("center", fixture.groups()));

        CoxMixedResult result = CoxMixedModel.fit(
            fixture.survival(), fixture.fixed(), List.of(frailty), null,
            mixedOptions(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(1, result.beta().length);
        assertTrue(result.beta()[0] > 0);
        assertTrue(result.hazardRatios()[0] > 1);
        assertTrue(result.standardErrors()[0] > 0);
        assertTrue(result.randomEffects("center").variance() > 0);
        assertEquals(12, result.randomEffects("center").modes().length);
        assertEquals(result.randomEffects("center"),
            result.ranef().get("center"));
        assertTrue(Double.isFinite(result.laplaceLogLikelihood()));
        assertTrue(!result.baselineHazard().isEmpty());
    }

    @Test
    void pedigreeFrailtyUsesRelationshipPrecisionAndNamedModes() {
        Fixture fixture = clusteredFixture();
        List<PedigreeIndividual> individuals = new ArrayList<>();
        individuals.add(PedigreeIndividual.founder("g0"));
        individuals.add(PedigreeIndividual.founder("g1"));
        individuals.add(new PedigreeIndividual("g2", "g0", "g1"));
        for (int group = 3; group < 12; group++)
            individuals.add(PedigreeIndividual.founder("g" + group));
        Pedigree pedigree = Pedigree.of(individuals);

        CoxPedigreeResult result = CoxPedigreeFrailty.fit(
            fixture.survival(), fixture.fixed(), fixture.groups(), pedigree,
            null, mixedOptions(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(12, result.ranef().size());
        assertEquals(result.ranef().get("g0"), result.frailty("g0"));
        assertTrue(result.frailtyVariance() > 0);
        assertTrue(result.hazardRatios()[0] > 1);
        assertTrue(Double.isFinite(result.pValues()[0]));
    }

    private static CoxMixedOptions mixedOptions() {
        return new CoxMixedOptions(CoxOptions.defaults(),
            new double[] {0.5}, 25, 0.02, 1e-6, 100);
    }

    private static Fixture clusteredFixture() {
        int groups = 12;
        int perGroup = 6;
        int rows = groups * perGroup;
        double[] time = new double[rows];
        boolean[] event = new boolean[rows];
        double[][] fixed = new double[rows][1];
        List<String> labels = new ArrayList<>();
        Random random = new Random(913);
        for (int group = 0; group < groups; group++) {
            double frailty = (group % 4 - 1.5) * 0.45;
            for (int member = 0; member < perGroup; member++) {
                int row = group * perGroup + member;
                double covariate = (member - 2.5) / 2.5;
                fixed[row][0] = covariate;
                double rate = Math.exp(0.75 * covariate + frailty);
                double failure = -Math.log(Math.max(1e-12,
                    random.nextDouble())) / rate;
                double censor = 1.7 + 0.15 * (member % 2);
                event[row] = failure <= censor;
                time[row] = event[row] ? failure : censor;
                labels.add("g" + group);
            }
        }
        return new Fixture(CoxSurvivalData.rightCensored(time, event),
            fixed, labels);
    }

    private record Fixture(
        CoxSurvivalData survival, double[][] fixed, List<String> groups) { }
}
