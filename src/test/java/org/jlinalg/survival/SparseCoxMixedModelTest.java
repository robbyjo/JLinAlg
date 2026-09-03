/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.junit.jupiter.api.Test;

class SparseCoxMixedModelTest {
    @Test
    void vanishingFrailtiesRecoverFixedCoxAndPreparedStateIsReusable() {
        double[] time = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        boolean[] event = {
            true, false, true, true, false, true,
            true, false, true, true, false, true};
        double[][] first = new double[time.length][1];
        double[][] second = new double[time.length][1];
        List<String> ids = new ArrayList<>();
        List<PedigreeIndividual> founders = new ArrayList<>();
        List<String> groups = new ArrayList<>();
        for (int row = 0; row < time.length; row++) {
            first[row][0] = (row % 5) - 2.0;
            second[row][0] = ((row * 3) % 7) - 3.0;
            String id = "id" + row;
            ids.add(id);
            founders.add(new PedigreeIndividual(id, null, null));
            groups.add(row < 6 ? "a" : "b");
        }
        CoxSurvivalData survival = CoxSurvivalData.rightCensored(time, event);
        PedigreeRandomEffectTerm genetic =
            PedigreeRandomEffectTerm.ofUninbred("animal", ids, founders);
        RandomEffectTerm batch = RandomEffectTerm.randomIntercept(
            "batch", groups);
        CoxMixedOptions options = new CoxMixedOptions(CoxOptions.defaults(),
            new double[] {0.1, 0.1}, 20, 1e-4, 1e-8, 1e4);

        CoxResult fixedFirst = CoxRegression.fit(survival, first, null,
            CoxOptions.defaults(), BackendPolicy.CPU);
        CoxResult fixedSecond = CoxRegression.fit(survival, second, null,
            CoxOptions.defaults(), BackendPolicy.CPU);
        try (SparseCoxMixedModel.Prepared prepared =
                SparseCoxMixedModel.prepare(survival, genetic,
                    List.of(batch), options, BackendPolicy.PREFERRED)) {
            CoxMixedResult sparseFirst = prepared.fitAtVariances(
                first, null, 1e-8, 1e-8);
            CoxMixedResult sparseSecond = prepared.fitAtVariances(
                second, null, 1e-8, 1e-8);
            assertTrue(sparseFirst.converged());
            assertTrue(sparseSecond.converged());
            assertEquals(fixedFirst.beta()[0], sparseFirst.beta()[0], 2e-6);
            assertEquals(fixedSecond.beta()[0], sparseSecond.beta()[0], 2e-6);
            assertEquals(fixedFirst.standardErrors()[0],
                sparseFirst.standardErrors()[0], 2e-6);
            assertEquals(fixedSecond.standardErrors()[0],
                sparseSecond.standardErrors()[0], 2e-6);
        }
    }
}
