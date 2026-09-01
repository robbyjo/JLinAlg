/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

class ExtendedMendelianRandomizationTest {
    @Test
    void steigerAndRobustMethodsRecoverForwardCommonEffect() {
        List<HarmonizedInstrument> instruments = instruments(false);
        SteigerResult steiger = SteigerFiltering.analyze(instruments,
            new double[] {10_000, 10_000, 10_000, 10_000, 10_000},
            new double[] {10_000, 10_000, 10_000, 10_000, 10_000});
        assertTrue(steiger.aggregateDirectionCorrect());
        assertEquals(5, steiger.retainedVariants().size());

        MrRapsResult raps = RobustMendelianRandomization.raps(instruments);
        assertEquals(0.5, raps.estimate().estimate(), 1e-6);
        assertEquals(0.5, ContaminationMixture.fit(instruments, 501)
            .estimate().estimate(), 0.02);
    }

    @Test
    void pressoFindsLargePleiotropicOutlier() {
        MrPressoResult result = MrPresso.analyze(instruments(true), 0.05);
        assertTrue(result.outlierVariants().contains("outlier"));
        assertEquals(0.5, result.outlierCorrectedEstimate().estimate(), 1e-10);
    }

    @Test
    void multivariableMrRecoversTwoDirectEffects() {
        List<MultivariableInstrument> values = new ArrayList<>();
        double[][] effects = {
            {0.10, 0.02}, {0.02, 0.12}, {0.15, -0.03},
            {-0.04, 0.16}, {0.12, 0.10}, {-0.08, 0.07}
        };
        for (int index = 0; index < effects.length; index++) {
            double outcome = 0.7 * effects[index][0] - 0.4 * effects[index][1];
            values.add(new MultivariableInstrument("v" + index,
                effects[index], new double[] {0.01, 0.01}, outcome, 0.01));
        }
        MultivariableMrResult result =
            MultivariableMendelianRandomization.ivw(values, List.of("x1", "x2"));
        assertEquals(0.7, result.beta()[0], 1e-12);
        assertEquals(-0.4, result.beta()[1], 1e-12);
    }

    @Test
    void overlapAwareReducesToErrorsInVariablesAtZeroCovariance() {
        List<HarmonizedInstrument> values = instruments(false);
        OverlapAwareMrResult result = OverlapAwareMendelianRandomization.ivw(
            values, new double[values.size()]);
        assertEquals(0.5, result.estimate().estimate(), 1e-9);
        assertTrue(Math.abs(WinnerCurseCorrection.correct(0.2, 0.03, 5.0)) < 0.2);
    }

    private static List<HarmonizedInstrument> instruments(boolean outlier) {
        List<HarmonizedInstrument> values = new ArrayList<>();
        for (int index = 1; index <= 5; index++) {
            double x = 0.05 * index;
            values.add(new HarmonizedInstrument("v" + index, "A", "C",
                x, 0.01, 0.5 * x, 0.01, 0.2, 0.2, false, false));
        }
        if (outlier) {
            values.add(new HarmonizedInstrument("outlier", "A", "C",
                0.2, 0.01, 1.0, 0.01, 0.2, 0.2, false, false));
        }
        return values;
    }
}
