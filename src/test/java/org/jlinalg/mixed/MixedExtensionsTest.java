/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class MixedExtensionsTest {
    @Test
    void sparseCorrelatedBlockProducesFiniteFit() {
        SparseCorrelatedRandomEffectBlock block = SparseCorrelatedRandomEffectBlock.of(
            "site", List.of("a", "a", "b", "b"), List.of("(Intercept)", "x"),
            new double[][] {{1, 0}, {1, 1}, {1, 0}, {1, 1}},
            new double[][] {{1, .2}, {.2, 1}});
        SparseLinearMixedModelResult result = SparseCorrelatedLinearMixedModel.fit(
            new double[] {1, 2, 1.5, 2.5},
            new double[] {1, 0, 1, 1, 1, 2, 1, 3}, 4, 2,
            List.of(block), RemlOptions.defaults(), BackendPolicy.CPU);
        assertTrue(result.converged() || Double.isFinite(result.logLikelihood()));
        assertTrue(Double.isFinite(result.beta()[0]));
    }

    @Test
    void profileIntervalFindsBothCrossings() {
        ProfileLikelihoodInterval interval = MixedModelProfile.interval(
            value -> -0.5 * value * value, 0.0, 0.0, 0.95,
            -4.0, 4.0, 32);
        assertTrue(interval.lowerFound());
        assertTrue(interval.upperFound());
        assertTrue(interval.lower() < 0.0 && interval.upper() > 0.0);
    }
}
