/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.penalized;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

final class SelectionAwareInferenceTest {
    @Test
    void deterministicSplitProducesHeldOutInference() {
        double[] y = new double[40]; double[][] x = new double[40][2]; for (int row = 0; row < y.length; row++) { x[row][0] = row; x[row][1] = Math.sin(row); y[row] = 1.0 + 2.0 * x[row][0] + 0.1 * x[row][1]; }
        SelectionAwareInferenceResult result = SelectionAwarePenalizedInference.fit(y, x, 0.001, ElasticNetOptions.builder().alpha(1.0).build(), 0.6, 1e-8, org.jlinalg.compute.BackendPolicy.PREFERRED);
        assertTrue(result.selectionObservations() > 0); assertTrue(result.inferenceObservations() > 0); assertTrue(Double.isFinite(result.inferenceFit().associationStatistics().pValues()[0]));
    }
}
