/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;

/** Sparse-equation Gaussian mixed model for grouped correlated blocks. */
public final class SparseCorrelatedLinearMixedModel {
    private SparseCorrelatedLinearMixedModel() { }

    /**
     * Fits correlated blocks without forming an observation-scale covariance.
     * Each block contributes one estimated variance scale; its supplied
     * covariance shape is repeated independently over groups.
     */
    public static SparseLinearMixedModelResult fit(double[] response,
            double[] fixedEffects, int rows, int columns,
            List<SparseCorrelatedRandomEffectBlock> blocks,
            RemlOptions options, BackendPolicy backendPolicy) {
        if (blocks == null || blocks.isEmpty())
            throw new IllegalArgumentException("at least one correlated block is required");
        List<RandomEffectTerm> terms = new ArrayList<>(blocks.size());
        List<SparsePrecisionMatrix> precision = new ArrayList<>(blocks.size());
        for (SparseCorrelatedRandomEffectBlock block : blocks) {
            if (block == null || block.term().observations() != rows)
                throw new IllegalArgumentException("correlated block rows must match response");
            terms.add(block.term()); precision.add(block.precision());
        }
        return SparseLinearMixedModel.fitWithPrecision(response, fixedEffects,
            rows, columns, terms, precision, options, backendPolicy);
    }
}
