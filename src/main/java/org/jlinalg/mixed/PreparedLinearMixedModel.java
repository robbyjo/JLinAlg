/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;

/**
 * Reusable Gaussian mixed-model structure for lme4-style response refits.
 * Fixed and random designs are validated and retained once.
 */
public final class PreparedLinearMixedModel {
    private final double[] fixedEffects;
    private final int observations;
    private final int fixedColumns;
    private final List<RandomEffectTerm> randomEffects;
    private final RemlOptions options;
    private final BackendPolicy backendPolicy;

    public PreparedLinearMixedModel(
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (fixedEffects == null || fixedEffects.length == 0)
            throw new IllegalArgumentException("fixed-effect design is required");
        this.observations = fixedEffects.length;
        this.fixedColumns = fixedEffects[0].length;
        this.fixedEffects = MatrixOps.rowMajor(fixedEffects, observations);
        this.randomEffects = List.copyOf(randomEffects);
        for (RandomEffectTerm term : this.randomEffects)
            if (term == null || term.observations() != observations)
                throw new IllegalArgumentException(
                    "random-effect rows must equal fixed-effect rows");
        this.options = Objects.requireNonNull(options, "options");
        this.backendPolicy = Objects.requireNonNull(backendPolicy, "backendPolicy");
    }

    public static PreparedLinearMixedModel of(
            double[][] fixedEffects, List<RandomEffectTerm> randomEffects) {
        return new PreparedLinearMixedModel(fixedEffects, randomEffects,
            RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits the retained model structure to a response. */
    public LinearMixedModelResult fit(double[] response) {
        validateResponse(response);
        return LinearMixedModel.fit(response, fixedEffects, observations,
            fixedColumns, randomEffects, options, backendPolicy);
    }

    /**
     * Refits a new response, warm-starting at a previous fit's variance
     * components while retaining every other optimization control.
     */
    public LinearMixedModelResult refit(
            LinearMixedModelResult previous, double[] response) {
        Objects.requireNonNull(previous, "previous");
        validateResponse(response);
        RemlOptions warm = options.toBuilder()
            .initialVariances(previous.reml().varianceComponents()).build();
        return LinearMixedModel.fit(response, fixedEffects, observations,
            fixedColumns, randomEffects, warm, backendPolicy);
    }

    public int observations() { return observations; }
    public int fixedEffectColumns() { return fixedColumns; }
    public List<RandomEffectTerm> randomEffects() { return randomEffects; }

    private void validateResponse(double[] response) {
        if (response == null || response.length != observations)
            throw new IllegalArgumentException(
                "response length must equal prepared observations");
        MatrixOps.requireFinite(response, "response");
    }
}
