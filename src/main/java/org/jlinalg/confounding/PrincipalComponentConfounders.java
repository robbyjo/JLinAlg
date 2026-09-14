/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** PCA sample-factor extraction for feature-by-sample matrices. */
public final class PrincipalComponentConfounders {
    private PrincipalComponentConfounders() { }

    /** PCA settings matching {@code prcomp(t(data))} orientation. */
    public record Options(int factors, boolean center, boolean scale) {
        public Options {
            if (factors < 1) throw new IllegalArgumentException("factors must be positive");
        }
        public static Options defaults(int factors) {
            return new Options(factors, true, false);
        }
    }

    /** Fits with the portable CPU backend. */
    public static LatentFactorResult fit(double[][] featuresBySamples, Options options) {
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            return fit(featuresBySamples, options, context.backend());
        }
    }

    /** Fits using the requested linear-algebra backend. */
    public static LatentFactorResult fit(double[][] featuresBySamples,
            Options options, ComputeBackend backend) {
        double[][] prepared = ConfounderMath.prepare(featuresBySamples,
            options.center(), options.scale());
        ConfounderMath.Svd factors = ConfounderMath.sampleSvd(
            prepared, options.factors(), backend);
        double[][] adjusted = ConfounderMath.removeFactors(
            featuresBySamples, factors.factors());
        return new LatentFactorResult(factors.factors(), factors.loadings(),
            factors.varianceExplained(), adjusted, 1, true, List.of());
    }
}
