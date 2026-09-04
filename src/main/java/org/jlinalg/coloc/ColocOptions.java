/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.coloc;

/** Priors and overlap controls for SuSiE colocalization. */
public record ColocOptions(
        double trait1Prior,
        double trait2Prior,
        double sharedPrior,
        double minimumPosteriorOverlap,
        boolean trimByPosterior,
        double[] trait1PriorWeights,
        double[] trait2PriorWeights) {

    public ColocOptions {
        if (!(trait1Prior > 0.0 && trait1Prior < 1.0)
                || !(trait2Prior > 0.0 && trait2Prior < 1.0)
                || !(sharedPrior > 0.0 && sharedPrior < 1.0)
                || sharedPrior > trait1Prior || sharedPrior > trait2Prior
                || !(minimumPosteriorOverlap >= 0.0
                    && minimumPosteriorOverlap <= 1.0)) {
            throw new IllegalArgumentException("invalid colocalization options");
        }
        trait1PriorWeights = copyWeights(trait1PriorWeights, "trait 1");
        trait2PriorWeights = copyWeights(trait2PriorWeights, "trait 2");
    }

    /** Defaults from {@code coloc::coloc.susie}. */
    public static ColocOptions defaults() {
        return new ColocOptions(
            1e-4, 1e-4, 5e-6, 0.5, true, null, null);
    }

    @Override public double[] trait1PriorWeights() {
        return trait1PriorWeights == null ? null : trait1PriorWeights.clone();
    }

    @Override public double[] trait2PriorWeights() {
        return trait2PriorWeights == null ? null : trait2PriorWeights.clone();
    }

    double[] rawTrait1PriorWeights() { return trait1PriorWeights; }
    double[] rawTrait2PriorWeights() { return trait2PriorWeights; }

    private static double[] copyWeights(double[] weights, String label) {
        if (weights == null) return null;
        double[] result = weights.clone();
        if (result.length == 0) {
            throw new IllegalArgumentException(label + " prior weights are empty");
        }
        for (double value : result) {
            if (!(value > 0.0) || !Double.isFinite(value)) {
                throw new IllegalArgumentException(
                    label + " prior weights must be finite and positive");
            }
        }
        return result;
    }
}
