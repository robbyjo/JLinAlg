/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.coloc;

import java.util.Arrays;
import java.util.List;

/** Pairwise signal and variant-level posteriors from SuSiE colocalization. */
public final class ColocSusieResult {
    private final List<String> commonVariants;
    private final List<ColocSignalPair> signalPairs;
    private final double[] sharedVariantPosterior;
    private final ColocOptions options;
    private final int skippedSignalPairs;

    ColocSusieResult(
            List<String> commonVariants,
            List<ColocSignalPair> signalPairs,
            double[] sharedVariantPosterior,
            ColocOptions options,
            int skippedSignalPairs) {
        this.commonVariants = List.copyOf(commonVariants);
        this.signalPairs = List.copyOf(signalPairs);
        this.sharedVariantPosterior = sharedVariantPosterior.clone();
        this.options = options;
        this.skippedSignalPairs = skippedSignalPairs;
    }

    public List<String> commonVariants() { return commonVariants; }
    public List<ColocSignalPair> signalPairs() { return signalPairs; }

    /**
     * Row-major signal-pair by common-variant posterior conditional on H4.
     */
    public double[] sharedVariantPosterior() {
        return sharedVariantPosterior.clone();
    }

    /** Returns the H4-conditional variant posterior for one signal pair. */
    public double[] sharedVariantPosterior(int signalPair) {
        if (signalPair < 0 || signalPair >= signalPairs.size()) {
            throw new IndexOutOfBoundsException("signal pair index is invalid");
        }
        int start = signalPair * commonVariants.size();
        return Arrays.copyOfRange(sharedVariantPosterior,
            start, start + commonVariants.size());
    }

    public ColocOptions options() { return options; }
    public int skippedSignalPairs() { return skippedSignalPairs; }
}
