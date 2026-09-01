/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import jdistlib.Beta;

/** Conventional frequency-dependent variant-set weights. */
public final class VariantWeights {
    private VariantWeights() { }

    /**
     * Returns a kernel-column weight whose square is a Beta(a,b) density.
     * Thus the resulting SKAT kernel uses the Beta density on its diagonal.
     */
    public static double betaKernel(
            double minorAlleleFrequency, double shape1, double shape2) {
        validate(minorAlleleFrequency, shape1, shape2);
        return Math.sqrt(Beta.density(
            minorAlleleFrequency, shape1, shape2, false));
    }

    /** Returns the direct Beta(a,b) burden coefficient. */
    public static double betaBurden(
            double minorAlleleFrequency, double shape1, double shape2) {
        validate(minorAlleleFrequency, shape1, shape2);
        return Beta.density(minorAlleleFrequency, shape1, shape2, false);
    }

    private static void validate(double maf, double shape1, double shape2) {
        if (!(maf > 0 && maf <= 0.5) || !(shape1 > 0) || !(shape2 > 0)
                || !Double.isFinite(maf) || !Double.isFinite(shape1)
                || !Double.isFinite(shape2))
            throw new IllegalArgumentException(
                "MAF must be in (0,0.5] and beta shapes must be positive");
    }
}
