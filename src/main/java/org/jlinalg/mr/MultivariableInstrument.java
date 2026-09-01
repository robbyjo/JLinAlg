/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** Summary associations for one instrument and multiple exposures. */
public record MultivariableInstrument(
        String variantId,
        double[] exposureEffects,
        double[] exposureStandardErrors,
        double outcomeEffect,
        double outcomeStandardError) {
    public MultivariableInstrument {
        if (variantId == null || variantId.isBlank()
                || exposureEffects == null || exposureStandardErrors == null) {
            throw new IllegalArgumentException("variant and exposure associations are required");
        }
        exposureEffects = exposureEffects.clone();
        exposureStandardErrors = exposureStandardErrors.clone();
    }
    @Override public double[] exposureEffects() { return exposureEffects.clone(); }
    @Override public double[] exposureStandardErrors() {
        return exposureStandardErrors.clone();
    }
}
