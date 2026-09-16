/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** Summary associations for one instrument, multiple exposures, and multiple outcomes. */
public record MultivariateInstrument(
        String variantId,
        double[] exposureEffects,
        double[] exposureStandardErrors,
        double[] outcomeEffects,
        double[] outcomeStandardErrors) {
    public MultivariateInstrument {
        if (variantId == null || variantId.isBlank()
                || exposureEffects == null || exposureStandardErrors == null
                || outcomeEffects == null || outcomeStandardErrors == null) {
            throw new IllegalArgumentException(
                "variant, exposure associations, and outcome associations are required");
        }
        exposureEffects = exposureEffects.clone();
        exposureStandardErrors = exposureStandardErrors.clone();
        outcomeEffects = outcomeEffects.clone();
        outcomeStandardErrors = outcomeStandardErrors.clone();
    }

    @Override public double[] exposureEffects() { return exposureEffects.clone(); }
    @Override public double[] exposureStandardErrors() {
        return exposureStandardErrors.clone();
    }
    @Override public double[] outcomeEffects() { return outcomeEffects.clone(); }
    @Override public double[] outcomeStandardErrors() {
        return outcomeStandardErrors.clone();
    }
}
