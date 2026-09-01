/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Exposure and outcome associations aligned to the same effect allele. */
public record HarmonizedInstrument(
        String variantId,
        String effectAllele,
        String otherAllele,
        double exposureEffect,
        double exposureStandardError,
        double outcomeEffect,
        double outcomeStandardError,
        double exposureEffectAlleleFrequency,
        double outcomeEffectAlleleFrequency,
        boolean outcomeEffectFlipped,
        boolean strandComplemented) {
}
