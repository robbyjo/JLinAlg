/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;

class AlleleHarmonizerTest {
    @Test
    void alignsSwapsAndFrequencyResolvedPalindromes() {
        List<SummaryAssociation> exposure = List.of(
            association("rs1", "a", "c", 0.10, 0.01, 0.20),
            association("rs2", "A", "T", 0.20, 0.02, 0.10),
            association("rs3", "C", "G", 0.30, 0.03, 0.50),
            association("rs4", "G", "A", 0.40, 0.04, 0.25),
            association("rs5", "A", "C", 0.00, 0.01, 0.20),
            association("rs6", "A", "C", 0.10, 0.01, 0.20),
            association("rs7", "A", "C", 0.10, 0.01, 0.20),
            association("rs7", "A", "C", 0.11, 0.01, 0.20),
            association("rs8", "A", "C", 0.10, 0.01, 0.20));
        List<SummaryAssociation> outcome = List.of(
            association("rs1", "C", "A", -0.05, 0.01, 0.80),
            association("rs2", "A", "T", -0.06, 0.01, 0.90),
            association("rs3", "C", "G", 0.07, 0.01, 0.50),
            association("rs4", "G", "C", 0.08, 0.01, 0.25),
            association("rs5", "A", "C", 0.01, 0.01, 0.20),
            association("rs7", "A", "C", 0.02, 0.01, 0.20),
            association("rs8", "A", "C", 0.02, 0.01, 0.20),
            association("rs8", "A", "C", 0.03, 0.01, 0.20));

        HarmonizationResult result = AlleleHarmonizer.harmonize(exposure, outcome);

        assertEquals(2, result.instruments().size());
        HarmonizedInstrument first = result.instruments().get(0);
        assertEquals("A", first.effectAllele());
        assertEquals(0.05, first.outcomeEffect(), 1e-15);
        assertEquals(0.20, first.outcomeEffectAlleleFrequency(), 1e-15);
        assertTrue(first.outcomeEffectFlipped());
        HarmonizedInstrument palindrome = result.instruments().get(1);
        assertEquals(0.06, palindrome.outcomeEffect(), 1e-15);
        assertEquals(0.10, palindrome.outcomeEffectAlleleFrequency(), 1e-15);
        assertTrue(palindrome.outcomeEffectFlipped());

        Map<String, HarmonizationExclusionReason> reasons = result.exclusions()
            .stream().collect(Collectors.toMap(
                HarmonizationExclusion::variantId,
                HarmonizationExclusion::reason));
        assertEquals(HarmonizationExclusionReason.PALINDROMIC_AMBIGUOUS,
            reasons.get("rs3"));
        assertEquals(HarmonizationExclusionReason.ALLELE_MISMATCH,
            reasons.get("rs4"));
        assertEquals(HarmonizationExclusionReason.ZERO_EXPOSURE_EFFECT,
            reasons.get("rs5"));
        assertEquals(HarmonizationExclusionReason.MISSING_OUTCOME,
            reasons.get("rs6"));
        assertEquals(HarmonizationExclusionReason.DUPLICATE_EXPOSURE,
            reasons.get("rs7"));
        assertEquals(HarmonizationExclusionReason.DUPLICATE_OUTCOME,
            reasons.get("rs8"));
    }

    @Test
    void validatesAssociationFields() {
        assertThrows(IllegalArgumentException.class,
            () -> new SummaryAssociation("rs1", "A", "A", 0.1, 0.01));
        assertThrows(IllegalArgumentException.class,
            () -> new SummaryAssociation("rs1", "I", "D", 0.1, 0.01));
        assertThrows(IllegalArgumentException.class,
            () -> new SummaryAssociation("rs1", "A", "C", 0.1, 0.0));
        assertThrows(IllegalArgumentException.class,
            () -> association("rs1", "A", "C", 0.1, 0.01, 1.1));
    }

    private static SummaryAssociation association(
            String id, String effectAllele, String otherAllele,
            double effect, double standardError, double frequency) {
        return new SummaryAssociation(id, effectAllele, otherAllele,
            effect, standardError, frequency);
    }
}
