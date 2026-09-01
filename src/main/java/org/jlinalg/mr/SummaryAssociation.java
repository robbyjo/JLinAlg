/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.Locale;

/** A biallelic SNP-trait association from summary statistics. */
public record SummaryAssociation(
        String variantId,
        String effectAllele,
        String otherAllele,
        double effect,
        double standardError,
        double effectAlleleFrequency) {

    /** Creates an association without effect-allele frequency. */
    public SummaryAssociation(
            String variantId,
            String effectAllele,
            String otherAllele,
            double effect,
            double standardError) {
        this(variantId, effectAllele, otherAllele, effect, standardError,
            Double.NaN);
    }

    /** Normalizes alleles and validates numerical fields. */
    public SummaryAssociation {
        if (variantId == null || variantId.isBlank()) {
            throw new IllegalArgumentException("variantId must not be null or blank");
        }
        effectAllele = normalizeAllele(effectAllele, "effectAllele");
        otherAllele = normalizeAllele(otherAllele, "otherAllele");
        if (effectAllele.equals(otherAllele)) {
            throw new IllegalArgumentException("effect and other allele must differ");
        }
        if (!Double.isFinite(effect)) {
            throw new IllegalArgumentException("effect must be finite");
        }
        if (!(standardError > 0.0) || !Double.isFinite(standardError)) {
            throw new IllegalArgumentException(
                "standardError must be finite and positive");
        }
        if (!Double.isNaN(effectAlleleFrequency)
                && (!(effectAlleleFrequency >= 0.0)
                    || !(effectAlleleFrequency <= 1.0)
                    || !Double.isFinite(effectAlleleFrequency))) {
            throw new IllegalArgumentException(
                "effectAlleleFrequency must be NaN or lie in [0, 1]");
        }
    }

    /** Returns whether an effect-allele frequency was supplied. */
    public boolean hasEffectAlleleFrequency() {
        return !Double.isNaN(effectAlleleFrequency);
    }

    private static String normalizeAllele(String allele, String name) {
        if (allele == null) {
            throw new IllegalArgumentException(name + " is required");
        }
        String normalized = allele.toUpperCase(Locale.ROOT);
        if (normalized.length() != 1 || "ACGT".indexOf(normalized.charAt(0)) < 0) {
            throw new IllegalArgumentException(name + " must be one of A, C, G, or T");
        }
        return normalized;
    }
}
