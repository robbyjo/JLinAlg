/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.Objects;

/** One additive alternate-allele dosage vector with genomic metadata. */
public record VariantRecord(
        String id,
        String chromosome,
        long position,
        String referenceAllele,
        String alternateAllele,
        double[] dosages,
        double imputationQuality) {
    public VariantRecord {
        id = require(id, "variant ID");
        chromosome = chromosome == null ? "" : chromosome.trim();
        referenceAllele = referenceAllele == null ? "" : referenceAllele.trim();
        alternateAllele = alternateAllele == null ? "" : alternateAllele.trim();
        if (position < 0) throw new IllegalArgumentException("position cannot be negative");
        Objects.requireNonNull(dosages, "dosages");
        if (dosages.length == 0)
            throw new IllegalArgumentException("dosages cannot be empty");
        if (!Double.isNaN(imputationQuality)
                && (!Double.isFinite(imputationQuality)
                    || imputationQuality < 0 || imputationQuality > 1)) {
            throw new IllegalArgumentException(
                "imputation quality must be NaN or in [0,1]");
        }
        dosages = dosages.clone();
    }

    @Override public double[] dosages() { return dosages.clone(); }

    /** Internal package view used by block pipeline code without another copy. */
    double[] dosagesView() { return dosages; }

    public boolean genomic() {
        return !chromosome.isEmpty() && position > 0
            && !referenceAllele.isEmpty() && !alternateAllele.isEmpty();
    }

    private static String require(String value, String name) {
        if (value == null || value.isBlank())
            throw new IllegalArgumentException(name + " is required");
        return value.trim();
    }
}
