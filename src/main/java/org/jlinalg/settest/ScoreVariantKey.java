/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import java.util.Objects;

/** Allele-specific key used by conditional score summaries.
 *
 * <p>Keys are always forward-strand {@code CHROM:POS:REF:ALT}. Alignment may
 * exchange REF and ALT, but never guesses strand complements or performs
 * normalization/liftover.</p>
 */
public record ScoreVariantKey(String chromosome, long position,
        String reference, String alternate) {
    public ScoreVariantKey {
        if (chromosome == null || !chromosome.matches("[A-Za-z0-9_.-]+")
                || position < 1 || reference == null || alternate == null
                || !reference.matches("[ACGT]+")
                || !alternate.matches("[ACGT]+")
                || reference.equals(alternate))
            throw new IllegalArgumentException(
                "score variant key must be CHROM:POS:REF:ALT with normalized uppercase alleles");
    }

    /** Parse an exact, unpadded forward-strand allele key. */
    public static ScoreVariantKey parse(String value) {
        if (value == null || value.isBlank() || !value.equals(value.trim()))
            throw new IllegalArgumentException("score variant key is blank or padded");
        String[] fields = value.split(":", -1);
        if (fields.length != 4)
            throw new IllegalArgumentException(
                "score variant key must be CHROM:POS:REF:ALT: " + value);
        try {
            long position = Long.parseLong(fields[1]);
            if (!fields[1].equals(Long.toString(position)))
                throw new NumberFormatException("position is not canonical");
            return new ScoreVariantKey(fields[0], position, fields[2], fields[3]);
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                "score variant position is invalid: " + value, exception);
        }
    }

    /** The same variant with the dosage/effect allele reversed. */
    public ScoreVariantKey swapped() {
        return new ScoreVariantKey(chromosome, position, alternate, reference);
    }

    /** +1 for an exact orientation and -1 for a REF/ALT exchange. */
    public int alignmentSign(ScoreVariantKey other) {
        Objects.requireNonNull(other, "other");
        if (equals(other)) return 1;
        if (swapped().equals(other)) return -1;
        throw new IllegalArgumentException("allele mismatch: " + this
            + " versus " + other + "; strand complements are not inferred");
    }

    /** Orientation-independent identity used only to detect/alignment-match swaps. */
    public String unorientedIdentity() {
        String first = reference.compareTo(alternate) <= 0 ? reference : alternate;
        String second = reference.compareTo(alternate) <= 0 ? alternate : reference;
        return chromosome + ':' + position + ':' + first + ':' + second;
    }

    @Override public String toString() {
        return chromosome + ":" + position + ":" + reference + ":" + alternate;
    }
}
