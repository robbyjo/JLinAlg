/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Aligned-cohort additive-dosage quality statistics. */
public record VariantStatistics(
        int samples,
        int calledSamples,
        int missingSamples,
        double missingRate,
        double alternateAlleleCount,
        double alternateAlleleFrequency,
        double minorAlleleCount,
        double minorAlleleFrequency,
        double dosageMean,
        double dosageVariance) {

    /** Computes diploid additive statistics from finite dosages in [0,2]. */
    public static VariantStatistics of(VariantRecord variant) {
        double[] values = variant.dosagesView();
        double sum = 0;
        double sumSquares = 0;
        int called = 0;
        for (double dosage : values) {
            if (!Double.isFinite(dosage)) continue;
            if (dosage < -1e-10 || dosage > 2 + 1e-10)
                throw new IllegalArgumentException(
                    "variant '" + variant.id()
                    + "' contains an additive dosage outside [0,2]: " + dosage);
            double bounded = Math.max(0, Math.min(2, dosage));
            sum += bounded;
            sumSquares += bounded * bounded;
            called++;
        }
        int missing = values.length - called;
        if (called == 0) {
            return new VariantStatistics(values.length, 0, missing, 1,
                Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                Double.NaN, Double.NaN);
        }
        double alleleNumber = 2.0 * called;
        double frequency = sum / alleleNumber;
        double maf = Math.min(frequency, 1 - frequency);
        double mean = sum / called;
        double variance = called > 1
            ? Math.max(0, (sumSquares - sum * sum / called) / (called - 1))
            : 0;
        return new VariantStatistics(values.length, called, missing,
            missing / (double) values.length, sum, frequency,
            alleleNumber * maf, maf, mean, variance);
    }
}
