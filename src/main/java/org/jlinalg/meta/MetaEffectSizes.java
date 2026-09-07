/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/** First-order effect-size constructions used by meta-analysis workflows. */
public final class MetaEffectSizes {
    private MetaEffectSizes() { }

    public static MetaEffectSize meanDifference(
            double mean1, double mean2, double standardDeviation1,
            double standardDeviation2, int sampleSize1, int sampleSize2) {
        requireSamples(standardDeviation1, sampleSize1);
        requireSamples(standardDeviation2, sampleSize2);
        return new MetaEffectSize("mean-difference", mean1 - mean2,
            Math.sqrt(standardDeviation1 * standardDeviation1 / sampleSize1
                + standardDeviation2 * standardDeviation2 / sampleSize2));
    }

    /** Hedges g with the small-sample J correction. */
    public static MetaEffectSize standardizedMeanDifference(
            double mean1, double mean2, double standardDeviation1,
            double standardDeviation2, int sampleSize1, int sampleSize2) {
        requireSamples(standardDeviation1, sampleSize1);
        requireSamples(standardDeviation2, sampleSize2);
        int degrees = sampleSize1 + sampleSize2 - 2;
        if (degrees < 2) throw new IllegalArgumentException(
            "standardized mean difference needs at least two residual degrees of freedom");
        double pooled = Math.sqrt(((sampleSize1 - 1.0) * standardDeviation1 * standardDeviation1
            + (sampleSize2 - 1.0) * standardDeviation2 * standardDeviation2) / degrees);
        double d = (mean1 - mean2) / pooled;
        double j = 1.0 - 3.0 / (4.0 * degrees - 1.0);
        double g = j * d;
        double variance = (sampleSize1 + sampleSize2) / (double) (sampleSize1 * sampleSize2)
            + d * d / (2.0 * degrees);
        return new MetaEffectSize("hedges-g", g, Math.sqrt(j * j * variance));
    }

    public static MetaEffectSize logOddsRatio(
            int exposedCases, int exposedNonCases,
            int controlCases, int controlNonCases) {
        requireCell(exposedCases); requireCell(exposedNonCases);
        requireCell(controlCases); requireCell(controlNonCases);
        double effect = Math.log((double) exposedCases * controlNonCases
            / ((double) exposedNonCases * controlCases));
        double variance = 1.0 / exposedCases + 1.0 / exposedNonCases
            + 1.0 / controlCases + 1.0 / controlNonCases;
        return new MetaEffectSize("log-odds-ratio", effect, Math.sqrt(variance));
    }

    public static MetaEffectSize logRiskRatio(
            int exposedCases, int exposedTotal,
            int controlCases, int controlTotal) {
        requireCell(exposedCases); requireCell(controlCases);
        if (exposedTotal <= exposedCases || controlTotal <= controlCases)
            throw new IllegalArgumentException("risk totals must exceed case counts");
        double effect = Math.log((double) exposedCases * controlTotal
            / ((double) controlCases * exposedTotal));
        double variance = 1.0 / exposedCases - 1.0 / exposedTotal
            + 1.0 / controlCases - 1.0 / controlTotal;
        return new MetaEffectSize("log-risk-ratio", effect, Math.sqrt(variance));
    }

    public static MetaEffectSize fisherZCorrelation(
            double correlation, int sampleSize) {
        if (!(correlation > -1.0 && correlation < 1.0)
                || !Double.isFinite(correlation) || sampleSize < 4)
            throw new IllegalArgumentException(
                "correlation must be in (-1,1) and sample size at least four");
        return new MetaEffectSize("fisher-z-correlation",
            0.5 * Math.log((1.0 + correlation) / (1.0 - correlation)),
            1.0 / Math.sqrt(sampleSize - 3.0));
    }

    private static void requireSamples(double standardDeviation, int sampleSize) {
        if (!(standardDeviation > 0.0) || !Double.isFinite(standardDeviation)
                || sampleSize < 2)
            throw new IllegalArgumentException("standard deviations must be positive and samples at least two");
    }

    private static void requireCell(int value) {
        if (value <= 0) throw new IllegalArgumentException(
            "2-by-2 table cells must be positive for an uncorrected effect size");
    }
}
