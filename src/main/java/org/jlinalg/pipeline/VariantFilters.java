/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.ArrayList;
import java.util.List;

/** Cohort-aligned variant filtering with auditable exclusion reasons. */
public final class VariantFilters {
    private VariantFilters() { }

    public static VariantFilterResult evaluate(
            VariantRecord variant, VariantFilterOptions options) {
        if (variant == null || options == null)
            throw new IllegalArgumentException("variant and options are required");
        VariantStatistics statistics = VariantStatistics.of(variant);
        List<VariantFilterReason> reasons = new ArrayList<>();
        if (statistics.calledSamples() == 0) {
            reasons.add(VariantFilterReason.NO_CALLED_SAMPLES);
            return new VariantFilterResult(variant, statistics, reasons);
        }
        if (statistics.missingRate() > options.maximumMissingRate())
            reasons.add(VariantFilterReason.TOO_MANY_MISSING);
        if (statistics.minorAlleleFrequency() < options.minimumMaf())
            reasons.add(VariantFilterReason.BELOW_MINIMUM_MAF);
        if (statistics.minorAlleleFrequency() > options.maximumMaf())
            reasons.add(VariantFilterReason.ABOVE_MAXIMUM_MAF);
        if (statistics.minorAlleleCount() < options.minimumMac())
            reasons.add(VariantFilterReason.BELOW_MINIMUM_MAC);
        if (statistics.minorAlleleCount() > options.maximumMac())
            reasons.add(VariantFilterReason.ABOVE_MAXIMUM_MAC);
        if (!Double.isNaN(options.minimumImputationQuality())
                && (Double.isNaN(variant.imputationQuality())
                    || variant.imputationQuality()
                        < options.minimumImputationQuality()))
            reasons.add(VariantFilterReason.BELOW_IMPUTATION_QUALITY);
        if (options.excludeMonomorphic()
                && !(statistics.dosageVariance() > 1e-14))
            reasons.add(VariantFilterReason.MONOMORPHIC);
        return new VariantFilterResult(variant, statistics, reasons);
    }
}
