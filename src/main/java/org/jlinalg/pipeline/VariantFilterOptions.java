/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Explicit aligned-cohort MAF, MAC, missingness, and quality thresholds. */
public record VariantFilterOptions(
        double minimumMaf,
        double maximumMaf,
        double minimumMac,
        double maximumMac,
        double maximumMissingRate,
        double minimumImputationQuality,
        boolean excludeMonomorphic) {
    public VariantFilterOptions {
        if (!Double.isFinite(minimumMaf) || !Double.isFinite(maximumMaf)
                || minimumMaf < 0 || maximumMaf > 0.5
                || minimumMaf > maximumMaf)
            throw new IllegalArgumentException(
                "MAF bounds must satisfy 0 <= minimum <= maximum <= 0.5");
        if (!Double.isFinite(minimumMac) || minimumMac < 0
                || Double.isNaN(maximumMac) || maximumMac < minimumMac)
            throw new IllegalArgumentException(
                "MAC bounds must satisfy 0 <= minimum <= maximum");
        if (!Double.isFinite(maximumMissingRate)
                || maximumMissingRate < 0 || maximumMissingRate > 1)
            throw new IllegalArgumentException(
                "maximum missing rate must be in [0,1]");
        if (!Double.isNaN(minimumImputationQuality)
                && (!Double.isFinite(minimumImputationQuality)
                    || minimumImputationQuality < 0
                    || minimumImputationQuality > 1))
            throw new IllegalArgumentException(
                "minimum imputation quality must be NaN or in [0,1]");
    }

    /** No implicit rare-variant cutoff; only monomorphic rows are excluded. */
    public static VariantFilterOptions defaults() {
        return new VariantFilterOptions(0, 0.5, 0,
            Double.POSITIVE_INFINITY, 1, Double.NaN, true);
    }

    public static Builder builder() { return new Builder(); }

    public static final class Builder {
        private double minimumMaf;
        private double maximumMaf = 0.5;
        private double minimumMac;
        private double maximumMac = Double.POSITIVE_INFINITY;
        private double maximumMissingRate = 1;
        private double minimumImputationQuality = Double.NaN;
        private boolean excludeMonomorphic = true;

        private Builder() { }
        public Builder minimumMaf(double value) { minimumMaf = value; return this; }
        public Builder maximumMaf(double value) { maximumMaf = value; return this; }
        public Builder minimumMac(double value) { minimumMac = value; return this; }
        public Builder maximumMac(double value) { maximumMac = value; return this; }
        public Builder maximumMissingRate(double value) {
            maximumMissingRate = value; return this;
        }
        public Builder minimumImputationQuality(double value) {
            minimumImputationQuality = value; return this;
        }
        public Builder excludeMonomorphic(boolean value) {
            excludeMonomorphic = value; return this;
        }
        public VariantFilterOptions build() {
            return new VariantFilterOptions(minimumMaf, maximumMaf,
                minimumMac, maximumMac, maximumMissingRate,
                minimumImputationQuality, excludeMonomorphic);
        }
    }
}
