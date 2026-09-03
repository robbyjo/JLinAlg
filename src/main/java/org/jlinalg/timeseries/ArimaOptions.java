/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import java.util.Objects;

/** Immutable controls for conditional Gaussian ARIMA fitting. */
public final class ArimaOptions {
    private final SeasonalArimaOrder seasonalOrder;
    private final boolean includeMean;
    private final boolean includeDrift;
    private final int optimizationStarts;
    private final int maximumFunctionEvaluations;
    private final double optimizationTolerance;

    private ArimaOptions(Builder builder) {
        seasonalOrder = builder.seasonalOrder;
        includeMean = builder.includeMean;
        includeDrift = builder.includeDrift;
        optimizationStarts = builder.optimizationStarts;
        maximumFunctionEvaluations = builder.maximumFunctionEvaluations;
        optimizationTolerance = builder.optimizationTolerance;
    }

    public static ArimaOptions defaults() { return builder().build(); }
    public static Builder builder() { return new Builder(); }

    public SeasonalArimaOrder seasonalOrder() { return seasonalOrder; }
    public boolean includeMean() { return includeMean; }
    public boolean includeDrift() { return includeDrift; }
    public int optimizationStarts() { return optimizationStarts; }
    public int maximumFunctionEvaluations() { return maximumFunctionEvaluations; }
    public double optimizationTolerance() { return optimizationTolerance; }

    /** Builder for {@link ArimaOptions}. */
    public static final class Builder {
        private SeasonalArimaOrder seasonalOrder = SeasonalArimaOrder.none();
        private boolean includeMean = true;
        private boolean includeDrift;
        private int optimizationStarts = 1;
        private int maximumFunctionEvaluations = 5_000;
        private double optimizationTolerance = 1e-8;

        private Builder() { }

        public Builder seasonalOrder(SeasonalArimaOrder value) {
            seasonalOrder = Objects.requireNonNull(value, "value");
            return this;
        }

        public Builder includeMean(boolean value) {
            includeMean = value;
            return this;
        }

        public Builder includeDrift(boolean value) {
            includeDrift = value;
            return this;
        }

        /** Number of deterministic optimizer starts; one is fastest. */
        public Builder optimizationStarts(int value) {
            optimizationStarts = value;
            return this;
        }

        public Builder maximumFunctionEvaluations(int value) {
            maximumFunctionEvaluations = value;
            return this;
        }

        public Builder optimizationTolerance(double value) {
            optimizationTolerance = value;
            return this;
        }

        public ArimaOptions build() {
            if (optimizationStarts < 1 || optimizationStarts > 20) {
                throw new IllegalArgumentException(
                    "optimizationStarts must be between 1 and 20");
            }
            if (maximumFunctionEvaluations < 20) {
                throw new IllegalArgumentException(
                    "maximumFunctionEvaluations must be at least 20");
            }
            if (!(optimizationTolerance > 0.0)
                    || !Double.isFinite(optimizationTolerance)) {
                throw new IllegalArgumentException(
                    "optimizationTolerance must be finite and positive");
            }
            return new ArimaOptions(this);
        }
    }
}
