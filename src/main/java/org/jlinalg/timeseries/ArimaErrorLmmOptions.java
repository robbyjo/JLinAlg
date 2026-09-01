/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import java.util.Objects;
import org.jlinalg.reml.RemlOptions;

/** Outer ARMA-correlation and inner REML controls for time-series LMMs. */
public final class ArimaErrorLmmOptions {
    private final SeasonalArimaOrder seasonalOrder;
    private final RemlOptions remlOptions;
    private final int maximumFunctionEvaluations;
    private final double optimizationTolerance;

    private ArimaErrorLmmOptions(Builder builder) {
        seasonalOrder = builder.seasonalOrder;
        remlOptions = builder.remlOptions;
        maximumFunctionEvaluations = builder.maximumFunctionEvaluations;
        optimizationTolerance = builder.optimizationTolerance;
    }

    public static ArimaErrorLmmOptions defaults() { return builder().build(); }
    public static Builder builder() { return new Builder(); }

    public SeasonalArimaOrder seasonalOrder() { return seasonalOrder; }
    public RemlOptions remlOptions() { return remlOptions; }
    public int maximumFunctionEvaluations() { return maximumFunctionEvaluations; }
    public double optimizationTolerance() { return optimizationTolerance; }

    /** Builder for {@link ArimaErrorLmmOptions}. */
    public static final class Builder {
        private SeasonalArimaOrder seasonalOrder = SeasonalArimaOrder.none();
        private RemlOptions remlOptions = RemlOptions.defaults();
        private int maximumFunctionEvaluations = 250;
        private double optimizationTolerance = 1e-4;

        private Builder() { }

        public Builder seasonalOrder(SeasonalArimaOrder value) {
            seasonalOrder = Objects.requireNonNull(value, "value");
            return this;
        }

        public Builder remlOptions(RemlOptions value) {
            remlOptions = Objects.requireNonNull(value, "value");
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

        public ArimaErrorLmmOptions build() {
            if (maximumFunctionEvaluations < 20) {
                throw new IllegalArgumentException(
                    "maximumFunctionEvaluations must be at least 20");
            }
            if (!(optimizationTolerance > 0.0)
                    || !Double.isFinite(optimizationTolerance)) {
                throw new IllegalArgumentException(
                    "optimizationTolerance must be finite and positive");
            }
            return new ArimaErrorLmmOptions(this);
        }
    }
}
