/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Objects;

/** Immutable controls shared by meta-analysis and meta-regression. */
public final class MetaAnalysisOptions {
    private final MetaAnalysisMethod method;
    private final TauSquaredEstimator tauSquaredEstimator;
    private final MetaInferenceMethod inferenceMethod;
    private final double confidenceLevel;
    private final int maximumIterations;
    private final double tolerance;

    private MetaAnalysisOptions(Builder builder) {
        method = builder.method;
        tauSquaredEstimator = builder.tauSquaredEstimator;
        inferenceMethod = builder.inferenceMethod;
        confidenceLevel = builder.confidenceLevel;
        maximumIterations = builder.maximumIterations;
        tolerance = builder.tolerance;
    }

    public static MetaAnalysisOptions fixedEffect() {
        return builder().method(MetaAnalysisMethod.FIXED_EFFECT).build();
    }

    public static MetaAnalysisOptions randomEffects() { return builder().build(); }
    public static Builder builder() { return new Builder(); }

    public MetaAnalysisMethod method() { return method; }
    public TauSquaredEstimator tauSquaredEstimator() { return tauSquaredEstimator; }
    public MetaInferenceMethod inferenceMethod() { return inferenceMethod; }
    public double confidenceLevel() { return confidenceLevel; }
    public int maximumIterations() { return maximumIterations; }
    public double tolerance() { return tolerance; }

    /** Builder with random-effects REML and Wald-z defaults. */
    public static final class Builder {
        private MetaAnalysisMethod method = MetaAnalysisMethod.RANDOM_EFFECT;
        private TauSquaredEstimator tauSquaredEstimator = TauSquaredEstimator.REML;
        private MetaInferenceMethod inferenceMethod = MetaInferenceMethod.NORMAL;
        private double confidenceLevel = 0.95;
        private int maximumIterations = 200;
        private double tolerance = 1e-10;

        private Builder() { }
        public Builder method(MetaAnalysisMethod value) {
            method = Objects.requireNonNull(value, "value"); return this;
        }
        public Builder tauSquaredEstimator(TauSquaredEstimator value) {
            tauSquaredEstimator = Objects.requireNonNull(value, "value"); return this;
        }
        public Builder inferenceMethod(MetaInferenceMethod value) {
            inferenceMethod = Objects.requireNonNull(value, "value"); return this;
        }
        public Builder confidenceLevel(double value) {
            confidenceLevel = value; return this;
        }
        public Builder maximumIterations(int value) {
            maximumIterations = value; return this;
        }
        public Builder tolerance(double value) { tolerance = value; return this; }

        public MetaAnalysisOptions build() {
            Objects.requireNonNull(method, "method");
            Objects.requireNonNull(tauSquaredEstimator, "tauSquaredEstimator");
            Objects.requireNonNull(inferenceMethod, "inferenceMethod");
            if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0))
                throw new IllegalArgumentException("confidence level must be in (0,1)");
            if (maximumIterations < 10)
                throw new IllegalArgumentException("maximum iterations must be at least 10");
            if (!(tolerance > 0.0) || !Double.isFinite(tolerance))
                throw new IllegalArgumentException("tolerance must be finite and positive");
            return new MetaAnalysisOptions(this);
        }
    }
}
