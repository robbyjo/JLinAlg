/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import java.util.Objects;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.reml.RemlOptions;

/** Immutable controls for outer PQL iterations and inner model fits. */
public final class GlmmPqlOptions {
    private final int maximumIterations;
    private final double relativeTolerance;
    private final GlmOptions initialGlmOptions;
    private final RemlOptions remlOptions;

    private GlmmPqlOptions(Builder builder) {
        maximumIterations = builder.maximumIterations;
        relativeTolerance = builder.relativeTolerance;
        initialGlmOptions = builder.initialGlmOptions;
        remlOptions = builder.remlOptions;
    }

    public static GlmmPqlOptions defaults() { return builder().build(); }
    public static Builder builder() { return new Builder(); }
    public int maximumIterations() { return maximumIterations; }
    public double relativeTolerance() { return relativeTolerance; }
    public GlmOptions initialGlmOptions() { return initialGlmOptions; }
    public RemlOptions remlOptions() { return remlOptions; }

    /** Builder for PQL controls. */
    public static final class Builder {
        private int maximumIterations = 30;
        private double relativeTolerance = 1e-6;
        private GlmOptions initialGlmOptions = GlmOptions.defaults();
        private RemlOptions remlOptions = RemlOptions.defaults();

        private Builder() { }

        public Builder maximumIterations(int value) {
            maximumIterations = value;
            return this;
        }
        public Builder relativeTolerance(double value) {
            relativeTolerance = value;
            return this;
        }
        public Builder initialGlmOptions(GlmOptions value) {
            initialGlmOptions = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder remlOptions(RemlOptions value) {
            remlOptions = Objects.requireNonNull(value, "value");
            return this;
        }
        public GlmmPqlOptions build() {
            if (maximumIterations < 1) {
                throw new IllegalArgumentException("maximumIterations must be positive");
            }
            if (!(relativeTolerance > 0.0) || !Double.isFinite(relativeTolerance)) {
                throw new IllegalArgumentException(
                    "relativeTolerance must be finite and positive");
            }
            return new GlmmPqlOptions(this);
        }
    }
}
