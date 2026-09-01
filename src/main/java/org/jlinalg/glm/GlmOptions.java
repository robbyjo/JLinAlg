/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

import java.util.Objects;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.ols.RankDeficiencyStrategy;

/** Immutable IRLS and inference controls. */
public final class GlmOptions {
    private final int maximumIterations;
    private final double relativeTolerance;
    private final double confidenceLevel;
    private final double[] initialCoefficients;
    private final RankDeficiencyStrategy rankDeficiencyStrategy;
    private final DispersionMode dispersionMode;
    private final MissingDataPolicy missingDataPolicy;

    private GlmOptions(Builder builder) {
        maximumIterations = builder.maximumIterations;
        relativeTolerance = builder.relativeTolerance;
        confidenceLevel = builder.confidenceLevel;
        initialCoefficients = builder.initialCoefficients == null
            ? null : builder.initialCoefficients.clone();
        rankDeficiencyStrategy = builder.rankDeficiencyStrategy;
        dispersionMode = builder.dispersionMode;
        missingDataPolicy = builder.missingDataPolicy;
    }

    public static GlmOptions defaults() { return builder().build(); }
    public static Builder builder() { return new Builder(); }
    public int maximumIterations() { return maximumIterations; }
    public double relativeTolerance() { return relativeTolerance; }
    public double confidenceLevel() { return confidenceLevel; }
    public double[] initialCoefficients() {
        return initialCoefficients == null ? null : initialCoefficients.clone();
    }
    public RankDeficiencyStrategy rankDeficiencyStrategy() {
        return rankDeficiencyStrategy;
    }
    public DispersionMode dispersionMode() { return dispersionMode; }
    public MissingDataPolicy missingDataPolicy() { return missingDataPolicy; }

    /** Builder for GLM controls. */
    public static final class Builder {
        private int maximumIterations = 100;
        private double relativeTolerance = 1e-9;
        private double confidenceLevel = 0.95;
        private double[] initialCoefficients;
        private RankDeficiencyStrategy rankDeficiencyStrategy =
            RankDeficiencyStrategy.ERROR;
        private DispersionMode dispersionMode = DispersionMode.FAMILY_DEFAULT;
        private MissingDataPolicy missingDataPolicy = MissingDataPolicy.ERROR;

        private Builder() { }

        public Builder maximumIterations(int value) {
            maximumIterations = value;
            return this;
        }
        public Builder relativeTolerance(double value) {
            relativeTolerance = value;
            return this;
        }
        public Builder confidenceLevel(double value) {
            confidenceLevel = value;
            return this;
        }
        public Builder initialCoefficients(double... values) {
            initialCoefficients = MatrixOps.finiteCopy(values, "initialCoefficients");
            return this;
        }
        public Builder rankDeficiencyStrategy(RankDeficiencyStrategy value) {
            rankDeficiencyStrategy = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder dispersionMode(DispersionMode value) {
            dispersionMode = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder missingDataPolicy(MissingDataPolicy value) {
            missingDataPolicy = Objects.requireNonNull(value, "value");
            return this;
        }
        public GlmOptions build() {
            if (maximumIterations < 1) {
                throw new IllegalArgumentException("maximumIterations must be positive");
            }
            if (!(relativeTolerance > 0.0) || !Double.isFinite(relativeTolerance)) {
                throw new IllegalArgumentException(
                    "relativeTolerance must be finite and positive");
            }
            if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
                throw new IllegalArgumentException(
                    "confidenceLevel must be strictly between zero and one");
            }
            return new GlmOptions(this);
        }
    }
}
