/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import org.jlinalg.internal.MatrixOps;

/** Immutable preprocessing and convergence controls for elastic-net fitting. */
public final class ElasticNetOptions {
    private final double alpha;
    private final boolean fitIntercept;
    private final boolean standardize;
    private final int maximumIterations;
    private final double relativeTolerance;
    private final double[] observationWeights;
    private final double[] penaltyFactors;

    private ElasticNetOptions(Builder builder) {
        alpha = builder.alpha;
        fitIntercept = builder.fitIntercept;
        standardize = builder.standardize;
        maximumIterations = builder.maximumIterations;
        relativeTolerance = builder.relativeTolerance;
        observationWeights = builder.observationWeights == null
            ? null : builder.observationWeights.clone();
        penaltyFactors = builder.penaltyFactors == null
            ? null : builder.penaltyFactors.clone();
    }

    public static ElasticNetOptions defaults() { return builder().build(); }
    public static Builder builder() { return new Builder(); }

    /** Elastic-net mixing parameter: zero is ridge and one is LASSO. */
    public double alpha() { return alpha; }
    public boolean fitIntercept() { return fitIntercept; }
    public boolean standardize() { return standardize; }
    public int maximumIterations() { return maximumIterations; }
    public double relativeTolerance() { return relativeTolerance; }
    public double[] observationWeights() {
        return observationWeights == null ? null : observationWeights.clone();
    }
    public double[] penaltyFactors() {
        return penaltyFactors == null ? null : penaltyFactors.clone();
    }

    /** Builder for {@link ElasticNetOptions}. */
    public static final class Builder {
        private double alpha = 0.5;
        private boolean fitIntercept = true;
        private boolean standardize = true;
        private int maximumIterations = 100_000;
        private double relativeTolerance = 1e-8;
        private double[] observationWeights;
        private double[] penaltyFactors;

        private Builder() { }

        public Builder alpha(double value) {
            alpha = value;
            return this;
        }

        public Builder fitIntercept(boolean value) {
            fitIntercept = value;
            return this;
        }

        public Builder standardize(boolean value) {
            standardize = value;
            return this;
        }

        public Builder maximumIterations(int value) {
            maximumIterations = value;
            return this;
        }

        public Builder relativeTolerance(double value) {
            relativeTolerance = value;
            return this;
        }

        /** Positive observation weights; fitting rescales them to sum to N. */
        public Builder observationWeights(double... values) {
            observationWeights = MatrixOps.finiteCopy(values, "observationWeights");
            return this;
        }

        /** Positive per-coefficient multipliers for both L1 and L2 penalties. */
        public Builder penaltyFactors(double... values) {
            penaltyFactors = MatrixOps.finiteCopy(values, "penaltyFactors");
            return this;
        }

        public ElasticNetOptions build() {
            if (!Double.isFinite(alpha) || alpha < 0.0 || alpha > 1.0) {
                throw new IllegalArgumentException("alpha must lie in [0, 1]");
            }
            if (maximumIterations < 1) {
                throw new IllegalArgumentException(
                    "maximumIterations must be positive");
            }
            if (!Double.isFinite(relativeTolerance)
                    || !(relativeTolerance > 0.0)) {
                throw new IllegalArgumentException(
                    "relativeTolerance must be finite and positive");
            }
            validatePositive(observationWeights, "observationWeights");
            validatePositive(penaltyFactors, "penaltyFactors");
            return new ElasticNetOptions(this);
        }

        private static void validatePositive(double[] values, String name) {
            if (values == null) {
                return;
            }
            for (double value : values) {
                if (!(value > 0.0)) {
                    throw new IllegalArgumentException(
                        name + " must contain only positive values");
                }
            }
        }
    }
}
