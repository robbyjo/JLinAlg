/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.reml;

import java.util.Objects;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;

/** Immutable controls for Fisher-scoring REML optimization. */
public final class RemlOptions {
    private final double[] initialVariances;
    private final int maximumIterations;
    private final double relativeTolerance;
    private final double scoreTolerance;
    private final double minimumVariance;
    private final double maximumVariance;
    private final double maximumLogVarianceStep;
    private final DegreesOfFreedomMethod degreesOfFreedomMethod;
    private final VarianceEstimation varianceEstimation;

    private RemlOptions(Builder builder) {
        this.initialVariances = builder.initialVariances == null
            ? null : builder.initialVariances.clone();
        this.maximumIterations = builder.maximumIterations;
        this.relativeTolerance = builder.relativeTolerance;
        this.scoreTolerance = builder.scoreTolerance;
        this.minimumVariance = builder.minimumVariance;
        this.maximumVariance = builder.maximumVariance;
        this.maximumLogVarianceStep = builder.maximumLogVarianceStep;
        this.degreesOfFreedomMethod = builder.degreesOfFreedomMethod;
        this.varianceEstimation = builder.varianceEstimation;
    }

    /** Returns default optimization controls. */
    public static RemlOptions defaults() {
        return builder().build();
    }

    /** Returns a new options builder. */
    public static Builder builder() {
        return new Builder();
    }

    /** Returns a builder initialized from these controls. */
    public Builder toBuilder() { return new Builder(this); }

    public double[] initialVariances() {
        return initialVariances == null ? null : initialVariances.clone();
    }
    public int maximumIterations() { return maximumIterations; }
    public double relativeTolerance() { return relativeTolerance; }
    public double scoreTolerance() { return scoreTolerance; }
    public double minimumVariance() { return minimumVariance; }
    public double maximumVariance() { return maximumVariance; }
    public double maximumLogVarianceStep() { return maximumLogVarianceStep; }
    public DegreesOfFreedomMethod degreesOfFreedomMethod() {
        return degreesOfFreedomMethod;
    }
    public VarianceEstimation varianceEstimation() { return varianceEstimation; }

    /** Builder for {@link RemlOptions}. */
    public static final class Builder {
        private double[] initialVariances;
        private int maximumIterations = 100;
        private double relativeTolerance = 1e-9;
        private double scoreTolerance = 1e-6;
        private double minimumVariance = 1e-10;
        private double maximumVariance = 1e10;
        private double maximumLogVarianceStep = 1.5;
        private DegreesOfFreedomMethod degreesOfFreedomMethod =
            DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION;
        private VarianceEstimation varianceEstimation = VarianceEstimation.REML;

        private Builder() {
        }

        private Builder(RemlOptions options) {
            initialVariances = options.initialVariances();
            maximumIterations = options.maximumIterations();
            relativeTolerance = options.relativeTolerance();
            scoreTolerance = options.scoreTolerance();
            minimumVariance = options.minimumVariance();
            maximumVariance = options.maximumVariance();
            maximumLogVarianceStep = options.maximumLogVarianceStep();
            degreesOfFreedomMethod = options.degreesOfFreedomMethod();
            varianceEstimation = options.varianceEstimation();
        }

        public Builder initialVariances(double... values) {
            this.initialVariances = MatrixOps.finiteCopy(values, "initialVariances");
            return this;
        }

        public Builder maximumIterations(int value) {
            this.maximumIterations = value;
            return this;
        }

        public Builder relativeTolerance(double value) {
            this.relativeTolerance = value;
            return this;
        }

        public Builder scoreTolerance(double value) {
            this.scoreTolerance = value;
            return this;
        }

        public Builder varianceBounds(double minimum, double maximum) {
            this.minimumVariance = minimum;
            this.maximumVariance = maximum;
            return this;
        }

        public Builder maximumLogVarianceStep(double value) {
            this.maximumLogVarianceStep = value;
            return this;
        }

        /** Selects the denominator-DF calculation for fixed-effect t tests. */
        public Builder degreesOfFreedomMethod(DegreesOfFreedomMethod value) {
            if (value != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION
                    && value != DegreesOfFreedomMethod.SATTERTHWAITE
                    && value != DegreesOfFreedomMethod.KENWARD_ROGER) {
                throw new IllegalArgumentException(
                    "REML supports residual-approximation, Satterthwaite, "
                        + "or Kenward-Roger DF");
            }
            this.degreesOfFreedomMethod = Objects.requireNonNull(value, "value");
            return this;
        }

        /** Selects restricted or profile maximum likelihood. */
        public Builder varianceEstimation(VarianceEstimation value) {
            this.varianceEstimation = Objects.requireNonNull(value, "value");
            return this;
        }

        public RemlOptions build() {
            if (maximumIterations < 1) {
                throw new IllegalArgumentException("maximumIterations must be positive");
            }
            if (!(relativeTolerance > 0.0) || !Double.isFinite(relativeTolerance)) {
                throw new IllegalArgumentException("relativeTolerance must be finite and positive");
            }
            if (!(scoreTolerance > 0.0) || !Double.isFinite(scoreTolerance)) {
                throw new IllegalArgumentException("scoreTolerance must be finite and positive");
            }
            if (!(minimumVariance > 0.0)
                    || !(maximumVariance > minimumVariance)
                    || !Double.isFinite(maximumVariance)) {
                throw new IllegalArgumentException("variance bounds are invalid");
            }
            if (!(maximumLogVarianceStep > 0.0)
                    || !Double.isFinite(maximumLogVarianceStep)) {
                throw new IllegalArgumentException(
                    "maximumLogVarianceStep must be finite and positive");
            }
            if (initialVariances != null) {
                for (double value : initialVariances) {
                    if (!(value > 0.0)) {
                        throw new IllegalArgumentException(
                            "initial variances must be strictly positive");
                    }
                }
            }
            Objects.requireNonNull(degreesOfFreedomMethod, "degreesOfFreedomMethod");
            Objects.requireNonNull(varianceEstimation, "varianceEstimation");
            return new RemlOptions(this);
        }
    }
}
