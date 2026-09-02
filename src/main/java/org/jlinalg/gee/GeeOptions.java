/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.Objects;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.model.MissingDataPolicy;

/** Immutable fitting, association, covariance, and adjustment controls. */
public final class GeeOptions {
    private final int maximumIterations;
    private final double relativeTolerance;
    private final double confidenceLevel;
    private final GeeCorrelation correlation;
    private final GeeCovariance covariance;
    private final GeeMethod method;
    private final GeeAssociation association;
    private final int dependenceOrder;
    private final double[] initialCoefficients;
    private final double[] fixedAssociation;
    private final int fixedAssociationDimension;
    private final double[] correlationDesign;
    private final int correlationDesignRows;
    private final int correlationDesignColumns;
    private final double[] scaleDesign;
    private final int scaleDesignRows;
    private final int scaleDesignColumns;
    private final boolean fixedDispersion;
    private final double dispersion;
    private final double jeffreysPower;
    private final double oddsRatioContinuityCorrection;
    private final MissingDataPolicy missingDataPolicy;

    private GeeOptions(Builder builder) {
        maximumIterations = builder.maximumIterations;
        relativeTolerance = builder.relativeTolerance;
        confidenceLevel = builder.confidenceLevel;
        correlation = builder.correlation;
        covariance = builder.covariance;
        method = builder.method;
        association = builder.association;
        dependenceOrder = builder.dependenceOrder;
        initialCoefficients = copy(builder.initialCoefficients);
        fixedAssociation = copy(builder.fixedAssociation);
        fixedAssociationDimension = builder.fixedAssociationDimension;
        correlationDesign = copy(builder.correlationDesign);
        correlationDesignRows = builder.correlationDesignRows;
        correlationDesignColumns = builder.correlationDesignColumns;
        scaleDesign = copy(builder.scaleDesign);
        scaleDesignRows = builder.scaleDesignRows;
        scaleDesignColumns = builder.scaleDesignColumns;
        fixedDispersion = builder.fixedDispersion;
        dispersion = builder.dispersion;
        jeffreysPower = builder.jeffreysPower;
        oddsRatioContinuityCorrection = builder.oddsRatioContinuityCorrection;
        missingDataPolicy = builder.missingDataPolicy;
    }

    public static GeeOptions defaults() { return builder().build(); }
    public static Builder builder() { return new Builder(); }
    public int maximumIterations() { return maximumIterations; }
    public double relativeTolerance() { return relativeTolerance; }
    public double confidenceLevel() { return confidenceLevel; }
    public GeeCorrelation correlation() { return correlation; }
    public GeeCovariance covariance() { return covariance; }
    public GeeMethod method() { return method; }
    public GeeAssociation association() { return association; }
    public int dependenceOrder() { return dependenceOrder; }
    public double[] initialCoefficients() { return copy(initialCoefficients); }
    public double[] fixedAssociation() { return copy(fixedAssociation); }
    public int fixedAssociationDimension() { return fixedAssociationDimension; }
    public double[] correlationDesign() { return copy(correlationDesign); }
    public int correlationDesignRows() { return correlationDesignRows; }
    public int correlationDesignColumns() { return correlationDesignColumns; }
    public double[] scaleDesign() { return copy(scaleDesign); }
    public int scaleDesignRows() { return scaleDesignRows; }
    public int scaleDesignColumns() { return scaleDesignColumns; }
    public boolean fixedDispersion() { return fixedDispersion; }
    public double dispersion() { return dispersion; }
    public double jeffreysPower() { return jeffreysPower; }
    public double oddsRatioContinuityCorrection() {
        return oddsRatioContinuityCorrection;
    }
    public MissingDataPolicy missingDataPolicy() { return missingDataPolicy; }

    private static double[] copy(double[] values) {
        return values == null ? null : values.clone();
    }

    /** Builder for GEE controls. */
    public static final class Builder {
        private int maximumIterations = 100;
        private double relativeTolerance = 1e-8;
        private double confidenceLevel = 0.95;
        private GeeCorrelation correlation = GeeCorrelation.INDEPENDENCE;
        private GeeCovariance covariance = GeeCovariance.ROBUST;
        private GeeMethod method = GeeMethod.ORDINARY;
        private GeeAssociation association = GeeAssociation.CORRELATION;
        private int dependenceOrder = 1;
        private double[] initialCoefficients;
        private double[] fixedAssociation;
        private int fixedAssociationDimension;
        private double[] correlationDesign;
        private int correlationDesignRows;
        private int correlationDesignColumns;
        private double[] scaleDesign;
        private int scaleDesignRows;
        private int scaleDesignColumns;
        private boolean fixedDispersion;
        private double dispersion = 1.0;
        private double jeffreysPower = 0.5;
        private double oddsRatioContinuityCorrection = 0.5;
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
        public Builder correlation(GeeCorrelation value) {
            correlation = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder covariance(GeeCovariance value) {
            covariance = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder method(GeeMethod value) {
            method = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder association(GeeAssociation value) {
            association = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder dependenceOrder(int value) {
            dependenceOrder = value;
            return this;
        }
        public Builder initialCoefficients(double... values) {
            initialCoefficients = MatrixOps.finiteCopy(values, "initialCoefficients");
            return this;
        }
        /** Supplies a symmetric correlation or odds-ratio matrix. */
        public Builder fixedAssociation(double[][] values) {
            fixedAssociation = square(values, "fixedAssociation");
            fixedAssociationDimension = values.length;
            return this;
        }
        /**
         * Supplies one row per lower-triangular wave pair and one column per
         * association parameter, following row order (1,0), (2,0), (2,1), ... .
         */
        public Builder correlationDesign(double[][] values) {
            correlationDesignRows = rows(values, "correlationDesign");
            correlationDesignColumns = values[0].length;
            correlationDesign = MatrixOps.rowMajor(values, values.length);
            return this;
        }
        /** Supplies an observation-level log-dispersion design matrix. */
        public Builder scaleDesign(double[][] values) {
            scaleDesignRows = rows(values, "scaleDesign");
            scaleDesignColumns = values[0].length;
            scaleDesign = MatrixOps.rowMajor(values, values.length);
            return this;
        }
        public Builder fixedDispersion(double value) {
            fixedDispersion = true;
            dispersion = value;
            return this;
        }
        public Builder estimateDispersion() {
            fixedDispersion = false;
            return this;
        }
        public Builder jeffreysPower(double value) {
            jeffreysPower = value;
            return this;
        }
        public Builder oddsRatioContinuityCorrection(double value) {
            oddsRatioContinuityCorrection = value;
            return this;
        }
        public Builder missingDataPolicy(MissingDataPolicy value) {
            missingDataPolicy = Objects.requireNonNull(value, "value");
            return this;
        }

        public GeeOptions build() {
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
            if (dependenceOrder < 1) {
                throw new IllegalArgumentException("dependenceOrder must be positive");
            }
            if (!(dispersion > 0.0) || !Double.isFinite(dispersion)) {
                throw new IllegalArgumentException("dispersion must be finite and positive");
            }
            if (!(jeffreysPower > 0.0) || !Double.isFinite(jeffreysPower)) {
                throw new IllegalArgumentException(
                    "jeffreysPower must be finite and positive");
            }
            if (!(oddsRatioContinuityCorrection > 0.0)
                    || !Double.isFinite(oddsRatioContinuityCorrection)) {
                throw new IllegalArgumentException(
                    "odds-ratio continuity correction must be finite and positive");
            }
            if (correlation == GeeCorrelation.FIXED && fixedAssociation == null) {
                throw new IllegalArgumentException(
                    "fixed correlation requires fixedAssociation");
            }
            if (correlation == GeeCorrelation.USER_DEFINED
                    && correlationDesign == null) {
                throw new IllegalArgumentException(
                    "user-defined correlation requires correlationDesign");
            }
            return new GeeOptions(this);
        }

        private static int rows(double[][] values, String name) {
            if (values == null || values.length == 0 || values[0] == null
                    || values[0].length == 0) {
                throw new IllegalArgumentException(name + " must not be empty");
            }
            return values.length;
        }

        private static double[] square(double[][] values, String name) {
            int dimension = rows(values, name);
            for (double[] row : values) {
                if (row == null || row.length != dimension) {
                    throw new IllegalArgumentException(name + " must be square");
                }
            }
            return MatrixOps.rowMajor(values, dimension);
        }
    }
}
