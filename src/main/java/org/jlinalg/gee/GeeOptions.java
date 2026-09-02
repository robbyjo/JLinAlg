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
    private final double associationTolerance;
    private final double scaleTolerance;
    private final double scoreTolerance;
    private final double confidenceLevel;
    private final GeeCorrelation correlation;
    private final GeeCovariance covariance;
    private final GeeMethod method;
    private final GeeAssociation association;
    private final GeeInference inference;
    private final GeeParameterLink associationLink;
    private final GeeParameterLink scaleLink;
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
    private final double associationDamping;
    private final double fayGraubardBound;
    private final int parallelism;
    private final int parallelThreshold;
    private final boolean exactClusterDeletion;
    private final boolean positiveDefiniteProjection;
    private final MissingDataPolicy missingDataPolicy;

    private GeeOptions(Builder builder) {
        maximumIterations = builder.maximumIterations;
        relativeTolerance = builder.relativeTolerance;
        associationTolerance = Double.isNaN(builder.associationTolerance)
            ? relativeTolerance : builder.associationTolerance;
        scaleTolerance = Double.isNaN(builder.scaleTolerance)
            ? relativeTolerance : builder.scaleTolerance;
        scoreTolerance = Double.isNaN(builder.scoreTolerance)
            ? relativeTolerance : builder.scoreTolerance;
        confidenceLevel = builder.confidenceLevel;
        correlation = builder.correlation;
        covariance = builder.covariance;
        method = builder.method;
        association = builder.association;
        inference = builder.inference;
        associationLink = builder.associationLink;
        scaleLink = builder.scaleLink;
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
        associationDamping = builder.associationDamping;
        fayGraubardBound = builder.fayGraubardBound;
        parallelism = builder.parallelism;
        parallelThreshold = builder.parallelThreshold;
        exactClusterDeletion = builder.exactClusterDeletion;
        positiveDefiniteProjection = builder.positiveDefiniteProjection;
        missingDataPolicy = builder.missingDataPolicy;
    }

    public static GeeOptions defaults() { return builder().build(); }
    public static Builder builder() { return new Builder(); }
    public static Builder builder(GeeOptions source) { return new Builder(source); }
    public Builder toBuilder() { return builder(this); }
    public int maximumIterations() { return maximumIterations; }
    public double relativeTolerance() { return relativeTolerance; }
    public double associationTolerance() { return associationTolerance; }
    public double scaleTolerance() { return scaleTolerance; }
    public double scoreTolerance() { return scoreTolerance; }
    public double confidenceLevel() { return confidenceLevel; }
    public GeeCorrelation correlation() { return correlation; }
    public GeeCovariance covariance() { return covariance; }
    public GeeMethod method() { return method; }
    public GeeAssociation association() { return association; }
    public GeeInference inference() { return inference; }
    public GeeParameterLink associationLink() { return associationLink; }
    public GeeParameterLink scaleLink() { return scaleLink; }
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
    public double associationDamping() { return associationDamping; }
    public double fayGraubardBound() { return fayGraubardBound; }
    public int parallelism() { return parallelism; }
    public int parallelThreshold() { return parallelThreshold; }
    public boolean exactClusterDeletion() { return exactClusterDeletion; }
    public boolean positiveDefiniteProjection() { return positiveDefiniteProjection; }
    public MissingDataPolicy missingDataPolicy() { return missingDataPolicy; }

    private static double[] copy(double[] values) {
        return values == null ? null : values.clone();
    }

    /** Builder for GEE controls. */
    public static final class Builder {
        private int maximumIterations = 100;
        private double relativeTolerance = 1e-8;
        private double associationTolerance = Double.NaN;
        private double scaleTolerance = Double.NaN;
        private double scoreTolerance = Double.NaN;
        private double confidenceLevel = 0.95;
        private GeeCorrelation correlation = GeeCorrelation.INDEPENDENCE;
        private GeeCovariance covariance = GeeCovariance.ROBUST;
        private GeeMethod method = GeeMethod.ORDINARY;
        private GeeAssociation association = GeeAssociation.CORRELATION;
        private GeeInference inference = GeeInference.ASYMPTOTIC;
        private GeeParameterLink associationLink = GeeParameterLink.IDENTITY;
        private GeeParameterLink scaleLink = GeeParameterLink.LOG;
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
        private double associationDamping = 0.5;
        private double fayGraubardBound = 0.75;
        private int parallelism = 1;
        private int parallelThreshold = 128;
        private boolean exactClusterDeletion;
        private boolean positiveDefiniteProjection = true;
        private MissingDataPolicy missingDataPolicy = MissingDataPolicy.ERROR;

        private Builder() { }

        private Builder(GeeOptions source) {
            maximumIterations = source.maximumIterations;
            relativeTolerance = source.relativeTolerance;
            associationTolerance = source.associationTolerance;
            scaleTolerance = source.scaleTolerance;
            scoreTolerance = source.scoreTolerance;
            confidenceLevel = source.confidenceLevel;
            correlation = source.correlation;
            covariance = source.covariance;
            method = source.method;
            association = source.association;
            inference = source.inference;
            associationLink = source.associationLink;
            scaleLink = source.scaleLink;
            dependenceOrder = source.dependenceOrder;
            initialCoefficients = copy(source.initialCoefficients);
            fixedAssociation = copy(source.fixedAssociation);
            fixedAssociationDimension = source.fixedAssociationDimension;
            correlationDesign = copy(source.correlationDesign);
            correlationDesignRows = source.correlationDesignRows;
            correlationDesignColumns = source.correlationDesignColumns;
            scaleDesign = copy(source.scaleDesign);
            scaleDesignRows = source.scaleDesignRows;
            scaleDesignColumns = source.scaleDesignColumns;
            fixedDispersion = source.fixedDispersion;
            dispersion = source.dispersion;
            jeffreysPower = source.jeffreysPower;
            oddsRatioContinuityCorrection = source.oddsRatioContinuityCorrection;
            associationDamping = source.associationDamping;
            fayGraubardBound = source.fayGraubardBound;
            parallelism = source.parallelism;
            parallelThreshold = source.parallelThreshold;
            exactClusterDeletion = source.exactClusterDeletion;
            positiveDefiniteProjection = source.positiveDefiniteProjection;
            missingDataPolicy = source.missingDataPolicy;
        }

        public Builder maximumIterations(int value) {
            maximumIterations = value;
            return this;
        }
        public Builder relativeTolerance(double value) {
            relativeTolerance = value;
            return this;
        }
        public Builder associationTolerance(double value) {
            associationTolerance = value;
            return this;
        }
        public Builder scaleTolerance(double value) {
            scaleTolerance = value;
            return this;
        }
        public Builder scoreTolerance(double value) {
            scoreTolerance = value;
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
        public Builder inference(GeeInference value) {
            inference = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder associationLink(GeeParameterLink value) {
            associationLink = Objects.requireNonNull(value, "value");
            return this;
        }
        public Builder scaleLink(GeeParameterLink value) {
            scaleLink = Objects.requireNonNull(value, "value");
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
        public Builder associationDamping(double value) {
            associationDamping = value;
            return this;
        }
        public Builder fayGraubardBound(double value) {
            fayGraubardBound = value;
            return this;
        }
        public Builder parallelism(int value) {
            parallelism = value;
            return this;
        }
        public Builder parallelThreshold(int value) {
            parallelThreshold = value;
            return this;
        }
        public Builder exactClusterDeletion(boolean value) {
            exactClusterDeletion = value;
            return this;
        }
        public Builder positiveDefiniteProjection(boolean value) {
            positiveDefiniteProjection = value;
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
            validateOptionalTolerance(associationTolerance, "associationTolerance");
            validateOptionalTolerance(scaleTolerance, "scaleTolerance");
            validateOptionalTolerance(scoreTolerance, "scoreTolerance");
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
            if (!(associationDamping > 0.0 && associationDamping <= 1.0)) {
                throw new IllegalArgumentException(
                    "associationDamping must be in (0, 1]");
            }
            if (!(fayGraubardBound > 0.0 && fayGraubardBound < 1.0)) {
                throw new IllegalArgumentException(
                    "fayGraubardBound must be in (0, 1)");
            }
            if (parallelism < 1 || parallelThreshold < 1) {
                throw new IllegalArgumentException(
                    "parallelism and parallelThreshold must be positive");
            }
            if (scaleLink == GeeParameterLink.FISHER_Z) {
                throw new IllegalArgumentException(
                    "Fisher-z is not a valid positive scale link");
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

        private static void validateOptionalTolerance(double value, String name) {
            if (!Double.isNaN(value)
                    && (!(value > 0.0) || !Double.isFinite(value))) {
                throw new IllegalArgumentException(name + " must be finite and positive");
            }
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
