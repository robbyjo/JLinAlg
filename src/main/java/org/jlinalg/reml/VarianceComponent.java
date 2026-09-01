/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.reml;

import java.util.List;
import java.util.Objects;
import org.jlinalg.internal.MatrixOps;

/** A named symmetric covariance basis matrix used in a REML model. */
public final class VarianceComponent {
    private final String name;
    private final int dimension;
    private final double[] covariance;

    /** Creates a component from a square row-major covariance basis matrix. */
    public VarianceComponent(String name, int dimension, double[] covariance) {
        this.name = Objects.requireNonNull(name, "name");
        if (name.isBlank()) {
            throw new IllegalArgumentException("component name must not be blank");
        }
        if (dimension < 1 || covariance == null
                || covariance.length != dimension * dimension) {
            throw new IllegalArgumentException("covariance dimensions are invalid");
        }
        this.dimension = dimension;
        this.covariance = MatrixOps.finiteCopy(covariance, "covariance");
        validateSymmetry(this.covariance, dimension);
    }

    /** Returns an identity component, normally used for residual variance. */
    public static VarianceComponent identity(String name, int dimension) {
        return new VarianceComponent(name, dimension, MatrixOps.identity(dimension));
    }

    /** Creates a diagonal covariance basis, for example inverse precision weights. */
    public static VarianceComponent diagonal(String name, double[] diagonal) {
        if (diagonal == null || diagonal.length == 0) {
            throw new IllegalArgumentException("diagonal is required");
        }
        double[] covariance = new double[diagonal.length * diagonal.length];
        for (int index = 0; index < diagonal.length; index++) {
            if (!(diagonal[index] > 0.0) || !Double.isFinite(diagonal[index])) {
                throw new IllegalArgumentException(
                    "covariance diagonal must be finite and positive");
            }
            covariance[index * diagonal.length + index] = diagonal[index];
        }
        return new VarianceComponent(name, diagonal.length, covariance);
    }

    /**
     * Creates an lme4-style random-intercept covariance basis. Observations
     * with equal, non-null group labels have covariance one.
     */
    public static VarianceComponent randomIntercept(
            String name, List<?> groups) {
        validateGroups(groups);
        int observations = groups.size();
        double[] covariance = new double[observations * observations];
        for (int row = 0; row < observations; row++) {
            for (int column = 0; column <= row; column++) {
                double value = groups.get(row).equals(groups.get(column))
                    ? 1.0 : 0.0;
                covariance[row * observations + column] = value;
                covariance[column * observations + row] = value;
            }
        }
        return new VarianceComponent(name, observations, covariance);
    }

    /**
     * Creates an uncorrelated random-slope covariance basis equivalent to
     * {@code (0 + covariate | group)} with a scalar slope variance.
     */
    public static VarianceComponent randomSlope(
            String name, List<?> groups, double[] covariate) {
        validateGroups(groups);
        if (covariate == null || covariate.length != groups.size()) {
            throw new IllegalArgumentException(
                "one random-slope covariate value is required per observation");
        }
        double[] values = MatrixOps.finiteCopy(covariate, "covariate");
        int observations = groups.size();
        double[] covariance = new double[observations * observations];
        for (int row = 0; row < observations; row++) {
            for (int column = 0; column <= row; column++) {
                double value = groups.get(row).equals(groups.get(column))
                    ? values[row] * values[column] : 0.0;
                covariance[row * observations + column] = value;
                covariance[column * observations + row] = value;
            }
        }
        return new VarianceComponent(name, observations, covariance);
    }

    /** Creates {@code K = Z Z'} for independent equal-variance random coefficients. */
    public static VarianceComponent fromRandomEffectDesign(
            String name, double[][] randomEffectDesign) {
        if (randomEffectDesign == null || randomEffectDesign.length == 0
                || randomEffectDesign[0] == null
                || randomEffectDesign[0].length == 0) {
            throw new IllegalArgumentException(
                "random-effect design must be a nonempty rectangular matrix");
        }
        int observations = randomEffectDesign.length;
        int coefficients = randomEffectDesign[0].length;
        double[] design = MatrixOps.rowMajor(randomEffectDesign, observations);
        double[] covariance = new double[observations * observations];
        for (int row = 0; row < observations; row++) {
            for (int column = 0; column <= row; column++) {
                double value = 0.0;
                for (int coefficient = 0; coefficient < coefficients; coefficient++) {
                    value += design[row * coefficients + coefficient]
                        * design[column * coefficients + coefficient];
                }
                covariance[row * observations + column] = value;
                covariance[column * observations + row] = value;
            }
        }
        return new VarianceComponent(name, observations, covariance);
    }

    public String name() { return name; }
    public int dimension() { return dimension; }
    public double[] covariance() { return covariance.clone(); }

    double[] covarianceView() {
        return covariance;
    }

    private static void validateGroups(List<?> groups) {
        if (groups == null || groups.isEmpty()) {
            throw new IllegalArgumentException("group labels are required");
        }
        for (Object group : groups) {
            if (group == null) {
                throw new IllegalArgumentException("group labels must not be null");
            }
        }
    }

    private static void validateSymmetry(double[] matrix, int dimension) {
        double maximum = 0.0;
        for (double value : matrix) {
            maximum = Math.max(maximum, Math.abs(value));
        }
        double tolerance = 1e-12 * Math.max(1.0, maximum);
        for (int row = 0; row < dimension; row++) {
            for (int column = row + 1; column < dimension; column++) {
                if (Math.abs(matrix[row * dimension + column]
                        - matrix[column * dimension + row]) > tolerance) {
                    throw new IllegalArgumentException(
                        "covariance component must be symmetric");
                }
            }
        }
    }
}
