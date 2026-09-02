/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.Objects;
import org.jlinalg.internal.MatrixOps;

/**
 * A univariate penalized B-spline smooth with an equally spaced knot sequence.
 *
 * <p>The coefficient penalty is the squared norm of an integer-order
 * difference operator. The basis and penalty factor are retained separately so
 * a fitter can use penalized least squares or the equivalent mixed model.</p>
 */
public final class PSplineTerm {
    private final String name;
    private final double[] covariate;
    private final int basisDimension;
    private final int degree;
    private final int differenceOrder;
    private final double lowerBoundary;
    private final double upperBoundary;
    private final double[] knots;
    private final double[] design;
    private final double[] penaltyFactor;

    private PSplineTerm(
            String name,
            double[] covariate,
            int basisDimension,
            int degree,
            int differenceOrder,
            double lowerBoundary,
            double upperBoundary) {
        this.name = Objects.requireNonNull(name, "name");
        if (name.isBlank()) {
            throw new IllegalArgumentException("smooth term name must not be blank");
        }
        this.covariate = MatrixOps.finiteCopy(covariate, "covariate");
        if (covariate.length < 2) {
            throw new IllegalArgumentException(
                "a smooth term requires at least two observations");
        }
        if (degree < 1 || degree > 5) {
            throw new IllegalArgumentException("spline degree must be in [1, 5]");
        }
        if (basisDimension <= degree
                || differenceOrder < 1
                || differenceOrder >= basisDimension) {
            throw new IllegalArgumentException(
                "basis dimension and difference order are incompatible");
        }
        if (!Double.isFinite(lowerBoundary) || !Double.isFinite(upperBoundary)
                || !(upperBoundary > lowerBoundary)) {
            throw new IllegalArgumentException(
                "spline boundaries must be finite and strictly increasing");
        }
        for (double value : covariate) {
            if (value < lowerBoundary || value > upperBoundary) {
                throw new IllegalArgumentException(
                    "covariate lies outside the requested spline boundaries");
            }
        }
        this.basisDimension = basisDimension;
        this.degree = degree;
        this.differenceOrder = differenceOrder;
        this.lowerBoundary = lowerBoundary;
        this.upperBoundary = upperBoundary;
        this.knots = knotSequence(
            basisDimension, degree, lowerBoundary, upperBoundary);
        this.design = basis(covariate);
        this.penaltyFactor = differenceOperator(
            basisDimension, differenceOrder);
    }

    /** Creates a cubic P-spline with a second-difference penalty. */
    public static PSplineTerm of(
            String name, double[] covariate, int basisDimension) {
        return of(name, covariate, basisDimension, 3, 2);
    }

    /** Creates a P-spline whose boundaries are the observed covariate range. */
    public static PSplineTerm of(
            String name,
            double[] covariate,
            int basisDimension,
            int degree,
            int differenceOrder) {
        double[] values = MatrixOps.finiteCopy(covariate, "covariate");
        double minimum = Double.POSITIVE_INFINITY;
        double maximum = Double.NEGATIVE_INFINITY;
        for (double value : values) {
            minimum = Math.min(minimum, value);
            maximum = Math.max(maximum, value);
        }
        return new PSplineTerm(name, values, basisDimension, degree,
            differenceOrder, minimum, maximum);
    }

    /** Creates a P-spline with explicit boundaries, useful for prediction. */
    public static PSplineTerm of(
            String name,
            double[] covariate,
            int basisDimension,
            int degree,
            int differenceOrder,
            double lowerBoundary,
            double upperBoundary) {
        return new PSplineTerm(name, covariate, basisDimension, degree,
            differenceOrder, lowerBoundary, upperBoundary);
    }

    public String name() { return name; }
    public int observations() { return covariate.length; }
    public int basisDimension() { return basisDimension; }
    public int degree() { return degree; }
    public int differenceOrder() { return differenceOrder; }
    public double lowerBoundary() { return lowerBoundary; }
    public double upperBoundary() { return upperBoundary; }
    public double[] covariate() { return covariate.clone(); }
    public double[] knots() { return knots.clone(); }
    public double[] design() { return design.clone(); }
    public double[] penaltyFactor() { return penaltyFactor.clone(); }

    /** Evaluates the fitted knot sequence at new covariate values. */
    public double[] basis(double[] values) {
        double[] checked = MatrixOps.finiteCopy(values, "values");
        double[] result = new double[checked.length * basisDimension];
        for (int row = 0; row < checked.length; row++) {
            double value = checked[row];
            if (value < lowerBoundary || value > upperBoundary) {
                throw new IllegalArgumentException(
                    "prediction value lies outside the spline boundaries");
            }
            evaluateBasis(value, result, row * basisDimension);
        }
        return result;
    }

    double[] designView() { return design; }
    double[] penaltyFactorView() { return penaltyFactor; }

    private void evaluateBasis(double value, double[] output, int offset) {
        if (value == upperBoundary) {
            output[offset + basisDimension - 1] = 1.0;
            return;
        }
        int initialCount = knots.length - 1;
        double[] previous = new double[initialCount];
        for (int index = 0; index < initialCount; index++) {
            if (value >= knots[index] && value < knots[index + 1]) {
                previous[index] = 1.0;
            }
        }
        for (int order = 1; order <= degree; order++) {
            int count = initialCount - order;
            double[] next = new double[count];
            for (int index = 0; index < count; index++) {
                double leftDenominator = knots[index + order] - knots[index];
                double rightDenominator =
                    knots[index + order + 1] - knots[index + 1];
                double left = leftDenominator == 0.0 ? 0.0
                    : (value - knots[index]) / leftDenominator * previous[index];
                double right = rightDenominator == 0.0 ? 0.0
                    : (knots[index + order + 1] - value) / rightDenominator
                        * previous[index + 1];
                next[index] = left + right;
            }
            previous = next;
        }
        System.arraycopy(previous, 0, output, offset, basisDimension);
    }

    private static double[] knotSequence(
            int dimension, int degree, double lower, double upper) {
        double[] result = new double[dimension + degree + 1];
        for (int index = 0; index <= degree; index++) {
            result[index] = lower;
            result[result.length - 1 - index] = upper;
        }
        int interior = dimension - degree - 1;
        for (int index = 1; index <= interior; index++) {
            result[degree + index] = lower
                + (upper - lower) * index / (interior + 1.0);
        }
        return result;
    }

    private static double[] differenceOperator(int dimension, int order) {
        int rows = dimension - order;
        double[] result = new double[rows * dimension];
        for (int row = 0; row < rows; row++) {
            for (int index = 0; index <= order; index++) {
                double sign = ((order - index) & 1) == 0 ? 1.0 : -1.0;
                result[row * dimension + row + index] =
                    sign * binomial(order, index);
            }
        }
        return result;
    }

    private static double binomial(int n, int k) {
        long value = 1L;
        for (int index = 1; index <= k; index++) {
            value = value * (n - index + 1L) / index;
        }
        return value;
    }
}
