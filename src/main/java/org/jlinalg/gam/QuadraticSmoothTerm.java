/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import org.jlinalg.internal.MatrixOps;

/** A compiled smooth basis with one or more quadratic coefficient penalties. */
public final class QuadraticSmoothTerm {
    private final String name;
    private final int observations;
    private final int columns;
    private final double[] design;
    private final List<double[]> penalties;

    /** Creates a validated row-major basis and symmetric penalty collection. */
    public QuadraticSmoothTerm(
            String name,
            int observations,
            int columns,
            double[] design,
            List<double[]> penalties) {
        this.name = Objects.requireNonNull(name, "name");
        if (name.isBlank() || observations < 2 || columns < 1
                || design == null || design.length != observations * columns
                || penalties == null || penalties.isEmpty()) {
            throw new IllegalArgumentException(
                "smooth name, basis dimensions, and penalties are required");
        }
        this.observations = observations;
        this.columns = columns;
        this.design = MatrixOps.finiteCopy(design, "smooth design");
        List<double[]> copied = new ArrayList<>(penalties.size());
        for (double[] penalty : penalties) {
            if (penalty == null || penalty.length != columns * columns) {
                throw new IllegalArgumentException(
                    "each penalty must be square with one row per basis column");
            }
            double[] value = MatrixOps.finiteCopy(penalty, "smooth penalty");
            validateSymmetry(value, columns);
            copied.add(value);
        }
        this.penalties = List.copyOf(copied);
    }

    /** Converts a P-spline difference factor to its quadratic penalty. */
    public static QuadraticSmoothTerm from(PSplineTerm term) {
        Objects.requireNonNull(term, "term");
        int columns = term.basisDimension();
        int rows = columns - term.differenceOrder();
        double[] factor = term.penaltyFactorView();
        double[] penalty = new double[columns * columns];
        for (int first = 0; first < columns; first++) {
            for (int second = 0; second <= first; second++) {
                double value = 0.0;
                for (int row = 0; row < rows; row++) {
                    value += factor[row * columns + first]
                        * factor[row * columns + second];
                }
                penalty[first * columns + second] = value;
                penalty[second * columns + first] = value;
            }
        }
        return new QuadraticSmoothTerm(term.name(), term.observations(),
            columns, term.designView(), List.of(penalty));
    }

    public String name() { return name; }
    public int observations() { return observations; }
    public int columns() { return columns; }
    public int penaltyCount() { return penalties.size(); }
    public double[] design() { return design.clone(); }
    public List<double[]> penalties() {
        List<double[]> result = new ArrayList<>(penalties.size());
        for (double[] penalty : penalties) result.add(penalty.clone());
        return List.copyOf(result);
    }

    double[] designView() { return design; }
    List<double[]> penaltyViews() { return penalties; }

    private static void validateSymmetry(double[] matrix, int dimension) {
        double maximum = 0.0;
        for (double value : matrix) maximum = Math.max(maximum, Math.abs(value));
        double tolerance = 1e-11 * Math.max(1.0, maximum);
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < row; column++) {
                if (Math.abs(matrix[row * dimension + column]
                        - matrix[column * dimension + row]) > tolerance) {
                    throw new IllegalArgumentException(
                        "smooth penalties must be symmetric");
                }
            }
        }
    }
}
