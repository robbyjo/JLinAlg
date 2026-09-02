/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

/** Fitted values and smoothing metadata for one additive smooth. */
public final class SmoothTermEstimate {
    private final PSplineTerm term;
    private final double[] fixedTransform;
    private final double[] fixedMeans;
    private final double[] fixedCoefficients;
    private final double[] randomTransform;
    private final double[] randomMeans;
    private final double[] randomCoefficients;
    private final double[] fittedValues;
    private final double smoothingParameter;
    private final double effectiveDegreesOfFreedom;

    SmoothTermEstimate(
            PSplineTerm term,
            double[] fixedTransform,
            double[] fixedMeans,
            double[] fixedCoefficients,
            double[] randomTransform,
            double[] randomMeans,
            double[] randomCoefficients,
            double[] fittedValues,
            double smoothingParameter,
            double effectiveDegreesOfFreedom) {
        this.term = term;
        this.fixedTransform = fixedTransform.clone();
        this.fixedMeans = fixedMeans.clone();
        this.fixedCoefficients = fixedCoefficients.clone();
        this.randomTransform = randomTransform.clone();
        this.randomMeans = randomMeans.clone();
        this.randomCoefficients = randomCoefficients.clone();
        this.fittedValues = fittedValues.clone();
        this.smoothingParameter = smoothingParameter;
        this.effectiveDegreesOfFreedom = effectiveDegreesOfFreedom;
    }

    public String name() { return term.name(); }
    public double smoothingParameter() { return smoothingParameter; }
    public double effectiveDegreesOfFreedom() {
        return effectiveDegreesOfFreedom;
    }
    public double[] fittedValues() { return fittedValues.clone(); }

    /** Predicts this centered smooth at new covariate values. */
    public double[] predict(double[] covariate) {
        double[] basis = term.basis(covariate);
        int rows = covariate.length;
        int basisColumns = term.basisDimension();
        double[] result = new double[rows];
        addContribution(result, basis, rows, basisColumns,
            fixedTransform, fixedMeans, fixedCoefficients);
        addContribution(result, basis, rows, basisColumns,
            randomTransform, randomMeans, randomCoefficients);
        return result;
    }

    private static void addContribution(
            double[] result,
            double[] basis,
            int rows,
            int basisColumns,
            double[] transform,
            double[] means,
            double[] coefficients) {
        int columns = coefficients.length;
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                double value = -means[column];
                for (int shared = 0; shared < basisColumns; shared++) {
                    value += basis[row * basisColumns + shared]
                        * transform[shared * columns + column];
                }
                result[row] += value * coefficients[column];
            }
        }
    }
}
