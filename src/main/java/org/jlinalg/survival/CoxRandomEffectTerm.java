/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;

/** Random-effect design and unit-variance Gaussian precision for Cox frailty. */
public final class CoxRandomEffectTerm {
    private final String name;
    private final int observations;
    private final int coefficients;
    private final double[] design;
    private final double[] precision;
    private final List<String> coefficientNames;

    public CoxRandomEffectTerm(
            String name, double[][] design, double[] precision,
            List<String> coefficientNames) {
        if (design == null || design.length == 0 || design[0] == null)
            throw new IllegalArgumentException("frailty design is required");
        this.name = Objects.requireNonNull(name, "name");
        observations = design.length;
        coefficients = design[0].length;
        this.design = MatrixOps.rowMajor(design, observations);
        if (precision == null || precision.length != coefficients * coefficients
                || coefficientNames == null
                || coefficientNames.size() != coefficients)
            throw new IllegalArgumentException(
                "frailty precision and coefficient names have invalid dimensions");
        this.precision = MatrixOps.finiteCopy(precision, "frailty precision");
        this.coefficientNames = List.copyOf(coefficientNames);
        validatePrecision();
    }

    private CoxRandomEffectTerm(
            String name, int observations, int coefficients,
            double[] design, double[] precision,
            List<String> coefficientNames) {
        this.name = name;
        this.observations = observations;
        this.coefficients = coefficients;
        this.design = design;
        this.precision = precision;
        this.coefficientNames = coefficientNames;
        validatePrecision();
    }

    public static CoxRandomEffectTerm independent(RandomEffectTerm term) {
        if (term == null)
            throw new IllegalArgumentException("random-effect term is required");
        return new CoxRandomEffectTerm(term.name(), term.observations(),
            term.coefficients(), term.design(),
            MatrixOps.identity(term.coefficients()), term.coefficientNames());
    }

    /**
     * Creates a correlated Gaussian term from a coefficient covariance.
     * A small relative diagonal regularization supports empirical GRMs that
     * are singular because of duplicate samples or finite marker rank.
     */
    public static CoxRandomEffectTerm fromCovariance(
            String name,
            double[][] design,
            double[] covariance,
            List<String> coefficientNames,
            double relativeDiagonalRegularization,
            BackendPolicy backendPolicy) {
        if (design == null || design.length == 0 || design[0] == null
                || covariance == null || coefficientNames == null
                || coefficientNames.isEmpty() || backendPolicy == null
                || !Double.isFinite(relativeDiagonalRegularization)
                || relativeDiagonalRegularization <= 0)
            throw new IllegalArgumentException(
                "covariance term inputs and positive regularization are required");
        int dimension = coefficientNames.size();
        if (covariance.length != dimension * dimension
                || design[0].length != dimension)
            throw new IllegalArgumentException(
                "covariance dimensions must match random coefficients");
        double[] regularized = MatrixOps.finiteCopy(
            covariance, "random-effect covariance");
        double diagonalScale = 0;
        double maximum = 0;
        for (double value : regularized)
            maximum = Math.max(maximum, Math.abs(value));
        double symmetryTolerance = 1e-10 * Math.max(1, maximum);
        for (int index = 0; index < dimension; index++) {
            if (!(regularized[index * dimension + index] > 0))
                throw new IllegalArgumentException(
                    "random-effect covariance diagonal must be positive");
            diagonalScale += regularized[index * dimension + index];
            for (int column = 0; column < index; column++)
                if (Math.abs(regularized[index * dimension + column]
                        - regularized[column * dimension + index])
                        > symmetryTolerance)
                    throw new IllegalArgumentException(
                        "random-effect covariance must be symmetric");
        }
        diagonalScale /= dimension;
        for (int index = 0; index < dimension; index++)
            regularized[index * dimension + index] +=
                relativeDiagonalRegularization * diagonalScale;
        double[] precision;
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            precision = context.backend().dpotrf(regularized, dimension)
                .solve(MatrixOps.identity(dimension), dimension);
        } catch (IllegalArgumentException | IllegalStateException
                | ArithmeticException exception) {
            throw new IllegalArgumentException(
                "regularized random-effect covariance is not positive definite",
                exception);
        }
        for (int row = 0; row < dimension; row++)
            for (int column = 0; column < row; column++) {
                double value = 0.5 * (precision[row * dimension + column]
                    + precision[column * dimension + row]);
                precision[row * dimension + column] = value;
                precision[column * dimension + row] = value;
            }
        return new CoxRandomEffectTerm(
            name, design, precision, coefficientNames);
    }

    public String name() { return name; }
    public int observations() { return observations; }
    public int coefficients() { return coefficients; }
    public double[] design() { return design.clone(); }
    public double[] precision() { return precision.clone(); }
    public List<String> coefficientNames() { return coefficientNames; }
    double[] designView() { return design; }
    double[] precisionView() { return precision; }

    private void validatePrecision() {
        if (name == null || name.isBlank() || coefficients < 1)
            throw new IllegalArgumentException("frailty term name/dimensions are invalid");
        for (int row = 0; row < coefficients; row++) {
            if (!(precision[row * coefficients + row] > 0))
                throw new IllegalArgumentException(
                    "frailty precision diagonal must be positive");
            for (int column = 0; column < row; column++) {
                double left = precision[row * coefficients + column];
                double right = precision[column * coefficients + row];
                if (Math.abs(left - right) > 1e-10
                        * Math.max(1, Math.max(Math.abs(left), Math.abs(right))))
                    throw new IllegalArgumentException(
                        "frailty precision must be symmetric");
            }
        }
    }
}
