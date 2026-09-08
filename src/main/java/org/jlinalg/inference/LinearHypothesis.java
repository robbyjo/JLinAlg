/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

import jdistlib.F;
import jdistlib.Gamma;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Wald tests for one or more linear contrasts of model coefficients.
 * Hypothesis rows are normalized internally. Exposed estimates/covariance retain
 * original contrast units and can underflow/overflow when those units are not
 * representable, even when the unit-invariant test statistic is representable. */
public final class LinearHypothesis {
    private LinearHypothesis() { }

    /** Tests {@code C beta = 0} using a finite-denominator F statistic. */
    public static ContrastTestResult fTest(
            double[] beta,
            double[] covariance,
            double[][] contrast,
            double denominatorDegreesOfFreedom) {
        return test(beta, covariance, contrast,
            denominatorDegreesOfFreedom, BackendPolicy.PREFERRED);
    }

    /** Tests {@code C beta = 0} using a finite-denominator F statistic. */
    public static ContrastTestResult test(
            double[] beta,
            double[] covariance,
            double[][] contrast,
            double denominatorDegreesOfFreedom,
            BackendPolicy backendPolicy) {
        if (!(denominatorDegreesOfFreedom > 0.0)) {
            throw new IllegalArgumentException(
                "denominator degrees of freedom must be positive");
        }
        Components values = components(beta, covariance, contrast, backendPolicy);
        double statistic = values.quadratic() / values.rows();
        double pValue = F.cumulative(statistic, values.rows(),
            denominatorDegreesOfFreedom, false, false);
        return new ContrastTestResult(values.estimates(), values.covariance(),
            values.rows(), denominatorDegreesOfFreedom,
            statistic, pValue, StatisticDistribution.F);
    }

    /** Tests {@code C beta = 0} using an asymptotic chi-square statistic. */
    public static ContrastTestResult chiSquareTest(
            double[] beta,
            double[] covariance,
            double[][] contrast) {
        Components values = components(
            beta, covariance, contrast, BackendPolicy.PREFERRED);
        double statistic = values.quadratic();
        double pValue = Gamma.cumulative(statistic,
            values.rows() / 2.0, 2.0, false, false);
        return new ContrastTestResult(values.estimates(), values.covariance(),
            values.rows(), Double.POSITIVE_INFINITY,
            statistic, pValue, StatisticDistribution.CHI_SQUARE);
    }

    private static Components components(
            double[] beta,
            double[] covariance,
            double[][] contrast,
            BackendPolicy backendPolicy) {
        if (beta == null || covariance == null
                || beta.length == 0 || covariance.length != (long) beta.length * beta.length
                || contrast == null || contrast.length == 0) {
            throw new IllegalArgumentException("contrast dimensions are invalid");
        }
        MatrixOps.requireFinite(beta, "coefficients");
        MatrixOps.requireFinite(covariance, "covariance");
        for (int i = 0; i < beta.length; i++) {
            if (covariance[i * beta.length + i] < 0.0) {
                throw new IllegalArgumentException("covariance diagonal must be nonnegative");
            }
            for (int j = 0; j < i; j++) {
                double a = covariance[i * beta.length + j], b = covariance[j * beta.length + i];
                double scale = Math.max(Math.max(Math.abs(a), Math.abs(b)),
                    Math.sqrt(covariance[i * beta.length + i])
                        * Math.sqrt(covariance[j * beta.length + j]));
                if (Math.abs(a - b) > 1e-12 * scale) {
                    throw new IllegalArgumentException("covariance must be symmetric");
                }
            }
        }
        int rows = contrast.length;
        double[] matrix = MatrixOps.rowMajor(contrast, rows);
        if (contrast[0].length != beta.length) {
            throw new IllegalArgumentException(
                "contrast columns must equal coefficient count");
        }
        double[] rowScales = new double[rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < beta.length; column++) {
                rowScales[row] = Math.max(rowScales[row], Math.abs(matrix[row * beta.length + column]));
            }
            if (rowScales[row] == 0) throw new IllegalArgumentException("contrast rows must be nonzero");
            for (int column = 0; column < beta.length; column++) matrix[row * beta.length + column] /= rowScales[row];
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            requireCovariance(covariance, beta.length, backend);
            double[] estimates = MatrixOps.multiply(
                backend, matrix, rows, beta.length, beta);
            double[] temporary = MatrixOps.multiply(
                backend, matrix, rows, beta.length,
                covariance, beta.length);
            double[] contrastCovariance = new double[rows * rows];
            jdistlib.accelerator.MatrixTranspose transpose =
                jdistlib.accelerator.MatrixTranspose.TRANSPOSE;
            backend.dgemm(jdistlib.accelerator.MatrixTranspose.NONE, transpose,
                rows, rows, beta.length, 1.0,
                temporary, matrix, 0.0, contrastCovariance);
            double[] solved = backend.dpotrf(
                contrastCovariance, rows).solve(estimates);
            double quadratic = backend.ddot(rows,
                estimates, 0, 1, solved, 0, 1);
            for (int row = 0; row < rows; row++) {
                estimates[row] *= rowScales[row];
                for (int column = 0; column < rows; column++) {
                    contrastCovariance[row * rows + column] = scaledProduct(
                        contrastCovariance[row * rows + column], rowScales[row], rowScales[column]);
                }
            }
            return new Components(estimates, contrastCovariance, rows, quadratic);
        } catch (IllegalArgumentException | IllegalStateException exception) {
            throw new IllegalArgumentException(
                "contrasts must be linearly independent and estimable", exception);
        }
    }

    private static double scaledProduct(double value, double first, double second) {
        if (value == 0.0 || !Double.isFinite(value)) return value;
        int a = Math.getExponent(value), b = Math.getExponent(first), c = Math.getExponent(second);
        return Math.scalb(Math.scalb(value, -a) * Math.scalb(first, -b) * Math.scalb(second, -c), a + b + c);
    }

    private static void requireCovariance(double[] covariance, int dimension, ComputeBackend backend) {
        double[] standardized = new double[covariance.length];
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                double first = Math.sqrt(covariance[row * dimension + row]);
                double second = Math.sqrt(covariance[column * dimension + column]);
                double value = covariance[row * dimension + column];
                if (first == 0 || second == 0) {
                    if (value != 0) throw new IllegalArgumentException("zero covariance diagonal requires a zero row");
                } else {
                    standardized[row * dimension + column] = value / first / second;
                    if (!Double.isFinite(standardized[row * dimension + column])) {
                        throw new IllegalArgumentException("invalid covariance scaling");
                    }
                }
            }
        }
        for (double eigenvalue : backend.dsyev(standardized, dimension).eigenvalues()) {
            if (eigenvalue < -1e-10 * dimension) throw new IllegalArgumentException("covariance must be positive semidefinite");
        }
    }

    private record Components(
            double[] estimates, double[] covariance,
            int rows, double quadratic) { }
}
