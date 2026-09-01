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

/** Fast Wald tests for one or more linear contrasts of model coefficients. */
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
                || covariance.length != beta.length * beta.length
                || contrast == null || contrast.length == 0) {
            throw new IllegalArgumentException("contrast dimensions are invalid");
        }
        int rows = contrast.length;
        double[] matrix = MatrixOps.rowMajor(contrast, rows);
        if (contrast[0].length != beta.length) {
            throw new IllegalArgumentException(
                "contrast columns must equal coefficient count");
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
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
            return new Components(estimates, contrastCovariance, rows, quadratic);
        } catch (IllegalArgumentException | IllegalStateException exception) {
            throw new IllegalArgumentException(
                "contrasts must be linearly independent and estimable", exception);
        }
    }

    private record Components(
            double[] estimates, double[] covariance,
            int rows, double quadratic) { }
}
