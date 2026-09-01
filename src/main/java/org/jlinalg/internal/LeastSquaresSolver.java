/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.internal;

import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.PivotedQrFactor;
import jdistlib.accelerator.SingularValueDecomposition;

/** Shared pivoted-QR and minimum-norm SVD least-squares solver. */
public final class LeastSquaresSolver {
    private LeastSquaresSolver() {
    }

    /** Solves a row-major least-squares problem and returns unscaled covariance. */
    public static Solution solve(
            double[] design,
            double[] response,
            int rows,
            int columns,
            boolean allowMinimumNorm,
            ComputeBackend backend) {
        if (rows >= columns) {
            PivotedQrFactor qr = backend.dgeqp3(design, rows, columns);
            if (qr.rank() == columns) {
                return new Solution(
                    qr.solveLeastSquares(response),
                    qrCovariance(qr), columns, false, 0.0);
            }
            if (!allowMinimumNorm) {
                throw new IllegalArgumentException(
                    "design matrix is rank deficient: rank " + qr.rank()
                    + " < " + columns);
            }
        } else if (!allowMinimumNorm) {
            throw new IllegalArgumentException(
                "design matrix has more columns than rows and cannot have full column rank");
        }
        return svdSolution(design, response, rows, columns, backend);
    }

    private static double[] qrCovariance(PivotedQrFactor qr) {
        int columns = qr.columns();
        double[] packed = qr.packed();
        double[] inverse = new double[columns * columns];
        for (int rightSide = 0; rightSide < columns; rightSide++) {
            for (int row = columns - 1; row >= 0; row--) {
                double value = row == rightSide ? 1.0 : 0.0;
                for (int column = row + 1; column < columns; column++) {
                    value -= packed[row * columns + column]
                        * inverse[column * columns + rightSide];
                }
                inverse[row * columns + rightSide] =
                    value / packed[row * columns + row];
            }
        }

        double[] permutedCovariance = new double[columns * columns];
        for (int row = 0; row < columns; row++) {
            for (int column = 0; column < columns; column++) {
                double value = 0.0;
                for (int shared = 0; shared < columns; shared++) {
                    value += inverse[row * columns + shared]
                        * inverse[column * columns + shared];
                }
                permutedCovariance[row * columns + column] = value;
            }
        }

        int[] pivot = qr.pivot();
        double[] covariance = new double[columns * columns];
        for (int row = 0; row < columns; row++) {
            for (int column = 0; column < columns; column++) {
                covariance[pivot[row] * columns + pivot[column]] =
                    permutedCovariance[row * columns + column];
            }
        }
        return covariance;
    }

    private static Solution svdSolution(
            double[] design,
            double[] response,
            int rows,
            int columns,
            ComputeBackend backend) {
        SingularValueDecomposition svd = backend.dgesvd(design, rows, columns);
        double[] singularValues = svd.singularValues();
        double tolerance = Math.max(rows, columns) * Math.ulp(1.0)
            * (singularValues.length == 0 ? 0.0 : singularValues[0]);
        int rank = 0;
        while (rank < singularValues.length && singularValues[rank] > tolerance) {
            rank++;
        }
        if (rank == 0) {
            throw new IllegalArgumentException("design matrix has zero numerical rank");
        }

        int components = svd.components();
        double[] left = svd.leftSingularVectors();
        double[] rightTransposed = svd.rightSingularVectorsTransposed();
        double[] scaledProjection = new double[rank];
        for (int component = 0; component < rank; component++) {
            double projection = 0.0;
            for (int row = 0; row < rows; row++) {
                projection += left[row * components + component] * response[row];
            }
            scaledProjection[component] = projection / singularValues[component];
        }

        double[] coefficients = new double[columns];
        for (int column = 0; column < columns; column++) {
            for (int component = 0; component < rank; component++) {
                coefficients[column] += rightTransposed[component * columns + column]
                    * scaledProjection[component];
            }
        }

        double[] covariance = new double[columns * columns];
        for (int row = 0; row < columns; row++) {
            for (int column = 0; column < columns; column++) {
                double value = 0.0;
                for (int component = 0; component < rank; component++) {
                    value += rightTransposed[component * columns + row]
                        * rightTransposed[component * columns + column]
                        / (singularValues[component] * singularValues[component]);
                }
                covariance[row * columns + column] = value;
            }
        }
        return new Solution(coefficients, covariance, rank, true, tolerance);
    }

    /** Least-squares coefficients, unscaled covariance, and rank metadata. */
    public record Solution(
            double[] coefficients,
            double[] unscaledCovariance,
            int rank,
            boolean minimumNorm,
            double tolerance) {
    }
}
