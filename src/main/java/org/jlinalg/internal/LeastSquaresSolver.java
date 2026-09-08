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
        int equilibratedRank = -1;
        double[] referenceDesign = null;
        if (rows >= columns) {
            double[] scales = new double[columns];
            for (int row = 0; row < rows; row++) {
                for (int column = 0; column < columns; column++) {
                    scales[column] = Math.max(scales[column],
                        Math.abs(design[row * columns + column]));
                }
            }
            boolean needsScaling = false;
            double maximumScale = 0.0;
            double minimumScale = Double.POSITIVE_INFINITY;
            for (double scale : scales) {
                needsScaling |= scale != 0.0 && (scale < 1e-6 || scale > 1e6);
                maximumScale = Math.max(maximumScale, scale);
                if (scale > 0.0) minimumScale = Math.min(minimumScale, scale);
            }
            needsScaling |= maximumScale / minimumScale > 1e6;
            double[] equilibrated = needsScaling ? design.clone() : design;
            for (int column = 0; column < columns; column++) {
                if (!needsScaling || scales[column] == 0.0) scales[column] = 1.0;
                if (needsScaling) {
                    for (int row = 0; row < rows; row++) {
                        equilibrated[row * columns + column] /= scales[column];
                    }
                }
            }
            if (needsScaling) referenceDesign = equilibrated;
            PivotedQrFactor qr = backend.dgeqp3(equilibrated, rows, columns);
            equilibratedRank = qr.rank();
            if (qr.rank() == columns) {
                double[] coefficients = qr.solveLeastSquares(response);
                double[] covariance = qrCovariance(qr);
                for (int row = 0; row < columns; row++) {
                    coefficients[row] /= scales[row];
                    for (int column = 0; column < columns; column++) {
                        covariance[row * columns + column] =
                            covariance[row * columns + column] / scales[row] / scales[column];
                    }
                }
                return new Solution(coefficients, covariance, columns, false, 0.0);
            }
            if (!allowMinimumNorm) {
                throw new IllegalArgumentException(
                    "design matrix is rank deficient: rank " + qr.rank()
                    + " < " + columns);
            }
        } else if (!allowMinimumNorm) {
            throw new IllegalArgumentException(
                "design matrix has more columns than rows and cannot have full column rank");
        } else {
            // Wide systems need scale-aware rank identification as well.
            referenceDesign = design.clone();
            for (int column = 0; column < columns; column++) {
                double scale = 0.0;
                for (int row = 0; row < rows; row++) {
                    scale = Math.max(scale, Math.abs(design[row * columns + column]));
                }
                if (scale > 0.0) for (int row = 0; row < rows; row++) {
                    referenceDesign[row * columns + column] /= scale;
                }
            }
        }
        // Equilibrating a deficient design changes the minimizing norm. Keep
        // the requested minimum norm in the caller's original coordinates.
        Solution reference = referenceDesign == null ? null
            : svdSolution(referenceDesign, response, rows, columns, equilibratedRank, backend);
        if (reference != null) equilibratedRank = reference.rank();
        Solution original = svdSolution(design, response, rows, columns, equilibratedRank, backend);
        if (reference != null) validateProjection(design, referenceDesign, response,
            rows, columns, original, reference);
        return original;
    }

    private static void validateProjection(double[] design, double[] referenceDesign,
            double[] response, int rows, int columns, Solution original, Solution reference) {
        double responseScale = Double.MIN_NORMAL;
        for (double value : response) responseScale = Math.max(responseScale, Math.abs(value));
        for (int row = 0; row < rows; row++) {
            double predicted = 0.0, expected = 0.0, leverage = 0.0, referenceLeverage = 0.0;
            for (int j = 0; j < columns; j++) {
                double x = design[row * columns + j], z = referenceDesign[row * columns + j];
                predicted += x * original.coefficients()[j];
                expected += z * reference.coefficients()[j];
                for (int k = 0; k < columns; k++) {
                    leverage += x * original.unscaledCovariance()[j * columns + k]
                        * design[row * columns + k];
                    referenceLeverage += z * reference.unscaledCovariance()[j * columns + k]
                        * referenceDesign[row * columns + k];
                }
            }
            if (!Double.isFinite(predicted) || !Double.isFinite(leverage)
                    || Math.abs(predicted / responseScale - expected / responseScale) > 1e-8
                    || Math.abs(leverage - referenceLeverage) > 1e-7 * (1.0 + Math.abs(referenceLeverage))) {
                throw new IllegalArgumentException("original-coordinate minimum-norm SVD is inaccurate "
                    + "at these column scales; reparameterize the deficient design explicitly");
            }
        }
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
            int equilibratedRank,
            ComputeBackend backend) {
        SingularValueDecomposition svd;
        try {
            svd = backend.dgesvd(design, rows, columns);
        } catch (IllegalArgumentException | IllegalStateException failure) {
            throw new IllegalArgumentException("original-coordinate minimum-norm SVD could not "
                + "resolve the design accurately; reparameterize the deficient design explicitly", failure);
        }
        double[] singularValues = svd.singularValues();
        double tolerance = Math.max(rows, columns) * Math.ulp(1.0)
            * (singularValues.length == 0 ? 0.0 : singularValues[0]);
        int rank = 0;
        while (rank < singularValues.length && singularValues[rank] > tolerance) {
            rank++;
        }
        if (equilibratedRank >= 0) {
            rank = equilibratedRank;
            if (rank > 0 && !(singularValues[rank - 1] > 0.0)) {
                throw new IllegalArgumentException(
                    "original-coordinate SVD cannot resolve the equilibrated design rank");
            }
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
        double[] projection = new double[columns * columns];
        for (int row = 0; row < columns; row++) {
            for (int column = 0; column < columns; column++) {
                for (int component = 0; component < rank; component++) {
                    projection[row * columns + column] +=
                        rightTransposed[component * columns + row]
                            * rightTransposed[component * columns + column];
                }
            }
        }
        return new Solution(coefficients, covariance, rank, true, tolerance, projection);
    }

    /** Whether a coordinate is uniquely estimable, not just fixed by minimum norm. */
    public static boolean estimableCoordinate(double[] projection, int dimension, int column) {
        return projection == null || Math.abs(1.0 - projection[column * dimension + column]) <= 1e-8;
    }

    /** Rejects contrasts outside the fitted design's row space. Null means full rank. */
    public static void requireEstimable(double[][] contrast, double[] projection, int dimension) {
        if (projection == null) return;
        if (contrast == null) throw new IllegalArgumentException("contrast is required");
        for (double[] row : contrast) {
            if (row == null || row.length != dimension) {
                throw new IllegalArgumentException("contrast columns must equal coefficient count");
            }
            double magnitude = 0.0;
            for (double value : row) magnitude = Math.max(magnitude, Math.abs(value));
            if (magnitude == 0.0) continue;
            for (int column = 0; column < dimension; column++) {
                double projected = 0.0;
                for (int shared = 0; shared < dimension; shared++) {
                    projected += row[shared] / magnitude * projection[shared * dimension + column];
                }
                if (!Double.isFinite(projected)
                        || Math.abs(projected - row[column] / magnitude) > 1e-8) {
                    throw new IllegalArgumentException("contrast is not estimable in the fitted design");
                }
            }
        }
    }

    /** Least-squares coefficients, unscaled covariance, and rank metadata. */
    public record Solution(
            double[] coefficients,
            double[] unscaledCovariance,
            int rank,
            boolean minimumNorm,
            double tolerance,
            double[] rowSpaceProjection) {
        /** Compatibility constructor for solutions constructed by callers. */
        public Solution(double[] coefficients, double[] unscaledCovariance,
                int rank, boolean minimumNorm, double tolerance) {
            this(coefficients, unscaledCovariance, rank, minimumNorm, tolerance, null);
        }
    }
}
