/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.internal.MatrixOps;

final class CoxMath {
    private CoxMath() { }

    static double[] solvePositive(
            ComputeBackend backend, double[] information,
            int dimension, double[] rightSide, double baseRidge) {
        return factor(backend, information, dimension, baseRidge)
            .solve(rightSide);
    }

    static double[] inversePositive(
            ComputeBackend backend, double[] information,
            int dimension, double baseRidge) {
        return factor(backend, information, dimension, baseRidge)
            .solve(MatrixOps.identity(dimension), dimension);
    }

    static CholeskyFactor factor(
            ComputeBackend backend, double[] matrix,
            int dimension, double baseRidge) {
        double scale = 0;
        for (int index = 0; index < dimension; index++)
            scale = Math.max(scale,
                Math.abs(matrix[index * dimension + index]));
        scale = Math.max(1, scale);
        for (int attempt = 0; attempt < 10; attempt++) {
            double[] regularized = matrix.clone();
            double ridge = attempt == 0 ? 0
                : Math.max(baseRidge, 1e-12)
                    * scale * Math.pow(10, attempt - 1);
            for (int index = 0; index < dimension; index++)
                regularized[index * dimension + index] += ridge;
            try {
                return backend.dpotrf(regularized, dimension);
            } catch (IllegalArgumentException | ArithmeticException exception) {
                if (attempt == 9) throw new IllegalArgumentException(
                    "Cox information matrix is not positive definite", exception);
            }
        }
        throw new IllegalStateException("unreachable Cox factorization path");
    }

    static double maximumAbsolute(double[] values) {
        double result = 0;
        for (double value : values) result = Math.max(result, Math.abs(value));
        return result;
    }
}
