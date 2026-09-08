/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

/** Numerical validation of caller-supplied genetic covariance matrices. */
public final class GeneticCovarianceValidation {
    private GeneticCovarianceValidation() { }

    /**
     * Checks finite symmetry and positive semidefiniteness, allowing singular
     * matrices and zero rows. Uses diagonally scaled, pivoted Cholesky with
     * roundoff tolerance {@code 1e-10}. Does not repair or mutate the input.
     * Work is O(n cubed), storage O(n squared); use only for external matrices,
     * not for cross-products already known to be positive semidefinite.
     */
    public static void requirePositiveSemidefinite(double[] matrix, int n) {
        if (n < 1 || matrix == null || (long) n * n != matrix.length)
            throw new IllegalArgumentException("covariance dimensions are invalid");
        double[] scale = new double[n];
        for (int i = 0; i < n; i++) {
            double diagonal = matrix[i * n + i];
            if (!Double.isFinite(diagonal) || diagonal < 0)
                throw new IllegalArgumentException("covariance diagonal must be finite and nonnegative");
            scale[i] = Math.sqrt(diagonal);
        }
        double[] a = new double[matrix.length];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
            double x = matrix[i * n + j], y = matrix[j * n + i];
            if (!Double.isFinite(x) || !Double.isFinite(y))
                throw new IllegalArgumentException("covariance must be finite");
            if (scale[i] == 0 || scale[j] == 0) {
                if (x != 0 || y != 0)
                    throw new IllegalArgumentException("a zero variance requires a zero covariance row");
            } else {
                x = x / scale[i] / scale[j];
                y = y / scale[i] / scale[j];
                if (Math.abs(x - y) > 1e-10 || Math.abs(x) > 1 + 1e-10)
                    throw new IllegalArgumentException("covariance must be symmetric positive semidefinite");
                a[i * n + j] = a[j * n + i] = 0.5 * x + 0.5 * y;
            }
        }
        int[] order = new int[n];
        for (int i = 0; i < n; i++) order[i] = i;
        for (int k = 0; k < n; k++) {
            int best = k;
            for (int i = k + 1; i < n; i++)
                if (a[order[i] * n + order[i]] > a[order[best] * n + order[best]]) best = i;
            int swap = order[k]; order[k] = order[best]; order[best] = swap;
            int pivot = order[k];
            double diagonal = a[pivot * n + pivot];
            if (diagonal <= 1e-10) {
                for (int i = k; i < n; i++) for (int j = k; j <= i; j++)
                    if (Math.abs(a[order[i] * n + order[j]]) > 1e-10)
                        throw new IllegalArgumentException("covariance is not positive semidefinite");
                return;
            }
            for (int i = k + 1; i < n; i++) {
                int row = order[i];
                double factor = a[row * n + pivot] / diagonal;
                for (int j = k + 1; j <= i; j++) {
                    int column = order[j];
                    double value = Math.fma(-factor, a[column * n + pivot], a[row * n + column]);
                    a[row * n + column] = a[column * n + row] = value;
                }
            }
        }
    }
}
