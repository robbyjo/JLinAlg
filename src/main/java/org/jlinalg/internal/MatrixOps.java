/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.internal;

import java.util.Arrays;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;

/** Internal row-major matrix helpers. */
public final class MatrixOps {
    private MatrixOps() {
    }

    /** Validates finite response and row-major design data. */
    public static void validateModelData(
            double[] response, double[] design, int rows, int columns) {
        if (rows < 1 || columns < 1) {
            throw new IllegalArgumentException("rows and columns must be positive");
        }
        if (response == null || response.length != rows) {
            throw new IllegalArgumentException("response length must equal rows");
        }
        if (design == null || design.length != rows * columns) {
            throw new IllegalArgumentException(
                "design length must equal rows * columns");
        }
        requireFinite(response, "response");
        requireFinite(design, "design");
    }

    /** Converts a rectangular two-dimensional matrix to row-major storage. */
    public static double[] rowMajor(double[][] matrix, int expectedRows) {
        double[] result = rowMajorUnchecked(matrix, expectedRows);
        requireFinite(result, "matrix");
        return result;
    }

    /** Converts a rectangular matrix while preserving non-finite values. */
    public static double[] rowMajorUnchecked(double[][] matrix, int expectedRows) {
        if (matrix == null || matrix.length != expectedRows || expectedRows < 1) {
            throw new IllegalArgumentException("matrix row count is invalid");
        }
        if (matrix[0] == null || matrix[0].length < 1) {
            throw new IllegalArgumentException("matrix must have at least one column");
        }
        int columns = matrix[0].length;
        double[] result = new double[expectedRows * columns];
        for (int row = 0; row < expectedRows; row++) {
            if (matrix[row] == null || matrix[row].length != columns) {
                throw new IllegalArgumentException("matrix must be rectangular");
            }
            System.arraycopy(matrix[row], 0, result, row * columns, columns);
        }
        return result;
    }

    /** Returns fitted values X beta. */
    public static double[] multiply(
            ComputeBackend backend, double[] matrix, int rows, int columns,
            double[] vector) {
        double[] result = new double[rows];
        backend.dgemv(MatrixTranspose.NONE, rows, columns,
            1.0, matrix, vector, 0.0, result);
        return result;
    }

    /** Returns A B for row-major matrices. */
    public static double[] multiply(
            ComputeBackend backend,
            double[] left, int leftRows, int shared,
            double[] right, int rightColumns) {
        double[] result = new double[leftRows * rightColumns];
        backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.NONE,
            leftRows, rightColumns, shared,
            1.0, left, right, 0.0, result);
        return result;
    }

    /** Returns A' B for row-major matrices with the same row count. */
    public static double[] transposeMultiply(
            ComputeBackend backend,
            double[] left, int rows, int leftColumns,
            double[] right, int rightColumns) {
        double[] result = new double[leftColumns * rightColumns];
        backend.dgemm(MatrixTranspose.TRANSPOSE, MatrixTranspose.NONE,
            leftColumns, rightColumns, rows,
            1.0, left, right, 0.0, result);
        return result;
    }

    /** Returns a square row-major identity matrix. */
    public static double[] identity(int dimension) {
        double[] result = new double[dimension * dimension];
        for (int index = 0; index < dimension; index++) {
            result[index * dimension + index] = 1.0;
        }
        return result;
    }

    /** Returns left - right. */
    public static double[] subtract(double[] left, double[] right) {
        if (left.length != right.length) {
            throw new IllegalArgumentException("array lengths must match");
        }
        double[] result = new double[left.length];
        for (int index = 0; index < left.length; index++) {
            result[index] = left[index] - right[index];
        }
        return result;
    }

    /** Returns a defensive copy after validating every element is finite. */
    public static double[] finiteCopy(double[] values, String name) {
        if (values == null) {
            throw new IllegalArgumentException(name + " is required");
        }
        requireFinite(values, name);
        return values.clone();
    }

    /** Ensures every array element is finite. */
    public static void requireFinite(double[] values, String name) {
        if (Arrays.stream(values).anyMatch(value -> !Double.isFinite(value))) {
            throw new IllegalArgumentException(name + " must contain only finite values");
        }
    }
}
