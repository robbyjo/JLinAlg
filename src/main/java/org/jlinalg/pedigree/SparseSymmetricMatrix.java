/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

/** Immutable symmetric matrix stored in full compressed-sparse-row form. */
public final class SparseSymmetricMatrix {
    private final int dimension;
    private final int[] rowPointers;
    private final int[] columnIndices;
    private final double[] values;

    SparseSymmetricMatrix(
            int dimension,
            int[] rowPointers,
            int[] columnIndices,
            double[] values) {
        this.dimension = dimension;
        this.rowPointers = rowPointers.clone();
        this.columnIndices = columnIndices.clone();
        this.values = values.clone();
    }

    public int dimension() { return dimension; }
    public int nonzeroCount() { return values.length; }
    public int[] rowPointers() { return rowPointers.clone(); }
    public int[] columnIndices() { return columnIndices.clone(); }
    public double[] values() { return values.clone(); }

    /** Returns one element, using binary search within its CSR row. */
    public double get(int row, int column) {
        checkIndex(row);
        checkIndex(column);
        int low = rowPointers[row];
        int high = rowPointers[row + 1] - 1;
        while (low <= high) {
            int middle = (low + high) >>> 1;
            int candidate = columnIndices[middle];
            if (candidate < column) {
                low = middle + 1;
            } else if (candidate > column) {
                high = middle - 1;
            } else {
                return values[middle];
            }
        }
        return 0.0;
    }

    /** Multiplies this matrix by a dense vector. */
    public double[] multiply(double[] vector) {
        if (vector == null || vector.length != dimension) {
            throw new IllegalArgumentException(
                "vector length must equal the sparse matrix dimension");
        }
        double[] result = new double[dimension];
        for (int row = 0; row < dimension; row++) {
            double value = 0.0;
            for (int index = rowPointers[row];
                    index < rowPointers[row + 1]; index++) {
                value += values[index] * vector[columnIndices[index]];
            }
            result[row] = value;
        }
        return result;
    }

    /** Materializes a row-major dense matrix for validation or small problems. */
    public double[] toDense() {
        double[] result = new double[dimension * dimension];
        for (int row = 0; row < dimension; row++) {
            for (int index = rowPointers[row];
                    index < rowPointers[row + 1]; index++) {
                result[row * dimension + columnIndices[index]] = values[index];
            }
        }
        return result;
    }

    private void checkIndex(int index) {
        if (index < 0 || index >= dimension) {
            throw new IndexOutOfBoundsException("matrix index: " + index);
        }
    }
}
