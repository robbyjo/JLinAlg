/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

/** Immutable zero-based full-CSR symmetric random-effect precision matrix. */
public final class SparsePrecisionMatrix {
    private final int dimension;
    private final int[] rowStarts;
    private final int[] columnIndices;
    private final double[] values;

    public SparsePrecisionMatrix(
            int dimension, int[] rowStarts,
            int[] columnIndices, double[] values) {
        if (dimension < 1 || rowStarts == null || columnIndices == null
                || values == null || rowStarts.length != dimension + 1
                || columnIndices.length != values.length
                || rowStarts[0] != 0
                || rowStarts[dimension] != values.length)
            throw new IllegalArgumentException(
                "sparse precision dimensions are invalid");
        for (int row = 0; row < dimension; row++) {
            if (rowStarts[row] > rowStarts[row + 1])
                throw new IllegalArgumentException(
                    "sparse precision row starts are not monotone");
            int previous = -1;
            for (int index = rowStarts[row];
                    index < rowStarts[row + 1]; index++) {
                int column = columnIndices[index];
                if (column < 0 || column >= dimension || column <= previous
                        || !Double.isFinite(values[index]))
                    throw new IllegalArgumentException(
                        "sparse precision entries are invalid");
                previous = column;
            }
        }
        this.dimension = dimension;
        this.rowStarts = rowStarts.clone();
        this.columnIndices = columnIndices.clone();
        this.values = values.clone();
    }

    public static SparsePrecisionMatrix identity(int dimension) {
        int[] starts = new int[dimension + 1];
        int[] columns = new int[dimension];
        double[] values = new double[dimension];
        for (int index = 0; index < dimension; index++) {
            starts[index] = index;
            columns[index] = index;
            values[index] = 1.0;
        }
        starts[dimension] = dimension;
        return new SparsePrecisionMatrix(
            dimension, starts, columns, values);
    }

    public int dimension() { return dimension; }
    public int nonzeroCount() { return values.length; }
    public int[] rowStarts() { return rowStarts.clone(); }
    public int[] columnIndices() { return columnIndices.clone(); }
    public double[] values() { return values.clone(); }
}
