/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.compute;

import java.util.Arrays;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.matrix.CsrMatrix;

/** Reusable SuiteSparse/CHOLMOD symbolic and numeric factorization. */
final class CholmodPreparedSparseCholesky implements PreparedSparseCholesky {
    private final int dimension;
    private final int[] rowStarts;
    private final int[] columnIndices;
    private final int[] permutation;
    private final int structuralNonzeros;
    private long handle;

    CholmodPreparedSparseCholesky(
            CsrMatrix matrix, MatrixTriangle triangle, SparseOrdering ordering) {
        if (matrix == null || triangle == null || ordering == null
                || matrix.rows() != matrix.columns())
            throw new IllegalArgumentException(
                "CHOLMOD requires a square matrix, triangle, and ordering");
        CholmodNative.requireAvailable();
        dimension = matrix.rows();
        rowStarts = matrix.rowStarts().clone();
        columnIndices = matrix.columnIndices().clone();
        structuralNonzeros = matrix.nonzeroCount();
        handle = CholmodNative.create(dimension, rowStarts, columnIndices,
            matrix.values(), triangle == MatrixTriangle.LOWER,
            ordering == SparseOrdering.NATURAL);
        if (handle == 0L) throw new IllegalStateException(
            "CHOLMOD did not return a factor handle");
        permutation = CholmodNative.permutation(handle);
    }

    @Override public int dimension() { return dimension; }
    @Override public int structuralNonzeroCount() { return structuralNonzeros; }

    @Override
    public synchronized int factorNonzeroCount() {
        requireOpen();
        return CholmodNative.factorNonzeroCount(handle);
    }

    @Override public int[] permutation() { return permutation.clone(); }

    @Override
    public synchronized double logDeterminant() {
        requireOpen();
        return CholmodNative.logDeterminant(handle);
    }

    @Override
    public synchronized void refactor(CsrMatrix matrix) {
        requireOpen();
        if (matrix == null || matrix.rows() != dimension
                || matrix.columns() != dimension
                || !Arrays.equals(rowStarts, matrix.rowStarts())
                || !Arrays.equals(columnIndices, matrix.columnIndices()))
            throw new IllegalArgumentException(
                "CHOLMOD refactor requires the original sparsity pattern");
        CholmodNative.refactor(handle, matrix.values());
    }

    @Override
    public synchronized void solveInPlace(
            double[] rightHandSide, int rightHandSides) {
        requireOpen();
        if (rightHandSide == null || rightHandSides < 1
                || rightHandSide.length != dimension * rightHandSides)
            throw new IllegalArgumentException(
                "CHOLMOD right-hand side dimensions are invalid");
        CholmodNative.solveInPlace(handle, rightHandSide, rightHandSides);
    }

    @Override
    public synchronized void close() {
        if (handle == 0L) return;
        CholmodNative.destroy(handle);
        handle = 0L;
    }

    private void requireOpen() {
        if (handle == 0L) throw new IllegalStateException(
            "CHOLMOD factor is closed");
    }
}
