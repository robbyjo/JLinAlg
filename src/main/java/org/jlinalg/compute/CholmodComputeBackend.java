/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.compute;

import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.ComputeCapabilities;
import jdistlib.accelerator.LuFactor;
import jdistlib.accelerator.MatrixTranspose;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PivotedQrFactor;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SingularValueDecomposition;
import jdistlib.accelerator.SparseCholeskyFactor;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.accelerator.SymmetricEigenDecomposition;
import jdistlib.accelerator.SymmetricIndefiniteFactor;
import jdistlib.matrix.CsrMatrix;

/** CHOLMOD sparse factorization with a JDistlib dense-operation delegate. */
final class CholmodComputeBackend implements ComputeBackend {
    private final ComputeBackend delegate;
    private final ComputeCapabilities capabilities;

    CholmodComputeBackend(ComputeBackend delegate) {
        if (delegate == null) throw new IllegalArgumentException(
            "CHOLMOD requires a dense backend delegate");
        CholmodNative.requireAvailable();
        this.delegate = delegate;
        ComputeCapabilities base = delegate.capabilities();
        capabilities = new ComputeCapabilities(
            id(), "CHOLMOD sparse + " + base.device(), true, false,
            base.globalMemoryBytes(), base.denseLinearAlgebra(), true,
            base.nativeFactorizations(), true, true, true,
            base.preparedDenseMatrices(), base.batchedLinearAlgebra());
    }

    @Override public String id() { return "cholmod+" + delegate.id(); }
    @Override public String selectedBackend() { return id(); }
    @Override public boolean automaticRouting() {
        return delegate.automaticRouting();
    }
    @Override public boolean available() { return CholmodNative.available(); }
    @Override public ComputeCapabilities capabilities() { return capabilities; }

    @Override
    public PreparedSparseCholesky prepareDcsrpotrf(
            CsrMatrix matrix, MatrixTriangle triangle, SparseOrdering ordering) {
        return new CholmodPreparedSparseCholesky(matrix, triangle, ordering);
    }

    @Override
    public SparseCholeskyFactor dcsrpotrf(
            CsrMatrix matrix, MatrixTriangle triangle, SparseOrdering ordering) {
        return delegate.dcsrpotrf(matrix, triangle, ordering);
    }

    @Override public double ddot(int n, double[] x, int xOffset, int xStride,
            double[] y, int yOffset, int yStride) {
        return delegate.ddot(n, x, xOffset, xStride, y, yOffset, yStride);
    }
    @Override public void dgemv(MatrixTranspose transpose, int rows, int columns,
            double alpha, double[] matrix, double[] x, double beta, double[] y) {
        delegate.dgemv(transpose, rows, columns, alpha, matrix, x, beta, y);
    }
    @Override public void dgemv(MatrixTranspose transpose, int rows, int columns,
            double alpha, double[] matrix, int matrixOffset, int leadingDimension,
            double[] x, int xOffset, int xStride, double beta, double[] y,
            int yOffset, int yStride) {
        delegate.dgemv(transpose, rows, columns, alpha, matrix, matrixOffset,
            leadingDimension, x, xOffset, xStride, beta, y, yOffset, yStride);
    }
    @Override public void dgemm(MatrixTranspose leftTranspose,
            MatrixTranspose rightTranspose, int rows, int columns, int shared,
            double alpha, double[] left, double[] right, double beta,
            double[] result) {
        delegate.dgemm(leftTranspose, rightTranspose, rows, columns, shared,
            alpha, left, right, beta, result);
    }
    @Override public void dgemm(MatrixTranspose leftTranspose,
            MatrixTranspose rightTranspose, int rows, int columns, int shared,
            double alpha, double[] left, int leftOffset, int leftLeading,
            double[] right, int rightOffset, int rightLeading, double beta,
            double[] result, int resultOffset, int resultLeading) {
        delegate.dgemm(leftTranspose, rightTranspose, rows, columns, shared,
            alpha, left, leftOffset, leftLeading, right, rightOffset,
            rightLeading, beta, result, resultOffset, resultLeading);
    }
    @Override public CholeskyFactor dpotrf(double[] matrix, int dimension) {
        return delegate.dpotrf(matrix, dimension);
    }
    @Override public LuFactor dgetrf(double[] matrix, int dimension) {
        return delegate.dgetrf(matrix, dimension);
    }
    @Override public SymmetricIndefiniteFactor dsytrf(
            double[] matrix, int dimension) {
        return delegate.dsytrf(matrix, dimension);
    }
    @Override public PivotedQrFactor dgeqp3(
            double[] matrix, int rows, int columns) {
        return delegate.dgeqp3(matrix, rows, columns);
    }
    @Override public SymmetricEigenDecomposition dsyev(
            double[] matrix, int dimension) {
        return delegate.dsyev(matrix, dimension);
    }
    @Override public SingularValueDecomposition dgesvd(
            double[] matrix, int rows, int columns) {
        return delegate.dgesvd(matrix, rows, columns);
    }
    @Override public void close() { }
}
