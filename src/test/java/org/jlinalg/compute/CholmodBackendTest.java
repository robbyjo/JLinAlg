/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.compute;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.matrix.CsrMatrix;
import org.junit.jupiter.api.Test;

final class CholmodBackendTest {
    @Test
    void factorsRefactorsAndSolvesOneBasedCsr() {
        BackendContext selected;
        try {
            selected = BackendContext.select(BackendPolicy.CHOLMOD);
        } catch (IllegalStateException | LinkageError unavailable) {
            assumeTrue(false, unavailable.getMessage());
            return;
        }
        try (BackendContext context = selected;
                PreparedSparseCholesky factor =
                    context.backend().prepareDcsrpotrf(
                        matrix(new double[] {4, 1, 3, 1, 2}),
                        MatrixTriangle.LOWER,
                        SparseOrdering.MINIMUM_DEGREE)) {
            assertTrue(context.provenance().selectedBackend()
                .startsWith("cholmod+"));
            assertEquals(Math.log(18.0), factor.logDeterminant(), 1e-12);
            double[] right = {1, 2, 3, 3, 2, 1};
            factor.solveInPlace(right, 2);
            assertSolutions(matrix(new double[] {4, 1, 3, 1, 2}),
                right, new double[] {1, 2, 3, 3, 2, 1}, 2);

            CsrMatrix refactored = matrix(new double[] {5, 1, 4, 1, 3});
            factor.refactor(refactored);
            assertEquals(Math.log(52.0), factor.logDeterminant(), 1e-12);
            double[] second = {2, -1, 4};
            factor.solveInPlace(second, 1);
            assertSolutions(refactored, second,
                new double[] {2, -1, 4}, 1);
        }
    }

    @Test
    void preferredSelectsCholmodWhenNativeLibraryIsPresent() {
        try (BackendContext context = BackendContext.preferred()) {
            if (!context.provenance().selectedBackend().startsWith("cholmod+"))
                assumeTrue(false, "CHOLMOD native library is not present");
            assertEquals(BackendPolicy.PREFERRED,
                context.provenance().requested());
        }
    }

    @Test
    void acceptsUpperTriangleCsr() {
        BackendContext selected;
        try {
            selected = BackendContext.select(BackendPolicy.CHOLMOD);
        } catch (IllegalStateException | LinkageError unavailable) {
            assumeTrue(false, unavailable.getMessage());
            return;
        }
        CsrMatrix upper = new CsrMatrix(3, 3,
            new double[] {4, 1, 3, 1, 2},
            new int[] {1, 2, 2, 3, 3},
            new int[] {1, 3, 5, 6});
        try (BackendContext context = selected;
                PreparedSparseCholesky factor =
                    context.backend().prepareDcsrpotrf(
                        upper, MatrixTriangle.UPPER,
                        SparseOrdering.NATURAL)) {
            assertEquals(Math.log(18.0), factor.logDeterminant(), 1e-12);
            double[] right = {1, 2, 3};
            factor.solveInPlace(right, 1);
            assertSolutions(upper, right, new double[] {1, 2, 3}, 1);
        }
    }

    private static CsrMatrix matrix(double[] values) {
        return new CsrMatrix(3, 3, values,
            new int[] {1, 1, 2, 2, 3},
            new int[] {1, 2, 4, 6});
    }

    private static void assertSolutions(
            CsrMatrix matrix, double[] solutions, double[] expected,
            int rightHandSides) {
        double[] dense = matrix.toDense();
        for (int row = 0; row < 3; row++) {
            for (int rhs = 0; rhs < rightHandSides; rhs++) {
                double actual = 0.0;
                for (int column = 0; column < 3; column++) {
                    double symmetric = dense[row * 3 + column];
                    if (symmetric == 0.0)
                        symmetric = dense[column * 3 + row];
                    actual += symmetric
                        * solutions[column * rightHandSides + rhs];
                }
                assertEquals(expected[row * rightHandSides + rhs],
                    actual, 1e-10);
            }
        }
    }
}
