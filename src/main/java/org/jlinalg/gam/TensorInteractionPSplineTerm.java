/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

/** Interaction-only tensor P-spline constrained against both marginal trends. */
public final class TensorInteractionPSplineTerm {
    private TensorInteractionPSplineTerm() { }

    /** Creates a {@code ti(x,z)} basis with two marginal penalties. */
    public static QuadraticSmoothTerm of(
            String name,
            double[] first,
            double[] second,
            int firstBasisDimension,
            int secondBasisDimension) {
        QuadraticSmoothTerm tensor = TensorProductPSplineTerm.of(name,
            first, second, firstBasisDimension, secondBasisDimension);
        double[] design = tensor.design();
        int rows = tensor.observations();
        int columns = tensor.columns();
        double[][] constraints = orthonormalConstraints(first, second);
        for (int column = 0; column < columns; column++) {
            for (double[] constraint : constraints) {
                double projection = 0.0;
                for (int row = 0; row < rows; row++) {
                    projection += constraint[row]
                        * design[row * columns + column];
                }
                for (int row = 0; row < rows; row++) {
                    design[row * columns + column] -=
                        projection * constraint[row];
                }
            }
        }
        return new QuadraticSmoothTerm(name, rows, columns,
            design, tensor.penalties());
    }

    private static double[][] orthonormalConstraints(
            double[] first, double[] second) {
        if (first == null || second == null || first.length != second.length
                || first.length < 3) {
            throw new IllegalArgumentException(
                "interaction covariates must have equal lengths of at least three");
        }
        double[][] candidates = new double[3][first.length];
        java.util.Arrays.fill(candidates[0], 1.0);
        System.arraycopy(first, 0, candidates[1], 0, first.length);
        System.arraycopy(second, 0, candidates[2], 0, second.length);
        java.util.ArrayList<double[]> accepted = new java.util.ArrayList<>();
        for (double[] candidate : candidates) {
            for (double[] basis : accepted) {
                double projection = dot(candidate, basis);
                for (int row = 0; row < candidate.length; row++) {
                    candidate[row] -= projection * basis[row];
                }
            }
            double norm = Math.sqrt(dot(candidate, candidate));
            if (norm > 1e-12 * Math.sqrt(candidate.length)) {
                for (int row = 0; row < candidate.length; row++) {
                    candidate[row] /= norm;
                }
                accepted.add(candidate);
            }
        }
        return accepted.toArray(double[][]::new);
    }

    private static double dot(double[] first, double[] second) {
        double result = 0.0;
        for (int index = 0; index < first.length; index++) {
            result += first[index] * second[index];
        }
        return result;
    }
}
