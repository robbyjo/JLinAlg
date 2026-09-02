/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;

/** Two-dimensional tensor-product P-spline with marginal penalties. */
public final class TensorProductPSplineTerm {
    private TensorProductPSplineTerm() { }

    /** Creates cubic marginal P-splines with second-difference penalties. */
    public static QuadraticSmoothTerm of(
            String name,
            double[] first,
            double[] second,
            int firstBasisDimension,
            int secondBasisDimension) {
        if (first == null || second == null || first.length != second.length) {
            throw new IllegalArgumentException(
                "tensor covariates must have equal non-null lengths");
        }
        PSplineTerm firstTerm = PSplineTerm.of(
            name + ".margin1", first, firstBasisDimension);
        PSplineTerm secondTerm = PSplineTerm.of(
            name + ".margin2", second, secondBasisDimension);
        int rows = first.length;
        int columns = firstBasisDimension * secondBasisDimension;
        double[] firstBasis = firstTerm.designView();
        double[] secondBasis = secondTerm.designView();
        double[] design = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            for (int firstColumn = 0;
                    firstColumn < firstBasisDimension; firstColumn++) {
                for (int secondColumn = 0;
                        secondColumn < secondBasisDimension; secondColumn++) {
                    int column = firstColumn * secondBasisDimension + secondColumn;
                    design[row * columns + column] =
                        firstBasis[row * firstBasisDimension + firstColumn]
                            * secondBasis[row * secondBasisDimension + secondColumn];
                }
            }
        }
        double[] firstPenalty = quadraticPenalty(firstTerm);
        double[] secondPenalty = quadraticPenalty(secondTerm);
        double[] tensorFirst = new double[columns * columns];
        double[] tensorSecond = new double[columns * columns];
        for (int firstRow = 0; firstRow < firstBasisDimension; firstRow++) {
            for (int secondRow = 0;
                    secondRow < secondBasisDimension; secondRow++) {
                int row = firstRow * secondBasisDimension + secondRow;
                for (int firstColumn = 0;
                        firstColumn < firstBasisDimension; firstColumn++) {
                    for (int secondColumn = 0;
                            secondColumn < secondBasisDimension; secondColumn++) {
                        int column = firstColumn * secondBasisDimension
                            + secondColumn;
                        if (secondRow == secondColumn) {
                            tensorFirst[row * columns + column] = firstPenalty[
                                firstRow * firstBasisDimension + firstColumn];
                        }
                        if (firstRow == firstColumn) {
                            tensorSecond[row * columns + column] = secondPenalty[
                                secondRow * secondBasisDimension + secondColumn];
                        }
                    }
                }
            }
        }
        return new QuadraticSmoothTerm(name, rows, columns, design,
            List.of(tensorFirst, tensorSecond));
    }

    private static double[] quadraticPenalty(PSplineTerm term) {
        int columns = term.basisDimension();
        int rows = columns - term.differenceOrder();
        double[] factor = term.penaltyFactorView();
        double[] result = new double[columns * columns];
        for (int first = 0; first < columns; first++) {
            for (int second = 0; second < columns; second++) {
                for (int row = 0; row < rows; row++) {
                    result[first * columns + second] +=
                        factor[row * columns + first]
                            * factor[row * columns + second];
                }
            }
        }
        return result;
    }
}
