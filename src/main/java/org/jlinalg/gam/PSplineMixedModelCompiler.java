/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import jdistlib.accelerator.ComputeBackend;

/** Shared P-spline null-space and random-effect model compilation. */
final class PSplineMixedModelCompiler {
    private PSplineMixedModelCompiler() { }

    static Compiled compile(
            double[] parametric,
            int rows,
            int parametricColumns,
            List<PSplineTerm> smoothTerms,
            ComputeBackend backend) {
        if (smoothTerms == null || smoothTerms.isEmpty()) {
            throw new IllegalArgumentException("at least one smooth term is required");
        }
        Set<String> names = new HashSet<>();
        List<Term> terms = new ArrayList<>(smoothTerms.size());
        int fixedColumns = parametricColumns;
        for (PSplineTerm smooth : smoothTerms) {
            if (smooth == null || smooth.observations() != rows
                    || !names.add(smooth.name())
                    || "residual".equals(smooth.name())) {
                throw new IllegalArgumentException(
                    "smooth terms must have unique names and matching rows");
            }
            Term term = decompose(smooth, backend);
            terms.add(term);
            fixedColumns += term.fixedColumns();
        }
        if (rows <= fixedColumns) {
            throw new IllegalArgumentException(
                "GAM requires more observations than unpenalized columns");
        }
        double[] fixed = new double[rows * fixedColumns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(parametric, row * parametricColumns,
                fixed, row * fixedColumns, parametricColumns);
        }
        int destination = parametricColumns;
        for (Term term : terms) {
            for (int row = 0; row < rows; row++) {
                System.arraycopy(term.fixedDesign(), row * term.fixedColumns(),
                    fixed, row * fixedColumns + destination, term.fixedColumns());
            }
            destination += term.fixedColumns();
        }
        return new Compiled(fixed, fixedColumns, List.copyOf(terms));
    }

    static double[] covarianceBasis(Term term) {
        int rows = term.term().observations();
        int columns = term.randomColumns();
        double[] design = term.randomDesign();
        double[] result = new double[rows * rows];
        for (int row = 0; row < rows; row++) {
            for (int other = 0; other <= row; other++) {
                double value = 0.0;
                for (int column = 0; column < columns; column++) {
                    value += design[row * columns + column]
                        * design[other * columns + column];
                }
                result[row * rows + other] = value;
                result[other * rows + row] = value;
            }
        }
        return result;
    }

    static double[] randomCoefficients(
            Term term, double[] projectedResidual, double variance) {
        int rows = term.term().observations();
        double[] result = new double[term.randomColumns()];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < term.randomColumns(); column++) {
                result[column] += variance
                    * term.randomDesign()[row * term.randomColumns() + column]
                    * projectedResidual[row];
            }
        }
        return result;
    }

    static double[] contribution(
            Term term,
            double[] fixedCoefficients,
            double[] randomCoefficients) {
        int rows = term.term().observations();
        double[] result = new double[rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < term.fixedColumns(); column++) {
                result[row] += term.fixedDesign()[row * term.fixedColumns() + column]
                    * fixedCoefficients[column];
            }
            for (int column = 0; column < term.randomColumns(); column++) {
                result[row] += term.randomDesign()[row * term.randomColumns() + column]
                    * randomCoefficients[column];
            }
        }
        return result;
    }

    static double randomEdf(Term term, double[] projection, double variance) {
        int rows = term.term().observations();
        double trace = 0.0;
        for (int column = 0; column < term.randomColumns(); column++) {
            for (int row = 0; row < rows; row++) {
                double left = term.randomDesign()[
                    row * term.randomColumns() + column];
                for (int other = 0; other < rows; other++) {
                    trace += left * projection[row * rows + other]
                        * term.randomDesign()[
                            other * term.randomColumns() + column];
                }
            }
        }
        return variance * trace;
    }

    private static Term decompose(
            PSplineTerm term, ComputeBackend backend) {
        int basisColumns = term.basisDimension();
        int penaltyRows = basisColumns - term.differenceOrder();
        double[] difference = term.penaltyFactorView();
        double[] gram = new double[penaltyRows * penaltyRows];
        for (int row = 0; row < penaltyRows; row++) {
            for (int column = 0; column <= row; column++) {
                double value = 0.0;
                for (int shared = 0; shared < basisColumns; shared++) {
                    value += difference[row * basisColumns + shared]
                        * difference[column * basisColumns + shared];
                }
                gram[row * penaltyRows + column] = value;
                gram[column * penaltyRows + row] = value;
            }
        }
        double[] solved = backend.dpotrf(gram, penaltyRows)
            .solve(difference, basisColumns);
        double[] randomTransform = new double[basisColumns * penaltyRows];
        for (int row = 0; row < basisColumns; row++) {
            for (int column = 0; column < penaltyRows; column++) {
                randomTransform[row * penaltyRows + column] =
                    solved[column * basisColumns + row];
            }
        }
        int nullColumns = term.differenceOrder();
        double[] fixedTransform = new double[basisColumns * nullColumns];
        for (int row = 0; row < basisColumns; row++) {
            double power = 1.0;
            for (int column = 0; column < nullColumns; column++) {
                fixedTransform[row * nullColumns + column] = power;
                power *= row;
            }
        }
        double[] rawFixed = multiply(term.designView(), term.observations(),
            basisColumns, fixedTransform, nullColumns);
        Centered centeredFixed = center(rawFixed, term.observations(), nullColumns);
        boolean[] keep = nonconstant(centeredFixed.matrix(),
            term.observations(), nullColumns);
        double[] reducedFixed = reduce(centeredFixed.matrix(),
            term.observations(), nullColumns, keep);
        double[] reducedTransform = reduce(fixedTransform,
            basisColumns, nullColumns, keep);
        double[] reducedMeans = reduce(centeredFixed.means(), keep);
        int kept = reducedMeans.length;
        double[] rawRandom = multiply(term.designView(), term.observations(),
            basisColumns, randomTransform, penaltyRows);
        Centered centeredRandom = center(rawRandom,
            term.observations(), penaltyRows);
        return new Term(term, reducedFixed, reducedTransform, reducedMeans, kept,
            centeredRandom.matrix(), randomTransform,
            centeredRandom.means(), penaltyRows);
    }

    private static Centered center(double[] matrix, int rows, int columns) {
        double[] result = matrix.clone();
        double[] means = new double[columns];
        for (int column = 0; column < columns; column++) {
            for (int row = 0; row < rows; row++) {
                means[column] += matrix[row * columns + column];
            }
            means[column] /= rows;
            for (int row = 0; row < rows; row++) {
                result[row * columns + column] -= means[column];
            }
        }
        return new Centered(result, means);
    }

    private static boolean[] nonconstant(
            double[] matrix, int rows, int columns) {
        boolean[] keep = new boolean[columns];
        double tolerance = 1e-12 * Math.sqrt(rows);
        for (int column = 0; column < columns; column++) {
            double norm = 0.0;
            for (int row = 0; row < rows; row++) {
                double value = matrix[row * columns + column];
                norm += value * value;
            }
            keep[column] = Math.sqrt(norm) > tolerance;
        }
        return keep;
    }

    private static double[] reduce(
            double[] matrix, int rows, int columns, boolean[] keep) {
        int kept = 0;
        for (boolean value : keep) if (value) kept++;
        double[] result = new double[rows * kept];
        for (int row = 0; row < rows; row++) {
            int destination = 0;
            for (int column = 0; column < columns; column++) {
                if (keep[column]) {
                    result[row * kept + destination++] =
                        matrix[row * columns + column];
                }
            }
        }
        return result;
    }

    private static double[] reduce(double[] vector, boolean[] keep) {
        int kept = 0;
        for (boolean value : keep) if (value) kept++;
        double[] result = new double[kept];
        int destination = 0;
        for (int index = 0; index < vector.length; index++) {
            if (keep[index]) result[destination++] = vector[index];
        }
        return result;
    }

    private static double[] multiply(
            double[] left, int rows, int shared,
            double[] right, int columns) {
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                for (int index = 0; index < shared; index++) {
                    result[row * columns + column] +=
                        left[row * shared + index]
                            * right[index * columns + column];
                }
            }
        }
        return result;
    }

    record Compiled(
            double[] fixedDesign,
            int fixedColumns,
            List<Term> terms) { }

    record Term(
            PSplineTerm term,
            double[] fixedDesign,
            double[] fixedTransform,
            double[] fixedMeans,
            int fixedColumns,
            double[] randomDesign,
            double[] randomTransform,
            double[] randomMeans,
            int randomColumns) { }

    private record Centered(double[] matrix, double[] means) { }
}
