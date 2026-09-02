/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.SymmetricEigenDecomposition;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Compiles arbitrary multi-penalty smooths into a whitened predictor. */
public final class QuadraticPenalizedPredictor {
    private QuadraticPenalizedPredictor() { }

    /**
     * Combines each term's penalties using its supplied smoothing parameters,
     * moves the joint null space into the fixed design, and whitens the
     * penalized space to a unit diagonal penalty.
     */
    public static PenalizedPredictor compile(
            double[][] parametricDesign,
            List<QuadraticSmoothTerm> smoothTerms,
            List<double[]> smoothingParameters,
            BackendPolicy backendPolicy) {
        if (parametricDesign == null || parametricDesign.length == 0
                || parametricDesign[0] == null || smoothTerms == null
                || smoothTerms.isEmpty() || smoothingParameters == null
                || smoothingParameters.size() != smoothTerms.size()
                || backendPolicy == null) {
            throw new IllegalArgumentException(
                "parametric design, smooths, smoothing parameters, and backend are required");
        }
        int rows = parametricDesign.length;
        int parametricColumns = parametricDesign[0].length;
        double[] parametric = MatrixOps.rowMajor(parametricDesign, rows);
        List<Decomposed> decomposed = new ArrayList<>(smoothTerms.size());
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            for (int index = 0; index < smoothTerms.size(); index++) {
                QuadraticSmoothTerm term = smoothTerms.get(index);
                if (term == null || term.observations() != rows) {
                    throw new IllegalArgumentException(
                        "smooth terms must match predictor observations");
                }
                double[] smoothing = smoothingParameters.get(index);
                if (smoothing == null
                        || smoothing.length != term.penaltyCount()) {
                    throw new IllegalArgumentException(
                        "one smoothing parameter is required per term penalty");
                }
                decomposed.add(decompose(
                    term, smoothing, context.backend()));
            }
        }
        int candidateFixedColumns = parametricColumns;
        int penalizedColumns = 0;
        for (Decomposed term : decomposed) {
            candidateFixedColumns += term.nullColumns();
            penalizedColumns += term.penalizedColumns();
        }
        double[] candidateFixed = candidateFixed(parametric, rows,
            parametricColumns, decomposed, candidateFixedColumns);
        IndependentColumns independent = independentColumns(candidateFixed,
            rows, candidateFixedColumns, parametricColumns);
        int fixedColumns = independent.columns();
        int columns = fixedColumns + penalizedColumns;
        if (rows <= fixedColumns) {
            throw new IllegalArgumentException(
                "predictor has too many unpenalized columns");
        }
        double[] design = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(independent.design(), row * fixedColumns,
                design, row * columns, fixedColumns);
        }
        int penalizedDestination = fixedColumns;
        for (Decomposed term : decomposed) {
            for (int row = 0; row < rows; row++) {
                System.arraycopy(term.penalizedDesign(),
                    row * term.penalizedColumns(),
                    design, row * columns + penalizedDestination,
                    term.penalizedColumns());
            }
            penalizedDestination += term.penalizedColumns();
        }
        double[] penalty = new double[columns];
        java.util.Arrays.fill(penalty, fixedColumns, columns, 1.0);
        return new PenalizedPredictor(
            design, penalty, rows, columns, parametricColumns);
    }

    private static double[] candidateFixed(
            double[] parametric,
            int rows,
            int parametricColumns,
            List<Decomposed> terms,
            int columns) {
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(parametric, row * parametricColumns,
                result, row * columns, parametricColumns);
        }
        int destination = parametricColumns;
        for (Decomposed term : terms) {
            for (int row = 0; row < rows; row++) {
                System.arraycopy(term.nullDesign(), row * term.nullColumns(),
                    result, row * columns + destination, term.nullColumns());
            }
            destination += term.nullColumns();
        }
        return result;
    }

    private static IndependentColumns independentColumns(
            double[] matrix,
            int rows,
            int columns,
            int requiredColumns) {
        List<double[]> accepted = new ArrayList<>();
        boolean[] keep = new boolean[columns];
        double maximumNorm = 0.0;
        for (int column = 0; column < columns; column++) {
            double[] residual = new double[rows];
            for (int row = 0; row < rows; row++) {
                residual[row] = matrix[row * columns + column];
            }
            maximumNorm = Math.max(maximumNorm, norm(residual));
            orthogonalize(residual, accepted);
            orthogonalize(residual, accepted);
            double residualNorm = norm(residual);
            boolean independent = residualNorm > 1e-10
                * Math.max(1.0, maximumNorm);
            if (!independent && column < requiredColumns) {
                throw new IllegalArgumentException(
                    "parametric design is rank deficient at column " + column);
            }
            if (independent) {
                keep[column] = true;
                for (int row = 0; row < rows; row++) {
                    residual[row] /= residualNorm;
                }
                accepted.add(residual);
            }
        }
        int kept = count(keep);
        return new IndependentColumns(reduce(matrix, rows, columns, keep), kept);
    }

    private static void orthogonalize(
            double[] residual, List<double[]> accepted) {
        for (double[] basis : accepted) {
            double projection = dot(residual, basis);
            for (int row = 0; row < residual.length; row++) {
                residual[row] -= projection * basis[row];
            }
        }
    }

    private static double dot(double[] first, double[] second) {
        double result = 0.0;
        for (int index = 0; index < first.length; index++) {
            result += first[index] * second[index];
        }
        return result;
    }

    private static double norm(double[] values) {
        return Math.sqrt(dot(values, values));
    }

    private record IndependentColumns(double[] design, int columns) { }

    private static Decomposed decompose(
            QuadraticSmoothTerm term,
            double[] smoothing,
            ComputeBackend backend) {
        int columns = term.columns();
        double[] combined = new double[columns * columns];
        for (int penalty = 0; penalty < smoothing.length; penalty++) {
            double scale = smoothing[penalty];
            if (!(scale > 0.0) || !Double.isFinite(scale)) {
                throw new IllegalArgumentException(
                    "smoothing parameters must be finite and positive");
            }
            double[] matrix = term.penaltyViews().get(penalty);
            for (int index = 0; index < combined.length; index++) {
                combined[index] += scale * matrix[index];
            }
        }
        SymmetricEigenDecomposition eigen = backend.dsyev(combined, columns);
        double[] values = eigen.eigenvalues();
        double[] vectors = eigen.eigenvectors();
        double maximum = 0.0;
        for (double value : values) maximum = Math.max(maximum, Math.abs(value));
        double tolerance = 1e-10 * Math.max(1.0, maximum);
        int nullColumns = 0;
        for (double value : values) {
            if (value < -tolerance) {
                throw new IllegalArgumentException(
                    "combined smooth penalty is not positive semidefinite");
            }
            if (value <= tolerance) nullColumns++;
        }
        int penalizedColumns = columns - nullColumns;
        double[] rawNull = new double[term.observations() * nullColumns];
        double[] rawPenalized =
            new double[term.observations() * penalizedColumns];
        int nullDestination = 0;
        int penalizedDestination = 0;
        for (int component = 0; component < columns; component++) {
            boolean unpenalized = values[component] <= tolerance;
            double scale = unpenalized ? 1.0
                : 1.0 / Math.sqrt(values[component]);
            int destination = unpenalized
                ? nullDestination++ : penalizedDestination++;
            int destinationColumns = unpenalized
                ? nullColumns : penalizedColumns;
            double[] target = unpenalized ? rawNull : rawPenalized;
            for (int row = 0; row < term.observations(); row++) {
                double value = 0.0;
                for (int basis = 0; basis < columns; basis++) {
                    value += term.designView()[row * columns + basis]
                        * vectors[basis * columns + component];
                }
                target[row * destinationColumns + destination] = value * scale;
            }
        }
        double[] centeredNull = center(rawNull,
            term.observations(), nullColumns);
        boolean[] keep = nonconstant(centeredNull,
            term.observations(), nullColumns);
        double[] reducedNull = reduce(centeredNull,
            term.observations(), nullColumns, keep);
        double[] centeredPenalized = center(rawPenalized,
            term.observations(), penalizedColumns);
        return new Decomposed(reducedNull, count(keep),
            centeredPenalized, penalizedColumns);
    }

    private static double[] center(double[] matrix, int rows, int columns) {
        double[] result = matrix.clone();
        for (int column = 0; column < columns; column++) {
            double mean = 0.0;
            for (int row = 0; row < rows; row++) {
                mean += matrix[row * columns + column];
            }
            mean /= rows;
            for (int row = 0; row < rows; row++) {
                result[row * columns + column] -= mean;
            }
        }
        return result;
    }

    private static boolean[] nonconstant(
            double[] matrix, int rows, int columns) {
        boolean[] keep = new boolean[columns];
        for (int column = 0; column < columns; column++) {
            double norm = 0.0;
            for (int row = 0; row < rows; row++) {
                double value = matrix[row * columns + column];
                norm += value * value;
            }
            keep[column] = Math.sqrt(norm) > 1e-12 * Math.sqrt(rows);
        }
        return keep;
    }

    private static int count(boolean[] values) {
        int result = 0;
        for (boolean value : values) if (value) result++;
        return result;
    }

    private static double[] reduce(
            double[] matrix, int rows, int columns, boolean[] keep) {
        int kept = count(keep);
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

    private record Decomposed(
            double[] nullDesign,
            int nullColumns,
            double[] penalizedDesign,
            int penalizedColumns) { }
}
