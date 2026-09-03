/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.LinearMixedModel;
import org.jlinalg.mixed.LinearMixedModelResult;
import org.jlinalg.mixed.RandomEffectEstimates;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.reml.RemlOptions;

/** Generalized additive model fitters. */
public final class Gam {
    private Gam() { }

    /** Fits an exact Gaussian GAM with REML-selected smoothing parameters. */
    public static GamResult fitGaussian(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms) {
        return fitGaussian(response, parametricDesign, smoothTerms,
            RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits an exact Gaussian GAM with explicit REML and backend controls. */
    public static GamResult fitGaussian(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(parametricDesign, response.length);
        return fitGaussian(response, rowMajor, response.length,
            parametricDesign[0].length, smoothTerms, options, backendPolicy);
    }

    /** Contiguous row-major overload for allocation-sensitive callers. */
    public static GamResult fitGaussian(
            double[] response,
            double[] parametricDesign,
            int rows,
            int parametricColumns,
            List<PSplineTerm> smoothTerms,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (backendPolicy == null)
            throw new IllegalArgumentException("backendPolicy is required");
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return fitGaussian(response, parametricDesign, rows,
                parametricColumns, smoothTerms, options, context);
        }
    }

    static GamResult fitGaussian(
            double[] response,
            double[] parametricDesign,
            int rows,
            int parametricColumns,
            List<PSplineTerm> smoothTerms,
            RemlOptions options,
            BackendContext context) {
        MatrixOps.validateModelData(
            response, parametricDesign, rows, parametricColumns);
        validateSmoothTerms(smoothTerms, rows);
        if (options == null || context == null) {
            throw new IllegalArgumentException(
                "options and backend context are required");
        }

        List<DecomposedTerm> decomposed = new ArrayList<>(smoothTerms.size());
        for (PSplineTerm term : smoothTerms) {
            decomposed.add(decompose(term, context.backend()));
        }

        int fixedColumns = parametricColumns;
        for (DecomposedTerm term : decomposed) {
            fixedColumns += term.fixedColumns();
        }
        if (rows <= fixedColumns) {
            throw new IllegalArgumentException(
                "GAM requires more observations than unpenalized columns");
        }
        double[] fixedDesign = new double[rows * fixedColumns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(parametricDesign, row * parametricColumns,
                fixedDesign, row * fixedColumns, parametricColumns);
        }

        int destination = parametricColumns;
        List<RandomEffectTerm> randomTerms = new ArrayList<>(decomposed.size());
        for (DecomposedTerm term : decomposed) {
            for (int row = 0; row < rows; row++) {
                System.arraycopy(term.fixedDesign(), row * term.fixedColumns(),
                    fixedDesign, row * fixedColumns + destination,
                    term.fixedColumns());
            }
            destination += term.fixedColumns();
            List<String> coefficientNames = new ArrayList<>(term.randomColumns());
            for (int column = 0; column < term.randomColumns(); column++) {
                coefficientNames.add(term.term().name() + ".pen" + (column + 1));
            }
            randomTerms.add(RandomEffectTerm.of(
                term.term().name(), term.randomDesign(), rows,
                term.randomColumns(), coefficientNames));
        }

        LinearMixedModelResult mixed;
        if (options.degreesOfFreedomMethod()
                == DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION) {
            SparseLinearMixedModelResult sparse =
                SparseLinearMixedModel.fitWithBackend(
                response, fixedDesign, rows, fixedColumns,
                randomTerms, options, context.backend(), context.provenance());
            mixed = LinearMixedModelResult.fromSparse(
                sparse, response, fixedDesign, rows, fixedColumns);
        } else {
            mixed = LinearMixedModel.fit(response, fixedDesign, rows,
                fixedColumns, randomTerms, options,
                context.provenance().requested());
        }
        double[] allFixed = mixed.beta();
        double[] componentVariances = mixed.reml().varianceComponents();
        double residualVariance = componentVariances[componentVariances.length - 1];
        List<SmoothTermEstimate> estimates = new ArrayList<>(decomposed.size());
        double totalEdf = parametricColumns;
        int fixedStart = parametricColumns;
        for (int index = 0; index < decomposed.size(); index++) {
            DecomposedTerm term = decomposed.get(index);
            double[] fixedCoefficients = java.util.Arrays.copyOfRange(
                allFixed, fixedStart, fixedStart + term.fixedColumns());
            fixedStart += term.fixedColumns();
            RandomEffectEstimates random = mixed.randomEffects(term.term().name());
            double[] randomCoefficients = random.estimates();
            double termVariance = componentVariances[index];
            double edf = term.fixedColumns();
            for (double predictionErrorVariance
                    : random.predictionErrorVariances()) {
                edf += 1.0 - predictionErrorVariance / termVariance;
            }
            edf = Math.max(term.fixedColumns(),
                Math.min(term.fixedColumns() + term.randomColumns(), edf));
            double[] fitted = contribution(
                term.fixedDesign(), rows, term.fixedColumns(), fixedCoefficients,
                term.randomDesign(), term.randomColumns(), randomCoefficients);
            estimates.add(new SmoothTermEstimate(
                term.term(), term.fixedTransform(), term.fixedMeans(),
                fixedCoefficients, term.randomTransform(), term.randomMeans(),
                randomCoefficients, fitted,
                residualVariance / termVariance, edf));
            totalEdf += edf;
        }
        return new GamResult(mixed, parametricColumns, estimates, totalEdf);
    }

    private static DecomposedTerm decompose(
            PSplineTerm term, ComputeBackend backend) {
        int basisColumns = term.basisDimension();
        int penaltyRows = basisColumns - term.differenceOrder();
        double[] difference = term.penaltyFactorView();
        double[] differenceGram = new double[penaltyRows * penaltyRows];
        for (int row = 0; row < penaltyRows; row++) {
            for (int column = 0; column <= row; column++) {
                double value = 0.0;
                for (int shared = 0; shared < basisColumns; shared++) {
                    value += difference[row * basisColumns + shared]
                        * difference[column * basisColumns + shared];
                }
                differenceGram[row * penaltyRows + column] = value;
                differenceGram[column * penaltyRows + row] = value;
            }
        }
        CholeskyFactor factor = backend.dpotrf(differenceGram, penaltyRows);
        double[] inverseTimesDifference = factor.solve(difference, basisColumns);
        double[] randomTransform = new double[basisColumns * penaltyRows];
        for (int row = 0; row < basisColumns; row++) {
            for (int column = 0; column < penaltyRows; column++) {
                randomTransform[row * penaltyRows + column] =
                    inverseTimesDifference[column * basisColumns + row];
            }
        }

        int nullColumns = term.differenceOrder();
        double[] candidateFixedTransform = new double[basisColumns * nullColumns];
        for (int row = 0; row < basisColumns; row++) {
            double power = 1.0;
            for (int column = 0; column < nullColumns; column++) {
                candidateFixedTransform[row * nullColumns + column] = power;
                power *= row;
            }
        }
        double[] candidateFixed = multiply(
            term.designView(), term.observations(), basisColumns,
            candidateFixedTransform, nullColumns);
        CenteredMatrix centeredFixed = centerAndDropConstant(
            candidateFixed, term.observations(), nullColumns,
            candidateFixedTransform, basisColumns);
        double[] randomDesign = multiply(
            term.designView(), term.observations(), basisColumns,
            randomTransform, penaltyRows);
        CenteredMatrix centeredRandom = center(
            randomDesign, term.observations(), penaltyRows);
        return new DecomposedTerm(term,
            centeredFixed.matrix(), centeredFixed.transform(), centeredFixed.means(),
            centeredFixed.columns(), centeredRandom.matrix(), randomTransform,
            centeredRandom.means(), penaltyRows);
    }

    private static CenteredMatrix centerAndDropConstant(
            double[] matrix,
            int rows,
            int columns,
            double[] transform,
            int transformRows) {
        CenteredMatrix centered = center(matrix, rows, columns);
        boolean[] keep = new boolean[columns];
        int kept = 0;
        double tolerance = 1e-12 * Math.sqrt(rows);
        for (int column = 0; column < columns; column++) {
            double norm = 0.0;
            for (int row = 0; row < rows; row++) {
                double value = centered.matrix()[row * columns + column];
                norm += value * value;
            }
            keep[column] = Math.sqrt(norm) > tolerance;
            if (keep[column]) kept++;
        }
        double[] reduced = new double[rows * kept];
        double[] reducedTransform = new double[transformRows * kept];
        double[] reducedMeans = new double[kept];
        int destination = 0;
        for (int column = 0; column < columns; column++) {
            if (!keep[column]) continue;
            for (int row = 0; row < rows; row++) {
                reduced[row * kept + destination] =
                    centered.matrix()[row * columns + column];
            }
            for (int row = 0; row < transformRows; row++) {
                reducedTransform[row * kept + destination] =
                    transform[row * columns + column];
            }
            reducedMeans[destination] = centered.means()[column];
            destination++;
        }
        return new CenteredMatrix(
            reduced, reducedTransform, reducedMeans, kept);
    }

    private static CenteredMatrix center(
            double[] matrix, int rows, int columns) {
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
        return new CenteredMatrix(result, new double[0], means, columns);
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

    private static double[] contribution(
            double[] fixed,
            int rows,
            int fixedColumns,
            double[] fixedCoefficients,
            double[] random,
            int randomColumns,
            double[] randomCoefficients) {
        double[] result = new double[rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < fixedColumns; column++) {
                result[row] += fixed[row * fixedColumns + column]
                    * fixedCoefficients[column];
            }
            for (int column = 0; column < randomColumns; column++) {
                result[row] += random[row * randomColumns + column]
                    * randomCoefficients[column];
            }
        }
        return result;
    }

    private static void validateSmoothTerms(
            List<PSplineTerm> terms, int rows) {
        if (terms == null || terms.isEmpty()) {
            throw new IllegalArgumentException("at least one smooth term is required");
        }
        Set<String> names = new HashSet<>();
        for (PSplineTerm term : terms) {
            if (term == null || term.observations() != rows) {
                throw new IllegalArgumentException(
                    "smooth terms must match the response length");
            }
            if (!names.add(term.name()) || "residual".equals(term.name())) {
                throw new IllegalArgumentException(
                    "smooth term names must be unique and not 'residual'");
            }
        }
    }

    private record CenteredMatrix(
            double[] matrix,
            double[] transform,
            double[] means,
            int columns) {
    }

    private record DecomposedTerm(
            PSplineTerm term,
            double[] fixedDesign,
            double[] fixedTransform,
            double[] fixedMeans,
            int fixedColumns,
            double[] randomDesign,
            double[] randomTransform,
            double[] randomMeans,
            int randomColumns) {
    }
}
