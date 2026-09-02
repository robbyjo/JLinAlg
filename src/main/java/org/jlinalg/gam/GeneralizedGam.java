/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmPql;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.glmm.GlmmPqlResult;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.VarianceComponent;

/** Generalized additive models fitted by PQL with REML smoothing updates. */
public final class GeneralizedGam {
    private static final double MINIMUM_WORKING_WEIGHT = 1e-12;
    private static final double MAXIMUM_WORKING_WEIGHT = 1e150;

    private GeneralizedGam() { }

    /** Fits with unit prior weights, zero offset, and default controls. */
    public static GeneralizedGamResult fit(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            GlmFamily family) {
        return fit(response, parametricDesign, smoothTerms, family,
            null, null, GlmmPqlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a fixed-dispersion generalized GAM with explicit controls. */
    public static GeneralizedGamResult fit(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(parametricDesign, response.length);
        return fit(response, rowMajor, response.length,
            parametricDesign[0].length, smoothTerms, family,
            priorWeights, offset, options, backendPolicy);
    }

    /** Contiguous row-major overload for allocation-sensitive callers. */
    public static GeneralizedGamResult fit(
            double[] response,
            double[] parametricDesign,
            int rows,
            int parametricColumns,
            List<PSplineTerm> smoothTerms,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(
            response, parametricDesign, rows, parametricColumns);
        validateSmoothTerms(smoothTerms, rows);
        if (family == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "family, options, and backendPolicy are required");
        }
        double[] weights = weights(priorWeights, rows);
        double[] offsets = offset == null ? new double[rows]
            : MatrixOps.finiteCopy(offset, "offset");
        if (offsets.length != rows) {
            throw new IllegalArgumentException("offset length must equal rows");
        }

        List<DecomposedTerm> terms = new ArrayList<>(smoothTerms.size());
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            for (PSplineTerm term : smoothTerms) {
                terms.add(decompose(term, context.backend()));
            }
        }
        ModelMatrices model = assemble(parametricDesign, rows,
            parametricColumns, terms);
        List<VarianceComponent> components = new ArrayList<>(terms.size());
        for (DecomposedTerm term : terms) {
            components.add(new VarianceComponent(term.term().name(), rows,
                covarianceBasis(term.randomDesign(), rows,
                    term.randomColumns())));
        }
        GlmmPqlResult fitted = GlmmPql.fit(response, model.fixedDesign(),
            rows, model.fixedColumns(), family, components, weights, offsets,
            options, backendPolicy);
        return estimates(response, rows, parametricColumns, family,
            weights, offsets, terms, model, components, fitted, backendPolicy);
    }

    private static GeneralizedGamResult estimates(
            double[] response,
            int rows,
            int parametricColumns,
            GlmFamily family,
            double[] weights,
            double[] offsets,
            List<DecomposedTerm> terms,
            ModelMatrices model,
            List<VarianceComponent> components,
            GlmmPqlResult fitted,
            BackendPolicy backendPolicy) {
        double[] predictor = fitted.linearPredictor();
        double[] means = fitted.fittedMeans();
        double[] workingResponse = new double[rows];
        double[] covariance = new double[rows * rows];
        double[] variances = fitted.varianceComponents();
        for (int row = 0; row < rows; row++) {
            double derivative = family.meanDerivative(predictor[row]);
            double variance = family.variance(means[row]);
            double workingWeight = clamp(weights[row]
                * derivative * derivative / variance,
                MINIMUM_WORKING_WEIGHT, MAXIMUM_WORKING_WEIGHT);
            workingResponse[row] = predictor[row]
                + (response[row] - means[row]) / derivative - offsets[row];
            covariance[row * rows + row] = 1.0 / workingWeight;
        }
        for (int component = 0; component < components.size(); component++) {
            double[] basis = components.get(component).covariance();
            for (int index = 0; index < covariance.length; index++) {
                covariance[index] += variances[component] * basis[index];
            }
        }

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            CholeskyFactor factor = backend.dpotrf(covariance, rows);
            double[] beta = fitted.beta();
            double[] fixedPredictor = MatrixOps.multiply(backend,
                model.fixedDesign(), rows, model.fixedColumns(), beta);
            double[] workingResidual = MatrixOps.subtract(
                workingResponse, fixedPredictor);
            double[] projected = factor.solve(workingResidual);

            double[] inverseCovariance = factor.solve(
                MatrixOps.identity(rows), rows);
            double[] inverseCovarianceFixed = factor.solve(
                model.fixedDesign(), model.fixedColumns());
            double[] temporary = MatrixOps.multiply(backend,
                inverseCovarianceFixed, rows, model.fixedColumns(),
                fitted.fixedEffectCovariance(), model.fixedColumns());
            double[] correction = new double[rows * rows];
            backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
                rows, rows, model.fixedColumns(), 1.0,
                temporary, inverseCovarianceFixed, 0.0, correction);
            double[] projection = MatrixOps.subtract(
                inverseCovariance, correction);

            List<SmoothTermEstimate> estimates = new ArrayList<>(terms.size());
            int fixedStart = parametricColumns;
            double totalEdf = parametricColumns;
            for (int index = 0; index < terms.size(); index++) {
                DecomposedTerm term = terms.get(index);
                double[] fixedCoefficients = Arrays.copyOfRange(beta,
                    fixedStart, fixedStart + term.fixedColumns());
                fixedStart += term.fixedColumns();
                double termVariance = variances[index];
                double[] randomCoefficients = transposeMultiply(
                    term.randomDesign(), rows, term.randomColumns(),
                    projected, termVariance);
                double[] fittedContribution = contribution(
                    term.fixedDesign(), rows, term.fixedColumns(),
                    fixedCoefficients, term.randomDesign(),
                    term.randomColumns(), randomCoefficients);
                double edf = term.fixedColumns() + termVariance
                    * traceQuadratic(term.randomDesign(), rows,
                        term.randomColumns(), projection);
                edf = Math.max(term.fixedColumns(),
                    Math.min(term.fixedColumns() + term.randomColumns(), edf));
                estimates.add(new SmoothTermEstimate(term.term(),
                    term.fixedTransform(), term.fixedMeans(), fixedCoefficients,
                    term.randomTransform(), term.randomMeans(),
                    randomCoefficients, fittedContribution,
                    1.0 / termVariance, edf));
                totalEdf += edf;
            }
            return new GeneralizedGamResult(
                family, fitted, parametricColumns, estimates, totalEdf);
        }
    }

    private static ModelMatrices assemble(
            double[] parametric,
            int rows,
            int parametricColumns,
            List<DecomposedTerm> terms) {
        int fixedColumns = parametricColumns;
        for (DecomposedTerm term : terms) fixedColumns += term.fixedColumns();
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
        for (DecomposedTerm term : terms) {
            for (int row = 0; row < rows; row++) {
                System.arraycopy(term.fixedDesign(), row * term.fixedColumns(),
                    fixed, row * fixedColumns + destination, term.fixedColumns());
            }
            destination += term.fixedColumns();
        }
        return new ModelMatrices(fixed, fixedColumns);
    }

    private static DecomposedTerm decompose(
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
        double[] inverseTimesDifference = backend.dpotrf(gram, penaltyRows)
            .solve(difference, basisColumns);
        double[] randomTransform = new double[basisColumns * penaltyRows];
        for (int row = 0; row < basisColumns; row++) {
            for (int column = 0; column < penaltyRows; column++) {
                randomTransform[row * penaltyRows + column] =
                    inverseTimesDifference[column * basisColumns + row];
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
        CenteredMatrix centeredFixed = centerAndDropConstant(rawFixed,
            term.observations(), nullColumns, fixedTransform, basisColumns);
        double[] rawRandom = multiply(term.designView(), term.observations(),
            basisColumns, randomTransform, penaltyRows);
        CenteredMatrix centeredRandom = center(rawRandom,
            term.observations(), penaltyRows);
        return new DecomposedTerm(term,
            centeredFixed.matrix(), centeredFixed.transform(),
            centeredFixed.means(), centeredFixed.columns(),
            centeredRandom.matrix(), randomTransform,
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
            reducedMeans[destination++] = centered.means()[column];
        }
        return new CenteredMatrix(
            reduced, reducedTransform, reducedMeans, kept);
    }

    private static CenteredMatrix center(
            double[] matrix, int rows, int columns) {
        double[] centered = matrix.clone();
        double[] means = new double[columns];
        for (int column = 0; column < columns; column++) {
            for (int row = 0; row < rows; row++) {
                means[column] += matrix[row * columns + column];
            }
            means[column] /= rows;
            for (int row = 0; row < rows; row++) {
                centered[row * columns + column] -= means[column];
            }
        }
        return new CenteredMatrix(centered, new double[0], means, columns);
    }

    private static double[] covarianceBasis(
            double[] design, int rows, int columns) {
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

    private static double traceQuadratic(
            double[] design, int rows, int columns, double[] projection) {
        double result = 0.0;
        for (int column = 0; column < columns; column++) {
            for (int row = 0; row < rows; row++) {
                double left = design[row * columns + column];
                for (int other = 0; other < rows; other++) {
                    result += left * projection[row * rows + other]
                        * design[other * columns + column];
                }
            }
        }
        return result;
    }

    private static double[] transposeMultiply(
            double[] design,
            int rows,
            int columns,
            double[] vector,
            double scale) {
        double[] result = new double[columns];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                result[column] += scale * design[row * columns + column]
                    * vector[row];
            }
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

    private static double[] weights(double[] supplied, int rows) {
        if (supplied == null) {
            double[] result = new double[rows];
            Arrays.fill(result, 1.0);
            return result;
        }
        if (supplied.length != rows) {
            throw new IllegalArgumentException("prior weight length must equal rows");
        }
        double[] result = MatrixOps.finiteCopy(supplied, "priorWeights");
        for (double value : result) {
            if (!(value > 0.0)) {
                throw new IllegalArgumentException(
                    "prior weights must be strictly positive");
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

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
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

    private record ModelMatrices(double[] fixedDesign, int fixedColumns) { }
}
