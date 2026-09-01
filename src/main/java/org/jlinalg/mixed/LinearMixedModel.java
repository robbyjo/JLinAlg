/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.Reml;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.RemlResult;
import org.jlinalg.reml.VarianceComponent;

/** Exact Gaussian REML for one or more independent random-effect terms. */
public final class LinearMixedModel {
    private static final String RESIDUAL_COMPONENT = "residual";

    private LinearMixedModel() { }

    /** Fits with default REML controls and preferred acceleration. */
    public static LinearMixedModelResult fit(
            double[] response,
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects) {
        return fit(response, fixedEffects, randomEffects,
            RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits from conventional fixed- and random-effect design matrices. */
    public static LinearMixedModelResult fit(
            double[] response,
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, rowMajor, response.length, fixedEffects[0].length,
            randomEffects, options, backendPolicy);
    }

    /** Fits from a contiguous row-major fixed-effect design matrix. */
    public static LinearMixedModelResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        return fitWithResidualCorrelation(response, fixedEffects, rows, columns,
            randomEffects, MatrixOps.identity(rows), options, backendPolicy);
    }

    /**
     * Fits a mixed model whose residual variance multiplies a caller-supplied
     * positive-definite correlation matrix.
     */
    public static LinearMixedModelResult fitWithResidualCorrelation(
            double[] response,
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects,
            double[] residualCorrelation,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(fixedEffects, response.length);
        return fitWithResidualCorrelation(response, rowMajor, response.length,
            fixedEffects[0].length, randomEffects, residualCorrelation,
            options, backendPolicy);
    }

    /** Contiguous-row-major overload with structured residual correlation. */
    public static LinearMixedModelResult fitWithResidualCorrelation(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<RandomEffectTerm> randomEffects,
            double[] residualCorrelation,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        validateTerms(randomEffects, rows);
        if (options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "options and backendPolicy are required");
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            List<VarianceComponent> components = new ArrayList<>(
                randomEffects.size() + 1);
            List<double[]> bases = new ArrayList<>(randomEffects.size());
            for (RandomEffectTerm term : randomEffects) {
                double[] basis = covarianceBasis(term, backend);
                bases.add(basis);
                components.add(new VarianceComponent(
                    term.name(), rows, basis));
            }
            VarianceComponent residual = new VarianceComponent(
                RESIDUAL_COMPONENT, rows, residualCorrelation);
            validateCorrelationDiagonal(residual.covariance(), rows);
            components.add(residual);
            RemlResult reml = Reml.fitWithKnownCovariance(
                response, fixedEffects, rows, columns, components,
                new double[rows * rows], options, context);
            return predictions(response, fixedEffects, rows, columns,
                randomEffects, bases, residual.covariance(), reml, backend);
        }
    }

    private static LinearMixedModelResult predictions(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<RandomEffectTerm> terms,
            List<double[]> bases,
            double[] residualCorrelation,
            RemlResult reml,
            ComputeBackend backend) {
        if (isIdentity(residualCorrelation, rows)) {
            return hendersonPredictions(response, fixedEffects, rows, columns,
                terms, reml, backend);
        }
        double[] variances = reml.varianceComponents();
        double[] covariance = new double[rows * rows];
        for (int term = 0; term < terms.size(); term++) {
            for (int index = 0; index < covariance.length; index++) {
                covariance[index] += variances[term] * bases.get(term)[index];
            }
        }
        double residualVariance = variances[terms.size()];
        for (int index = 0; index < covariance.length; index++) {
            covariance[index] += residualVariance * residualCorrelation[index];
        }
        CholeskyFactor factor = backend.dpotrf(covariance, rows);
        double[] marginalResidual = reml.residuals();
        double[] projectedResidual = factor.solve(marginalResidual);

        double[] inverseCovariance = factor.solve(MatrixOps.identity(rows), rows);
        double[] inverseCovarianceFixed = factor.solve(fixedEffects, columns);
        double[] temporary = MatrixOps.multiply(
            backend, inverseCovarianceFixed, rows, columns,
            reml.fixedEffectCovariance(), columns);
        double[] fixedCorrection = new double[rows * rows];
        backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
            rows, rows, columns, 1.0,
            temporary, inverseCovarianceFixed, 0.0, fixedCorrection);
        double[] projection = MatrixOps.subtract(
            inverseCovariance, fixedCorrection);

        double[] conditionalFitted = MatrixOps.multiply(
            backend, fixedEffects, rows, columns, reml.fixedEffects());
        List<RandomEffectEstimates> estimates = new ArrayList<>(terms.size());
        for (int termIndex = 0; termIndex < terms.size(); termIndex++) {
            RandomEffectTerm term = terms.get(termIndex);
            double variance = variances[termIndex];
            double[] values = transposeMultiply(
                term, projectedResidual, variance, backend);
            double[] termFitted = multiply(term, values, backend);
            for (int row = 0; row < rows; row++) {
                conditionalFitted[row] += termFitted[row];
            }

            double[] reductionDiagonal = projectionDiagonal(
                term, projection, backend);
            double[] predictionErrorVariances = new double[term.coefficients()];
            for (int coefficient = 0;
                    coefficient < term.coefficients(); coefficient++) {
                predictionErrorVariances[coefficient] = clamp(
                    variance - variance * variance
                        * reductionDiagonal[coefficient],
                    0.0, variance);
            }
            estimates.add(new RandomEffectEstimates(
                term.name(), term.coefficientNames(), variance,
                values, predictionErrorVariances));
        }
        return new LinearMixedModelResult(reml, estimates,
            conditionalFitted, MatrixOps.subtract(response, conditionalFitted));
    }

    private static LinearMixedModelResult hendersonPredictions(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<RandomEffectTerm> terms,
            RemlResult reml,
            ComputeBackend backend) {
        int randomColumns = 0;
        int[] starts = new int[terms.size()];
        for (int term = 0; term < terms.size(); term++) {
            starts[term] = randomColumns;
            randomColumns += terms.get(term).coefficients();
        }
        int dimension = columns + randomColumns;
        double[] variances = reml.varianceComponents();
        double inverseResidual = 1.0 / variances[terms.size()];
        double[] equations = new double[dimension * dimension];
        double[] rightSide = new double[dimension];
        int[] randomIndices = new int[randomColumns];
        double[] randomValues = new double[randomColumns];

        for (int row = 0; row < rows; row++) {
            int nonzeros = rowEntries(
                terms, starts, row, randomIndices, randomValues);
            double weightedResponse = inverseResidual * response[row];
            for (int fixed = 0; fixed < columns; fixed++) {
                double value = fixedEffects[row * columns + fixed];
                rightSide[fixed] += value * weightedResponse;
                for (int second = 0; second < columns; second++) {
                    equations[fixed * dimension + second] += inverseResidual
                        * value * fixedEffects[row * columns + second];
                }
                for (int random = 0; random < nonzeros; random++) {
                    int coefficient = columns + randomIndices[random];
                    double cross = inverseResidual * value * randomValues[random];
                    equations[fixed * dimension + coefficient] += cross;
                    equations[coefficient * dimension + fixed] += cross;
                }
            }
            for (int first = 0; first < nonzeros; first++) {
                int firstCoefficient = columns + randomIndices[first];
                rightSide[firstCoefficient] += randomValues[first] * weightedResponse;
                for (int second = 0; second < nonzeros; second++) {
                    int secondCoefficient = columns + randomIndices[second];
                    equations[firstCoefficient * dimension + secondCoefficient] +=
                        inverseResidual * randomValues[first] * randomValues[second];
                }
            }
        }
        for (int term = 0; term < terms.size(); term++) {
            double precision = 1.0 / variances[term];
            for (int coefficient = 0;
                    coefficient < terms.get(term).coefficients(); coefficient++) {
                int index = columns + starts[term] + coefficient;
                equations[index * dimension + index] += precision;
            }
        }

        CholeskyFactor factor = backend.dpotrf(equations, dimension);
        double[] solution = factor.solve(rightSide);
        double[] coefficientCovariance = factor.solve(
            MatrixOps.identity(dimension), dimension);
        double[] conditionalFitted = new double[rows];
        List<RandomEffectEstimates> estimates = new ArrayList<>(terms.size());
        for (int row = 0; row < rows; row++) {
            for (int fixed = 0; fixed < columns; fixed++) {
                conditionalFitted[row] += fixedEffects[row * columns + fixed]
                    * solution[fixed];
            }
        }
        for (int term = 0; term < terms.size(); term++) {
            RandomEffectTerm value = terms.get(term);
            double[] modes = new double[value.coefficients()];
            double[] pev = new double[value.coefficients()];
            for (int coefficient = 0; coefficient < value.coefficients(); coefficient++) {
                int index = columns + starts[term] + coefficient;
                modes[coefficient] = solution[index];
                pev[coefficient] = clamp(
                    coefficientCovariance[index * dimension + index],
                    0.0, variances[term]);
            }
            double[] contribution = multiply(value, modes, backend);
            for (int row = 0; row < rows; row++) {
                conditionalFitted[row] += contribution[row];
            }
            estimates.add(new RandomEffectEstimates(
                value.name(), value.coefficientNames(), variances[term], modes, pev));
        }
        return new LinearMixedModelResult(reml, estimates,
            conditionalFitted, MatrixOps.subtract(response, conditionalFitted));
    }

    private static int rowEntries(
            List<RandomEffectTerm> terms,
            int[] starts,
            int row,
            int[] indices,
            double[] values) {
        int count = 0;
        for (int termIndex = 0; termIndex < terms.size(); termIndex++) {
            RandomEffectTerm term = terms.get(termIndex);
            if (term.sparse()) {
                int[] pointers = term.rowPointersView();
                int[] columns = term.columnIndicesView();
                double[] sparse = term.sparseValuesView();
                for (int position = pointers[row]; position < pointers[row + 1]; position++) {
                    if (sparse[position] != 0.0) {
                        indices[count] = starts[termIndex] + columns[position];
                        values[count++] = sparse[position];
                    }
                }
            } else {
                double[] dense = term.denseDesignView();
                for (int column = 0; column < term.coefficients(); column++) {
                    double value = dense[row * term.coefficients() + column];
                    if (value != 0.0) {
                        indices[count] = starts[termIndex] + column;
                        values[count++] = value;
                    }
                }
            }
        }
        return count;
    }

    private static boolean isIdentity(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                double expected = row == column ? 1.0 : 0.0;
                if (matrix[row * dimension + column] != expected) return false;
            }
        }
        return true;
    }

    private static double[] covarianceBasis(
            RandomEffectTerm term, ComputeBackend backend) {
        double[] result = new double[term.observations() * term.observations()];
        if (!term.sparse()) {
            backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
                term.observations(), term.observations(), term.coefficients(),
                1.0, term.denseDesignView(), term.denseDesignView(), 0.0, result);
            return result;
        }
        int[] pointers = term.rowPointersView();
        int[] columns = term.columnIndicesView();
        double[] values = term.sparseValuesView();
        for (int row = 0; row < term.observations(); row++) {
            for (int column = 0; column <= row; column++) {
                int left = pointers[row];
                int right = pointers[column];
                double value = 0.0;
                while (left < pointers[row + 1]
                        && right < pointers[column + 1]) {
                    int leftColumn = columns[left];
                    int rightColumn = columns[right];
                    if (leftColumn < rightColumn) {
                        left++;
                    } else if (leftColumn > rightColumn) {
                        right++;
                    } else {
                        value += values[left++] * values[right++];
                    }
                }
                result[row * term.observations() + column] = value;
                result[column * term.observations() + row] = value;
            }
        }
        return result;
    }

    private static double[] transposeMultiply(
            RandomEffectTerm term,
            double[] vector,
            double scale,
            ComputeBackend backend) {
        double[] result = new double[term.coefficients()];
        if (!term.sparse()) {
            backend.dgemv(MatrixTranspose.TRANSPOSE,
                term.observations(), term.coefficients(), scale,
                term.denseDesignView(), vector, 0.0, result);
            return result;
        }
        int[] pointers = term.rowPointersView();
        int[] columns = term.columnIndicesView();
        double[] values = term.sparseValuesView();
        for (int row = 0; row < term.observations(); row++) {
            for (int index = pointers[row]; index < pointers[row + 1]; index++) {
                result[columns[index]] += scale * values[index] * vector[row];
            }
        }
        return result;
    }

    private static double[] multiply(
            RandomEffectTerm term,
            double[] vector,
            ComputeBackend backend) {
        if (!term.sparse()) {
            return MatrixOps.multiply(backend, term.denseDesignView(),
                term.observations(), term.coefficients(), vector);
        }
        double[] result = new double[term.observations()];
        int[] pointers = term.rowPointersView();
        int[] columns = term.columnIndicesView();
        double[] values = term.sparseValuesView();
        for (int row = 0; row < term.observations(); row++) {
            for (int index = pointers[row]; index < pointers[row + 1]; index++) {
                result[row] += values[index] * vector[columns[index]];
            }
        }
        return result;
    }

    private static double[] projectionDiagonal(
            RandomEffectTerm term,
            double[] projection,
            ComputeBackend backend) {
        if (!term.sparse()) {
            double[] projectedDesign = MatrixOps.multiply(
                backend, projection, term.observations(), term.observations(),
                term.denseDesignView(), term.coefficients());
            double[] result = new double[term.coefficients()];
            for (int coefficient = 0;
                    coefficient < term.coefficients(); coefficient++) {
                for (int row = 0; row < term.observations(); row++) {
                    result[coefficient] += term.denseDesignView()[
                        row * term.coefficients() + coefficient]
                        * projectedDesign[row * term.coefficients() + coefficient];
                }
            }
            return result;
        }
        double[] result = new double[term.coefficients()];
        int[] pointers = term.rowPointersView();
        int[] columns = term.columnIndicesView();
        double[] values = term.sparseValuesView();
        int observations = term.observations();
        for (int row = 0; row < observations; row++) {
            for (int left = pointers[row]; left < pointers[row + 1]; left++) {
                int coefficient = columns[left];
                for (int other = 0; other < observations; other++) {
                    for (int right = pointers[other];
                            right < pointers[other + 1]; right++) {
                        if (columns[right] == coefficient) {
                            result[coefficient] += values[left]
                                * projection[row * observations + other]
                                * values[right];
                            break;
                        }
                        if (columns[right] > coefficient) {
                            break;
                        }
                    }
                }
            }
        }
        return result;
    }

    private static void validateTerms(
            List<RandomEffectTerm> terms, int observations) {
        if (terms == null || terms.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one random-effect term is required");
        }
        Set<String> names = new HashSet<>();
        for (RandomEffectTerm term : terms) {
            if (term == null || term.observations() != observations) {
                throw new IllegalArgumentException(
                    "random-effect terms must match the observation count");
            }
            if (!names.add(term.name()) || RESIDUAL_COMPONENT.equals(term.name())) {
                throw new IllegalArgumentException(
                    "random-effect term names must be unique and not 'residual'");
            }
        }
    }

    private static void validateCorrelationDiagonal(
            double[] correlation, int dimension) {
        for (int index = 0; index < dimension; index++) {
            if (Math.abs(correlation[index * dimension + index] - 1.0) > 1e-10) {
                throw new IllegalArgumentException(
                    "residual correlation must have a unit diagonal");
            }
        }
    }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }
}
