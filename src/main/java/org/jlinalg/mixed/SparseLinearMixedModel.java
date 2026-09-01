/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.OptimizationResult;
import jdistlib.matrix.CsrMatrix;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceEstimation;

/** Sparse-equation REML/ML for independent Gaussian random-effect terms. */
public final class SparseLinearMixedModel {
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    private static final double INVALID_OBJECTIVE = Double.MAX_VALUE / 16.0;

    private SparseLinearMixedModel() { }

    public static SparseLinearMixedModelResult fit(
            double[] response,
            double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null)
            throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, fixed, response.length, fixedEffects[0].length,
            randomEffects, options, backendPolicy);
    }

    /** Fits without constructing observation-scale covariance matrices. */
    public static SparseLinearMixedModelResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        return fitWithPrecision(response, fixedEffects, rows, columns,
            randomEffects, null, options, backendPolicy);
    }

    /** Fits terms with caller-supplied sparse coefficient precision bases. */
    public static SparseLinearMixedModelResult fitWithPrecision(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        validate(randomEffects, rows, options, backendPolicy);
        if (options.degreesOfFreedomMethod()
                != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION)
            throw new IllegalArgumentException(
                "sparse LMM currently supports residual-approximation DF; "
                    + "use the dense reference fitter for Satterthwaite or Kenward-Roger");
        CombinedDesign design = combine(randomEffects, rows);
        List<SparsePrecisionMatrix> precisions = precisions(
            randomEffects, precisionBases);
        SparsePattern pattern = crossProductPattern(
            design, randomEffects, precisions);
        double[] initial = initialLogVariances(
            response, randomEffects.size() + 1, options);
        double[] lower = new double[initial.length];
        double[] upper = new double[initial.length];
        Arrays.fill(lower, Math.log(options.minimumVariance()));
        Arrays.fill(upper, Math.log(options.maximumVariance()));
        int maximumEvaluations = Math.max(100,
            options.maximumIterations() * 10);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] precisionLogDeterminants = precisionLogDeterminants(
                precisions, backend);
            CsrMatrix initialMatrix = pattern.matrix(exp(initial));
            try (PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                    initialMatrix, MatrixTriangle.LOWER,
                    SparseOrdering.MINIMUM_DEGREE)) {
                Objective objective = new Objective(response, fixedEffects,
                    rows, columns, randomEffects, design, pattern, factor,
                    precisionLogDeterminants,
                    options.varianceEstimation(), backend);
                OptimizationResult optimized = Bobyqa.bobyqa(initial.clone(),
                    lower, upper, parameters -> objective.value(parameters),
                    Math.min(2 * initial.length + 1,
                        (initial.length + 1) * (initial.length + 2) / 2),
                    Math.min(1.0, options.maximumLogVarianceStep()),
                    Math.max(1e-7, options.relativeTolerance()),
                    maximumEvaluations,
                    true);
                Evaluation fitted = objective.evaluate(optimized.mX);
                double[] variances = exp(optimized.mX);
                double[] standardErrors = new double[columns];
                for (int column = 0; column < columns; column++)
                    standardErrors[column] = Math.sqrt(Math.max(0.0,
                        fitted.fixedCovariance()[column * columns + column]));
                double degrees = rows - columns - 1.0;
                if (!(degrees > 0.0))
                    throw new IllegalArgumentException(
                        "sparse LMM requires positive denominator DF");
                AssociationStatistics association = AssociationStatistics.studentT(
                    fitted.beta(), standardErrors, degrees,
                    DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
                List<RandomEffectEstimates> estimates = estimates(
                    randomEffects, design.termStarts(), variances,
                    fitted.randomModes());
                List<String> names = new ArrayList<>(randomEffects.size() + 1);
                for (RandomEffectTerm term : randomEffects) names.add(term.name());
                names.add("residual");
                return new SparseLinearMixedModelResult(names, variances,
                    association, fitted.fixedCovariance(), estimates,
                    fitted.conditionalFitted(), fitted.conditionalResiduals(),
                    fitted.logLikelihood(), options.varianceEstimation(),
                    optimized.numFunctionCalls,
                    optimized.numFunctionCalls < maximumEvaluations,
                    design.columns(), pattern.values().length,
                    factor.factorNonzeroCount(), context.provenance());
            }
        }
    }

    private static final class Objective {
        private final double[] response;
        private final double[] fixed;
        private final int rows;
        private final int columns;
        private final List<RandomEffectTerm> terms;
        private final CombinedDesign design;
        private final SparsePattern pattern;
        private final PreparedSparseCholesky factor;
        private final double[] precisionLogDeterminants;
        private final VarianceEstimation estimation;
        private final ComputeBackend backend;

        Objective(double[] response, double[] fixed, int rows, int columns,
                List<RandomEffectTerm> terms, CombinedDesign design,
                SparsePattern pattern, PreparedSparseCholesky factor,
                double[] precisionLogDeterminants,
                VarianceEstimation estimation, ComputeBackend backend) {
            this.response = response;
            this.fixed = fixed;
            this.rows = rows;
            this.columns = columns;
            this.terms = terms;
            this.design = design;
            this.pattern = pattern;
            this.factor = factor;
            this.precisionLogDeterminants = precisionLogDeterminants;
            this.estimation = estimation;
            this.backend = backend;
        }

        double value(double[] parameters) {
            try {
                Evaluation result = evaluate(parameters);
                return Double.isFinite(result.logLikelihood())
                    ? -result.logLikelihood() : INVALID_OBJECTIVE;
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return INVALID_OBJECTIVE;
            }
        }

        Evaluation evaluate(double[] logVariances) {
            double[] variances = exp(logVariances);
            factor.refactor(pattern.matrix(variances));
            double residualVariance = variances[variances.length - 1];
            double[] inverseY = inverseMarginal(response, residualVariance);
            double[] inverseX = new double[rows * columns];
            double[] source = new double[rows];
            for (int column = 0; column < columns; column++) {
                for (int row = 0; row < rows; row++)
                    source[row] = fixed[row * columns + column];
                double[] solved = inverseMarginal(source, residualVariance);
                for (int row = 0; row < rows; row++)
                    inverseX[row * columns + column] = solved[row];
            }
            double[] information = MatrixOps.transposeMultiply(backend,
                fixed, rows, columns, inverseX, columns);
            symmetrize(information, columns);
            CholeskyFactor fixedFactor = backend.dpotrf(information, columns);
            double[] rightSide = new double[columns];
            backend.dgemv(jdistlib.accelerator.MatrixTranspose.TRANSPOSE,
                rows, columns, 1.0, fixed, inverseY, 0.0, rightSide);
            double[] beta = fixedFactor.solve(rightSide);
            double quadratic = backend.ddot(rows, response, 0, 1,
                inverseY, 0, 1) - backend.ddot(columns, rightSide, 0, 1,
                beta, 0, 1);
            double logDeterminant = rows * Math.log(residualVariance)
                + factor.logDeterminant();
            for (int term = 0; term < terms.size(); term++)
                logDeterminant += terms.get(term).coefficients()
                    * Math.log(variances[term])
                    - precisionLogDeterminants[term];
            boolean restricted = estimation == VarianceEstimation.REML;
            double logLikelihood = -0.5 * (
                (restricted ? rows - columns : rows) * LOG_TWO_PI
                    + logDeterminant
                    + (restricted ? fixedFactor.logDeterminant() : 0.0)
                    + Math.max(0.0, quadratic));
            double[] fixedCovariance = fixedFactor.solve(
                MatrixOps.identity(columns), columns);
            symmetrize(fixedCovariance, columns);
            double[] residual = MatrixOps.subtract(response,
                MatrixOps.multiply(backend, fixed, rows, columns, beta));
            double[] randomRightSide = transposeMultiply(design, residual);
            for (int index = 0; index < randomRightSide.length; index++)
                randomRightSide[index] /= residualVariance;
            double[] randomModes = factor.solve(randomRightSide);
            double[] conditionalFitted = MatrixOps.multiply(
                backend, fixed, rows, columns, beta);
            double[] randomFitted = multiply(design, randomModes);
            for (int row = 0; row < rows; row++)
                conditionalFitted[row] += randomFitted[row];
            return new Evaluation(beta, fixedCovariance, randomModes,
                conditionalFitted,
                MatrixOps.subtract(response, conditionalFitted),
                logLikelihood);
        }

        private double[] inverseMarginal(
                double[] vector, double residualVariance) {
            double[] rightSide = transposeMultiply(design, vector);
            for (int index = 0; index < rightSide.length; index++)
                rightSide[index] /= residualVariance;
            double[] random = factor.solve(rightSide);
            double[] fitted = multiply(design, random);
            double[] result = new double[rows];
            for (int row = 0; row < rows; row++)
                result[row] = (vector[row] - fitted[row]) / residualVariance;
            return result;
        }
    }

    private static SparsePattern crossProductPattern(
            CombinedDesign design, List<RandomEffectTerm> terms,
            List<SparsePrecisionMatrix> precisions) {
        int dimension = design.columns();
        TreeMap<Long, Double> cross = new TreeMap<>();
        for (int row = 0; row < design.rows(); row++) {
            for (int left = design.rowStarts()[row];
                    left < design.rowStarts()[row + 1]; left++) {
                int leftColumn = design.columnIndices()[left];
                for (int right = design.rowStarts()[row]; right <= left; right++) {
                    int rightColumn = design.columnIndices()[right];
                    int matrixRow = Math.max(leftColumn, rightColumn);
                    int matrixColumn = Math.min(leftColumn, rightColumn);
                    long key = (long) matrixRow * dimension + matrixColumn;
                    cross.merge(key, design.values()[left]
                        * design.values()[right], Double::sum);
                }
            }
        }
        TreeMap<Long, PrecisionEntry> precision = new TreeMap<>();
        for (int term = 0; term < terms.size(); term++) {
            SparsePrecisionMatrix basis = precisions.get(term);
            int[] starts = basis.rowStarts();
            int[] columns = basis.columnIndices();
            double[] values = basis.values();
            int offset = design.termStarts()[term];
            for (int row = 0; row < basis.dimension(); row++)
                for (int index = starts[row]; index < starts[row + 1]; index++) {
                    int column = columns[index];
                    if (row < column) continue;
                    long key = (long) (offset + row) * dimension
                        + offset + column;
                    precision.put(key, new PrecisionEntry(term, values[index]));
                    cross.putIfAbsent(key, 0.0);
                }
        }
        int[] rowStarts = new int[dimension + 1];
        int[] columnIndices = new int[cross.size()];
        double[] values = new double[cross.size()];
        int[] precisionTerms = new int[cross.size()];
        double[] precisionValues = new double[cross.size()];
        Arrays.fill(precisionTerms, -1);
        int position = 0;
        int currentRow = 0;
        for (Map.Entry<Long, Double> entry : cross.entrySet()) {
            int row = (int) (entry.getKey() / dimension);
            int column = (int) (entry.getKey() % dimension);
            while (currentRow < row) rowStarts[++currentRow] = position;
            columnIndices[position] = column;
            values[position] = entry.getValue();
            PrecisionEntry precisionEntry = precision.get(entry.getKey());
            if (precisionEntry != null) {
                precisionTerms[position] = precisionEntry.term();
                precisionValues[position] = precisionEntry.value();
            }
            position++;
        }
        while (currentRow < dimension) rowStarts[++currentRow] = position;
        return new SparsePattern(dimension, rowStarts, columnIndices,
            values, precisionTerms, precisionValues);
    }

    private static List<SparsePrecisionMatrix> precisions(
            List<RandomEffectTerm> terms,
            List<SparsePrecisionMatrix> supplied) {
        if (supplied == null) {
            List<SparsePrecisionMatrix> result = new ArrayList<>(terms.size());
            for (RandomEffectTerm term : terms)
                result.add(SparsePrecisionMatrix.identity(term.coefficients()));
            return List.copyOf(result);
        }
        if (supplied.size() != terms.size())
            throw new IllegalArgumentException(
                "one precision basis is required per random-effect term");
        for (int term = 0; term < terms.size(); term++)
            if (supplied.get(term) == null
                    || supplied.get(term).dimension()
                        != terms.get(term).coefficients())
                throw new IllegalArgumentException(
                    "precision dimensions must match random coefficients");
        return List.copyOf(supplied);
    }

    private static double[] precisionLogDeterminants(
            List<SparsePrecisionMatrix> precisions, ComputeBackend backend) {
        double[] result = new double[precisions.size()];
        for (int term = 0; term < precisions.size(); term++) {
            SparsePrecisionMatrix value = precisions.get(term);
            int[] rows = value.rowStarts();
            int[] columns = value.columnIndices();
            double[] numeric = value.values();
            int lowerCount = 0;
            for (int row = 0; row < value.dimension(); row++)
                for (int index = rows[row]; index < rows[row + 1]; index++)
                    if (columns[index] <= row) lowerCount++;
            int[] lowerRows = new int[value.dimension() + 1];
            int[] lowerColumns = new int[lowerCount];
            double[] lowerValues = new double[lowerCount];
            int position = 0;
            for (int row = 0; row < value.dimension(); row++) {
                lowerRows[row] = position + 1;
                for (int index = rows[row]; index < rows[row + 1]; index++)
                    if (columns[index] <= row) {
                        lowerColumns[position] = columns[index] + 1;
                        lowerValues[position++] = numeric[index];
                    }
            }
            lowerRows[value.dimension()] = position + 1;
            CsrMatrix matrix = new CsrMatrix(value.dimension(),
                value.dimension(), lowerValues, lowerColumns, lowerRows);
            result[term] = backend.dcsrpotrf(matrix, MatrixTriangle.LOWER,
                SparseOrdering.MINIMUM_DEGREE).logDeterminant();
        }
        return result;
    }

    private static CombinedDesign combine(
            List<RandomEffectTerm> terms, int rows) {
        int[] starts = new int[terms.size()];
        int columns = 0;
        int nonzeros = 0;
        for (int term = 0; term < terms.size(); term++) {
            starts[term] = columns;
            columns += terms.get(term).coefficients();
            nonzeros += terms.get(term).nonzeroCount();
        }
        int[] rowStarts = new int[rows + 1];
        int[] columnIndices = new int[nonzeros];
        double[] values = new double[nonzeros];
        TermData[] data = new TermData[terms.size()];
        for (int term = 0; term < terms.size(); term++)
            data[term] = termData(terms.get(term));
        int position = 0;
        for (int row = 0; row < rows; row++) {
            rowStarts[row] = position;
            for (int term = 0; term < terms.size(); term++) {
                TermData value = data[term];
                for (int index = value.rowStarts()[row];
                        index < value.rowStarts()[row + 1]; index++) {
                    if (value.values()[index] != 0.0) {
                        columnIndices[position] = starts[term]
                            + value.columnIndices()[index];
                        values[position++] = value.values()[index];
                    }
                }
            }
        }
        rowStarts[rows] = position;
        return new CombinedDesign(rows, columns, rowStarts,
            Arrays.copyOf(columnIndices, position),
            Arrays.copyOf(values, position), starts);
    }

    private static TermData termData(RandomEffectTerm term) {
        if (term.sparse())
            return new TermData(term.rowPointers(), term.columnIndices(),
                term.sparseValues());
        double[] dense = term.design();
        int[] starts = new int[term.observations() + 1];
        int count = 0;
        for (double value : dense) if (value != 0.0) count++;
        int[] columns = new int[count];
        double[] values = new double[count];
        int position = 0;
        for (int row = 0; row < term.observations(); row++) {
            starts[row] = position;
            for (int column = 0; column < term.coefficients(); column++) {
                double value = dense[row * term.coefficients() + column];
                if (value != 0.0) {
                    columns[position] = column;
                    values[position++] = value;
                }
            }
        }
        starts[term.observations()] = position;
        return new TermData(starts, columns, values);
    }

    private static double[] transposeMultiply(
            CombinedDesign design, double[] vector) {
        double[] result = new double[design.columns()];
        for (int row = 0; row < design.rows(); row++)
            for (int index = design.rowStarts()[row];
                    index < design.rowStarts()[row + 1]; index++)
                result[design.columnIndices()[index]] +=
                    design.values()[index] * vector[row];
        return result;
    }

    private static double[] multiply(
            CombinedDesign design, double[] vector) {
        double[] result = new double[design.rows()];
        for (int row = 0; row < design.rows(); row++)
            for (int index = design.rowStarts()[row];
                    index < design.rowStarts()[row + 1]; index++)
                result[row] += design.values()[index]
                    * vector[design.columnIndices()[index]];
        return result;
    }

    private static List<RandomEffectEstimates> estimates(
            List<RandomEffectTerm> terms, int[] starts,
            double[] variances, double[] modes) {
        List<RandomEffectEstimates> result = new ArrayList<>(terms.size());
        for (int term = 0; term < terms.size(); term++) {
            RandomEffectTerm value = terms.get(term);
            double[] selected = Arrays.copyOfRange(modes, starts[term],
                starts[term] + value.coefficients());
            double[] unavailablePev = new double[value.coefficients()];
            Arrays.fill(unavailablePev, Double.NaN);
            result.add(new RandomEffectEstimates(value.name(),
                value.coefficientNames(), variances[term], selected,
                unavailablePev));
        }
        return List.copyOf(result);
    }

    private static double[] initialLogVariances(
            double[] response, int count, RemlOptions options) {
        double[] supplied = options.initialVariances();
        if (supplied != null && supplied.length != count)
            throw new IllegalArgumentException(
                "initial variance count must equal random terms plus residual");
        double[] result = new double[count];
        if (supplied != null) {
            for (int index = 0; index < count; index++)
                result[index] = Math.log(supplied[index]);
            return result;
        }
        double mean = Arrays.stream(response).average().orElse(0.0);
        double variance = 0.0;
        for (double value : response) variance += (value - mean) * (value - mean);
        variance = Math.max(options.minimumVariance(),
            variance / Math.max(1, response.length - 1) / count);
        Arrays.fill(result, Math.log(variance));
        return result;
    }

    private static void validate(
            List<RandomEffectTerm> terms, int rows,
            RemlOptions options, BackendPolicy backendPolicy) {
        if (terms == null || terms.isEmpty() || options == null
                || backendPolicy == null)
            throw new IllegalArgumentException(
                "random effects, options, and backend are required");
        for (RandomEffectTerm term : terms)
            if (term == null || term.observations() != rows)
                throw new IllegalArgumentException(
                    "random-effect rows must equal observations");
    }

    private static double[] exp(double[] values) {
        double[] result = new double[values.length];
        for (int index = 0; index < values.length; index++)
            result[index] = Math.exp(values[index]);
        return result;
    }

    private static void symmetrize(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++)
            for (int column = 0; column < row; column++) {
                double value = 0.5 * (matrix[row * dimension + column]
                    + matrix[column * dimension + row]);
                matrix[row * dimension + column] = value;
                matrix[column * dimension + row] = value;
            }
    }

    private record TermData(int[] rowStarts, int[] columnIndices,
                            double[] values) { }
    private record CombinedDesign(int rows, int columns, int[] rowStarts,
                                  int[] columnIndices, double[] values,
                                  int[] termStarts) { }
    private record SparsePattern(int dimension, int[] rowStarts,
                                 int[] columnIndices, double[] values,
                                 int[] precisionTerms,
                                 double[] precisionValues) {
        CsrMatrix matrix(double[] variances) {
            double residual = variances[variances.length - 1];
            double[] numeric = values.clone();
            for (int index = 0; index < numeric.length; index++) {
                numeric[index] /= residual;
                if (precisionTerms[index] >= 0)
                    numeric[index] += precisionValues[index]
                        / variances[precisionTerms[index]];
            }
            int[] oneBasedColumns = columnIndices.clone();
            int[] oneBasedRows = rowStarts.clone();
            for (int index = 0; index < oneBasedColumns.length; index++)
                oneBasedColumns[index]++;
            for (int index = 0; index < oneBasedRows.length; index++)
                oneBasedRows[index]++;
            return new CsrMatrix(dimension, dimension, numeric,
                oneBasedColumns, oneBasedRows);
        }
    }
    private record PrecisionEntry(int term, double value) { }
    private record Evaluation(double[] beta, double[] fixedCovariance,
                              double[] randomModes,
                              double[] conditionalFitted,
                              double[] conditionalResiduals,
                              double logLikelihood) { }
}
