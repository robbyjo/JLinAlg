/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.matrix.CsrMatrix;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;

/**
 * First-order Laplace GLMM using sparse grouped designs and coefficient-space
 * precision matrices. Observation-scale covariance matrices are never formed.
 */
public final class SparseGlmmLaplace {
    private static final double MINIMUM_WORKING_WEIGHT = 1e-12;

    private SparseGlmmLaplace() { }

    /** Fits a grouped sparse GLMM with default controls and backend selection. */
    public static GlmmLaplaceResult fit(
            double[] response, double[][] fixedEffects,
            GlmFamily family, List<RandomEffectTerm> randomEffects) {
        if (response == null || fixedEffects == null
                || fixedEffects.length == 0 || fixedEffects[0] == null)
            throw new IllegalArgumentException(
                "response and fixed effects are required");
        return fit(response, MatrixOps.rowMajor(fixedEffects, response.length),
            response.length, fixedEffects[0].length, family, randomEffects,
            null, null, GlmmLaplaceOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    public static GlmmLaplaceResult fit(
            double[] response, double[] fixedEffects, int rows, int columns,
            GlmFamily family, List<RandomEffectTerm> randomEffects,
            double[] priorWeights, double[] offset,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        return fitWithPrecision(response, fixedEffects, rows, columns, family,
            randomEffects, null, priorWeights, offset, options, backendPolicy);
    }

    public static GlmmLaplaceResult fitWithPrecision(
            double[] response, double[] fixedEffects, int rows, int columns,
            GlmFamily family, List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            double[] priorWeights, double[] offset,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        try (Prepared prepared = prepareWithPrecision(rows, family,
                randomEffects, precisionBases, options, backendPolicy)) {
            return prepared.fit(response, fixedEffects, columns,
                priorWeights, offset);
        }
    }

    /** Prepares scan-owned sparse structure and one numerical factor per worker. */
    public static Prepared prepare(
            int rows, GlmFamily family,
            List<RandomEffectTerm> randomEffects,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        return prepareWithPrecision(rows, family, randomEffects, null,
            options, backendPolicy);
    }

    /** Prepares a scan with caller-supplied random-coefficient precisions. */
    public static Prepared prepareWithPrecision(
            int rows, GlmFamily family,
            List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        validate(rows, family, randomEffects, options, backendPolicy);
        return new Prepared(rows, family, randomEffects,
            precisions(randomEffects, precisionBases), options, backendPolicy);
    }

    /** Reusable scan state. Close it after all worker tasks have completed. */
    public static final class Prepared implements AutoCloseable {
        private final int rows;
        private final GlmFamily family;
        private final List<RandomEffectTerm> terms;
        private final GlmmLaplaceOptions options;
        private final CombinedDesign design;
        private final SparsePattern pattern;
        private final PrecisionData precisionData;
        private final BackendContext context;
        private final ComputeBackend backend;
        private final ConcurrentLinkedQueue<PreparedSparseCholesky> factors =
            new ConcurrentLinkedQueue<>();
        private final ThreadLocal<PreparedSparseCholesky> localFactor;
        private volatile double[] sharedLogVariances;
        private volatile boolean closed;

        private Prepared(
                int rows, GlmFamily family,
                List<RandomEffectTerm> terms,
                List<SparsePrecisionMatrix> precisions,
                GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
            this.rows = rows;
            this.family = family;
            this.terms = List.copyOf(terms);
            this.options = options;
            BackendContext selectedContext =
                BackendContext.select(backendPolicy);
            ComputeBackend selectedBackend = selectedContext.backend();
            CombinedDesign combinedDesign;
            PrecisionData preparedPrecision;
            SparsePattern preparedPattern;
            try {
                combinedDesign = combine(this.terms, rows);
                preparedPrecision = precisionData(
                    this.terms, precisions, combinedDesign, selectedBackend);
                preparedPattern = pattern(
                    combinedDesign, preparedPrecision);
            } catch (RuntimeException | Error failure) {
                selectedContext.close();
                throw failure;
            }
            context = selectedContext;
            backend = selectedBackend;
            design = combinedDesign;
            precisionData = preparedPrecision;
            pattern = preparedPattern;
            double[] initial = initial(options, terms.size());
            sharedLogVariances = initial;
            localFactor = ThreadLocal.withInitial(() -> {
                PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                    pattern.matrix(new double[rows], design, precisionData,
                        initial), MatrixTriangle.LOWER,
                    SparseOrdering.MINIMUM_DEGREE);
                factors.add(factor);
                return factor;
            });
        }

        public GlmmLaplaceResult fit(
                double[] response, double[] fixedEffects, int columns) {
            return fit(response, fixedEffects, columns, null, null);
        }

        public GlmmLaplaceResult fit(
                double[] response, double[] fixedEffects, int columns,
                double[] priorWeights, double[] offset) {
            if (closed) throw new IllegalStateException(
                "prepared sparse GLMM is closed");
            MatrixOps.validateModelData(response, fixedEffects, rows, columns);
            double[] weights = weights(priorWeights, rows);
            double[] offsets = offsets(offset, rows);
            for (int row = 0; row < rows; row++)
                family.validateResponse(response[row], weights[row]);
            double[] logVariances = sharedLogVariances.clone();
            PreparedSparseCholesky factor = localFactor.get();
            Mode best = mode(response, fixedEffects, columns, weights, offsets,
                logVariances, factor, null);
            int outerIterations = 0;
            double step = options.initialLogVarianceStep();
            boolean outerConverged = false;
            for (int sweep = 1;
                    sweep <= options.maximumOuterIterations(); sweep++) {
                outerIterations = sweep;
                boolean improved = false;
                for (int component = 0;
                        component < logVariances.length; component++) {
                    double original = logVariances[component];
                    double selected = original;
                    Mode coordinateBest = best;
                    for (double direction : new double[] {-1.0, 1.0}) {
                        double trial = clamp(original + direction * step);
                        if (trial == original) continue;
                        logVariances[component] = trial;
                        Mode candidate = mode(response, fixedEffects, columns,
                            weights, offsets, logVariances, factor, best);
                        if (candidate.laplaceLogLikelihood()
                                > coordinateBest.laplaceLogLikelihood()) {
                            coordinateBest = candidate;
                            selected = trial;
                        }
                    }
                    logVariances[component] = selected;
                    if (coordinateBest.laplaceLogLikelihood()
                            > best.laplaceLogLikelihood()
                            + options.relativeTolerance()
                            * (1.0 + Math.abs(
                                best.laplaceLogLikelihood()))) {
                        best = coordinateBest;
                        improved = true;
                    }
                }
                if (!improved) step *= 0.5;
                if (step <= options.relativeTolerance() * 10.0) {
                    outerConverged = true;
                    break;
                }
            }
            best = mode(response, fixedEffects, columns, weights, offsets,
                logVariances, factor, best);
            return result(best, columns, logVariances, outerIterations,
                outerConverged && best.converged());
        }

        /** Sets deterministic variance starts shared by subsequent workers. */
        public void warmStart(double... varianceComponents) {
            if (closed) throw new IllegalStateException(
                "prepared sparse GLMM is closed");
            if (varianceComponents == null
                    || varianceComponents.length != terms.size())
                throw new IllegalArgumentException(
                    "one variance is required per random-effect term");
            double[] values = new double[varianceComponents.length];
            for (int index = 0; index < values.length; index++) {
                double value = varianceComponents[index];
                if (!(value > 0.0) || !Double.isFinite(value))
                    throw new IllegalArgumentException(
                        "warm-start variances must be finite and positive");
                values[index] = Math.log(value);
            }
            sharedLogVariances = values;
        }

        private double clamp(double value) {
            return Math.max(Math.log(options.minimumVariance()),
                Math.min(Math.log(options.maximumVariance()), value));
        }

        private Mode mode(
                double[] response, double[] fixed, int fixedColumns,
                double[] priorWeights, double[] offsets,
                double[] logVariances, PreparedSparseCholesky factor,
                Mode start) {
            int randomColumns = design.columns();
            double[] beta = start == null
                ? new double[fixedColumns] : start.beta().clone();
            double[] random = start == null
                ? new double[randomColumns] : start.random().clone();
            if (start == null)
                initializeIntercept(response, fixed, fixedColumns,
                    offsets, beta);
            double[] means = new double[rows];
            double[] linear = new double[rows];
            double[] workingWeights = new double[rows];
            boolean converged = false;
            int iterations = 0;
            Schur solved = null;
            for (int iteration = 1;
                    iteration <= options.maximumModeIterations(); iteration++) {
                iterations = iteration;
                double[] workingResponse = working(
                    response, fixed, fixedColumns, beta, random,
                    priorWeights, offsets, means, linear, workingWeights);
                solved = solveWorking(fixed, fixedColumns, workingResponse,
                    workingWeights, logVariances, factor, false);
                double change = Math.max(relativeChange(beta, solved.beta()),
                    relativeChange(random, solved.random()));
                beta = solved.beta();
                random = solved.random();
                if (change <= options.relativeTolerance()) {
                    converged = true;
                    break;
                }
            }
            double[] workingResponse = working(response, fixed, fixedColumns,
                beta, random, priorWeights, offsets, means, linear,
                workingWeights);
            solved = solveWorking(fixed, fixedColumns, workingResponse,
                workingWeights, logVariances, factor, true);
            beta = solved.beta();
            random = solved.random();
            double conditional = 0.0;
            for (int row = 0; row < rows; row++) {
                double eta = offsets[row]
                    + fixedValue(fixed, row, fixedColumns, beta)
                    + design.rowProduct(row, random);
                linear[row] = eta;
                means[row] = family.inverseLink(eta);
                conditional += family.logLikelihood(response[row], means[row],
                    priorWeights[row], 1.0);
            }
            if (!Double.isFinite(conditional))
                throw new IllegalArgumentException(
                    "Laplace fitting requires a finite family likelihood");
            conditional -= 0.5 * precisionData.quadratic(
                random, logVariances);
            double laplace = conditional
                + 0.5 * precisionData.logDeterminant(logVariances)
                - 0.5 * factor.logDeterminant();
            return new Mode(beta, random, solved.fixedCovariance(),
                linear, means, laplace, iterations, converged);
        }

        private double[] working(
                double[] response, double[] fixed, int fixedColumns,
                double[] beta, double[] random, double[] priorWeights,
                double[] offsets, double[] means, double[] linear,
                double[] workingWeights) {
            double[] workingResponse = new double[rows];
            for (int row = 0; row < rows; row++) {
                double eta = offsets[row]
                    + fixedValue(fixed, row, fixedColumns, beta)
                    + design.rowProduct(row, random);
                linear[row] = eta;
                double mean = family.inverseLink(eta);
                means[row] = mean;
                double derivative = family.meanDerivative(eta);
                double variance = family.variance(mean);
                double weight = priorWeights[row] * derivative * derivative
                    / variance;
                workingWeights[row] = Math.max(
                    MINIMUM_WORKING_WEIGHT, weight);
                workingResponse[row] = eta - offsets[row]
                    + (response[row] - mean) / derivative;
            }
            return workingResponse;
        }

        private Schur solveWorking(
                double[] fixed, int fixedColumns, double[] workingResponse,
                double[] workingWeights, double[] logVariances,
                PreparedSparseCholesky factor, boolean covariance) {
            factor.refactor(pattern.matrix(workingWeights, design,
                precisionData, logVariances));
            int randomColumns = design.columns();
            double[] fixedInformation = new double[fixedColumns * fixedColumns];
            double[] cross = new double[fixedColumns * randomColumns];
            double[] fixedRight = new double[fixedColumns];
            double[] randomRight = new double[randomColumns];
            for (int row = 0; row < rows; row++) {
                double weight = workingWeights[row];
                double z = workingResponse[row];
                for (int left = 0; left < fixedColumns; left++) {
                    double xleft = fixed[row * fixedColumns + left];
                    fixedRight[left] += weight * xleft * z;
                    for (int right = 0; right <= left; right++)
                        fixedInformation[left * fixedColumns + right] +=
                            weight * xleft
                            * fixed[row * fixedColumns + right];
                    for (int index = design.rowStarts()[row];
                            index < design.rowStarts()[row + 1]; index++)
                        cross[left * randomColumns
                            + design.columnIndices()[index]] += weight * xleft
                            * design.values()[index];
                }
                for (int index = design.rowStarts()[row];
                        index < design.rowStarts()[row + 1]; index++)
                    randomRight[design.columnIndices()[index]] += weight
                        * design.values()[index] * z;
            }
            symmetrizeLower(fixedInformation, fixedColumns);
            double[] randomSolves = new double[
                randomColumns * (fixedColumns + 1)];
            for (int row = 0; row < randomColumns; row++) {
                randomSolves[row * (fixedColumns + 1)] = randomRight[row];
                for (int column = 0; column < fixedColumns; column++)
                    randomSolves[row * (fixedColumns + 1) + column + 1] =
                        cross[column * randomColumns + row];
            }
            factor.solveInPlace(randomSolves, fixedColumns + 1);
            double[] schur = fixedInformation.clone();
            double[] schurRight = fixedRight.clone();
            for (int left = 0; left < fixedColumns; left++)
                for (int randomColumn = 0;
                        randomColumn < randomColumns; randomColumn++) {
                    double value = cross[left * randomColumns + randomColumn];
                    schurRight[left] -= value
                        * randomSolves[randomColumn * (fixedColumns + 1)];
                    for (int right = 0; right < fixedColumns; right++)
                        schur[left * fixedColumns + right] -= value
                            * randomSolves[randomColumn * (fixedColumns + 1)
                                + right + 1];
                }
            symmetrize(schur, fixedColumns);
            CholeskyFactor fixedFactor = backend.dpotrf(schur, fixedColumns);
            double[] beta = fixedFactor.solve(schurRight);
            double[] random = new double[randomColumns];
            for (int row = 0; row < randomColumns; row++) {
                double value = randomSolves[row * (fixedColumns + 1)];
                for (int column = 0; column < fixedColumns; column++)
                    value -= randomSolves[row * (fixedColumns + 1)
                        + column + 1] * beta[column];
                random[row] = value;
            }
            double[] fixedCovariance = covariance
                ? fixedFactor.solve(MatrixOps.identity(fixedColumns),
                    fixedColumns) : null;
            return new Schur(beta, random, fixedCovariance);
        }

        private GlmmLaplaceResult result(
                Mode mode, int fixedColumns, double[] logVariances,
                int outerIterations, boolean converged) {
            double[] standardErrors = new double[fixedColumns];
            for (int column = 0; column < fixedColumns; column++)
                standardErrors[column] = Math.sqrt(Math.max(0.0,
                    mode.fixedCovariance()[column * fixedColumns + column]));
            double[] variances = new double[logVariances.length];
            for (int index = 0; index < variances.length; index++)
                variances[index] = Math.exp(logVariances[index]);
            List<String> names = new ArrayList<>(terms.size());
            Map<String, double[]> predictors = new LinkedHashMap<>();
            for (int term = 0; term < terms.size(); term++) {
                RandomEffectTerm value = terms.get(term);
                names.add(value.name());
                double[] coefficients = Arrays.copyOfRange(mode.random(),
                    design.termStarts()[term],
                    design.termStarts()[term] + value.coefficients());
                predictors.put(value.name(), multiply(value, coefficients));
            }
            return new GlmmLaplaceResult(family.name(), names, variances,
                AssociationStatistics.normal(mode.beta(), standardErrors),
                mode.fixedCovariance(), predictors, mode.linear(), mode.means(),
                mode.laplaceLogLikelihood(), outerIterations,
                mode.iterations(), converged);
        }

        private void initializeIntercept(
                double[] response, double[] fixed, int columns,
                double[] offsets, double[] beta) {
            int intercept = -1;
            for (int column = 0; column < columns; column++) {
                boolean allOne = true;
                for (int row = 0; row < rows; row++)
                    if (Math.abs(fixed[row * columns + column] - 1.0)
                            > 1e-12) {
                        allOne = false;
                        break;
                    }
                if (allOne) {
                    intercept = column;
                    break;
                }
            }
            if (intercept < 0) return;
            double total = 0.0;
            for (double value : response) total += value;
            double mean = Math.max(1e-7,
                Math.min(1.0 - 1e-7, total / rows));
            double offsetMean = 0.0;
            for (double value : offsets) offsetMean += value;
            beta[intercept] = family.link(mean) - offsetMean / rows;
        }

        @Override public void close() {
            if (closed) return;
            closed = true;
            localFactor.remove();
            for (PreparedSparseCholesky factor : factors) factor.close();
            factors.clear();
            context.close();
        }
    }

    private static SparsePattern pattern(
            CombinedDesign design, PrecisionData precision) {
        int dimension = design.columns();
        TreeMap<Long, Boolean> keys = new TreeMap<>();
        for (int row = 0; row < design.rows(); row++)
            for (int left = design.rowStarts()[row];
                    left < design.rowStarts()[row + 1]; left++)
                for (int right = design.rowStarts()[row];
                        right <= left; right++) {
                    int r = Math.max(design.columnIndices()[left],
                        design.columnIndices()[right]);
                    int c = Math.min(design.columnIndices()[left],
                        design.columnIndices()[right]);
                    keys.put((long) r * dimension + c, Boolean.TRUE);
                }
        for (int index = 0; index < precision.rows().length; index++)
            keys.put((long) precision.rows()[index] * dimension
                + precision.columns()[index], Boolean.TRUE);
        int[] rowStarts = new int[dimension + 1];
        int[] columns = new int[keys.size()];
        int position = 0;
        int currentRow = 0;
        for (long key : keys.keySet()) {
            int row = (int) (key / dimension);
            while (currentRow < row) rowStarts[++currentRow] = position;
            columns[position++] = (int) (key % dimension);
        }
        while (currentRow < dimension) rowStarts[++currentRow] = position;
        return new SparsePattern(dimension, rowStarts, columns);
    }

    private static PrecisionData precisionData(
            List<RandomEffectTerm> terms,
            List<SparsePrecisionMatrix> precisions,
            CombinedDesign design, ComputeBackend backend) {
        List<Integer> rows = new ArrayList<>();
        List<Integer> columns = new ArrayList<>();
        List<Double> values = new ArrayList<>();
        List<Integer> termIndices = new ArrayList<>();
        double[] logDeterminants = new double[terms.size()];
        for (int term = 0; term < terms.size(); term++) {
            SparsePrecisionMatrix precision = precisions.get(term);
            int[] starts = precision.rowStarts();
            int[] cols = precision.columnIndices();
            double[] numeric = precision.values();
            int lowerCount = 0;
            for (int row = 0; row < precision.dimension(); row++)
                for (int index = starts[row]; index < starts[row + 1]; index++)
                    if (cols[index] <= row) lowerCount++;
            int[] lowerRows = new int[precision.dimension() + 1];
            int[] lowerColumns = new int[lowerCount];
            double[] lowerValues = new double[lowerCount];
            int lowerPosition = 0;
            for (int row = 0; row < precision.dimension(); row++) {
                lowerRows[row] = lowerPosition + 1;
                for (int index = starts[row]; index < starts[row + 1]; index++) {
                    int column = cols[index];
                    if (column <= row) {
                        lowerColumns[lowerPosition] = column + 1;
                        lowerValues[lowerPosition++] = numeric[index];
                        rows.add(design.termStarts()[term] + row);
                        columns.add(design.termStarts()[term] + column);
                        values.add(numeric[index]);
                        termIndices.add(term);
                    }
                }
            }
            lowerRows[precision.dimension()] = lowerPosition + 1;
            CsrMatrix lower = new CsrMatrix(precision.dimension(),
                precision.dimension(), lowerValues, lowerColumns, lowerRows);
            logDeterminants[term] = backend.dcsrpotrf(lower,
                MatrixTriangle.LOWER,
                SparseOrdering.MINIMUM_DEGREE).logDeterminant();
        }
        return new PrecisionData(toIntArray(rows), toIntArray(columns),
            toDoubleArray(values), toIntArray(termIndices), logDeterminants,
            terms.stream().mapToInt(RandomEffectTerm::coefficients).toArray());
    }

    private static CombinedDesign combine(
            List<RandomEffectTerm> terms, int rows) {
        int[] termStarts = new int[terms.size()];
        int columns = 0;
        int nonzeros = 0;
        for (int term = 0; term < terms.size(); term++) {
            termStarts[term] = columns;
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
                    double entry = value.values()[index];
                    if (entry != 0.0) {
                        columnIndices[position] = termStarts[term]
                            + value.columnIndices()[index];
                        values[position++] = entry;
                    }
                }
            }
        }
        rowStarts[rows] = position;
        return new CombinedDesign(rows, columns, rowStarts,
            Arrays.copyOf(columnIndices, position),
            Arrays.copyOf(values, position), termStarts);
    }

    private static TermData termData(RandomEffectTerm term) {
        if (term.sparse())
            return new TermData(term.rowPointers(), term.columnIndices(),
                term.sparseValues());
        double[] dense = term.design();
        int count = 0;
        for (double value : dense) if (value != 0.0) count++;
        int[] starts = new int[term.observations() + 1];
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

    private static double[] multiply(
            RandomEffectTerm term, double[] coefficients) {
        double[] result = new double[term.observations()];
        if (term.sparse()) {
            int[] starts = term.rowPointers();
            int[] columns = term.columnIndices();
            double[] values = term.sparseValues();
            for (int row = 0; row < result.length; row++)
                for (int index = starts[row]; index < starts[row + 1]; index++)
                    result[row] += values[index] * coefficients[columns[index]];
        } else {
            double[] design = term.design();
            for (int row = 0; row < result.length; row++)
                for (int column = 0; column < coefficients.length; column++)
                    result[row] += design[row * coefficients.length + column]
                        * coefficients[column];
        }
        return result;
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
                "one precision matrix is required per random-effect term");
        for (int index = 0; index < terms.size(); index++)
            if (supplied.get(index) == null
                    || supplied.get(index).dimension()
                        != terms.get(index).coefficients())
                throw new IllegalArgumentException(
                    "precision dimensions must match random coefficients");
        return List.copyOf(supplied);
    }

    private static void validate(
            int rows, GlmFamily family,
            List<RandomEffectTerm> terms,
            GlmmLaplaceOptions options, BackendPolicy policy) {
        if (rows < 1 || family == null || terms == null || terms.isEmpty()
                || options == null || policy == null)
            throw new IllegalArgumentException(
                "rows, family, random effects, controls, and backend are required");
        java.util.HashSet<String> names = new java.util.HashSet<>();
        for (RandomEffectTerm term : terms)
            if (term == null || term.observations() != rows
                    || !names.add(term.name()))
                throw new IllegalArgumentException(
                    "random effects need unique names and matching rows");
        double[] supplied = options.initialVariances();
        if (supplied != null && supplied.length != terms.size())
            throw new IllegalArgumentException(
                "one initial variance is required per random-effect term");
    }

    private static double[] initial(
            GlmmLaplaceOptions options, int components) {
        double[] supplied = options.initialVariances();
        double[] result = new double[components];
        if (supplied == null) return result;
        for (int index = 0; index < components; index++)
            result[index] = Math.log(supplied[index]);
        return result;
    }

    private static double[] weights(double[] supplied, int rows) {
        if (supplied == null) {
            double[] result = new double[rows];
            Arrays.fill(result, 1.0);
            return result;
        }
        if (supplied.length != rows)
            throw new IllegalArgumentException(
                "prior-weight length must match observations");
        double[] result = supplied.clone();
        for (double value : result)
            if (!(value > 0.0) || !Double.isFinite(value))
                throw new IllegalArgumentException(
                    "prior weights must be finite and positive");
        return result;
    }

    private static double[] offsets(double[] supplied, int rows) {
        if (supplied == null) return new double[rows];
        if (supplied.length != rows)
            throw new IllegalArgumentException(
                "offset length must match observations");
        return MatrixOps.finiteCopy(supplied, "offset");
    }

    private static double fixedValue(
            double[] fixed, int row, int columns, double[] beta) {
        double value = 0.0;
        for (int column = 0; column < columns; column++)
            value += fixed[row * columns + column] * beta[column];
        return value;
    }

    private static double relativeChange(double[] previous, double[] next) {
        double maximum = 0.0;
        for (int index = 0; index < previous.length; index++)
            maximum = Math.max(maximum, Math.abs(next[index] - previous[index])
                / (1.0 + Math.abs(previous[index])));
        return maximum;
    }

    private static void symmetrizeLower(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++)
            for (int column = 0; column < row; column++)
                matrix[column * dimension + row] =
                    matrix[row * dimension + column];
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

    private static int[] toIntArray(List<Integer> values) {
        int[] result = new int[values.size()];
        for (int index = 0; index < result.length; index++)
            result[index] = values.get(index);
        return result;
    }

    private static double[] toDoubleArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int index = 0; index < result.length; index++)
            result[index] = values.get(index);
        return result;
    }

    private record TermData(
            int[] rowStarts, int[] columnIndices, double[] values) { }

    private record CombinedDesign(
            int rows, int columns, int[] rowStarts,
            int[] columnIndices, double[] values, int[] termStarts) {
        double rowProduct(int row, double[] coefficients) {
            double result = 0.0;
            for (int index = rowStarts[row]; index < rowStarts[row + 1]; index++)
                result += values[index] * coefficients[columnIndices[index]];
            return result;
        }
    }

    private record SparsePattern(
            int dimension, int[] rowStarts, int[] columnIndices) {
        CsrMatrix matrix(
                double[] weights, CombinedDesign design,
                PrecisionData precision, double[] logVariances) {
            double[] numeric = new double[columnIndices.length];
            for (int observation = 0;
                    observation < design.rows(); observation++) {
                double weight = weights[observation];
                for (int left = design.rowStarts()[observation];
                        left < design.rowStarts()[observation + 1]; left++)
                    for (int right = design.rowStarts()[observation];
                            right <= left; right++) {
                        int row = Math.max(design.columnIndices()[left],
                            design.columnIndices()[right]);
                        int column = Math.min(design.columnIndices()[left],
                            design.columnIndices()[right]);
                        int position = Arrays.binarySearch(columnIndices,
                            rowStarts[row], rowStarts[row + 1], column);
                        numeric[position] += weight * design.values()[left]
                            * design.values()[right];
                    }
            }
            for (int index = 0; index < precision.values().length; index++) {
                int row = precision.rows()[index];
                int position = Arrays.binarySearch(columnIndices,
                    rowStarts[row], rowStarts[row + 1],
                    precision.columns()[index]);
                numeric[position] += precision.values()[index]
                    / Math.exp(logVariances[precision.terms()[index]]);
            }
            int[] csrRows = rowStarts.clone();
            int[] csrColumns = columnIndices.clone();
            for (int index = 0; index < csrRows.length; index++) csrRows[index]++;
            for (int index = 0; index < csrColumns.length; index++)
                csrColumns[index]++;
            return new CsrMatrix(dimension, dimension, numeric,
                csrColumns, csrRows);
        }
    }

    private record PrecisionData(
            int[] rows, int[] columns, double[] values, int[] terms,
            double[] logDeterminants, int[] dimensions) {
        double logDeterminant(double[] logVariances) {
            double result = 0.0;
            for (int term = 0; term < logDeterminants.length; term++)
                result += logDeterminants[term]
                    - dimensions[term] * logVariances[term];
            return result;
        }

        double quadratic(double[] random, double[] logVariances) {
            double result = 0.0;
            for (int index = 0; index < values.length; index++) {
                int row = rows[index];
                int column = columns[index];
                double contribution = values[index] * random[row]
                    * random[column]
                    / Math.exp(logVariances[terms[index]]);
                result += row == column ? contribution : 2.0 * contribution;
            }
            return result;
        }
    }

    private record Schur(
            double[] beta, double[] random, double[] fixedCovariance) { }

    private record Mode(
            double[] beta, double[] random, double[] fixedCovariance,
            double[] linear, double[] means, double laplaceLogLikelihood,
            int iterations, boolean converged) { }
}
