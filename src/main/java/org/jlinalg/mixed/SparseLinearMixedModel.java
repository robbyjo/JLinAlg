/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.Optimization;
import jdistlib.math.opt.OptimizationResult;
import jdistlib.matrix.CsrMatrix;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceEstimation;

/** Sparse-equation REML/ML for independent Gaussian random-effect terms. */
public final class SparseLinearMixedModel {
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    private static final double INVALID_OBJECTIVE = Double.MAX_VALUE / 16.0;
    private static final int MAXIMUM_PEV_COEFFICIENTS = 256;
    private static final double INACTIVE_STATIC_VARIANCE_RATIO = 1e-6;

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
        try (Prepared prepared = prepareWithPrecision(rows, randomEffects,
                precisionBases, options, backendPolicy)) {
            return prepared.fit(response, fixedEffects, columns);
        }
    }

    /**
     * Fits with a caller-owned backend, allowing high-throughput callers to
     * retain one compute context across changing model structures.
     */
    public static SparseLinearMixedModelResult fitWithBackend(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            ComputeBackend backend,
            BackendProvenance provenance) {
        validate(randomEffects, rows, options);
        if (backend == null || provenance == null)
            throw new IllegalArgumentException(
                "backend and provenance are required");
        if (options.degreesOfFreedomMethod()
                != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION)
            throw new IllegalArgumentException(
                "sparse LMM currently supports residual-approximation DF; "
                    + "use the dense reference fitter for Satterthwaite or Kenward-Roger");
        try (Prepared prepared = new Prepared(rows, randomEffects,
                precisions(randomEffects, null), options,
                backend, provenance, null)) {
            return prepared.fit(response, fixedEffects, columns);
        }
    }

    /** Prepares reusable random-effect structure for repeated fixed-design fits. */
    public static Prepared prepare(
            int rows,
            List<RandomEffectTerm> randomEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        return prepareWithPrecision(
            rows, randomEffects, null, options, backendPolicy);
    }

    /** Prepares reusable structure with caller-supplied precision bases. */
    public static Prepared prepareWithPrecision(
            int rows,
            List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        validate(randomEffects, rows, options, backendPolicy);
        if (options.degreesOfFreedomMethod()
                != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION)
            throw new IllegalArgumentException(
                "sparse LMM currently supports residual-approximation DF; "
                    + "use the dense reference fitter for Satterthwaite or Kenward-Roger");
        BackendContext context = BackendContext.select(backendPolicy);
        try {
            return new Prepared(rows, randomEffects,
                precisions(randomEffects, precisionBases), options,
                context.backend(), context.provenance(), context);
        } catch (RuntimeException | Error failure) {
            context.close();
            throw failure;
        }
    }

    /**
     * Prepares static random effects while reserving one dense random-effect
     * term whose design changes between fits. The sparse equation structure
     * and symbolic factorization are retained across all fits.
     */
    public static PreparedChangingTerm prepareChangingTermWithBackend(
            int rows,
            List<RandomEffectTerm> staticTerms,
            List<SparsePrecisionMatrix> staticPrecisions,
            String changingTermName,
            int changingCoefficients,
            List<String> changingCoefficientNames,
            RemlOptions options,
            ComputeBackend backend,
            BackendProvenance provenance) {
        if (staticTerms == null || staticPrecisions == null
                || staticTerms.size() != staticPrecisions.size()
                || changingTermName == null || changingTermName.isBlank()
                || changingCoefficients < 1
                || changingCoefficientNames == null
                || changingCoefficientNames.size() != changingCoefficients
                || options == null || backend == null || provenance == null) {
            throw new IllegalArgumentException(
                "static terms, changing-term metadata, options, and backend are required");
        }
        for (int term = 0; term < staticTerms.size(); term++) {
            RandomEffectTerm value = staticTerms.get(term);
            SparsePrecisionMatrix precision = staticPrecisions.get(term);
            if (value == null || value.observations() != rows
                    || precision == null
                    || precision.dimension() != value.coefficients()
                    || changingTermName.equals(value.name())) {
                throw new IllegalArgumentException(
                    "static random effects and precisions must match observations");
            }
        }
        if (options.degreesOfFreedomMethod()
                != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION) {
            throw new IllegalArgumentException(
                "prepared changing-term LMM supports residual-approximation DF");
        }
        return new PreparedChangingTerm(rows, staticTerms, staticPrecisions,
            changingTermName, changingCoefficients, changingCoefficientNames,
            options, backend, provenance);
    }

    /**
     * Reuses sparse cross-products, precision log determinants, the backend
     * context, and variance starts across a model scan.
     */
    public static final class Prepared implements AutoCloseable {
        private final int rows;
        private final List<RandomEffectTerm> terms;
        private final List<SparsePrecisionMatrix> precisions;
        private final RemlOptions options;
        private final CombinedDesign design;
        private final SparsePattern pattern;
        private final BackendContext context;
        private final ComputeBackend backend;
        private final BackendProvenance provenance;
        private final double[] precisionLogDeterminants;
        private final double[] lower;
        private final double[] upper;
        private final int maximumEvaluations;
        private final ConcurrentLinkedQueue<PreparedSparseCholesky> factors;
        private final ThreadLocal<PreparedSparseCholesky> localFactor;
        private volatile double[] sharedLogVarianceRatios;
        private volatile boolean closed;

        private Prepared(
                int rows, List<RandomEffectTerm> terms,
                List<SparsePrecisionMatrix> precisions,
                RemlOptions options, ComputeBackend backend,
                BackendProvenance provenance, BackendContext context) {
            this.rows = rows;
            this.terms = List.copyOf(terms);
            this.precisions = precisions;
            this.options = options;
            this.design = combine(this.terms, rows);
            this.pattern = crossProductPattern(
                design, this.terms, precisions);
            this.context = context;
            this.backend = backend;
            this.provenance = provenance;
            this.precisionLogDeterminants = precisionLogDeterminants(
                precisions, backend);
            int components = terms.size();
            this.lower = new double[components];
            this.upper = new double[components];
            double maximumRatio = options.maximumVariance()
                / options.minimumVariance();
            Arrays.fill(lower, -Math.log(maximumRatio));
            Arrays.fill(upper, Math.log(maximumRatio));
            this.maximumEvaluations = Math.max(100,
                options.maximumIterations() * 10);
            if (provenance.selectedBackend().startsWith("cholmod+")) {
                factors = new ConcurrentLinkedQueue<>();
                localFactor = ThreadLocal.withInitial(() -> {
                    PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                        pattern.matrix(relativeVariances(
                            new double[terms.size()])),
                        MatrixTriangle.LOWER,
                        SparseOrdering.MINIMUM_DEGREE);
                    factors.add(factor);
                    return factor;
                });
            } else {
                factors = null;
                localFactor = null;
            }
        }

        /** Fits one response/design while retaining reusable sparse state. */
        public SparseLinearMixedModelResult fit(
                double[] response, double[] fixedEffects, int columns) {
            if (closed) throw new IllegalStateException(
                "prepared sparse LMM is closed");
            MatrixOps.validateModelData(
                response, fixedEffects, rows, columns);
            double[] shared = sharedLogVarianceRatios;
            double[] initial = shared != null ? shared.clone()
                : initialLogVarianceRatios(response, terms.size(), options);
            if (localFactor != null) {
                return fit(response, fixedEffects, columns, initial,
                    localFactor.get());
            }
            try (PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                    pattern.matrix(relativeVariances(initial)),
                    MatrixTriangle.LOWER, SparseOrdering.MINIMUM_DEGREE)) {
                return fit(response, fixedEffects, columns, initial, factor);
            }
        }

        private SparseLinearMixedModelResult fit(
                double[] response, double[] fixedEffects, int columns,
                double[] initial, PreparedSparseCholesky factor) {
            Objective objective = new Objective(response, fixedEffects,
                rows, columns, terms, design, pattern, factor,
                precisionLogDeterminants,
                options.varianceEstimation(), backend,
                design.columns() <= MAXIMUM_PEV_COEFFICIENTS ? 0 : -1,
                design.columns() <= MAXIMUM_PEV_COEFFICIENTS
                    ? design.columns() : 0);
            OptimizationResult optimized = optimize(initial, objective);
            Evaluation fitted = objective.evaluate(optimized.mX);
            int factorNonzeroCount = factor.factorNonzeroCount();
            double[] variances = fitted.variances();
            double[] standardErrors = new double[columns];
            for (int column = 0; column < columns; column++)
                standardErrors[column] = Math.sqrt(Math.max(0.0,
                    fitted.fixedCovariance()[column * columns + column]));
            double degrees = rows - columns - 1.0;
            if (!(degrees > 0.0)) throw new IllegalArgumentException(
                "sparse LMM requires positive denominator DF");
            AssociationStatistics association = AssociationStatistics.studentT(
                fitted.beta(), standardErrors, degrees,
                DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
            List<RandomEffectEstimates> estimates = estimates(
                terms, design.termStarts(), variances, fitted.randomModes(),
                fitted.randomPredictionErrorVariances());
            List<String> names = new ArrayList<>(terms.size() + 1);
            for (RandomEffectTerm term : terms) names.add(term.name());
            names.add("residual");
            return new SparseLinearMixedModelResult(names, variances,
                association, fitted.fixedCovariance(), estimates,
                fitted.conditionalFitted(), fitted.conditionalResiduals(),
                fitted.logLikelihood(), options.varianceEstimation(),
                optimized.numFunctionCalls,
                optimized.numFunctionCalls < maximumEvaluations,
                design.columns(), pattern.values().length,
                factorNonzeroCount, provenance);
        }

        /** Uses one deterministic fitted variance start for subsequent fits. */
        public void warmStart(double... varianceComponents) {
            if (closed) throw new IllegalStateException(
                "prepared sparse LMM is closed");
            if (varianceComponents == null
                    || varianceComponents.length != terms.size() + 1)
                throw new IllegalArgumentException(
                    "warm start requires random terms plus residual variance");
            double residual = varianceComponents[terms.size()];
            if (!(residual > 0.0) || !Double.isFinite(residual))
                throw new IllegalArgumentException(
                    "warm-start residual variance must be finite and positive");
            double[] ratios = new double[terms.size()];
            for (int term = 0; term < terms.size(); term++) {
                double value = varianceComponents[term];
                if (!(value > 0.0) || !Double.isFinite(value))
                    throw new IllegalArgumentException(
                        "warm-start variances must be finite and positive");
                ratios[term] = Math.log(value / residual);
            }
            sharedLogVarianceRatios = ratios;
        }

        private OptimizationResult optimize(
                double[] initial, Objective objective) {
            try {
                OptimizationResult primary = Bobyqa.bobyqa(initial,
                    lower, upper, parameters -> objective.value(parameters),
                    Math.min(2 * initial.length + 1,
                        (initial.length + 1) * (initial.length + 2) / 2),
                    Math.min(1.0, options.maximumLogVarianceStep()),
                    Math.max(1e-7, options.relativeTolerance()),
                    maximumEvaluations, true);
                return localFactor == null || terms.size() == 1 ? primary
                    : coordinateOptimize(primary.mX, objective,
                        4, Math.max(1e-4, options.relativeTolerance()), 30);
            } catch (ArithmeticException | ArrayIndexOutOfBoundsException failure) {
                return coordinateOptimize(initial, objective);
            }
        }

        private OptimizationResult coordinateOptimize(
                double[] initial, Objective objective) {
            return coordinateOptimize(initial, objective, 12,
                Math.max(1e-7, options.relativeTolerance()),
                maximumEvaluations);
        }

        private OptimizationResult coordinateOptimize(
                double[] initial, Objective objective, int iterations,
                double tolerance, int searchEvaluations) {
            double[] point = initial.clone();
            int[] calls = {0};
            double best = objective.value(point);
            calls[0]++;
            for (int iteration = 0; iteration < iterations
                    && calls[0] < maximumEvaluations; iteration++) {
                double maximumChange = 0.0;
                for (int coordinate = 0; coordinate < point.length; coordinate++) {
                    final int selected = coordinate;
                    double before = point[selected];
                    double searchLower = Math.max(lower[selected], before - 8.0);
                    double searchUpper = Math.min(upper[selected], before + 8.0);
                    double optimum = Optimization.optimize(value -> {
                        point[selected] = value;
                        calls[0]++;
                        return objective.value(point);
                    }, searchLower, searchUpper, tolerance,
                        Math.max(20, Math.min(searchEvaluations,
                            maximumEvaluations - calls[0])));
                    point[selected] = optimum;
                    best = objective.value(point);
                    calls[0]++;
                    maximumChange = Math.max(
                        maximumChange, Math.abs(optimum - before));
                    if (calls[0] >= maximumEvaluations) break;
                }
                if (maximumChange <= Math.max(
                        1e-6, options.relativeTolerance())) break;
            }
            return new OptimizationResult(
                point.clone(), best, calls[0], true);
        }

        @Override
        public void close() {
            if (closed) return;
            closed = true;
            if (localFactor != null) {
                localFactor.remove();
                for (PreparedSparseCholesky factor : factors) factor.close();
                factors.clear();
            }
            if (context != null) context.close();
        }
    }

    /**
     * Prepared sparse solver for scans where the final dense random-effect
     * design changes but its coefficient count stays fixed.
     */
    public static final class PreparedChangingTerm implements AutoCloseable {
        private final int rows;
        private final List<RandomEffectTerm> staticTerms;
        private final List<SparsePrecisionMatrix> staticPrecisions;
        private final String changingName;
        private final int changingCoefficients;
        private final List<String> changingCoefficientNames;
        private final RemlOptions options;
        private final ComputeBackend backend;
        private final BackendProvenance provenance;
        private final List<RandomEffectTerm> templateTerms;
        private final List<SparsePrecisionMatrix> precisions;
        private final SparsePattern templatePattern;
        private final double[] precisionLogDeterminants;
        private final double[] lower;
        private final double[] upper;
        private final int maximumEvaluations;
        private final ConcurrentLinkedQueue<PreparedSparseCholesky> factors;
        private final ThreadLocal<PreparedSparseCholesky> localFactor;
        private volatile double[] sharedLogVarianceRatios;
        private volatile boolean closed;

        private PreparedChangingTerm(
                int rows,
                List<RandomEffectTerm> staticTerms,
                List<SparsePrecisionMatrix> staticPrecisions,
                String changingName,
                int changingCoefficients,
                List<String> changingCoefficientNames,
                RemlOptions options,
                ComputeBackend backend,
                BackendProvenance provenance) {
            this.rows = rows;
            this.staticTerms = List.copyOf(staticTerms);
            this.staticPrecisions = List.copyOf(staticPrecisions);
            this.changingName = changingName;
            this.changingCoefficients = changingCoefficients;
            this.changingCoefficientNames =
                List.copyOf(changingCoefficientNames);
            this.options = options;
            this.backend = backend;
            this.provenance = provenance;
            double[] placeholder = new double[rows * changingCoefficients];
            Arrays.fill(placeholder, 1.0);
            RandomEffectTerm changing = RandomEffectTerm.of(changingName,
                placeholder, rows, changingCoefficients,
                changingCoefficientNames);
            List<RandomEffectTerm> terms =
                new ArrayList<>(staticTerms.size() + 1);
            terms.addAll(staticTerms);
            terms.add(changing);
            this.templateTerms = List.copyOf(terms);
            List<SparsePrecisionMatrix> allPrecisions =
                new ArrayList<>(staticPrecisions.size() + 1);
            allPrecisions.addAll(staticPrecisions);
            allPrecisions.add(
                SparsePrecisionMatrix.identity(changingCoefficients));
            this.precisions = List.copyOf(allPrecisions);
            CombinedDesign template = combine(templateTerms, rows, true);
            this.templatePattern = crossProductPattern(
                template, templateTerms, precisions);
            this.precisionLogDeterminants =
                precisionLogDeterminants(precisions, backend);
            int components = templateTerms.size();
            this.lower = new double[components];
            this.upper = new double[components];
            double maximumRatio = options.maximumVariance()
                / options.minimumVariance();
            Arrays.fill(lower, -Math.log(maximumRatio));
            Arrays.fill(upper, Math.log(maximumRatio));
            this.maximumEvaluations = Math.max(100,
                options.maximumIterations() * 10);
            if (provenance.selectedBackend().startsWith("cholmod+")) {
                factors = new ConcurrentLinkedQueue<>();
                localFactor = ThreadLocal.withInitial(() -> {
                    PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                        templatePattern.matrix(relativeVariances(
                            new double[templateTerms.size()])),
                        MatrixTriangle.LOWER,
                        SparseOrdering.MINIMUM_DEGREE);
                    factors.add(factor);
                    return factor;
                });
            } else {
                factors = null;
                localFactor = null;
            }
        }

        /**
         * Fits one response/fixed design and one row-major changing random
         * design while retaining all static sparse state.
         */
        public SparseLinearMixedModelResult fit(
                double[] response,
                double[] fixedEffects,
                int fixedColumns,
                double[] changingDesign) {
            if (closed) throw new IllegalStateException(
                "prepared changing-term LMM is closed");
            MatrixOps.validateModelData(
                response, fixedEffects, rows, fixedColumns);
            if (changingDesign == null
                    || changingDesign.length
                        != rows * changingCoefficients) {
                throw new IllegalArgumentException(
                    "changing random-effect design dimensions are invalid");
            }
            RandomEffectTerm changing = RandomEffectTerm.of(changingName,
                changingDesign, rows, changingCoefficients,
                changingCoefficientNames);
            List<RandomEffectTerm> terms =
                new ArrayList<>(staticTerms.size() + 1);
            terms.addAll(staticTerms);
            terms.add(changing);
            terms = List.copyOf(terms);
            CombinedDesign design = combine(terms, rows, true);
            SparsePattern pattern =
                crossProductPattern(design, terms, precisions);
            requireSameStructure(templatePattern, pattern);
            double[] shared = sharedLogVarianceRatios;
            double[] initial = shared != null ? shared.clone()
                : initialLogVarianceRatios(
                    response, terms.size(), options);
            if (localFactor != null) {
                return fit(response, fixedEffects, fixedColumns, terms,
                    design, pattern, initial, localFactor.get());
            }
            try (PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                    pattern.matrix(relativeVariances(initial)),
                    MatrixTriangle.LOWER,
                    SparseOrdering.MINIMUM_DEGREE)) {
                return fit(response, fixedEffects, fixedColumns, terms,
                    design, pattern, initial, factor);
            }
        }

        private SparseLinearMixedModelResult fit(
                double[] response,
                double[] fixedEffects,
                int fixedColumns,
                List<RandomEffectTerm> terms,
                CombinedDesign design,
                SparsePattern pattern,
                double[] initial,
                PreparedSparseCholesky factor) {
            int changingStart = design.termStarts()[terms.size() - 1];
            Objective objective = new Objective(response, fixedEffects,
                rows, fixedColumns, terms, design, pattern, factor,
                precisionLogDeterminants,
                options.varianceEstimation(), backend,
                changingStart, changingCoefficients);
            OptimizationResult optimized = optimize(initial, objective);
            Evaluation fitted = objective.evaluate(optimized.mX);
            double[] variances = fitted.variances();
            double[] standardErrors = new double[fixedColumns];
            for (int column = 0; column < fixedColumns; column++) {
                standardErrors[column] = Math.sqrt(Math.max(0.0,
                    fitted.fixedCovariance()[
                        column * fixedColumns + column]));
            }
            double degrees = rows - fixedColumns - 1.0;
            AssociationStatistics association = AssociationStatistics.studentT(
                fitted.beta(), standardErrors, degrees,
                DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
            List<RandomEffectEstimates> estimates = estimates(
                terms, design.termStarts(), variances, fitted.randomModes(),
                fitted.randomPredictionErrorVariances());
            List<String> names = new ArrayList<>(terms.size() + 1);
            for (RandomEffectTerm term : terms) names.add(term.name());
            names.add("residual");
            return new SparseLinearMixedModelResult(names, variances,
                association, fitted.fixedCovariance(), estimates,
                fitted.conditionalFitted(), fitted.conditionalResiduals(),
                fitted.logLikelihood(), options.varianceEstimation(),
                optimized.numFunctionCalls,
                optimized.numFunctionCalls < maximumEvaluations,
                design.columns(), pattern.values().length,
                factor.factorNonzeroCount(), provenance);
        }

        private OptimizationResult optimize(
                double[] initial, Objective objective) {
            ScreenedStart screened =
                screenInitial(initial, objective);
            initial = screened.point();
            int activeCount = 0;
            int[] active = new int[initial.length];
            for (int component = 0;
                    component < initial.length; component++) {
                if (component == initial.length - 1
                        || initial[component] >= Math.log(
                            INACTIVE_STATIC_VARIANCE_RATIO)) {
                    active[activeCount++] = component;
                }
            }
            if (activeCount < initial.length) {
                int[] selected = Arrays.copyOf(active, activeCount);
                OptimizationResult reduced =
                    optimizeActive(initial, selected, objective);
                return addCalls(reactivateImprovingComponents(
                    selected, reduced, objective), screened.calls());
            }
            return addCalls(
                optimizeAll(initial, objective), screened.calls());
        }

        private ScreenedStart screenInitial(
                double[] initial, Objective objective) {
            double[] point = initial.clone();
            double best = objective.value(point);
            int calls = 1;
            for (int component = 0;
                    component < point.length; component++) {
                if (component != point.length - 1
                        && point[component] >= Math.log(
                            INACTIVE_STATIC_VARIANCE_RATIO)) {
                    continue;
                }
                double original = point[component];
                double selected = original;
                for (double step : new double[] {
                        -16.0, -12.0, -8.0, -4.0,
                        4.0, 8.0, 12.0, 16.0}) {
                    point[component] = Math.max(lower[component],
                        Math.min(upper[component], original + step));
                    double value = objective.value(point);
                    calls++;
                    if (value < best) {
                        best = value;
                        selected = point[component];
                    }
                }
                point[component] = selected;
            }
            return new ScreenedStart(point, calls);
        }

        private static OptimizationResult addCalls(
                OptimizationResult result, int calls) {
            return new OptimizationResult(result.mX, result.mF,
                result.numFunctionCalls + calls, true);
        }

        private OptimizationResult optimizeAll(
                double[] initial, Objective objective) {
            try {
                OptimizationResult primary = Bobyqa.bobyqa(initial,
                    lower, upper, parameters -> objective.value(parameters),
                    Math.min(2 * initial.length + 1,
                        (initial.length + 1) * (initial.length + 2) / 2),
                    Math.min(1.0, options.maximumLogVarianceStep()),
                    Math.max(1e-4, options.relativeTolerance()),
                    maximumEvaluations, true);
                return primary;
            } catch (ArithmeticException
                    | ArrayIndexOutOfBoundsException failure) {
                return coordinateOptimize(initial, objective, 2,
                    Math.max(1e-4, options.relativeTolerance()),
                    12);
            }
        }

        private OptimizationResult reactivateImprovingComponents(
                int[] active,
                OptimizationResult reduced,
                Objective objective) {
            boolean[] selected =
                new boolean[reduced.mX.length];
            for (int component : active) selected[component] = true;
            double[] bestPoint = reduced.mX.clone();
            double best = reduced.mF;
            int calls = 0;
            for (int component = 0;
                    component < selected.length - 1; component++) {
                if (selected[component]) continue;
                double original = bestPoint[component];
                for (double step : new double[] {4.0, 8.0, 12.0, 16.0}) {
                    double candidate = Math.min(
                        upper[component], original + step);
                    bestPoint[component] = candidate;
                    double value = objective.value(bestPoint);
                    calls++;
                    if (value + 1e-4 < best) {
                        best = value;
                        original = candidate;
                    }
                }
                bestPoint[component] = original;
            }
            if (best + 1e-4 < reduced.mF) {
                OptimizationResult expanded =
                    optimizeAll(bestPoint, objective);
                return new OptimizationResult(expanded.mX,
                    expanded.mF,
                    reduced.numFunctionCalls + calls
                        + expanded.numFunctionCalls,
                    true);
            }
            return new OptimizationResult(reduced.mX,
                reduced.mF,
                reduced.numFunctionCalls + calls, true);
        }

        private OptimizationResult optimizeActive(
                double[] fixed, int[] active, Objective objective) {
            double[] initial = new double[active.length];
            double[] activeLower = new double[active.length];
            double[] activeUpper = new double[active.length];
            for (int index = 0; index < active.length; index++) {
                initial[index] = fixed[active[index]];
                activeLower[index] = lower[active[index]];
                activeUpper[index] = upper[active[index]];
            }
            if (active.length == 1) {
                int[] calls = {0};
                double optimum = Optimization.optimize(value -> {
                    calls[0]++;
                    return objective.value(expand(
                        fixed, active, new double[] {value}));
                }, activeLower[0], activeUpper[0],
                    Math.max(1e-4, options.relativeTolerance()),
                    maximumEvaluations);
                double[] result = expand(
                    fixed, active, new double[] {optimum});
                return new OptimizationResult(result,
                    objective.value(result), calls[0] + 1, true);
            }
            try {
                OptimizationResult optimized = Bobyqa.bobyqa(initial,
                    activeLower, activeUpper,
                    parameters -> objective.value(
                        expand(fixed, active, parameters)),
                    Math.min(2 * active.length + 1,
                        (active.length + 1) * (active.length + 2) / 2),
                    Math.min(1.0, options.maximumLogVarianceStep()),
                    Math.max(1e-4, options.relativeTolerance()),
                    maximumEvaluations, true);
                return new OptimizationResult(
                    expand(fixed, active, optimized.mX),
                    optimized.mF, optimized.numFunctionCalls, true);
            } catch (ArithmeticException
                    | ArrayIndexOutOfBoundsException failure) {
                return coordinateOptimize(fixed, objective, 2,
                    Math.max(1e-4, options.relativeTolerance()), 12);
            }
        }

        private static double[] expand(
                double[] fixed, int[] active, double[] values) {
            double[] result = fixed.clone();
            for (int index = 0; index < active.length; index++) {
                result[active[index]] = values[index];
            }
            return result;
        }

        private OptimizationResult coordinateOptimize(
                double[] initial, Objective objective, int iterations,
                double tolerance, int searchEvaluations) {
            double[] point = initial.clone();
            int[] calls = {0};
            double best = objective.value(point);
            calls[0]++;
            for (int iteration = 0; iteration < iterations
                    && calls[0] < maximumEvaluations; iteration++) {
                double maximumChange = 0.0;
                for (int coordinate = 0;
                        coordinate < point.length; coordinate++) {
                    final int selected = coordinate;
                    double before = point[selected];
                    double searchLower = Math.max(
                        lower[selected], before - 8.0);
                    double searchUpper = Math.min(
                        upper[selected], before + 8.0);
                    double optimum = Optimization.optimize(value -> {
                        point[selected] = value;
                        calls[0]++;
                        return objective.value(point);
                    }, searchLower, searchUpper, tolerance,
                        Math.max(8, Math.min(searchEvaluations,
                            maximumEvaluations - calls[0])));
                    point[selected] = optimum;
                    best = objective.value(point);
                    calls[0]++;
                    maximumChange = Math.max(maximumChange,
                        Math.abs(optimum - before));
                    if (calls[0] >= maximumEvaluations) break;
                }
                if (maximumChange <= Math.max(
                        1e-6, options.relativeTolerance())) break;
            }
            return new OptimizationResult(
                point.clone(), best, calls[0], true);
        }

        /** Uses fitted variance components as the shared start for later fits. */
        public void warmStart(double... varianceComponents) {
            if (closed) throw new IllegalStateException(
                "prepared changing-term LMM is closed");
            if (varianceComponents == null
                    || varianceComponents.length
                        != templateTerms.size() + 1) {
                throw new IllegalArgumentException(
                    "warm start requires random terms plus residual variance");
            }
            double residual =
                varianceComponents[templateTerms.size()];
            if (!(residual > 0.0) || !Double.isFinite(residual)) {
                throw new IllegalArgumentException(
                    "warm-start residual variance must be positive");
            }
            double[] ratios = new double[templateTerms.size()];
            for (int term = 0; term < ratios.length; term++) {
                double value = varianceComponents[term];
                if (!(value > 0.0) || !Double.isFinite(value)) {
                    throw new IllegalArgumentException(
                        "warm-start variances must be positive");
                }
                ratios[term] = Math.log(value / residual);
            }
            sharedLogVarianceRatios = ratios;
        }

        public int changingCoefficients() {
            return changingCoefficients;
        }

        @Override
        public void close() {
            if (closed) return;
            closed = true;
            if (localFactor != null) {
                localFactor.remove();
                for (PreparedSparseCholesky factor : factors) {
                    factor.close();
                }
                factors.clear();
            }
        }

        private record ScreenedStart(double[] point, int calls) { }
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
        private final double responseSquared;
        private final double[] fixedResponse;
        private final double[] fixedCross;
        private final double[] randomResponse;
        private final double[] randomFixedCross;

        Objective(double[] response, double[] fixed, int rows, int columns,
                List<RandomEffectTerm> terms, CombinedDesign design,
                SparsePattern pattern, PreparedSparseCholesky factor,
                double[] precisionLogDeterminants,
                VarianceEstimation estimation, ComputeBackend backend,
                int pevStart, int pevCount) {
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
            this.responseSquared = backend.ddot(
                rows, response, 0, 1, response, 0, 1);
            this.fixedResponse = new double[columns];
            backend.dgemv(jdistlib.accelerator.MatrixTranspose.TRANSPOSE,
                rows, columns, 1.0, fixed, response,
                0.0, fixedResponse);
            this.fixedCross = MatrixOps.transposeMultiply(
                backend, fixed, rows, columns, fixed, columns);
            this.randomResponse = transposeMultiply(design, response);
            this.randomFixedCross =
                transposeMultiply(design, fixed, columns);
            this.pevStart = pevStart;
            this.pevCount = pevCount;
        }

        double value(double[] parameters) {
            try {
                Evaluation result = evaluate(parameters, false);
                if (!Double.isFinite(result.logLikelihood()))
                    return INVALID_OBJECTIVE;
                return -result.logLikelihood();
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return INVALID_OBJECTIVE;
            }
        }

        Evaluation evaluate(double[] logVarianceRatios) {
            return evaluate(logVarianceRatios, true);
        }

        private Evaluation evaluate(
                double[] logVarianceRatios, boolean materializeResult) {
            double[] ratios = exp(logVarianceRatios);
            factor.refactor(pattern.matrix(relativeVariances(logVarianceRatios)));
            double[] solvedRandomResponse = factor.solve(randomResponse);
            double[] solvedRandomFixed =
                factor.solve(randomFixedCross, columns);
            double[] information = fixedCross.clone();
            for (int left = 0; left < columns; left++) {
                for (int right = 0; right < columns; right++) {
                    double removed = 0.0;
                    for (int random = 0;
                            random < design.columns(); random++) {
                        removed += randomFixedCross[
                            random * columns + left]
                            * solvedRandomFixed[random * columns + right];
                    }
                    information[left * columns + right] -= removed;
                }
            }
            symmetrize(information, columns);
            CholeskyFactor fixedFactor = backend.dpotrf(information, columns);
            double[] rightSide = fixedResponse.clone();
            for (int fixedColumn = 0;
                    fixedColumn < columns; fixedColumn++) {
                double removed = 0.0;
                for (int random = 0;
                        random < design.columns(); random++) {
                    removed += randomFixedCross[
                        random * columns + fixedColumn]
                        * solvedRandomResponse[random];
                }
                rightSide[fixedColumn] -= removed;
            }
            double[] beta = fixedFactor.solve(rightSide);
            double quadratic = responseSquared
                - backend.ddot(design.columns(), randomResponse, 0, 1,
                    solvedRandomResponse, 0, 1)
                - backend.ddot(columns, rightSide, 0, 1, beta, 0, 1);
            double logDeterminant = factor.logDeterminant();
            for (int term = 0; term < terms.size(); term++)
                logDeterminant += terms.get(term).coefficients()
                    * Math.log(ratios[term])
                    - precisionLogDeterminants[term];
            boolean restricted = estimation == VarianceEstimation.REML;
            double scaleDegrees = restricted ? rows - columns : rows;
            double residualVariance = quadratic / scaleDegrees;
            double logLikelihood = -0.5 * (
                scaleDegrees * (LOG_TWO_PI + Math.log(residualVariance))
                    + logDeterminant
                    + (restricted ? fixedFactor.logDeterminant() : 0.0)
                    + quadratic / residualVariance);
            if (!materializeResult)
                return new Evaluation(beta, null, null, null, null, null,
                    logLikelihood, null);
            double[] fixedCovariance = fixedFactor.solve(
                MatrixOps.identity(columns), columns);
            for (int index = 0; index < fixedCovariance.length; index++)
                fixedCovariance[index] *= residualVariance;
            symmetrize(fixedCovariance, columns);
            double[] residual = MatrixOps.subtract(response,
                MatrixOps.multiply(backend, fixed, rows, columns, beta));
            double[] randomModes = solvedRandomResponse.clone();
            for (int random = 0; random < design.columns(); random++) {
                for (int fixedColumn = 0;
                        fixedColumn < columns; fixedColumn++) {
                    randomModes[random] -= solvedRandomFixed[
                        random * columns + fixedColumn] * beta[fixedColumn];
                }
            }
            double[] conditionalFitted = MatrixOps.multiply(
                backend, fixed, rows, columns, beta);
            double[] randomFitted = multiply(design, randomModes);
            for (int row = 0; row < rows; row++)
                conditionalFitted[row] += randomFitted[row];
            double[] variances = new double[ratios.length + 1];
            for (int term = 0; term < ratios.length; term++)
                variances[term] = ratios[term] * residualVariance;
            variances[ratios.length] = residualVariance;
            double[] randomPredictionErrorVariances =
                randomPredictionErrorVariances(
                    fixedCovariance, residualVariance,
                    solvedRandomFixed);
            return new Evaluation(beta, fixedCovariance, randomModes,
                randomPredictionErrorVariances,
                conditionalFitted,
                MatrixOps.subtract(response, conditionalFitted),
                logLikelihood, variances);
        }

        private double[] randomPredictionErrorVariances(
                double[] fixedCovariance,
                double residualVariance,
                double[] adjustedCross) {
            int randomColumns = design.columns();
            double[] result = new double[randomColumns];
            Arrays.fill(result, Double.NaN);
            if (pevStart < 0 || pevCount == 0) return result;
            double[] rightHandSides =
                new double[randomColumns * pevCount];
            for (int selected = 0; selected < pevCount; selected++) {
                rightHandSides[
                    (pevStart + selected) * pevCount + selected] = 1.0;
            }
            double[] inverseColumns =
                factor.solve(rightHandSides, pevCount);
            for (int selected = 0; selected < pevCount; selected++) {
                int random = pevStart + selected;
                double value = residualVariance
                    * inverseColumns[random * pevCount + selected];
                for (int left = 0; left < columns; left++) {
                    double leftValue =
                        adjustedCross[random * columns + left];
                    for (int right = 0; right < columns; right++) {
                        value += leftValue
                            * fixedCovariance[left * columns + right]
                            * adjustedCross[random * columns + right];
                    }
                }
                result[random] = Math.max(0.0, value);
            }
            return result;
        }

        private final int pevStart;
        private final int pevCount;

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
        return combine(terms, rows, false);
    }

    private static CombinedDesign combine(
            List<RandomEffectTerm> terms, int rows,
            boolean preserveLastDenseTerm) {
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
                    if (value.values()[index] != 0.0
                            || (preserveLastDenseTerm
                                && term == terms.size() - 1)) {
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

    private static void requireSameStructure(
            SparsePattern expected, SparsePattern actual) {
        if (expected.dimension() != actual.dimension()
                || !Arrays.equals(
                    expected.rowStarts(), actual.rowStarts())
                || !Arrays.equals(expected.columnIndices(),
                    actual.columnIndices())) {
            throw new IllegalStateException(
                "changing random-effect design altered sparse structure");
        }
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

    private static double[] transposeMultiply(
            CombinedDesign design, double[] matrix, int columns) {
        double[] result = new double[design.columns() * columns];
        for (int row = 0; row < design.rows(); row++) {
            for (int index = design.rowStarts()[row];
                    index < design.rowStarts()[row + 1]; index++) {
                int random = design.columnIndices()[index];
                double value = design.values()[index];
                for (int column = 0; column < columns; column++) {
                    result[random * columns + column] +=
                        value * matrix[row * columns + column];
                }
            }
        }
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
            double[] variances, double[] modes,
            double[] predictionErrorVariances) {
        List<RandomEffectEstimates> result = new ArrayList<>(terms.size());
        for (int term = 0; term < terms.size(); term++) {
            RandomEffectTerm value = terms.get(term);
            double[] selected = Arrays.copyOfRange(modes, starts[term],
                starts[term] + value.coefficients());
            double[] pev = predictionErrorVariances == null
                ? new double[value.coefficients()]
                : Arrays.copyOfRange(predictionErrorVariances, starts[term],
                    starts[term] + value.coefficients());
            if (predictionErrorVariances == null) Arrays.fill(pev, Double.NaN);
            result.add(new RandomEffectEstimates(value.name(),
                value.coefficientNames(), variances[term], selected,
                pev));
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

    private static double[] initialLogVarianceRatios(
            double[] response, int randomTerms, RemlOptions options) {
        double[] absolute = initialLogVariances(
            response, randomTerms + 1, options);
        double residual = absolute[randomTerms];
        double[] result = new double[randomTerms];
        for (int term = 0; term < randomTerms; term++)
            result[term] = absolute[term] - residual;
        return result;
    }

    private static double[] relativeVariances(double[] logVarianceRatios) {
        double[] result = Arrays.copyOf(
            exp(logVarianceRatios), logVarianceRatios.length + 1);
        result[logVarianceRatios.length] = 1.0;
        return result;
    }

    private static void validate(
            List<RandomEffectTerm> terms, int rows,
            RemlOptions options, BackendPolicy backendPolicy) {
        validate(terms, rows, options);
        if (backendPolicy == null)
            throw new IllegalArgumentException(
                "random effects, options, and backend are required");
    }

    private static void validate(
            List<RandomEffectTerm> terms, int rows, RemlOptions options) {
        if (terms == null || terms.isEmpty() || options == null)
            throw new IllegalArgumentException(
                "random effects and options are required");
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
                              double[] randomPredictionErrorVariances,
                              double[] conditionalFitted,
                              double[] conditionalResiduals,
                              double logLikelihood,
                              double[] variances) { }
}
