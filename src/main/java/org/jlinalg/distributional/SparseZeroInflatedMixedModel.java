/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.random.RandomGenerator;
import java.util.random.RandomGeneratorFactory;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.OptimizationResult;
import jdistlib.matrix.CsrMatrix;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;

/**
 * Frequentist ZIP and NB2-ZINB mixed models with sparse Gaussian random
 * effects in either distributional process. Random coefficients are integrated
 * by first-order Laplace approximation; BOBYQA optimizes only the non-random
 * parameters.
 */
public final class SparseZeroInflatedMixedModel {
    private static final double INVALID_OBJECTIVE = 1e100;
    private static final double MINIMUM_PROBABILITY = 1e-12;

    private SparseZeroInflatedMixedModel() { }

    /** Fits a zero-inflated Poisson mixed model. */
    public static ZeroInflatedMixedResult fitPoisson(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            double[] countOffset,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, countFixed, countColumns,
            zeroFixed, zeroColumns, null, 0,
            countRandomEffects, precisionBases, List.of(), List.of(),
            List.of(), countOffset,
            CountFamily.POISSON, options, backendPolicy);
    }

    /** Fits a zero-inflated NB2 mixed model with a fixed-effects size model. */
    public static ZeroInflatedMixedResult fitNegativeBinomial(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            double[] dispersionFixed, int dispersionColumns,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            double[] countOffset,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, countFixed, countColumns,
            zeroFixed, zeroColumns, dispersionFixed, dispersionColumns,
            countRandomEffects, precisionBases, List.of(), List.of(),
            List.of(), countOffset,
            CountFamily.NEGATIVE_BINOMIAL, options, backendPolicy);
    }

    /** Fits a ZIP model with independent and correlated two-process effects. */
    public static ZeroInflatedMixedResult fitPoisson(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> countPrecisionBases,
            List<RandomEffectTerm> zeroRandomEffects,
            List<SparsePrecisionMatrix> zeroPrecisionBases,
            List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
            double[] countOffset,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, countFixed, countColumns,
            zeroFixed, zeroColumns, null, 0,
            countRandomEffects, countPrecisionBases,
            zeroRandomEffects, zeroPrecisionBases, correlatedEffects,
            countOffset, CountFamily.POISSON, options, backendPolicy);
    }

    /** Fits an NB2-ZINB model with independent and correlated effects. */
    public static ZeroInflatedMixedResult fitNegativeBinomial(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            double[] dispersionFixed, int dispersionColumns,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> countPrecisionBases,
            List<RandomEffectTerm> zeroRandomEffects,
            List<SparsePrecisionMatrix> zeroPrecisionBases,
            List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
            double[] countOffset,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, countFixed, countColumns,
            zeroFixed, zeroColumns, dispersionFixed, dispersionColumns,
            countRandomEffects, countPrecisionBases,
            zeroRandomEffects, zeroPrecisionBases, correlatedEffects,
            countOffset, CountFamily.NEGATIVE_BINOMIAL, options,
            backendPolicy);
    }

    /** Prepares a ZIP random structure for repeated, concurrent fits. */
    public static Prepared preparePoisson(
            int rows,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> countPrecisionBases,
            List<RandomEffectTerm> zeroRandomEffects,
            List<SparsePrecisionMatrix> zeroPrecisionBases,
            List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        return new Prepared(rows, CountFamily.POISSON, countRandomEffects,
            countPrecisionBases, zeroRandomEffects, zeroPrecisionBases,
            correlatedEffects, options, backendPolicy);
    }

    /** Prepares an NB2-ZINB random structure for repeated, concurrent fits. */
    public static Prepared prepareNegativeBinomial(
            int rows,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> countPrecisionBases,
            List<RandomEffectTerm> zeroRandomEffects,
            List<SparsePrecisionMatrix> zeroPrecisionBases,
            List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        return new Prepared(rows, CountFamily.NEGATIVE_BINOMIAL,
            countRandomEffects, countPrecisionBases, zeroRandomEffects,
            zeroPrecisionBases, correlatedEffects, options, backendPolicy);
    }

    private static ZeroInflatedMixedResult fit(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            double[] dispersionFixed, int dispersionColumns,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> countPrecisionBases,
            List<RandomEffectTerm> zeroRandomEffects,
            List<SparsePrecisionMatrix> zeroPrecisionBases,
            List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
            double[] countOffset,
            CountFamily family,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        validate(response, countFixed, countColumns, zeroFixed, zeroColumns,
            dispersionFixed, dispersionColumns, countRandomEffects,
            countPrecisionBases, zeroRandomEffects, zeroPrecisionBases,
            correlatedEffects, countOffset, family, options, backendPolicy);
        try (Prepared prepared = new Prepared(response.length, family,
                countRandomEffects, countPrecisionBases, zeroRandomEffects,
                zeroPrecisionBases, correlatedEffects, options, backendPolicy)) {
            return prepared.fit(response, countFixed, countColumns, zeroFixed,
                zeroColumns, dispersionFixed, dispersionColumns, countOffset);
        }
    }

    /** Reusable symbolic sparse analysis with one numerical factor per worker. */
    public static final class Prepared implements AutoCloseable {
        private final int rows;
        private final CountFamily family;
        private final ZeroInflatedMixedOptions options;
        private final RandomStructure structure;
        private final CombinedDesign design;
        private final PrecisionData precision;
        private final SparsePattern pattern;
        private final BackendContext context;
        private final ComputeBackend backend;
        private final ConcurrentLinkedQueue<PreparedSparseCholesky> factors =
            new ConcurrentLinkedQueue<>();
        private final ThreadLocal<PreparedSparseCholesky> localFactor;
        private volatile ExecutorService gradientExecutor;
        private volatile boolean closed;

        private Prepared(
                int rows, CountFamily family,
                List<RandomEffectTerm> countRandomEffects,
                List<SparsePrecisionMatrix> countPrecisionBases,
                List<RandomEffectTerm> zeroRandomEffects,
                List<SparsePrecisionMatrix> zeroPrecisionBases,
                List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
                ZeroInflatedMixedOptions options,
                BackendPolicy backendPolicy) {
            validateStructure(rows, countRandomEffects, countPrecisionBases,
                zeroRandomEffects, zeroPrecisionBases, correlatedEffects,
                options, backendPolicy);
            this.rows = rows;
            this.family = family;
            this.options = options;
            List<SparsePrecisionMatrix> countPrecisions =
                countPrecisionBases == null ? identities(countRandomEffects)
                    : List.copyOf(countPrecisionBases);
            List<SparsePrecisionMatrix> zeroPrecisions =
                zeroPrecisionBases == null ? identities(zeroRandomEffects)
                    : List.copyOf(zeroPrecisionBases);
            BackendContext selected = BackendContext.select(backendPolicy);
            ComputeBackend selectedBackend = selected.backend();
            RandomStructure preparedStructure;
            CombinedDesign preparedDesign;
            PrecisionData preparedPrecision;
            SparsePattern preparedPattern;
            try {
                preparedStructure = randomStructure(countRandomEffects,
                    countPrecisions, zeroRandomEffects, zeroPrecisions,
                    correlatedEffects);
                preparedDesign = combine(preparedStructure.terms(),
                    preparedStructure.countTerms(), rows);
                preparedPrecision = precisionData(preparedStructure,
                    preparedDesign, selectedBackend);
                preparedPattern = pattern(preparedDesign, preparedPrecision);
            } catch (RuntimeException | Error failure) {
                selected.close();
                throw failure;
            }
            context = selected;
            backend = selectedBackend;
            structure = preparedStructure;
            design = preparedDesign;
            precision = preparedPrecision;
            pattern = preparedPattern;
            localFactor = ThreadLocal.withInitial(() -> {
                double[] zeroWeights = new double[rows];
                double[] covariance = new double[
                    structure.varianceCount() + structure.correlationCount()];
                Arrays.fill(covariance, 0, structure.varianceCount(),
                    Math.log(0.5));
                PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                    pattern.matrix(zeroWeights, zeroWeights, zeroWeights,
                        design, precision, covariance, 0.0),
                    MatrixTriangle.LOWER, SparseOrdering.MINIMUM_DEGREE);
                factors.add(factor);
                return factor;
            });
        }

        /** Fits one response/design set using this prepared sparse structure. */
        public ZeroInflatedMixedResult fit(
                double[] response,
                double[] countFixed, int countColumns,
                double[] zeroFixed, int zeroColumns,
                double[] dispersionFixed, int dispersionColumns,
                double[] countOffset) {
            return fitInternal(response, countFixed, countColumns, zeroFixed,
                zeroColumns, dispersionFixed, dispersionColumns, countOffset,
                false);
        }

        /**
         * Fits and then computes an observed numerical Hessian of the marginal
         * Laplace objective. BOBYQA itself remains derivative-free.
         */
        public ZeroInflatedMixedResult fitWithInference(
                double[] response,
                double[] countFixed, int countColumns,
                double[] zeroFixed, int zeroColumns,
                double[] dispersionFixed, int dispersionColumns,
                double[] countOffset) {
            return fitInternal(response, countFixed, countColumns, zeroFixed,
                zeroColumns, dispersionFixed, dispersionColumns, countOffset,
                true);
        }

        /** Profiles one optimizer-scale parameter while refitting all others. */
        public ZeroInflatedProfile profile(
                ZeroInflatedMixedResult fitted,
                double[] response,
                double[] countFixed, int countColumns,
                double[] zeroFixed, int zeroColumns,
                double[] dispersionFixed, int dispersionColumns,
                double[] countOffset, int parameterIndex,
                double[] grid) {
            if (closed) throw new IllegalStateException("prepared model is closed");
            if (fitted == null || grid == null || grid.length == 0
                    || parameterIndex < 0
                    || parameterIndex >= fitted.outerParameterEstimates().length)
                throw new IllegalArgumentException(
                    "fitted result, parameter index, and profile grid are required");
            List<String> expectedNames = outerParameterNames(countColumns,
                zeroColumns, family == CountFamily.NEGATIVE_BINOMIAL
                    ? dispersionColumns : 0, structure);
            if (!expectedNames.equals(fitted.outerParameterNames()))
                throw new IllegalArgumentException(
                    "fitted result does not match this prepared structure");
            validateData(response, countFixed, countColumns, zeroFixed,
                zeroColumns, dispersionFixed, dispersionColumns, countOffset,
                family, rows);
            double[] offsets = countOffset == null ? new double[rows]
                : MatrixOps.finiteCopy(countOffset, "countOffset");
            double[] start = fitted.outerParameterEstimates();
            double[] lower = lower(start.length, countColumns, zeroColumns,
                dispersionColumns, structure.varianceCount(),
                structure.correlationCount(), family, options);
            double[] upper = upper(start.length, countColumns, zeroColumns,
                dispersionColumns, structure.varianceCount(),
                structure.correlationCount(), family, options);
            Objective objective = new Objective(response, countFixed,
                countColumns, zeroFixed, zeroColumns, dispersionFixed,
                dispersionColumns, offsets, family, design, precision,
                pattern, localFactor.get(), options);
            double[] logLikelihoods = new double[grid.length];
            boolean[] converged = new boolean[grid.length];
            for (int point = 0; point < grid.length; point++) {
                if (!Double.isFinite(grid[point]) || grid[point] < lower[parameterIndex]
                        || grid[point] > upper[parameterIndex])
                    throw new IllegalArgumentException(
                        "profile grid lies outside optimizer bounds");
                double[] reducedStart = remove(start, parameterIndex);
                double[] reducedLower = remove(lower, parameterIndex);
                double[] reducedUpper = remove(upper, parameterIndex);
                final double fixedValue = grid[point];
                java.util.function.ToDoubleFunction<double[]> function = reduced ->
                    objective.coldValue(insert(reduced, parameterIndex,
                        fixedValue));
                OptimizationResult optimized = optimizeReduced(reducedStart,
                    reducedLower, reducedUpper, function, options);
                double[] full = insert(optimized.mX == null
                    ? reducedStart : optimized.mX, parameterIndex, fixedValue);
                logLikelihoods[point] = -objective.coldValue(full);
                converged[point] = optimized.numFunctionCalls
                    < options.maximumOuterEvaluations();
            }
            return new ZeroInflatedProfile(
                fitted.outerParameterNames().get(parameterIndex), grid,
                logLikelihoods, converged);
        }

        /**
         * Runs a deterministic conditional parametric bootstrap. Failed fits
         * retain NaN rows; callers can parallelize calls on this prepared object.
         */
        public ZeroInflatedBootstrapResult parametricBootstrap(
                ZeroInflatedMixedResult fitted,
                double[] countFixed, int countColumns,
                double[] zeroFixed, int zeroColumns,
                double[] dispersionFixed, int dispersionColumns,
                double[] countOffset, int replicates, long seed) {
            if (fitted == null || replicates < 1)
                throw new IllegalArgumentException(
                    "fitted result and positive replicate count are required");
            List<String> expectedNames = outerParameterNames(countColumns,
                zeroColumns, family == CountFamily.NEGATIVE_BINOMIAL
                    ? dispersionColumns : 0, structure);
            if (!expectedNames.equals(fitted.outerParameterNames()))
                throw new IllegalArgumentException(
                    "fitted result does not match this prepared structure");
            int parameters = fitted.outerParameterNames().size();
            double[] estimates = new double[replicates * parameters];
            Arrays.fill(estimates, Double.NaN);
            int successes = 0;
            RandomGenerator random = RandomGeneratorFactory
                .of("L64X128MixRandom").create(seed);
            for (int replicate = 0; replicate < replicates; replicate++) {
                double[] simulated = simulate(fitted, random, family);
                try {
                    ZeroInflatedMixedResult refitted = fit(simulated,
                        countFixed, countColumns, zeroFixed, zeroColumns,
                        dispersionFixed, dispersionColumns, countOffset);
                    double[] values = refitted.outerParameterEstimates();
                    System.arraycopy(values, 0, estimates,
                        replicate * parameters, parameters);
                    successes++;
                } catch (IllegalArgumentException | IllegalStateException
                        exception) {
                    // A failed bootstrap fit is represented explicitly by NaN.
                }
            }
            return new ZeroInflatedBootstrapResult(
                fitted.outerParameterNames(), estimates, replicates, successes);
        }

        private ZeroInflatedMixedResult fitInternal(
                double[] response,
                double[] countFixed, int countColumns,
                double[] zeroFixed, int zeroColumns,
                double[] dispersionFixed, int dispersionColumns,
                double[] countOffset, boolean inference) {
            if (closed)
                throw new IllegalStateException("prepared model is closed");
            validateData(response, countFixed, countColumns, zeroFixed,
                zeroColumns, dispersionFixed, dispersionColumns, countOffset,
                family, rows);
            double[] offsets = countOffset == null ? new double[rows]
                : MatrixOps.finiteCopy(countOffset, "countOffset");
            double[] initial = initial(response, countFixed, countColumns,
                zeroFixed, zeroColumns, dispersionFixed, dispersionColumns,
                structure.varianceCount(), structure.correlationCount(),
                family, options);
            double[] lower = lower(initial.length, countColumns, zeroColumns,
                dispersionColumns, structure.varianceCount(),
                structure.correlationCount(), family, options);
            double[] upper = upper(initial.length, countColumns, zeroColumns,
                dispersionColumns, structure.varianceCount(),
                structure.correlationCount(), family, options);
            PreparedSparseCholesky factor = localFactor.get();
                Objective objective = new Objective(response,
                    countFixed, countColumns, zeroFixed, zeroColumns,
                    dispersionFixed, dispersionColumns, offsets, family,
                    design, precision, pattern, factor, options);
                int interpolationPoints = Math.min(2 * initial.length + 1,
                    (initial.length + 1) * (initial.length + 2) / 2);
                OptimizationResult optimized;
                int priorEvaluations = 0;
                boolean outerConverged;
                boolean derivativeFree = options.outerOptimizer()
                        == ZeroInflatedOuterOptimizer.BOBYQA
                    || (options.outerOptimizer()
                            == ZeroInflatedOuterOptimizer.AUTO
                        && design.columns() > 128);
                if (derivativeFree) {
                    optimized = bobyqa(initial, lower, upper, objective,
                        interpolationPoints, options);
                    outerConverged = optimized.numFunctionCalls
                        < options.maximumOuterEvaluations();
                } else {
                    ThreadLocal<Objective> workerObjectives =
                        ThreadLocal.withInitial(() ->
                            objective.fork(localFactor.get()));
                    optimized = boundedBfgs(initial, lower, upper, objective,
                        options, initial.length >= 8
                            && options.maximumGradientThreads() > 1
                                ? gradientExecutor() : null,
                        workerObjectives);
                    outerConverged = optimized.numFunctionCalls
                        < options.maximumOuterEvaluations();
                    if (options.outerOptimizer()
                            == ZeroInflatedOuterOptimizer.AUTO
                            && !outerConverged) {
                        priorEvaluations = optimized.numFunctionCalls;
                        optimized = bobyqa(optimized.mX, lower, upper,
                            objective, interpolationPoints, options);
                        outerConverged = optimized.numFunctionCalls
                            < options.maximumOuterEvaluations();
                    }
                }
                double[] parameters = optimized.mX == null
                    ? initial : optimized.mX;
                InferenceData inferenceData = inference
                    ? numericalInference(objective, parameters, lower, upper,
                        countColumns, zeroColumns,
                        family == CountFamily.NEGATIVE_BINOMIAL
                            ? dispersionColumns : 0)
                    : InferenceData.unavailable();
                Evaluation fitted = objective.evaluate(parameters);
                return result(fitted, parameters, countColumns, zeroColumns,
                    dispersionColumns, family, structure, design,
                    pattern, factor.factorNonzeroCount(),
                    priorEvaluations + optimized.numFunctionCalls,
                    outerConverged, inferenceData,
                    options);
        }

        public int randomCoefficientCount() { return design.columns(); }
        public int sparseEquationNonzeroCount() {
            return pattern.columnIndices().length;
        }
        /** Number of lazily created worker-local numerical factors. */
        public int numericFactorCount() { return factors.size(); }

        private synchronized ExecutorService gradientExecutor() {
            if (gradientExecutor == null) {
                int threads = Math.min(options.maximumGradientThreads(),
                    Runtime.getRuntime().availableProcessors());
                gradientExecutor = Executors.newFixedThreadPool(threads,
                    runnable -> {
                        Thread thread = new Thread(runnable,
                            "zero-inflated-gradient");
                        thread.setDaemon(true);
                        return thread;
                    });
            }
            return gradientExecutor;
        }

        @Override public void close() {
            if (closed) return;
            closed = true;
            if (gradientExecutor != null) gradientExecutor.shutdownNow();
            for (PreparedSparseCholesky factor : factors) factor.close();
            factors.clear();
            localFactor.remove();
            context.close();
        }
    }

    private static OptimizationResult coordinateOptimize(
            double[] initial, double[] lower, double[] upper,
            Objective objective, ZeroInflatedMixedOptions options) {
        double[] point = initial.clone();
        double best = objective.value(point);
        int calls = 1;
        double step = options.initialTrustRadius();
        while (calls < options.maximumOuterEvaluations()
                && step > options.relativeTolerance()) {
            boolean improved = false;
            for (int parameter = 0; parameter < point.length
                    && calls < options.maximumOuterEvaluations(); parameter++) {
                double original = point[parameter];
                for (int direction : new int[] {-1, 1}) {
                    point[parameter] = Math.max(lower[parameter],
                        Math.min(upper[parameter], original + direction * step));
                    double candidate = objective.value(point);
                    calls++;
                    if (candidate < best) {
                        best = candidate;
                        original = point[parameter];
                        improved = true;
                    } else {
                        point[parameter] = original;
                    }
                }
            }
            if (!improved) step *= 0.5;
        }
        return new OptimizationResult(point, best, calls,
            calls < options.maximumOuterEvaluations());
    }

    private static OptimizationResult bobyqa(
            double[] initial, double[] lower, double[] upper,
            Objective objective, int interpolationPoints,
            ZeroInflatedMixedOptions options) {
        try {
            return Bobyqa.bobyqa(initial, lower, upper,
                objective::value, interpolationPoints,
                options.initialTrustRadius(),
                Math.max(1e-7, options.relativeTolerance()),
                options.maximumOuterEvaluations(), true);
        } catch (ArithmeticException | ArrayIndexOutOfBoundsException
                exception) {
            return coordinateOptimize(initial, lower, upper,
                objective, options);
        }
    }

    /** Bounded BFGS using central differences of the full Laplace objective. */
    private static OptimizationResult boundedBfgs(
            double[] initial, double[] lower, double[] upper,
            Objective objective, ZeroInflatedMixedOptions options,
            ExecutorService gradientExecutor,
            ThreadLocal<Objective> workerObjectives) {
        int dimensions = initial.length;
        double[] point = initial.clone();
        ObjectivePoint center = objective.point(point, null);
        double value = center.value();
        double[] randomMode = center.random();
        int calls = 1;
        GradientEvaluation differentiated = numericalGradient(objective,
            point, value, randomMode, lower, upper,
            options.maximumOuterEvaluations() - calls, gradientExecutor,
            workerObjectives);
        calls += differentiated.calls();
        if (!differentiated.complete())
            return exhausted(point, value, options.maximumOuterEvaluations());
        double[] gradient = differentiated.gradient();
        double[] inverseHessian = MatrixOps.identity(dimensions);
        double gradientTolerance = Math.max(1e-5,
            Math.sqrt(options.relativeTolerance()));
        boolean converged = projectedGradientMaximum(point, gradient,
            lower, upper) <= gradientTolerance;
        while (!converged && calls < options.maximumOuterEvaluations()) {
            double[] direction = new double[dimensions];
            for (int row = 0; row < dimensions; row++)
                for (int column = 0; column < dimensions; column++)
                    direction[row] -= inverseHessian[row * dimensions + column]
                        * gradient[column];
            projectDirection(point, direction, lower, upper);
            double directionalDerivative = dot(gradient, direction);
            if (!(directionalDerivative < 0.0)) {
                for (int index = 0; index < dimensions; index++)
                    direction[index] = -gradient[index];
                projectDirection(point, direction, lower, upper);
                directionalDerivative = dot(gradient, direction);
                inverseHessian = MatrixOps.identity(dimensions);
            }
            double directionMaximum = maximumAbsolute(direction);
            if (!(directionMaximum > 0.0)) {
                converged = true;
                break;
            }
            double scale = Math.min(1.0,
                options.initialTrustRadius() / directionMaximum);
            double[] candidate = new double[dimensions];
            double candidateValue = Double.POSITIVE_INFINITY;
            double[] candidateRandom = randomMode;
            boolean accepted = false;
            while (scale >= 1e-10
                    && calls < options.maximumOuterEvaluations()) {
                for (int index = 0; index < dimensions; index++)
                    candidate[index] = Math.max(lower[index],
                        Math.min(upper[index], point[index]
                            + scale * direction[index]));
                ObjectivePoint candidatePoint = objective.point(candidate,
                    randomMode);
                candidateValue = candidatePoint.value();
                candidateRandom = candidatePoint.random();
                calls++;
                if (Double.isFinite(candidateValue)
                        && candidateValue <= value
                            + 1e-4 * scale * directionalDerivative) {
                    accepted = true;
                    break;
                }
                scale *= 0.5;
            }
            if (!accepted) break;
            differentiated = numericalGradient(objective, candidate,
                candidateValue, candidateRandom, lower, upper,
                options.maximumOuterEvaluations() - calls, gradientExecutor,
                workerObjectives);
            calls += differentiated.calls();
            if (!differentiated.complete()) {
                point = candidate.clone();
                value = candidateValue;
                randomMode = candidateRandom;
                break;
            }
            double[] nextGradient = differentiated.gradient();
            double[] displacement = new double[dimensions];
            double[] gradientChange = new double[dimensions];
            for (int index = 0; index < dimensions; index++) {
                displacement[index] = candidate[index] - point[index];
                gradientChange[index] = nextGradient[index] - gradient[index];
            }
            double curvature = dot(displacement, gradientChange);
            double curvatureFloor = 1e-12 * Math.sqrt(
                dot(displacement, displacement)
                    * dot(gradientChange, gradientChange));
            if (curvature > curvatureFloor) {
                double[] hessianChange = new double[dimensions];
                for (int row = 0; row < dimensions; row++)
                    for (int column = 0; column < dimensions; column++)
                        hessianChange[row] += inverseHessian[
                            row * dimensions + column]
                            * gradientChange[column];
                double changeQuadratic = dot(gradientChange, hessianChange);
                double coefficient = (curvature + changeQuadratic)
                    / (curvature * curvature);
                for (int row = 0; row < dimensions; row++)
                    for (int column = 0; column < dimensions; column++)
                        inverseHessian[row * dimensions + column] += coefficient
                            * displacement[row] * displacement[column]
                            - (hessianChange[row] * displacement[column]
                                + displacement[row] * hessianChange[column])
                                / curvature;
            } else inverseHessian = MatrixOps.identity(dimensions);
            double relativeStep = relativeChange(point, candidate);
            point = candidate.clone();
            value = candidateValue;
            randomMode = candidateRandom;
            gradient = nextGradient;
            converged = relativeStep <= options.relativeTolerance()
                || projectedGradientMaximum(point, gradient, lower, upper)
                    <= gradientTolerance;
        }
        return new OptimizationResult(point, value,
            converged ? calls : options.maximumOuterEvaluations(), true);
    }

    private static GradientEvaluation numericalGradient(
            Objective objective, double[] point, double center,
            double[] centerRandom, double[] lower, double[] upper,
            int availableCalls, ExecutorService executor,
            ThreadLocal<Objective> workerObjectives) {
        if (executor != null && availableCalls >= point.length)
            return parallelNumericalGradient(point, center, centerRandom,
                lower, upper, executor, workerObjectives);
        double[] gradient = new double[point.length];
        double[] trial = point.clone();
        int calls = 0;
        for (int parameter = 0; parameter < point.length; parameter++) {
            double step = 1e-5 * (1.0 + Math.abs(point[parameter]));
            double below = Math.max(lower[parameter], point[parameter] - step);
            double above = Math.min(upper[parameter], point[parameter] + step);
            if (above == point[parameter] && below == point[parameter]) continue;
            if (calls + 1 > availableCalls)
                return new GradientEvaluation(gradient, calls, false);
            if (above > point[parameter]) {
                trial[parameter] = above;
                double upperValue = objective.valueFrom(trial, centerRandom);
                gradient[parameter] = (upperValue - center)
                    / (above - point[parameter]);
            } else {
                trial[parameter] = below;
                double lowerValue = objective.valueFrom(trial, centerRandom);
                gradient[parameter] = (center - lowerValue)
                    / (point[parameter] - below);
            }
            calls++;
            trial[parameter] = point[parameter];
        }
        return new GradientEvaluation(gradient, calls, true);
    }

    private static GradientEvaluation parallelNumericalGradient(
            double[] point, double center, double[] centerRandom,
            double[] lower, double[] upper, ExecutorService executor,
            ThreadLocal<Objective> workerObjectives) {
        List<Future<Double>> futures = new ArrayList<>(point.length);
        for (int parameter = 0; parameter < point.length; parameter++) {
            final int index = parameter;
            futures.add(executor.submit(() -> gradientComponent(
                workerObjectives.get(), point, center, centerRandom,
                lower, upper, index)));
        }
        double[] gradient = new double[point.length];
        try {
            for (int parameter = 0; parameter < point.length; parameter++)
                gradient[parameter] = futures.get(parameter).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            futures.forEach(value -> value.cancel(true));
            return new GradientEvaluation(gradient, point.length, false);
        } catch (ExecutionException exception) {
            futures.forEach(value -> value.cancel(true));
            return new GradientEvaluation(gradient, point.length, false);
        }
        return new GradientEvaluation(gradient, point.length, true);
    }

    private static double gradientComponent(
            Objective objective, double[] point, double center,
            double[] centerRandom, double[] lower, double[] upper,
            int parameter) {
        double[] trial = point.clone();
        double step = 1e-5 * (1.0 + Math.abs(point[parameter]));
        double above = Math.min(upper[parameter], point[parameter] + step);
        if (above > point[parameter]) {
            trial[parameter] = above;
            return (objective.valueFrom(trial, centerRandom) - center)
                / (above - point[parameter]);
        }
        double below = Math.max(lower[parameter], point[parameter] - step);
        if (below < point[parameter]) {
            trial[parameter] = below;
            return (center - objective.valueFrom(trial, centerRandom))
                / (point[parameter] - below);
        }
        return 0.0;
    }

    private static OptimizationResult exhausted(
            double[] point, double value, int maximumCalls) {
        return new OptimizationResult(point, value, maximumCalls, true);
    }

    private static void projectDirection(
            double[] point, double[] direction,
            double[] lower, double[] upper) {
        for (int index = 0; index < point.length; index++) {
            if ((point[index] <= lower[index] && direction[index] < 0.0)
                    || (point[index] >= upper[index]
                        && direction[index] > 0.0))
                direction[index] = 0.0;
        }
    }

    private static double projectedGradientMaximum(
            double[] point, double[] gradient,
            double[] lower, double[] upper) {
        double result = 0.0;
        for (int index = 0; index < point.length; index++) {
            double value = gradient[index];
            if ((point[index] <= lower[index] && value > 0.0)
                    || (point[index] >= upper[index] && value < 0.0))
                value = 0.0;
            result = Math.max(result, Math.abs(value));
        }
        return result;
    }

    private static double dot(double[] first, double[] second) {
        double result = 0.0;
        for (int index = 0; index < first.length; index++)
            result += first[index] * second[index];
        return result;
    }

    private static OptimizationResult optimizeReduced(
            double[] initial, double[] lower, double[] upper,
            java.util.function.ToDoubleFunction<double[]> objective,
            ZeroInflatedMixedOptions options) {
        int interpolationPoints = Math.min(2 * initial.length + 1,
            (initial.length + 1) * (initial.length + 2) / 2);
        try {
            return Bobyqa.bobyqa(initial, lower, upper,
                objective::applyAsDouble,
                interpolationPoints, options.initialTrustRadius(),
                Math.max(1e-7, options.relativeTolerance()),
                options.maximumOuterEvaluations(), true);
        } catch (ArithmeticException | ArrayIndexOutOfBoundsException failure) {
            double[] point = initial.clone();
            double best = objective.applyAsDouble(point);
            int calls = 1;
            double step = options.initialTrustRadius();
            while (calls < options.maximumOuterEvaluations()
                    && step > options.relativeTolerance()) {
                boolean improved = false;
                for (int parameter = 0; parameter < point.length
                        && calls < options.maximumOuterEvaluations(); parameter++) {
                    double original = point[parameter];
                    for (int direction : new int[] {-1, 1}) {
                        point[parameter] = Math.max(lower[parameter],
                            Math.min(upper[parameter],
                                original + direction * step));
                        double candidate = objective.applyAsDouble(point);
                        calls++;
                        if (candidate < best) {
                            best = candidate;
                            original = point[parameter];
                            improved = true;
                        } else point[parameter] = original;
                    }
                }
                if (!improved) step *= 0.5;
            }
            return new OptimizationResult(point, best, calls,
                calls < options.maximumOuterEvaluations());
        }
    }

    private static double[] remove(double[] values, int removed) {
        double[] result = new double[values.length - 1];
        System.arraycopy(values, 0, result, 0, removed);
        System.arraycopy(values, removed + 1, result, removed,
            values.length - removed - 1);
        return result;
    }

    private static double[] insert(
            double[] values, int inserted, double value) {
        double[] result = new double[values.length + 1];
        System.arraycopy(values, 0, result, 0, inserted);
        result[inserted] = value;
        System.arraycopy(values, inserted, result, inserted + 1,
            values.length - inserted);
        return result;
    }

    private static double[] simulate(
            ZeroInflatedMixedResult fitted, RandomGenerator random,
            CountFamily family) {
        double[] means = fitted.conditionalCountMeans();
        double[] probabilities = fitted.structuralZeroProbabilities();
        double[] sizes = fitted.sizes();
        double[] result = new double[means.length];
        for (int row = 0; row < result.length; row++) {
            if (random.nextDouble() < probabilities[row]) {
                result[row] = 0.0;
            } else if (family == CountFamily.POISSON) {
                result[row] = poisson(random, means[row]);
            } else {
                double intensity = gamma(random, sizes[row],
                    means[row] / sizes[row]);
                result[row] = poisson(random, intensity);
            }
        }
        return result;
    }

    private static int poisson(RandomGenerator random, double mean) {
        if (mean < 12.0) {
            double threshold = Math.exp(-mean);
            double product = 1.0;
            int value = -1;
            do {
                value++;
                product *= random.nextDouble();
            } while (product > threshold);
            return value;
        }
        double squareRoot = Math.sqrt(2.0 * mean);
        double logMean = Math.log(mean);
        double constant = mean * logMean - SpecialFunctions.logGamma(mean + 1.0);
        while (true) {
            double tangent = Math.tan(Math.PI * random.nextDouble());
            int value = (int) Math.floor(squareRoot * tangent + mean);
            if (value < 0) continue;
            double accept = 0.9 * (1.0 + tangent * tangent)
                * Math.exp(value * logMean
                    - SpecialFunctions.logGamma(value + 1.0) - constant);
            if (random.nextDouble() <= accept) return value;
        }
    }

    private static double gamma(
            RandomGenerator random, double shape, double scale) {
        if (shape < 1.0)
            return gamma(random, shape + 1.0, scale)
                * Math.pow(random.nextDouble(), 1.0 / shape);
        double d = shape - 1.0 / 3.0;
        double c = 1.0 / Math.sqrt(9.0 * d);
        while (true) {
            double normal = random.nextGaussian();
            double factor = 1.0 + c * normal;
            if (factor <= 0.0) continue;
            factor = factor * factor * factor;
            double uniform = random.nextDouble();
            if (uniform < 1.0 - 0.0331 * normal * normal * normal * normal
                    || Math.log(uniform) < 0.5 * normal * normal
                        + d * (1.0 - factor + Math.log(factor)))
                return scale * d * factor;
        }
    }

    private static InferenceData numericalInference(
            Objective objective, double[] parameters,
            double[] lower, double[] upper,
            int countColumns, int zeroColumns, int dispersionColumns) {
        int dimensions = parameters.length;
        double[] steps = new double[dimensions];
        for (int index = 0; index < dimensions; index++) {
            double proposed = 2e-4 * (1.0 + Math.abs(parameters[index]));
            double room = Math.min(parameters[index] - lower[index],
                upper[index] - parameters[index]);
            steps[index] = Math.min(proposed, 0.45 * room);
            if (!(steps[index] > 1e-8)) return InferenceData.unavailable();
        }
        double center = objective.coldValue(parameters);
        double[] hessian = new double[dimensions * dimensions];
        for (int left = 0; left < dimensions; left++) {
            double[] plus = parameters.clone();
            double[] minus = parameters.clone();
            plus[left] += steps[left];
            minus[left] -= steps[left];
            double diagonal = (objective.coldValue(plus) - 2.0 * center
                + objective.coldValue(minus)) / (steps[left] * steps[left]);
            hessian[left * dimensions + left] = diagonal;
            for (int right = 0; right < left; right++) {
                double[] pp = parameters.clone();
                double[] pm = parameters.clone();
                double[] mp = parameters.clone();
                double[] mm = parameters.clone();
                pp[left] += steps[left];
                pp[right] += steps[right];
                pm[left] += steps[left];
                pm[right] -= steps[right];
                mp[left] -= steps[left];
                mp[right] += steps[right];
                mm[left] -= steps[left];
                mm[right] -= steps[right];
                double value = (objective.coldValue(pp)
                    - objective.coldValue(pm) - objective.coldValue(mp)
                    + objective.coldValue(mm))
                    / (4.0 * steps[left] * steps[right]);
                hessian[left * dimensions + right] = value;
                hessian[right * dimensions + left] = value;
            }
        }
        double[] inverse = inverseWithRidge(hessian, dimensions);
        if (inverse == null) return InferenceData.unavailable();
        int reported = countColumns + zeroColumns + dispersionColumns;
        double[] covariance = new double[reported * reported];
        for (int row = 0; row < reported; row++)
            System.arraycopy(inverse, row * dimensions, covariance,
                row * reported, reported);
        List<String> names = new ArrayList<>(reported);
        for (int index = 0; index < countColumns; index++)
            names.add("count[" + index + "]");
        for (int index = 0; index < zeroColumns; index++)
            names.add("zero[" + index + "]");
        for (int index = 0; index < dispersionColumns; index++)
            names.add("log-size[" + index + "]");
        return new InferenceData(List.copyOf(names), covariance);
    }

    private static double[] inverseWithRidge(double[] matrix, int dimension) {
        for (double ridge : new double[] {0.0, 1e-9, 1e-7, 1e-5}) {
            double[] augmented = new double[dimension * 2 * dimension];
            int width = 2 * dimension;
            for (int row = 0; row < dimension; row++) {
                System.arraycopy(matrix, row * dimension, augmented,
                    row * width, dimension);
                augmented[row * width + row] += ridge;
                augmented[row * width + dimension + row] = 1.0;
            }
            boolean valid = true;
            for (int pivot = 0; pivot < dimension; pivot++) {
                int selected = pivot;
                for (int row = pivot + 1; row < dimension; row++)
                    if (Math.abs(augmented[row * width + pivot])
                            > Math.abs(augmented[selected * width + pivot]))
                        selected = row;
                if (!(Math.abs(augmented[selected * width + pivot]) > 1e-12)) {
                    valid = false;
                    break;
                }
                if (selected != pivot)
                    for (int column = 0; column < width; column++) {
                        double value = augmented[pivot * width + column];
                        augmented[pivot * width + column] =
                            augmented[selected * width + column];
                        augmented[selected * width + column] = value;
                    }
                double scale = augmented[pivot * width + pivot];
                for (int column = 0; column < width; column++)
                    augmented[pivot * width + column] /= scale;
                for (int row = 0; row < dimension; row++) {
                    if (row == pivot) continue;
                    double multiple = augmented[row * width + pivot];
                    for (int column = 0; column < width; column++)
                        augmented[row * width + column] -= multiple
                            * augmented[pivot * width + column];
                }
            }
            if (valid) {
                double[] inverse = new double[dimension * dimension];
                for (int row = 0; row < dimension; row++)
                    System.arraycopy(augmented, row * width + dimension,
                        inverse, row * dimension, dimension);
                return inverse;
            }
        }
        return null;
    }

    private static List<String> diagnostics(
            double[] parameters, int varianceStart,
            RandomStructure structure, DataState data,
            ZeroInflatedMixedOptions options, boolean converged) {
        List<String> warnings = new ArrayList<>();
        if (!converged) warnings.add("optimizer did not satisfy all tolerances");
        for (int index = 0; index < structure.varianceCount(); index++) {
            double variance = Math.exp(parameters[varianceStart + index]);
            if (variance <= options.minimumVariance() * 1.05)
                warnings.add("variance component " + index
                    + " is on the lower boundary");
        }
        double[] covariance = Arrays.copyOfRange(parameters, varianceStart,
            parameters.length);
        for (CorrelatedLayout pair : structure.correlated())
            if (Math.abs(structure.correlation(
                    covariance, pair.correlation())) >= 0.95)
                warnings.add("correlation " + pair.name()
                    + " is near its guarded boundary");
        double minimum = Arrays.stream(data.zeroProbabilities()).min()
            .orElse(0.5);
        double maximum = Arrays.stream(data.zeroProbabilities()).max()
            .orElse(0.5);
        if (minimum <= 1e-8 || maximum >= 1.0 - 1e-8)
            warnings.add("structural-zero probabilities are near a boundary");
        return List.copyOf(warnings);
    }

    private static ZeroInflatedMixedResult result(
            Evaluation fitted, double[] parameters,
            int countColumns, int zeroColumns, int dispersionColumns,
            CountFamily family, RandomStructure structure,
            CombinedDesign design, SparsePattern pattern,
            int factorNonzeros, int evaluations, boolean outerConverged,
            InferenceData inference, ZeroInflatedMixedOptions options) {
        int zeroStart = countColumns;
        int dispersionStart = zeroStart + zeroColumns;
        int varianceStart = dispersionStart
            + (family == CountFamily.NEGATIVE_BINOMIAL
                ? dispersionColumns : 0);
        double[] count = Arrays.copyOfRange(parameters, 0, zeroStart);
        double[] zero = Arrays.copyOfRange(
            parameters, zeroStart, dispersionStart);
        double[] dispersion = family == CountFamily.NEGATIVE_BINOMIAL
            ? Arrays.copyOfRange(parameters, dispersionStart, varianceStart)
            : new double[0];
        double[] covarianceParameters = Arrays.copyOfRange(parameters,
            varianceStart, parameters.length);
        double[] variances = structure.variances(covarianceParameters);
        List<String> names = new ArrayList<>(structure.terms().size());
        Map<String, double[]> random = new LinkedHashMap<>();
        for (int term = 0; term < structure.terms().size(); term++) {
            RandomEffectTerm value = structure.terms().get(term);
            String name = structure.resultNames().get(term);
            names.add(name);
            random.put(name, Arrays.copyOfRange(fitted.random(),
                design.termStarts()[term],
                design.termStarts()[term] + value.coefficients()));
        }
        DataState data = fitted.data();
        List<String> warnings = diagnostics(parameters, varianceStart,
            structure, data, options, outerConverged && fitted.converged());
        return new ZeroInflatedMixedResult(
            family == CountFamily.POISSON
                ? "zero-inflated-poisson" : "zero-inflated-negative-binomial-2",
            count, zero, dispersion, names, variances,
            structure.correlations(covarianceParameters), random,
            data.means(), data.zeroProbabilities(), data.sizes(),
            data.responseMeans(), data.totalZeroProbabilities(),
            fitted.laplaceLogLikelihood(), evaluations, fitted.iterations(),
            outerConverged && fitted.converged(), design.columns(),
            pattern.columnIndices().length, factorNonzeros,
            inference.parameterNames(), inference.covariance(), warnings,
            outerParameterNames(countColumns, zeroColumns,
                family == CountFamily.NEGATIVE_BINOMIAL
                    ? dispersionColumns : 0, structure), parameters);
    }

    private static List<String> outerParameterNames(
            int countColumns, int zeroColumns, int dispersionColumns,
            RandomStructure structure) {
        List<String> names = new ArrayList<>();
        for (int index = 0; index < countColumns; index++)
            names.add("count[" + index + "]");
        for (int index = 0; index < zeroColumns; index++)
            names.add("zero[" + index + "]");
        for (int index = 0; index < dispersionColumns; index++)
            names.add("log-size[" + index + "]");
        String[] varianceNames = new String[structure.varianceCount()];
        for (int term = 0; term < structure.terms().size(); term++) {
            int variance = structure.termVarianceIndices()[term];
            if (varianceNames[variance] == null)
                varianceNames[variance] = "log-variance["
                    + structure.resultNames().get(term) + "]";
        }
        names.addAll(Arrays.asList(varianceNames));
        for (CorrelatedLayout pair : structure.correlated())
            names.add("fisher-z[" + pair.name() + "]");
        return List.copyOf(names);
    }

    private static final class Objective {
        private final double[] response;
        private final double[] countFixed;
        private final int countColumns;
        private final double[] zeroFixed;
        private final int zeroColumns;
        private final double[] dispersionFixed;
        private final int dispersionColumns;
        private final double[] offsets;
        private final CountFamily family;
        private final CombinedDesign design;
        private final PrecisionData precision;
        private final SparsePattern pattern;
        private final PreparedSparseCholesky factor;
        private final ZeroInflatedMixedOptions options;
        private final double[] responseLogFactorials;
        private final double[] means;
        private final double[] zeroProbabilities;
        private final double[] sizes;
        private final double[] responseMeans;
        private final double[] totalZeroProbabilities;
        private final double[] countScores;
        private final double[] zeroScores;
        private final double[] countCurvatures;
        private final double[] zeroCurvatures;
        private final double[] crossCurvatures;
        private double[] warmRandom;

        Objective(
                double[] response,
                double[] countFixed, int countColumns,
                double[] zeroFixed, int zeroColumns,
                double[] dispersionFixed, int dispersionColumns,
                double[] offsets, CountFamily family,
                CombinedDesign design, PrecisionData precision,
                SparsePattern pattern, PreparedSparseCholesky factor,
                ZeroInflatedMixedOptions options) {
            this.response = response;
            this.countFixed = countFixed;
            this.countColumns = countColumns;
            this.zeroFixed = zeroFixed;
            this.zeroColumns = zeroColumns;
            this.dispersionFixed = dispersionFixed;
            this.dispersionColumns = dispersionColumns;
            this.offsets = offsets;
            this.family = family;
            this.design = design;
            this.precision = precision;
            this.pattern = pattern;
            this.factor = factor;
            this.options = options;
            responseLogFactorials = new double[response.length];
            for (int row = 0; row < response.length; row++)
                responseLogFactorials[row] = SpecialFunctions.logGamma(
                    response[row] + 1.0);
            means = new double[response.length];
            zeroProbabilities = new double[response.length];
            sizes = family == CountFamily.NEGATIVE_BINOMIAL
                ? new double[response.length] : new double[0];
            responseMeans = new double[response.length];
            totalZeroProbabilities = new double[response.length];
            countScores = new double[response.length];
            zeroScores = new double[response.length];
            countCurvatures = new double[response.length];
            zeroCurvatures = new double[response.length];
            crossCurvatures = new double[response.length];
        }

        private Objective(Objective source, PreparedSparseCholesky factor) {
            response = source.response;
            countFixed = source.countFixed;
            countColumns = source.countColumns;
            zeroFixed = source.zeroFixed;
            zeroColumns = source.zeroColumns;
            dispersionFixed = source.dispersionFixed;
            dispersionColumns = source.dispersionColumns;
            offsets = source.offsets;
            family = source.family;
            design = source.design;
            precision = source.precision;
            pattern = source.pattern;
            this.factor = factor;
            options = source.options;
            responseLogFactorials = source.responseLogFactorials;
            means = new double[response.length];
            zeroProbabilities = new double[response.length];
            sizes = family == CountFamily.NEGATIVE_BINOMIAL
                ? new double[response.length] : new double[0];
            responseMeans = new double[response.length];
            totalZeroProbabilities = new double[response.length];
            countScores = new double[response.length];
            zeroScores = new double[response.length];
            countCurvatures = new double[response.length];
            zeroCurvatures = new double[response.length];
            crossCurvatures = new double[response.length];
        }

        Objective fork(PreparedSparseCholesky workerFactor) {
            return new Objective(this, workerFactor);
        }

        double value(double[] parameters) {
            try {
                Evaluation result = evaluate(parameters);
                return Double.isFinite(result.laplaceLogLikelihood())
                    ? -result.laplaceLogLikelihood() : INVALID_OBJECTIVE;
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return INVALID_OBJECTIVE;
            }
        }

        double coldValue(double[] parameters) {
            warmRandom = null;
            return value(parameters);
        }

        ObjectivePoint point(double[] parameters, double[] initialRandom) {
            warmRandom = initialRandom == null ? null : initialRandom.clone();
            try {
                Evaluation result = evaluate(parameters);
                double value = -result.laplaceLogLikelihood();
                return new ObjectivePoint(
                    Double.isFinite(value) ? value : INVALID_OBJECTIVE,
                    result.random().clone());
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return new ObjectivePoint(INVALID_OBJECTIVE,
                    initialRandom == null ? new double[design.columns()]
                        : initialRandom.clone());
            }
        }

        double valueFrom(double[] parameters, double[] initialRandom) {
            return point(parameters, initialRandom).value();
        }

        Evaluation evaluate(double[] parameters) {
            int varianceStart = countColumns + zeroColumns
                + (family == CountFamily.NEGATIVE_BINOMIAL
                    ? dispersionColumns : 0);
            double[] logVariances = Arrays.copyOfRange(
                parameters, varianceStart, parameters.length);
            double[] random = warmRandom == null
                ? new double[design.columns()] : warmRandom.clone();
            boolean converged = false;
            int iterations = 0;
            for (int iteration = 1;
                    iteration <= options.maximumModeIterations(); iteration++) {
                iterations = iteration;
                DataState data = data(parameters, random, false);
                double[] gradient = randomGradient(
                    data.countScores(), data.zeroScores(), random,
                    logVariances);
                double gradientMaximum = maximumAbsolute(gradient);
                double damping = 0.0;
                boolean factored = false;
                for (int attempt = 0; attempt < 12; attempt++) {
                    try {
                        factor.refactor(pattern.matrix(data.countCurvatures(),
                            data.zeroCurvatures(), data.crossCurvatures(),
                            design, precision, logVariances, damping));
                        factored = true;
                        break;
                    } catch (IllegalArgumentException exception) {
                        damping = damping == 0.0 ? 1e-8 : damping * 10.0;
                    }
                }
                if (!factored) {
                    throw new IllegalArgumentException(
                        "random-effect Hessian is not positive definite");
                }
                double[] step = gradient.clone();
                factor.solveInPlace(step, 1);
                double joint = data.logLikelihood()
                    - 0.5 * precision.quadratic(random, logVariances);
                double scale = 1.0;
                double[] candidate = null;
                for (int attempt = 0; attempt < 30; attempt++) {
                    double[] trial = random.clone();
                    for (int index = 0; index < trial.length; index++)
                        trial[index] += scale * step[index];
                    double trialJoint = logLikelihood(parameters, trial)
                        - 0.5 * precision.quadratic(trial, logVariances);
                    if (Double.isFinite(trialJoint) && trialJoint >= joint) {
                        candidate = trial;
                        break;
                    }
                    scale *= 0.5;
                }
                if (candidate == null) {
                    if (gradientMaximum <= Math.sqrt(
                            options.relativeTolerance())) {
                        converged = true;
                        break;
                    }
                    throw new IllegalArgumentException(
                        "random-effect Newton step could not improve the mode");
                }
                double change = relativeChange(random, candidate);
                random = candidate;
                if (change <= options.relativeTolerance()
                        && gradientMaximum
                            <= Math.sqrt(options.relativeTolerance())) {
                    converged = true;
                    break;
                }
            }
            DataState fitted = data(parameters, random, true);
            factor.refactor(pattern.matrix(fitted.countCurvatures(),
                fitted.zeroCurvatures(), fitted.crossCurvatures(), design,
                precision, logVariances, 0.0));
            double joint = fitted.logLikelihood()
                - 0.5 * precision.quadratic(random, logVariances);
            double laplace = joint
                + 0.5 * precision.logDeterminant(logVariances)
                - 0.5 * factor.logDeterminant();
            if (!Double.isFinite(laplace)) {
                throw new IllegalArgumentException(
                    "Laplace likelihood is not finite");
            }
            warmRandom = random.clone();
            return new Evaluation(random, fitted, laplace,
                iterations, converged);
        }

        private DataState data(
                double[] parameters, double[] random,
                boolean fittedOutputs) {
            int zeroStart = countColumns;
            int dispersionStart = zeroStart + zeroColumns;
            int varianceStart = dispersionStart
                + (family == CountFamily.NEGATIVE_BINOMIAL
                    ? dispersionColumns : 0);
            if (parameters.length != varianceStart
                    + precision.structure().varianceCount()
                    + precision.structure().correlationCount()) {
                throw new IllegalArgumentException(
                    "outer parameter dimensions do not match the model");
            }
            double logLikelihood = 0.0;
            RowLikelihoodWorkspace likelihood = new RowLikelihoodWorkspace();
            for (int row = 0; row < response.length; row++) {
                double countEta = offsets[row]
                    + fixedValue(countFixed, row, countColumns, parameters, 0)
                    + design.rowProduct(row, random, true);
                double zeroEta = fixedValue(zeroFixed, row, zeroColumns,
                    parameters, zeroStart)
                    + design.rowProduct(row, random, false);
                double mean = safeExp(countEta);
                double probability = logistic(zeroEta);
                double size = family == CountFamily.NEGATIVE_BINOMIAL
                    ? safeExp(fixedValue(dispersionFixed, row,
                        dispersionColumns, parameters, dispersionStart))
                    : Double.POSITIVE_INFINITY;
                fillRowLikelihood(response[row], mean, probability, size,
                    family, responseLogFactorials[row], fittedOutputs,
                    likelihood);
                logLikelihood += likelihood.logLikelihood;
                if (fittedOutputs) {
                    means[row] = mean;
                    zeroProbabilities[row] = probability;
                    if (sizes.length > 0) sizes[row] = size;
                    responseMeans[row] = (1.0 - probability) * mean;
                    totalZeroProbabilities[row] =
                        likelihood.totalZeroProbability;
                }
                countScores[row] = likelihood.countScore;
                zeroScores[row] = likelihood.zeroScore;
                countCurvatures[row] = likelihood.countCurvature;
                zeroCurvatures[row] = likelihood.zeroCurvature;
                crossCurvatures[row] = likelihood.crossCurvature;
            }
            return new DataState(logLikelihood, countScores, zeroScores,
                countCurvatures, zeroCurvatures, crossCurvatures, means,
                zeroProbabilities, sizes, responseMeans, totalZeroProbabilities);
        }

        /** Likelihood-only pass used by the Newton line search. */
        private double logLikelihood(double[] parameters, double[] random) {
            int zeroStart = countColumns;
            int dispersionStart = zeroStart + zeroColumns;
            double result = 0.0;
            for (int row = 0; row < response.length; row++) {
                double countEta = offsets[row]
                    + fixedValue(countFixed, row, countColumns, parameters, 0)
                    + design.rowProduct(row, random, true);
                double zeroEta = fixedValue(zeroFixed, row, zeroColumns,
                    parameters, zeroStart)
                    + design.rowProduct(row, random, false);
                double mean = safeExp(countEta);
                double probability = logistic(zeroEta);
                double size = family == CountFamily.NEGATIVE_BINOMIAL
                    ? safeExp(fixedValue(dispersionFixed, row,
                        dispersionColumns, parameters, dispersionStart))
                    : Double.POSITIVE_INFINITY;
                result += rowLogLikelihood(response[row], mean, probability,
                    size, family, responseLogFactorials[row]);
            }
            return result;
        }

        private double[] randomGradient(
                double[] countScores, double[] zeroScores,
                double[] random, double[] logVariances) {
            double[] gradient = new double[design.columns()];
            for (int row = 0; row < design.rows(); row++)
                for (int index = design.rowStarts()[row];
                        index < design.rowStarts()[row + 1]; index++)
                    gradient[design.columnIndices()[index]] +=
                        (design.countColumns()[design.columnIndices()[index]]
                            ? countScores[row] : zeroScores[row])
                        * design.values()[index];
            precision.subtractProduct(random, logVariances, gradient);
            return gradient;
        }
    }

    private static RowLikelihood rowLikelihood(
            double response, double mean, double zeroProbability,
            double size, CountFamily family) {
        RowLikelihoodWorkspace result = new RowLikelihoodWorkspace();
        fillRowLikelihood(response, mean, zeroProbability, size, family,
            SpecialFunctions.logGamma(response + 1.0), true, result);
        return new RowLikelihood(result.logLikelihood, result.countScore,
            result.zeroScore, result.countCurvature, result.zeroCurvature,
            result.crossCurvature, result.totalZeroProbability);
    }

    private static void fillRowLikelihood(
            double response, double mean, double zeroProbability,
            double size, CountFamily family, double responseLogFactorial,
            boolean totalZeroNeeded, RowLikelihoodWorkspace result) {
        double logCount;
        double logCountZero;
        double baseZeroScore;
        double baseZeroScoreDerivative;
        double positiveScore;
        double positiveCurvature;
        if (family == CountFamily.POISSON) {
            logCountZero = -mean;
            logCount = response * Math.log(mean) - mean
                - responseLogFactorial;
            baseZeroScore = -mean;
            baseZeroScoreDerivative = -mean;
            positiveScore = response - mean;
            positiveCurvature = mean;
        } else {
            double total = size + mean;
            logCountZero = size * (Math.log(size) - Math.log(total));
            logCount = SpecialFunctions.logGamma(response + size)
                - SpecialFunctions.logGamma(size)
                - responseLogFactorial
                + size * (Math.log(size) - Math.log(total))
                + response * (Math.log(mean) - Math.log(total));
            baseZeroScore = -size * mean / total;
            baseZeroScoreDerivative =
                -size * size * mean / (total * total);
            positiveScore = size * (response - mean) / total;
            positiveCurvature =
                size * mean * (size + response) / (total * total);
        }
        double logOneMinusPi = logOneMinusLogistic(zeroProbability);
        if (response > 0.0) {
            double totalZero = 0.0;
            if (totalZeroNeeded) {
                double logZeroMass = logAddExp(logLogistic(zeroProbability),
                    logOneMinusPi + logCountZero);
                totalZero = Math.exp(logZeroMass);
            }
            result.set(logOneMinusPi + logCount, positiveScore,
                -zeroProbability, positiveCurvature,
                zeroProbability * (1.0 - zeroProbability), 0.0, totalZero);
            return;
        }
        double logZeroMass = logAddExp(logLogistic(zeroProbability),
            logOneMinusPi + logCountZero);
        double totalZero = Math.exp(logZeroMass);
        double countPosterior = Math.exp(
            logOneMinusPi + logCountZero - logZeroMass);
        double structuralPosterior = 1.0 - countPosterior;
        double score = countPosterior * baseZeroScore;
        double derivative = countPosterior
            * ((1.0 - countPosterior) * baseZeroScore * baseZeroScore
                + baseZeroScoreDerivative);
        double zeroScore = structuralPosterior - zeroProbability;
        double zeroCurvature = zeroProbability * (1.0 - zeroProbability)
            - countPosterior * structuralPosterior;
        double crossCurvature = countPosterior * structuralPosterior
            * baseZeroScore;
        result.set(logZeroMass, score, zeroScore, -derivative,
            zeroCurvature, crossCurvature, totalZero);
    }

    private static double rowLogLikelihood(
            double response, double mean, double zeroProbability,
            double size, CountFamily family, double responseLogFactorial) {
        double logOneMinusPi = logOneMinusLogistic(zeroProbability);
        if (response > 0.0) {
            if (family == CountFamily.POISSON) {
                return logOneMinusPi + response * Math.log(mean) - mean
                    - responseLogFactorial;
            }
            double total = size + mean;
            return logOneMinusPi
                + SpecialFunctions.logGamma(response + size)
                - SpecialFunctions.logGamma(size) - responseLogFactorial
                + size * (Math.log(size) - Math.log(total))
                + response * (Math.log(mean) - Math.log(total));
        }
        double logCountZero = family == CountFamily.POISSON
            ? -mean : size * (Math.log(size) - Math.log(size + mean));
        return logAddExp(logLogistic(zeroProbability),
            logOneMinusPi + logCountZero);
    }

    static double[] countLikelihoodDerivatives(
            double response, double countPredictor,
            double zeroPredictor, double size,
            boolean negativeBinomial) {
        RowLikelihood value = rowLikelihood(response, safeExp(countPredictor),
            logistic(zeroPredictor), size,
            negativeBinomial ? CountFamily.NEGATIVE_BINOMIAL
                : CountFamily.POISSON);
        return new double[] {value.logLikelihood(), value.countScore(),
            value.countCurvature()};
    }

    static double[] likelihoodDerivatives(
            double response, double countPredictor,
            double zeroPredictor, double size,
            boolean negativeBinomial) {
        RowLikelihood value = rowLikelihood(response, safeExp(countPredictor),
            logistic(zeroPredictor), size,
            negativeBinomial ? CountFamily.NEGATIVE_BINOMIAL
                : CountFamily.POISSON);
        return new double[] {value.logLikelihood(), value.countScore(),
            value.zeroScore(), value.countCurvature(), value.zeroCurvature(),
            value.crossCurvature()};
    }

    private static double logLogistic(double probability) {
        return Math.log(Math.max(MINIMUM_PROBABILITY, probability));
    }

    private static double logOneMinusLogistic(double probability) {
        return Math.log(Math.max(MINIMUM_PROBABILITY, 1.0 - probability));
    }

    private static double logistic(double value) {
        if (value >= 0.0) {
            double exponential = Math.exp(-Math.min(40.0, value));
            return Math.max(MINIMUM_PROBABILITY,
                Math.min(1.0 - MINIMUM_PROBABILITY,
                    1.0 / (1.0 + exponential)));
        }
        double exponential = Math.exp(Math.max(-40.0, value));
        return Math.max(MINIMUM_PROBABILITY,
            Math.min(1.0 - MINIMUM_PROBABILITY,
                exponential / (1.0 + exponential)));
    }

    private static double safeExp(double value) {
        return Math.exp(Math.max(-40.0, Math.min(40.0, value)));
    }

    private static double logAddExp(double first, double second) {
        double maximum = Math.max(first, second);
        return maximum + Math.log(
            Math.exp(first - maximum) + Math.exp(second - maximum));
    }

    private static double fixedValue(
            double[] design, int row, int columns,
            double[] coefficients, int start) {
        double value = 0.0;
        for (int column = 0; column < columns; column++)
            value += design[row * columns + column]
                * coefficients[start + column];
        return value;
    }

    private static double maximumAbsolute(double[] values) {
        double result = 0.0;
        for (double value : values) result = Math.max(result, Math.abs(value));
        return result;
    }

    private static double relativeChange(double[] previous, double[] next) {
        double result = 0.0;
        for (int index = 0; index < previous.length; index++)
            result = Math.max(result, Math.abs(next[index] - previous[index])
                / (1.0 + Math.abs(previous[index])));
        return result;
    }

    private static double[] initial(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            double[] dispersionFixed, int dispersionColumns,
            int varianceCount, int correlationCount, CountFamily family,
            ZeroInflatedMixedOptions options) {
        int dispersionCount = family == CountFamily.NEGATIVE_BINOMIAL
            ? dispersionColumns : 0;
        double[] result = new double[
            countColumns + zeroColumns + dispersionCount + varianceCount
                + correlationCount];
        double mean = Arrays.stream(response).average().orElseThrow();
        int zeros = 0;
        double variance = 0.0;
        for (double value : response) {
            if (value == 0.0) zeros++;
            variance += (value - mean) * (value - mean);
        }
        variance /= Math.max(1, response.length - 1);
        double size = variance > mean
            ? mean * mean / Math.max(1e-8, variance - mean) : 100.0;
        size = Math.max(options.minimumSize(),
            Math.min(options.maximumSize(), size));
        double baseZero = family == CountFamily.POISSON
            ? Math.exp(-Math.max(1e-8, mean))
            : Math.pow(size / (size + Math.max(1e-8, mean)), size);
        double observedZero = zeros / (double) response.length;
        double inflation = (observedZero - baseZero)
            / Math.max(1e-8, 1.0 - baseZero);
        inflation = Math.max(0.01, Math.min(0.9, inflation));
        double countMean = Math.max(1e-8, mean / (1.0 - inflation));
        setIntercept(result, 0, countFixed, response.length, countColumns,
            Math.log(countMean));
        setIntercept(result, countColumns, zeroFixed, response.length,
            zeroColumns, Math.log(inflation / (1.0 - inflation)));
        int dispersionStart = countColumns + zeroColumns;
        if (family == CountFamily.NEGATIVE_BINOMIAL)
            setIntercept(result, dispersionStart, dispersionFixed,
                response.length, dispersionColumns, Math.log(size));
        int varianceStart = dispersionStart + dispersionCount;
        double[] supplied = options.initialVariances();
        for (int term = 0; term < varianceCount; term++) {
            double initialVariance = supplied == null ? 0.5 : supplied[term];
            result[varianceStart + term] = Math.log(initialVariance);
        }
        // Correlations use guarded Fisher-z parameters and start at independence.
        return result;
    }

    private static void setIntercept(
            double[] coefficients, int start, double[] design,
            int rows, int columns, double value) {
        for (int column = 0; column < columns; column++) {
            boolean intercept = true;
            for (int row = 0; row < rows; row++)
                if (Math.abs(design[row * columns + column] - 1.0) > 1e-12) {
                    intercept = false;
                    break;
                }
            if (intercept) {
                coefficients[start + column] = value;
                return;
            }
        }
    }

    private static double[] lower(
            int dimensions, int countColumns, int zeroColumns,
            int dispersionColumns, int varianceCount, int correlationCount,
            CountFamily family,
            ZeroInflatedMixedOptions options) {
        double[] result = new double[dimensions];
        Arrays.fill(result, -options.maximumAbsoluteCoefficient());
        int dispersionStart = countColumns + zeroColumns;
        int varianceStart = dispersionStart
            + (family == CountFamily.NEGATIVE_BINOMIAL
                ? dispersionColumns : 0);
        if (family == CountFamily.NEGATIVE_BINOMIAL
                && dispersionColumns == 1) {
            result[dispersionStart] = Math.log(options.minimumSize());
        }
        for (int term = 0; term < varianceCount; term++)
            result[varianceStart + term] = Math.log(options.minimumVariance());
        for (int correlation = 0; correlation < correlationCount; correlation++)
            result[varianceStart + varianceCount + correlation] = -4.0;
        return result;
    }

    private static double[] upper(
            int dimensions, int countColumns, int zeroColumns,
            int dispersionColumns, int varianceCount, int correlationCount,
            CountFamily family,
            ZeroInflatedMixedOptions options) {
        double[] result = new double[dimensions];
        Arrays.fill(result, options.maximumAbsoluteCoefficient());
        int dispersionStart = countColumns + zeroColumns;
        int varianceStart = dispersionStart
            + (family == CountFamily.NEGATIVE_BINOMIAL
                ? dispersionColumns : 0);
        if (family == CountFamily.NEGATIVE_BINOMIAL
                && dispersionColumns == 1) {
            result[dispersionStart] = Math.log(options.maximumSize());
        }
        for (int term = 0; term < varianceCount; term++)
            result[varianceStart + term] = Math.log(options.maximumVariance());
        for (int correlation = 0; correlation < correlationCount; correlation++)
            result[varianceStart + varianceCount + correlation] = 4.0;
        return result;
    }

    private static double[] varianceParameters(
            double[] parameters, int countColumns, int zeroColumns,
            int dispersionColumns, CountFamily family) {
        int start = countColumns + zeroColumns
            + (family == CountFamily.NEGATIVE_BINOMIAL
                ? dispersionColumns : 0);
        return Arrays.copyOfRange(parameters, start, parameters.length);
    }

    private static void validate(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            double[] dispersionFixed, int dispersionColumns,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> countPrecisionBases,
            List<RandomEffectTerm> zeroRandomEffects,
            List<SparsePrecisionMatrix> zeroPrecisionBases,
            List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
            double[] offset, CountFamily family,
            ZeroInflatedMixedOptions options, BackendPolicy backendPolicy) {
        if (response == null || response.length == 0
                || countColumns < 1 || zeroColumns < 1
                || countFixed == null
                || countFixed.length != response.length * countColumns
                || zeroFixed == null
                || zeroFixed.length != response.length * zeroColumns
                || countRandomEffects == null || zeroRandomEffects == null
                || correlatedEffects == null
                || countRandomEffects.isEmpty() && zeroRandomEffects.isEmpty()
                    && correlatedEffects.isEmpty()
                || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "response, fixed designs, random effects, controls, and backend are required");
        }
        MatrixOps.requireFinite(response, "response");
        MatrixOps.requireFinite(countFixed, "countFixed");
        MatrixOps.requireFinite(zeroFixed, "zeroFixed");
        for (double value : response)
            if (value < 0.0 || value != Math.rint(value))
                throw new IllegalArgumentException(
                    "zero-inflated responses must be nonnegative integers");
        if (offset != null && offset.length != response.length)
            throw new IllegalArgumentException(
                "count offset length must match observations");
        if (family == CountFamily.NEGATIVE_BINOMIAL) {
            if (dispersionColumns < 1 || dispersionFixed == null
                    || dispersionFixed.length
                        != response.length * dispersionColumns)
                throw new IllegalArgumentException(
                    "ZINB requires a matching fixed-effects size design");
            MatrixOps.requireFinite(dispersionFixed, "dispersionFixed");
        } else if (dispersionColumns != 0 || dispersionFixed != null) {
            throw new IllegalArgumentException(
                "ZIP does not accept a dispersion design");
        }
        java.util.HashSet<String> names = new java.util.HashSet<>();
        for (RandomEffectTerm term : countRandomEffects)
            if (term == null || term.observations() != response.length
                    || !names.add("count:" + term.name()))
                throw new IllegalArgumentException(
                    "random effects need unique names and matching rows");
        for (RandomEffectTerm term : zeroRandomEffects)
            if (term == null || term.observations() != response.length
                    || !names.add("zero:" + term.name()))
                throw new IllegalArgumentException(
                    "random effects need unique names and matching rows");
        for (CorrelatedZeroInflatedRandomEffect effect : correlatedEffects)
            if (effect == null
                    || effect.countTerm().observations() != response.length
                    || !names.add("correlated:" + effect.name()))
                throw new IllegalArgumentException(
                    "correlated effects need unique names and matching rows");
        validatePrecisions(countRandomEffects, countPrecisionBases);
        validatePrecisions(zeroRandomEffects, zeroPrecisionBases);
        int varianceCount = countRandomEffects.size() + zeroRandomEffects.size()
            + 2 * correlatedEffects.size();
        double[] supplied = options.initialVariances();
        if (supplied != null && supplied.length != varianceCount)
            throw new IllegalArgumentException(
                "one initial variance is required per independent process variance");
    }

    private static void validateData(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            double[] dispersionFixed, int dispersionColumns,
            double[] offset, CountFamily family, int expectedRows) {
        if (response == null || response.length != expectedRows
                || countColumns < 1 || zeroColumns < 1
                || countFixed == null
                || countFixed.length != expectedRows * countColumns
                || zeroFixed == null
                || zeroFixed.length != expectedRows * zeroColumns)
            throw new IllegalArgumentException(
                "response and fixed designs must match the prepared rows");
        MatrixOps.requireFinite(response, "response");
        MatrixOps.requireFinite(countFixed, "countFixed");
        MatrixOps.requireFinite(zeroFixed, "zeroFixed");
        for (double value : response)
            if (value < 0.0 || value != Math.rint(value))
                throw new IllegalArgumentException(
                    "zero-inflated responses must be nonnegative integers");
        if (offset != null && offset.length != expectedRows)
            throw new IllegalArgumentException(
                "count offset length must match observations");
        if (family == CountFamily.NEGATIVE_BINOMIAL) {
            if (dispersionColumns < 1 || dispersionFixed == null
                    || dispersionFixed.length
                        != expectedRows * dispersionColumns)
                throw new IllegalArgumentException(
                    "ZINB requires a matching fixed-effects size design");
            MatrixOps.requireFinite(dispersionFixed, "dispersionFixed");
        } else if (dispersionColumns != 0 || dispersionFixed != null) {
            throw new IllegalArgumentException(
                "ZIP does not accept a dispersion design");
        }
    }

    private static void validateStructure(
            int rows,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> countPrecisionBases,
            List<RandomEffectTerm> zeroRandomEffects,
            List<SparsePrecisionMatrix> zeroPrecisionBases,
            List<CorrelatedZeroInflatedRandomEffect> correlatedEffects,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        if (rows < 1 || countRandomEffects == null
                || zeroRandomEffects == null || correlatedEffects == null
                || countRandomEffects.isEmpty() && zeroRandomEffects.isEmpty()
                    && correlatedEffects.isEmpty()
                || options == null || backendPolicy == null)
            throw new IllegalArgumentException(
                "random effects, controls, and backend are required");
        for (RandomEffectTerm term : countRandomEffects)
            if (term == null || term.observations() != rows)
                throw new IllegalArgumentException(
                    "count random-effect rows must match prepared rows");
        for (RandomEffectTerm term : zeroRandomEffects)
            if (term == null || term.observations() != rows)
                throw new IllegalArgumentException(
                    "zero random-effect rows must match prepared rows");
        for (CorrelatedZeroInflatedRandomEffect effect : correlatedEffects)
            if (effect == null || effect.countTerm().observations() != rows)
                throw new IllegalArgumentException(
                    "correlated random-effect rows must match prepared rows");
        validatePrecisions(countRandomEffects, countPrecisionBases);
        validatePrecisions(zeroRandomEffects, zeroPrecisionBases);
        int varianceCount = countRandomEffects.size() + zeroRandomEffects.size()
            + 2 * correlatedEffects.size();
        double[] supplied = options.initialVariances();
        if (supplied != null && supplied.length != varianceCount)
            throw new IllegalArgumentException(
                "one initial variance is required per independent process variance");
    }

    private static void validatePrecisions(
            List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases) {
        if (precisionBases != null) {
            if (precisionBases.size() != randomEffects.size())
                throw new IllegalArgumentException(
                    "one precision is required per random-effect term");
            for (int term = 0; term < randomEffects.size(); term++)
                if (precisionBases.get(term) == null
                        || precisionBases.get(term).dimension()
                            != randomEffects.get(term).coefficients())
                    throw new IllegalArgumentException(
                        "precision dimensions must match random effects");
        }
    }

    private static List<SparsePrecisionMatrix> identities(
            List<RandomEffectTerm> terms) {
        List<SparsePrecisionMatrix> result = new ArrayList<>(terms.size());
        for (RandomEffectTerm term : terms)
            result.add(SparsePrecisionMatrix.identity(term.coefficients()));
        return List.copyOf(result);
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

    private static RandomStructure randomStructure(
            List<RandomEffectTerm> countTerms,
            List<SparsePrecisionMatrix> countPrecisions,
            List<RandomEffectTerm> zeroTerms,
            List<SparsePrecisionMatrix> zeroPrecisions,
            List<CorrelatedZeroInflatedRandomEffect> correlated) {
        List<RandomEffectTerm> terms = new ArrayList<>();
        List<String> resultNames = new ArrayList<>();
        terms.addAll(countTerms);
        countTerms.forEach(term -> resultNames.add(term.name()));
        int[] pairedCountTerms = new int[correlated.size()];
        for (int index = 0; index < correlated.size(); index++) {
            pairedCountTerms[index] = terms.size();
            terms.add(correlated.get(index).countTerm());
            resultNames.add(correlated.get(index).name() + ":count");
        }
        int countTermCount = terms.size();
        terms.addAll(zeroTerms);
        zeroTerms.forEach(term -> resultNames.add("zero:" + term.name()));
        int[] pairedZeroTerms = new int[correlated.size()];
        for (int index = 0; index < correlated.size(); index++) {
            pairedZeroTerms[index] = terms.size();
            terms.add(correlated.get(index).zeroTerm());
            resultNames.add(correlated.get(index).name() + ":zero");
        }
        int independentCount = countTerms.size() + zeroTerms.size();
        int[] termVariance = new int[terms.size()];
        for (int index = 0; index < countTerms.size(); index++)
            termVariance[index] = index;
        for (int index = 0; index < zeroTerms.size(); index++)
            termVariance[countTermCount + index] = countTerms.size() + index;
        List<CorrelatedLayout> layouts = new ArrayList<>();
        for (int index = 0; index < correlated.size(); index++) {
            int variance = independentCount + 2 * index;
            termVariance[pairedCountTerms[index]] = variance;
            termVariance[pairedZeroTerms[index]] = variance + 1;
            layouts.add(new CorrelatedLayout(correlated.get(index).name(),
                pairedCountTerms[index], pairedZeroTerms[index],
                correlated.get(index).coefficientPrecision(), variance,
                variance + 1, index));
        }
        List<Integer> independentTerms = new ArrayList<>();
        List<SparsePrecisionMatrix> independentPrecisions = new ArrayList<>();
        for (int index = 0; index < countTerms.size(); index++) {
            independentTerms.add(index);
            independentPrecisions.add(countPrecisions.get(index));
        }
        for (int index = 0; index < zeroTerms.size(); index++) {
            independentTerms.add(countTermCount + index);
            independentPrecisions.add(zeroPrecisions.get(index));
        }
        return new RandomStructure(List.copyOf(terms),
            List.copyOf(resultNames), countTermCount,
            toIntArray(independentTerms), List.copyOf(independentPrecisions),
            List.copyOf(layouts), termVariance,
            independentCount + 2 * correlated.size(), correlated.size());
    }

    private static PrecisionData precisionData(
            RandomStructure structure, CombinedDesign design,
            ComputeBackend backend) {
        List<Integer> rows = new ArrayList<>();
        List<Integer> columns = new ArrayList<>();
        List<Double> values = new ArrayList<>();
        List<Integer> scaleIndices = new ArrayList<>();
        List<Integer> kinds = new ArrayList<>();
        int independentCount = structure.independentTerms().length;
        double[] independentLogDeterminants = new double[independentCount];
        int[] independentDimensions = new int[independentCount];
        for (int independent = 0; independent < independentCount; independent++) {
            int term = structure.independentTerms()[independent];
            SparsePrecisionMatrix precision =
                structure.independentPrecisions().get(independent);
            independentLogDeterminants[independent] =
                precisionLogDeterminant(precision, backend);
            independentDimensions[independent] = precision.dimension();
            addLowerPrecision(precision, design.termStarts()[term],
                structure.termVarianceIndices()[term], 0, rows, columns,
                values, scaleIndices, kinds);
        }
        int pairCount = structure.correlated().size();
        double[] correlatedLogDeterminants = new double[pairCount];
        int[] correlatedDimensions = new int[pairCount];
        for (int pair = 0; pair < pairCount; pair++) {
            CorrelatedLayout layout = structure.correlated().get(pair);
            SparsePrecisionMatrix precision = layout.precision();
            correlatedLogDeterminants[pair] =
                precisionLogDeterminant(precision, backend);
            correlatedDimensions[pair] = precision.dimension();
            addLowerPrecision(precision,
                design.termStarts()[layout.countTerm()], pair, 1,
                rows, columns, values, scaleIndices, kinds);
            addLowerPrecision(precision,
                design.termStarts()[layout.zeroTerm()], pair, 2,
                rows, columns, values, scaleIndices, kinds);
            int[] starts = precision.rowStarts();
            int[] cols = precision.columnIndices();
            double[] numeric = precision.values();
            int zeroStart = design.termStarts()[layout.zeroTerm()];
            int countStart = design.termStarts()[layout.countTerm()];
            for (int row = 0; row < precision.dimension(); row++)
                for (int entry = starts[row]; entry < starts[row + 1]; entry++) {
                    rows.add(zeroStart + row);
                    columns.add(countStart + cols[entry]);
                    values.add(numeric[entry]);
                    scaleIndices.add(pair);
                    kinds.add(3);
                }
        }
        return new PrecisionData(toIntArray(rows), toIntArray(columns),
            toDoubleArray(values), toIntArray(scaleIndices),
            toIntArray(kinds), independentLogDeterminants,
            independentDimensions, correlatedLogDeterminants,
            correlatedDimensions, structure);
    }

    private static void addLowerPrecision(
            SparsePrecisionMatrix precision, int termStart,
            int scaleIndex, int kind,
            List<Integer> rows, List<Integer> columns, List<Double> values,
            List<Integer> scaleIndices, List<Integer> kinds) {
        int[] starts = precision.rowStarts();
        int[] cols = precision.columnIndices();
        double[] numeric = precision.values();
        for (int row = 0; row < precision.dimension(); row++)
            for (int entry = starts[row]; entry < starts[row + 1]; entry++)
                if (cols[entry] <= row) {
                    rows.add(termStart + row);
                    columns.add(termStart + cols[entry]);
                    values.add(numeric[entry]);
                    scaleIndices.add(scaleIndex);
                    kinds.add(kind);
                }
    }

    private static double precisionLogDeterminant(
            SparsePrecisionMatrix precision, ComputeBackend backend) {
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
            for (int index = starts[row]; index < starts[row + 1]; index++)
                if (cols[index] <= row) {
                    lowerColumns[lowerPosition] = cols[index] + 1;
                    lowerValues[lowerPosition++] = numeric[index];
                }
        }
        lowerRows[precision.dimension()] = lowerPosition + 1;
        CsrMatrix lower = new CsrMatrix(precision.dimension(),
            precision.dimension(), lowerValues, lowerColumns, lowerRows);
        return backend.dcsrpotrf(lower, MatrixTriangle.LOWER,
            SparseOrdering.MINIMUM_DEGREE).logDeterminant();
    }

    private static CombinedDesign combine(
            List<RandomEffectTerm> terms, int countTerms, int rows) {
        int[] termStarts = new int[terms.size()];
        int columns = 0;
        int nonzeros = 0;
        for (int term = 0; term < terms.size(); term++) {
            termStarts[term] = columns;
            columns += terms.get(term).coefficients();
            nonzeros += terms.get(term).nonzeroCount();
        }
        int[] rowStarts = new int[rows + 1];
        int[] countEnds = new int[rows];
        int[] columnIndices = new int[nonzeros];
        double[] values = new double[nonzeros];
        TermData[] data = new TermData[terms.size()];
        for (int term = 0; term < terms.size(); term++)
            data[term] = termData(terms.get(term));
        int position = 0;
        for (int row = 0; row < rows; row++) {
            rowStarts[row] = position;
            for (int term = 0; term < terms.size(); term++) {
                if (term == countTerms) countEnds[row] = position;
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
            if (countTerms == terms.size()) countEnds[row] = position;
        }
        rowStarts[rows] = position;
        boolean[] countColumns = new boolean[columns];
        int countColumnLimit = countTerms == terms.size()
            ? columns : termStarts[countTerms];
        Arrays.fill(countColumns, 0, countColumnLimit, true);
        return new CombinedDesign(rows, columns, rowStarts,
            Arrays.copyOf(columnIndices, position),
            Arrays.copyOf(values, position), termStarts, countColumns,
            countEnds);
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

    private enum CountFamily { POISSON, NEGATIVE_BINOMIAL }

    private record RowLikelihood(
            double logLikelihood, double countScore,
            double zeroScore, double countCurvature,
            double zeroCurvature, double crossCurvature,
            double totalZeroProbability) { }

    private static final class RowLikelihoodWorkspace {
        private double logLikelihood;
        private double countScore;
        private double zeroScore;
        private double countCurvature;
        private double zeroCurvature;
        private double crossCurvature;
        private double totalZeroProbability;

        private void set(
                double logLikelihood, double countScore, double zeroScore,
                double countCurvature, double zeroCurvature,
                double crossCurvature, double totalZeroProbability) {
            this.logLikelihood = logLikelihood;
            this.countScore = countScore;
            this.zeroScore = zeroScore;
            this.countCurvature = countCurvature;
            this.zeroCurvature = zeroCurvature;
            this.crossCurvature = crossCurvature;
            this.totalZeroProbability = totalZeroProbability;
        }
    }

    private record DataState(
            double logLikelihood, double[] countScores, double[] zeroScores,
            double[] countCurvatures, double[] zeroCurvatures,
            double[] crossCurvatures, double[] means,
            double[] zeroProbabilities, double[] sizes, double[] responseMeans,
            double[] totalZeroProbabilities) { }

    private record Evaluation(
            double[] random, DataState data, double laplaceLogLikelihood,
            int iterations, boolean converged) { }

    private record InferenceData(
            List<String> parameterNames, double[] covariance) {
        static InferenceData unavailable() {
            return new InferenceData(List.of(), new double[0]);
        }
    }

    private record GradientEvaluation(
            double[] gradient, int calls, boolean complete) { }

    private record ObjectivePoint(double value, double[] random) { }

    private record TermData(
            int[] rowStarts, int[] columnIndices, double[] values) { }

    private record CombinedDesign(
            int rows, int columns, int[] rowStarts,
            int[] columnIndices, double[] values, int[] termStarts,
            boolean[] countColumns, int[] countEnds) {
        double rowProduct(
                int row, double[] coefficients, boolean countProcess) {
            double result = 0.0;
            int start = countProcess ? rowStarts[row] : countEnds[row];
            int end = countProcess ? countEnds[row] : rowStarts[row + 1];
            for (int index = start; index < end; index++)
                result += values[index] * coefficients[columnIndices[index]];
            return result;
        }
    }

    private record SparsePattern(
            int dimension, int[] rowStarts, int[] columnIndices) {
        CsrMatrix matrix(
                double[] countWeights, double[] zeroWeights,
                double[] crossWeights, CombinedDesign design,
                PrecisionData precision, double[] logVariances,
                double damping) {
            double[] numeric = new double[columnIndices.length];
            for (int observation = 0;
                    observation < design.rows(); observation++) {
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
                        boolean leftCount = design.countColumns()[
                            design.columnIndices()[left]];
                        boolean rightCount = design.countColumns()[
                            design.columnIndices()[right]];
                        double weight = leftCount == rightCount
                            ? leftCount ? countWeights[observation]
                                : zeroWeights[observation]
                            : crossWeights[observation];
                        numeric[position] += weight * design.values()[left]
                            * design.values()[right];
                    }
            }
            for (int index = 0; index < precision.values().length; index++) {
                int row = precision.rows()[index];
                int position = Arrays.binarySearch(columnIndices,
                    rowStarts[row], rowStarts[row + 1],
                    precision.columns()[index]);
                numeric[position] += precision.scaledValue(
                    index, logVariances);
            }
            if (damping > 0.0) {
                for (int row = 0; row < dimension; row++) {
                    int position = Arrays.binarySearch(columnIndices,
                        rowStarts[row], rowStarts[row + 1], row);
                    numeric[position] += damping;
                }
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
            int[] rows, int[] columns, double[] values, int[] scaleIndices,
            int[] kinds, double[] independentLogDeterminants,
            int[] independentDimensions, double[] correlatedLogDeterminants,
            int[] correlatedDimensions, RandomStructure structure) {
        double scaledValue(int entry, double[] covarianceParameters) {
            int kind = kinds[entry];
            if (kind == 0)
                return values[entry]
                    / Math.exp(covarianceParameters[scaleIndices[entry]]);
            CorrelatedLayout pair = structure.correlated().get(
                scaleIndices[entry]);
            double countVariance = Math.exp(
                covarianceParameters[pair.countVariance()]);
            double zeroVariance = Math.exp(
                covarianceParameters[pair.zeroVariance()]);
            double correlation = structure.correlation(
                covarianceParameters, pair.correlation());
            double denominator = 1.0 - correlation * correlation;
            if (kind == 1)
                return values[entry] / (countVariance * denominator);
            if (kind == 2)
                return values[entry] / (zeroVariance * denominator);
            return -values[entry] * correlation
                / (Math.sqrt(countVariance * zeroVariance) * denominator);
        }

        double logDeterminant(double[] logVariances) {
            double result = 0.0;
            for (int term = 0;
                    term < independentLogDeterminants.length; term++)
                result += independentLogDeterminants[term]
                    - independentDimensions[term] * logVariances[term];
            for (int index = 0;
                    index < correlatedLogDeterminants.length; index++) {
                CorrelatedLayout pair = structure.correlated().get(index);
                double correlation = structure.correlation(
                    logVariances, pair.correlation());
                result += 2.0 * correlatedLogDeterminants[index]
                    - correlatedDimensions[index]
                        * (logVariances[pair.countVariance()]
                            + logVariances[pair.zeroVariance()]
                            + Math.log(1.0 - correlation * correlation));
            }
            return result;
        }

        double quadratic(double[] random, double[] logVariances) {
            double result = 0.0;
            for (int index = 0; index < values.length; index++) {
                int row = rows[index];
                int column = columns[index];
                double contribution = scaledValue(index, logVariances)
                    * random[row] * random[column];
                result += row == column ? contribution : 2.0 * contribution;
            }
            return result;
        }

        void subtractProduct(
                double[] random, double[] logVariances, double[] gradient) {
            for (int index = 0; index < values.length; index++) {
                int row = rows[index];
                int column = columns[index];
                double scaled = scaledValue(index, logVariances);
                gradient[row] -= scaled * random[column];
                if (row != column) gradient[column] -= scaled * random[row];
            }
        }
    }

    private record CorrelatedLayout(
            String name, int countTerm, int zeroTerm,
            SparsePrecisionMatrix precision, int countVariance,
            int zeroVariance, int correlation) { }

    private record RandomStructure(
            List<RandomEffectTerm> terms, List<String> resultNames,
            int countTerms, int[] independentTerms,
            List<SparsePrecisionMatrix> independentPrecisions,
            List<CorrelatedLayout> correlated, int[] termVarianceIndices,
            int varianceCount, int correlationCount) {
        double correlation(double[] parameters, int index) {
            return 0.99 * Math.tanh(parameters[varianceCount + index]);
        }

        double[] variances(double[] parameters) {
            double[] result = new double[terms.size()];
            for (int term = 0; term < terms.size(); term++)
                result[term] = Math.exp(parameters[termVarianceIndices[term]]);
            return result;
        }

        Map<String, Double> correlations(double[] parameters) {
            Map<String, Double> result = new LinkedHashMap<>();
            for (CorrelatedLayout pair : correlated)
                result.put(pair.name(), correlation(
                    parameters, pair.correlation()));
            return result;
        }
    }
}
