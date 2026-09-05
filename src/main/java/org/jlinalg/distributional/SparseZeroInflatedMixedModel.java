/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
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
 * Frequentist ZIP and NB2-ZINB mixed models with count-side sparse Gaussian
 * random effects. Random coefficients are integrated by first-order Laplace
 * approximation; BOBYQA optimizes only the non-random parameters.
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
            countRandomEffects, precisionBases, countOffset,
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
            countRandomEffects, precisionBases, countOffset,
            CountFamily.NEGATIVE_BINOMIAL, options, backendPolicy);
    }

    private static ZeroInflatedMixedResult fit(
            double[] response,
            double[] countFixed, int countColumns,
            double[] zeroFixed, int zeroColumns,
            double[] dispersionFixed, int dispersionColumns,
            List<RandomEffectTerm> countRandomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            double[] countOffset,
            CountFamily family,
            ZeroInflatedMixedOptions options,
            BackendPolicy backendPolicy) {
        validate(response, countFixed, countColumns, zeroFixed, zeroColumns,
            dispersionFixed, dispersionColumns, countRandomEffects,
            precisionBases, countOffset, family, options, backendPolicy);
        int rows = response.length;
        double[] offsets = countOffset == null
            ? new double[rows]
            : MatrixOps.finiteCopy(countOffset, "countOffset");
        List<SparsePrecisionMatrix> precisions = precisionBases == null
            ? identities(countRandomEffects) : List.copyOf(precisionBases);

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            CombinedDesign design = combine(countRandomEffects, rows);
            PrecisionData precision = precisionData(
                countRandomEffects, precisions, design, backend);
            SparsePattern pattern = pattern(design, precision);
            double[] initial = initial(response, countFixed, countColumns,
                zeroFixed, zeroColumns, dispersionFixed, dispersionColumns,
                countRandomEffects.size(), family, options);
            double[] lower = lower(initial.length, countColumns, zeroColumns,
                dispersionColumns, countRandomEffects.size(), family, options);
            double[] upper = upper(initial.length, countColumns, zeroColumns,
                dispersionColumns, countRandomEffects.size(), family, options);
            double[] initialWeights = new double[rows];
            try (PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                    pattern.matrix(initialWeights, design, precision,
                        varianceParameters(initial, countColumns, zeroColumns,
                            dispersionColumns, family), 0.0),
                    MatrixTriangle.LOWER, SparseOrdering.MINIMUM_DEGREE)) {
                Objective objective = new Objective(response,
                    countFixed, countColumns, zeroFixed, zeroColumns,
                    dispersionFixed, dispersionColumns, offsets, family,
                    design, precision, pattern, factor, options);
                int interpolationPoints = Math.min(2 * initial.length + 1,
                    (initial.length + 1) * (initial.length + 2) / 2);
                OptimizationResult optimized;
                try {
                    optimized = Bobyqa.bobyqa(initial, lower, upper,
                        objective::value, interpolationPoints,
                        options.initialTrustRadius(),
                        Math.max(1e-7, options.relativeTolerance()),
                        options.maximumOuterEvaluations(), true);
                } catch (ArithmeticException | ArrayIndexOutOfBoundsException
                        exception) {
                    optimized = coordinateOptimize(initial, lower, upper,
                        objective, options);
                }
                double[] parameters = optimized.mX == null
                    ? initial : optimized.mX;
                Evaluation fitted = objective.evaluate(parameters);
                return result(fitted, parameters, countColumns, zeroColumns,
                    dispersionColumns, family, countRandomEffects, design,
                    pattern, factor.factorNonzeroCount(),
                    optimized.numFunctionCalls,
                    optimized.numFunctionCalls
                        < options.maximumOuterEvaluations());
            }
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

    private static ZeroInflatedMixedResult result(
            Evaluation fitted, double[] parameters,
            int countColumns, int zeroColumns, int dispersionColumns,
            CountFamily family, List<RandomEffectTerm> terms,
            CombinedDesign design, SparsePattern pattern,
            int factorNonzeros, int evaluations, boolean outerConverged) {
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
        double[] variances = new double[terms.size()];
        List<String> names = new ArrayList<>(terms.size());
        Map<String, double[]> random = new LinkedHashMap<>();
        for (int term = 0; term < terms.size(); term++) {
            RandomEffectTerm value = terms.get(term);
            names.add(value.name());
            variances[term] = Math.exp(parameters[varianceStart + term]);
            random.put(value.name(), Arrays.copyOfRange(fitted.random(),
                design.termStarts()[term],
                design.termStarts()[term] + value.coefficients()));
        }
        DataState data = fitted.data();
        return new ZeroInflatedMixedResult(
            family == CountFamily.POISSON
                ? "zero-inflated-poisson" : "zero-inflated-negative-binomial-2",
            count, zero, dispersion, names, variances, random,
            data.means(), data.zeroProbabilities(), data.sizes(),
            data.responseMeans(), data.totalZeroProbabilities(),
            fitted.laplaceLogLikelihood(), evaluations, fitted.iterations(),
            outerConverged && fitted.converged(), design.columns(),
            pattern.columnIndices().length, factorNonzeros);
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
                DataState data = data(parameters, random, true);
                double[] gradient = randomGradient(
                    data.scores(), random, logVariances);
                double gradientMaximum = maximumAbsolute(gradient);
                double damping = 0.0;
                boolean factored = false;
                for (int attempt = 0; attempt < 12; attempt++) {
                    try {
                        factor.refactor(pattern.matrix(data.curvatures(), design,
                            precision, logVariances, damping));
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
                    DataState trialData = data(parameters, trial, false);
                    double trialJoint = trialData.logLikelihood()
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
            factor.refactor(pattern.matrix(fitted.curvatures(), design,
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
                double[] parameters, double[] random, boolean derivatives) {
            int zeroStart = countColumns;
            int dispersionStart = zeroStart + zeroColumns;
            int varianceStart = dispersionStart
                + (family == CountFamily.NEGATIVE_BINOMIAL
                    ? dispersionColumns : 0);
            if (parameters.length != varianceStart
                    + precision.logDeterminants().length) {
                throw new IllegalArgumentException(
                    "outer parameter dimensions do not match the model");
            }
            int rows = response.length;
            double[] means = new double[rows];
            double[] zeroProbabilities = new double[rows];
            double[] sizes = family == CountFamily.NEGATIVE_BINOMIAL
                ? new double[rows] : new double[0];
            double[] responseMeans = new double[rows];
            double[] totalZeroProbabilities = new double[rows];
            double[] scores = derivatives ? new double[rows] : null;
            double[] curvatures = derivatives ? new double[rows] : null;
            double logLikelihood = 0.0;
            for (int row = 0; row < rows; row++) {
                double countEta = offsets[row]
                    + fixedValue(countFixed, row, countColumns, parameters, 0)
                    + design.rowProduct(row, random);
                double zeroEta = fixedValue(zeroFixed, row, zeroColumns,
                    parameters, zeroStart);
                double mean = safeExp(countEta);
                double probability = logistic(zeroEta);
                double size = family == CountFamily.NEGATIVE_BINOMIAL
                    ? safeExp(fixedValue(dispersionFixed, row,
                        dispersionColumns, parameters, dispersionStart))
                    : Double.POSITIVE_INFINITY;
                RowLikelihood likelihood = rowLikelihood(
                    response[row], mean, probability, size, family);
                logLikelihood += likelihood.logLikelihood();
                means[row] = mean;
                zeroProbabilities[row] = probability;
                if (sizes.length > 0) sizes[row] = size;
                responseMeans[row] = (1.0 - probability) * mean;
                totalZeroProbabilities[row] = likelihood.totalZeroProbability();
                if (derivatives) {
                    scores[row] = likelihood.countScore();
                    curvatures[row] = likelihood.countCurvature();
                }
            }
            return new DataState(logLikelihood, scores, curvatures, means,
                zeroProbabilities, sizes, responseMeans,
                totalZeroProbabilities);
        }

        private double[] randomGradient(
                double[] scores, double[] random, double[] logVariances) {
            double[] gradient = new double[design.columns()];
            for (int row = 0; row < design.rows(); row++)
                for (int index = design.rowStarts()[row];
                        index < design.rowStarts()[row + 1]; index++)
                    gradient[design.columnIndices()[index]] += scores[row]
                        * design.values()[index];
            precision.subtractProduct(random, logVariances, gradient);
            return gradient;
        }
    }

    private static RowLikelihood rowLikelihood(
            double response, double mean, double zeroProbability,
            double size, CountFamily family) {
        double logCount;
        double logCountZero;
        double baseZeroScore;
        double baseZeroScoreDerivative;
        double positiveScore;
        double positiveCurvature;
        if (family == CountFamily.POISSON) {
            logCountZero = -mean;
            logCount = response * Math.log(mean) - mean
                - SpecialFunctions.logGamma(response + 1.0);
            baseZeroScore = -mean;
            baseZeroScoreDerivative = -mean;
            positiveScore = response - mean;
            positiveCurvature = mean;
        } else {
            double total = size + mean;
            logCountZero = size * (Math.log(size) - Math.log(total));
            logCount = SpecialFunctions.logGamma(response + size)
                - SpecialFunctions.logGamma(size)
                - SpecialFunctions.logGamma(response + 1.0)
                + size * (Math.log(size) - Math.log(total))
                + response * (Math.log(mean) - Math.log(total));
            baseZeroScore = -size * mean / total;
            baseZeroScoreDerivative =
                -size * size * mean / (total * total);
            positiveScore = size * (response - mean) / total;
            positiveCurvature =
                size * mean * (size + response) / (total * total);
        }
        double logPi = logLogistic(zeroProbability);
        double logOneMinusPi = logOneMinusLogistic(zeroProbability);
        double logZeroMass = logAddExp(
            logPi, logOneMinusPi + logCountZero);
        double totalZero = Math.exp(logZeroMass);
        if (response > 0.0) {
            return new RowLikelihood(logOneMinusPi + logCount,
                positiveScore, positiveCurvature, totalZero);
        }
        double countPosterior = Math.exp(
            logOneMinusPi + logCountZero - logZeroMass);
        double score = countPosterior * baseZeroScore;
        double derivative = countPosterior
            * ((1.0 - countPosterior) * baseZeroScore * baseZeroScore
                + baseZeroScoreDerivative);
        return new RowLikelihood(logZeroMass, score, -derivative, totalZero);
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
            int varianceCount, CountFamily family,
            ZeroInflatedMixedOptions options) {
        int dispersionCount = family == CountFamily.NEGATIVE_BINOMIAL
            ? dispersionColumns : 0;
        double[] result = new double[
            countColumns + zeroColumns + dispersionCount + varianceCount];
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
            int dispersionColumns, int varianceCount, CountFamily family,
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
        return result;
    }

    private static double[] upper(
            int dimensions, int countColumns, int zeroColumns,
            int dispersionColumns, int varianceCount, CountFamily family,
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
            List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            double[] offset, CountFamily family,
            ZeroInflatedMixedOptions options, BackendPolicy backendPolicy) {
        if (response == null || response.length == 0
                || countColumns < 1 || zeroColumns < 1
                || countFixed == null
                || countFixed.length != response.length * countColumns
                || zeroFixed == null
                || zeroFixed.length != response.length * zeroColumns
                || randomEffects == null || randomEffects.isEmpty()
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
        for (RandomEffectTerm term : randomEffects)
            if (term == null || term.observations() != response.length
                    || !names.add(term.name()))
                throw new IllegalArgumentException(
                    "random effects need unique names and matching rows");
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
        double[] supplied = options.initialVariances();
        if (supplied != null && supplied.length != randomEffects.size())
            throw new IllegalArgumentException(
                "one initial variance is required per random-effect term");
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
            double countCurvature, double totalZeroProbability) { }

    private record DataState(
            double logLikelihood, double[] scores, double[] curvatures,
            double[] means, double[] zeroProbabilities, double[] sizes,
            double[] responseMeans, double[] totalZeroProbabilities) { }

    private record Evaluation(
            double[] random, DataState data, double laplaceLogLikelihood,
            int iterations, boolean converged) { }

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
                PrecisionData precision, double[] logVariances,
                double damping) {
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

        void subtractProduct(
                double[] random, double[] logVariances, double[] gradient) {
            for (int index = 0; index < values.length; index++) {
                int row = rows[index];
                int column = columns[index];
                double scaled = values[index]
                    / Math.exp(logVariances[terms[index]]);
                gradient[row] -= scaled * random[column];
                if (row != column) gradient[column] -= scaled * random[row];
            }
        }
    }
}
