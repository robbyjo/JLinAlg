/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Objects;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import jdistlib.ChiSquare;
import jdistlib.Normal;
import jdistlib.T;

/**
 * Prepared, allocation-light batch of independent intercept-only
 * meta-analyses. Inputs are row-major: all study values for analysis zero,
 * followed by all study values for analysis one, and so on.
 */
public final class PreparedMetaAnalysisBatch {
    private static final double LOG_TWO = Math.log(2.0);
    private static final double LOG_TEN = Math.log(10.0);
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    private static final double GOLDEN_RATIO = (Math.sqrt(5.0) - 1.0) / 2.0;
    private static final int PARALLEL_CHUNK_SIZE = 2048;

    private final double[] effects;
    private final double[] variances;
    private final int analyses;
    private final int studies;

    /** Validates and copies row-major effects and sampling standard errors. */
    public PreparedMetaAnalysisBatch(
            double[] effects, double[] standardErrors,
            int analyses, int studies) {
        Objects.requireNonNull(effects, "effects");
        Objects.requireNonNull(standardErrors, "standardErrors");
        if (analyses < 1)
            throw new IllegalArgumentException("analyses must be positive");
        if (studies < 2)
            throw new IllegalArgumentException("at least two studies are required");
        int length;
        try {
            length = Math.multiplyExact(analyses, studies);
        } catch (ArithmeticException exception) {
            throw new IllegalArgumentException("batch dimensions overflow", exception);
        }
        if (effects.length != length || standardErrors.length != length)
            throw new IllegalArgumentException(
                "row-major inputs must have analyses * studies values");
        this.effects = effects.clone();
        this.variances = new double[length];
        for (int index = 0; index < length; index++) {
            double effect = this.effects[index];
            double standardError = standardErrors[index];
            if (!Double.isFinite(effect))
                throw new IllegalArgumentException("effect sizes must be finite");
            if (!(standardError > 0.0) || !Double.isFinite(standardError))
                throw new IllegalArgumentException(
                    "sampling standard errors must be finite and positive");
            variances[index] = standardError * standardError;
        }
        this.analyses = analyses;
        this.studies = studies;
    }

    public int analyses() { return analyses; }
    public int studies() { return studies; }

    /** Fits with one deterministic worker thread. */
    public MetaAnalysisBatchResult fit(MetaAnalysisOptions options) {
        return fit(options, 1);
    }

    /**
     * Fits every row, splitting independent row chunks across the requested
     * number of worker threads.
     */
    public MetaAnalysisBatchResult fit(
            MetaAnalysisOptions options, int parallelism) {
        Objects.requireNonNull(options, "options");
        if (parallelism < 1)
            throw new IllegalArgumentException("parallelism must be positive");
        Output output = new Output(analyses,
            options.method() == MetaAnalysisMethod.RANDOM_EFFECT);
        double degreesOfFreedom = studies - 1.0;
        double critical = MetaAnalysis.critical(options, degreesOfFreedom);
        if (parallelism == 1 || analyses <= PARALLEL_CHUNK_SIZE) {
            fitRange(0, analyses, options, critical, degreesOfFreedom, output);
        } else {
            int chunks = (analyses + PARALLEL_CHUNK_SIZE - 1)
                / PARALLEL_CHUNK_SIZE;
            ForkJoinPool pool = new ForkJoinPool(Math.min(parallelism, chunks));
            try {
                pool.submit(() -> IntStream.range(0, chunks).parallel()
                    .forEach(chunk -> {
                        int from = chunk * PARALLEL_CHUNK_SIZE;
                        int to = Math.min(analyses,
                            from + PARALLEL_CHUNK_SIZE);
                        fitRange(from, to, options, critical,
                            degreesOfFreedom, output);
                    })).get();
            } catch (InterruptedException exception) {
                Thread.currentThread().interrupt();
                throw new IllegalStateException(
                    "meta-analysis batch was interrupted", exception);
            } catch (ExecutionException exception) {
                Throwable cause = exception.getCause();
                if (cause instanceof RuntimeException runtime) throw runtime;
                throw new IllegalStateException(
                    "meta-analysis batch failed", cause);
            } finally {
                pool.shutdown();
            }
        }
        return output.result(options);
    }

    private void fitRange(
            int from, int to, MetaAnalysisOptions options,
            double critical, double degreesOfFreedom, Output output) {
        for (int analysis = from; analysis < to; analysis++)
            fitOne(analysis, options, critical, degreesOfFreedom, output);
    }

    private void fitOne(
            int analysis, MetaAnalysisOptions options,
            double critical, double degreesOfFreedom, Output output) {
        int offset = analysis * studies;
        Sums fixed = sums(offset, 0.0);
        double q = fixed.q();
        double tauSquared = estimateTauSquared(offset, fixed, options);
        Sums fitted = tauSquared == 0.0 ? fixed : sums(offset, tauSquared);
        double scale = MetaAnalysis.inferenceScale(
            options, fitted.q(), degreesOfFreedom);
        double standardError = Math.sqrt(scale / fitted.sumWeight());
        double beta = fitted.mean();
        double statistic = beta / standardError;
        Probability probability = probability(
            statistic, degreesOfFreedom, options.inferenceMethod());

        output.pooledEffectSizes[analysis] = beta;
        output.standardErrors[analysis] = standardError;
        output.statistics[analysis] = statistic;
        output.pValues[analysis] = probability.pValue();
        output.negativeLog10PValues[analysis] =
            probability.negativeLog10PValue();
        output.confidenceLower[analysis] = beta - critical * standardError;
        output.confidenceUpper[analysis] = beta + critical * standardError;
        if (options.method() == MetaAnalysisMethod.RANDOM_EFFECT) {
            double predictionStandardError = Math.sqrt(
                tauSquared + standardError * standardError);
            output.predictionLower[analysis] =
                beta - critical * predictionStandardError;
            output.predictionUpper[analysis] =
                beta + critical * predictionStandardError;
        }
        output.cochranQ[analysis] = q;
        output.cochranQPValues[analysis] = ChiSquare.cumulative(
            q, degreesOfFreedom, false, false);
        output.tauSquared[analysis] = tauSquared;
        output.iSquared[analysis] = q > 0.0
            ? 100.0 * Math.max(0.0, (q - degreesOfFreedom) / q) : 0.0;
        output.hSquared[analysis] = Math.max(1.0, q / degreesOfFreedom);
    }

    private double estimateTauSquared(
            int offset, Sums fixed, MetaAnalysisOptions options) {
        if (options.method() == MetaAnalysisMethod.FIXED_EFFECT) return 0.0;
        return switch (options.tauSquaredEstimator()) {
            case DERSIMONIAN_LAIRD -> dersimonianLaird(offset, fixed);
            case PAULE_MANDEL -> pauleMandel(offset, fixed, options);
            case REML -> reml(offset, options);
        };
    }

    private double dersimonianLaird(int offset, Sums fixed) {
        double sumSquaredWeight = 0.0;
        for (int study = 0; study < studies; study++) {
            double weight = 1.0 / variances[offset + study];
            sumSquaredWeight += weight * weight;
        }
        double denominator = fixed.sumWeight()
            - sumSquaredWeight / fixed.sumWeight();
        return denominator > 0.0
            ? Math.max(0.0, (fixed.q() - (studies - 1.0)) / denominator)
            : 0.0;
    }

    private double pauleMandel(
            int offset, Sums fixed, MetaAnalysisOptions options) {
        double target = studies - 1.0;
        if (fixed.q() <= target) return 0.0;
        double upper = startingUpper(offset);
        while (sums(offset, upper).q() > target && upper < 1e12)
            upper *= 4.0;
        double lower = 0.0;
        for (int iteration = 0;
                iteration < options.maximumIterations(); iteration++) {
            double middle = 0.5 * (lower + upper);
            if (sums(offset, middle).q() > target) lower = middle;
            else upper = middle;
            if (upper - lower <= options.tolerance()
                    * Math.max(1.0, upper)) break;
        }
        return 0.5 * (lower + upper);
    }

    private double reml(int offset, MetaAnalysisOptions options) {
        double upper = startingUpper(offset);
        double atZero = restrictedObjective(offset, 0.0);
        double previous = atZero;
        double atUpper = restrictedObjective(offset, upper);
        while (atUpper < previous && upper < 1e12) {
            previous = atUpper;
            upper *= 4.0;
            atUpper = restrictedObjective(offset, upper);
        }
        double left = 0.0;
        double right = upper;
        double first = right - GOLDEN_RATIO * (right - left);
        double second = left + GOLDEN_RATIO * (right - left);
        double firstValue = restrictedObjective(offset, first);
        double secondValue = restrictedObjective(offset, second);
        for (int iteration = 0;
                iteration < options.maximumIterations(); iteration++) {
            if (firstValue < secondValue) {
                right = second;
                second = first;
                secondValue = firstValue;
                first = right - GOLDEN_RATIO * (right - left);
                firstValue = restrictedObjective(offset, first);
            } else {
                left = first;
                first = second;
                firstValue = secondValue;
                second = left + GOLDEN_RATIO * (right - left);
                secondValue = restrictedObjective(offset, second);
            }
            if (right - left <= options.tolerance()
                    * Math.max(1.0, right)) break;
        }
        double candidate = 0.5 * (left + right);
        return atZero <= restrictedObjective(offset, candidate)
            ? 0.0 : candidate;
    }

    private double restrictedObjective(int offset, double tauSquared) {
        double sumLogVariance = 0.0;
        for (int study = 0; study < studies; study++) {
            int index = offset + study;
            sumLogVariance += Math.log(variances[index] + tauSquared);
        }
        Sums values = sums(offset, tauSquared);
        return sumLogVariance + Math.log(values.sumWeight()) + values.q()
            + (studies - 1.0) * LOG_TWO_PI;
    }

    private double startingUpper(int offset) {
        double mean = 0.0;
        for (int study = 0; study < studies; study++)
            mean += effects[offset + study];
        mean /= studies;
        double variance = 0.0;
        for (int study = 0; study < studies; study++) {
            double difference = effects[offset + study] - mean;
            variance += difference * difference;
        }
        return Math.max(1e-8, variance / (studies - 1.0));
    }

    private Sums sums(int offset, double tauSquared) {
        double sumWeight = 0.0;
        double sumWeightedEffect = 0.0;
        for (int study = 0; study < studies; study++) {
            int index = offset + study;
            double weight = 1.0 / (variances[index] + tauSquared);
            double effect = effects[index];
            sumWeight += weight;
            sumWeightedEffect += weight * effect;
        }
        double mean = sumWeightedEffect / sumWeight;
        double q = 0.0;
        for (int study = 0; study < studies; study++) {
            int index = offset + study;
            double weight = 1.0 / (variances[index] + tauSquared);
            double residual = effects[index] - mean;
            q += weight * residual * residual;
        }
        return new Sums(sumWeight, mean, q);
    }

    private static Probability probability(
            double statistic, double degreesOfFreedom,
            MetaInferenceMethod inferenceMethod) {
        double magnitude = Math.abs(statistic);
        double logUpper = inferenceMethod == MetaInferenceMethod.NORMAL
            ? Normal.cumulative(magnitude, 0.0, 1.0, false, true)
            : T.cumulative(magnitude, degreesOfFreedom, false, true);
        double logP = Math.min(0.0, LOG_TWO + logUpper);
        return new Probability(Math.exp(logP), -logP / LOG_TEN);
    }

    private record Sums(double sumWeight, double mean, double q) { }
    private record Probability(double pValue, double negativeLog10PValue) { }

    private static final class Output {
        private final double[] pooledEffectSizes;
        private final double[] standardErrors;
        private final double[] statistics;
        private final double[] pValues;
        private final double[] negativeLog10PValues;
        private final double[] confidenceLower;
        private final double[] confidenceUpper;
        private final double[] predictionLower;
        private final double[] predictionUpper;
        private final double[] cochranQ;
        private final double[] cochranQPValues;
        private final double[] tauSquared;
        private final double[] iSquared;
        private final double[] hSquared;

        Output(int analyses, boolean randomEffects) {
            pooledEffectSizes = new double[analyses];
            standardErrors = new double[analyses];
            statistics = new double[analyses];
            pValues = new double[analyses];
            negativeLog10PValues = new double[analyses];
            confidenceLower = new double[analyses];
            confidenceUpper = new double[analyses];
            predictionLower = new double[analyses];
            predictionUpper = new double[analyses];
            cochranQ = new double[analyses];
            cochranQPValues = new double[analyses];
            tauSquared = new double[analyses];
            iSquared = new double[analyses];
            hSquared = new double[analyses];
            if (!randomEffects) {
                java.util.Arrays.fill(predictionLower, Double.NaN);
                java.util.Arrays.fill(predictionUpper, Double.NaN);
            }
        }

        MetaAnalysisBatchResult result(MetaAnalysisOptions options) {
            return new MetaAnalysisBatchResult(options.method(),
                options.tauSquaredEstimator(), options.inferenceMethod(),
                pooledEffectSizes, standardErrors, statistics, pValues,
                negativeLog10PValues, confidenceLower, confidenceUpper,
                predictionLower, predictionUpper, cochranQ, cochranQPValues,
                tauSquared, iSquared, hSquared);
        }
    }
}
