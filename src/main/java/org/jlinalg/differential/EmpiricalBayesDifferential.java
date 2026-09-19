/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.Normal;
import jdistlib.T;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.OlsResult;

/**
 * Cross-feature moderated Gaussian/voom inference and a distinct
 * size-factor/negative-binomial estimator with dispersion shrinkage.
 */
public final class EmpiricalBayesDifferential {
    private static final double LOG_TWO = Math.log(2.0);
    private static final double MINIMUM_VARIANCE = 1e-12;

    private EmpiricalBayesDifferential() { }

    /** Fits limma-style moderated linear models to continuous feature rows. */
    public static DifferentialFit fitContinuous(double[][] features,
            double[][] design, double[] contrast) {
        Dimensions dimensions = validate(features, design, contrast, false);
        FeatureFit[] fits = new FeatureFit[features.length];
        double[] raw = new double[features.length];
        for (int feature = 0; feature < features.length; feature++) {
            fits[feature] = gaussian(features[feature], design, contrast, null);
            raw[feature] = fits[feature].variance;
        }
        Prior prior = variancePrior(raw, fits[0].degreesOfFreedom);
        return moderated("limma", fits, raw, prior,
            ones(dimensions.observations), false);
    }

    /**
     * Fits a voom-style log-count mean/variance trend, precision-weighted
     * linear models, and cross-feature moderated t statistics.
     */
    public static DifferentialFit fitVoom(double[][] counts,
            double[][] design, double[] contrast) {
        Dimensions dimensions = validate(counts, design, contrast, true);
        double[] sizeFactors = medianRatioSizeFactors(counts);
        double scale = geometricMean(Arrays.stream(counts)
            .mapToDouble(EmpiricalBayesDifferential::sum).toArray());
        double[][] logCpm = new double[counts.length][dimensions.observations];
        FeatureFit[] initial = new FeatureFit[counts.length];
        double[] average = new double[counts.length];
        double[] squareRootDeviation = new double[counts.length];
        for (int feature = 0; feature < counts.length; feature++) {
            for (int sample = 0; sample < dimensions.observations; sample++) {
                double library = Math.max(1.0, sizeFactors[sample] * scale);
                logCpm[feature][sample] = Math.log(
                    (counts[feature][sample] + 0.5) * 1e6 / (library + 1.0))
                    / LOG_TWO;
                average[feature] += logCpm[feature][sample]
                    / dimensions.observations;
            }
            initial[feature] = gaussian(logCpm[feature], design, contrast, null);
            squareRootDeviation[feature] = Math.sqrt(Math.sqrt(
                Math.max(initial[feature].variance, MINIMUM_VARIANCE)));
        }
        FeatureFit[] fits = new FeatureFit[counts.length];
        double[] raw = new double[counts.length];
        VoomTrend trendFit = new VoomTrend(average, squareRootDeviation,
            Arrays.stream(initial).flatMapToDouble(f -> Arrays.stream(f.fitted)).toArray());
        for (int feature = 0; feature < counts.length; feature++) {
            double[] weights = new double[dimensions.observations];
            for (int sample = 0; sample < weights.length; sample++) {
                double trend = trendFit.at(initial[feature].fitted[sample]);
                weights[sample] = 1.0 / Math.max(MINIMUM_VARIANCE,
                    trend * trend * trend * trend);
            }
            fits[feature] = gaussian(logCpm[feature], design, contrast, weights);
            raw[feature] = fits[feature].variance;
        }
        Prior prior = variancePrior(raw, fits[0].degreesOfFreedom);
        return moderated("voom", fits, raw, prior, sizeFactors, false);
    }

    /**
     * Fits median-ratio normalized negative-binomial log-link models with
     * empirical log-dispersion shrinkage and model-based Wald inference.
     */
    public static DifferentialFit fitNegativeBinomial(double[][] counts,
            double[][] design, double[] contrast) {
        Dimensions dimensions = validate(counts, design, contrast, true);
        double[] sizeFactors = medianRatioSizeFactors(counts);
        int residualDf = dimensions.observations - dimensions.parameters;
        double[] rawDispersion = new double[counts.length];
        double[] baseMean = new double[counts.length];
        for (int feature = 0; feature < counts.length; feature++) {
            double[] normalized = new double[dimensions.observations];
            for (int sample = 0; sample < normalized.length; sample++) {
                normalized[sample] = counts[feature][sample] / sizeFactors[sample];
                baseMean[feature] += normalized[sample] / normalized.length;
            }
            double variance = sampleVariance(normalized, baseMean[feature]);
            double poisson = 0.0;
            for (int sample = 0; sample < normalized.length; sample++)
                poisson += baseMean[feature] / sizeFactors[sample] / normalized.length;
            rawDispersion[feature] = Math.max(1e-8,
                (variance - poisson) / Math.max(MINIMUM_VARIANCE,
                    baseMean[feature] * baseMean[feature]));
        }
        Prior prior = logPrior(rawDispersion, residualDf);
        List<DifferentialResult> results = new ArrayList<>();
        for (int feature = 0; feature < counts.length; feature++) {
            double logDispersion = (residualDf * Math.log(rawDispersion[feature])
                + prior.degreesOfFreedom * Math.log(prior.variance))
                / (residualDf + prior.degreesOfFreedom);
            double dispersion = Math.exp(logDispersion);
            NegativeBinomialFit fit = negativeBinomial(counts[feature], design,
                contrast, sizeFactors, dispersion);
            double statistic = fit.converged ? fit.effect / fit.standardError : Double.NaN;
            double pValue = fit.converged ? Math.min(1.0, 2.0 * Normal.cumulative(
                Math.abs(statistic), 0.0, 1.0, false, false)) : Double.NaN;
            results.add(new DifferentialResult(feature, fit.effect,
                fit.effect / LOG_TWO, fit.standardError, statistic,
                fit.converged ? Double.POSITIVE_INFINITY : Double.NaN, pValue, baseMean[feature],
                rawDispersion[feature], dispersion, dispersion, Double.NaN,
                fit.converged, fit.iterations));
        }
        return new DifferentialFit("negative-binomial", results,
            prior.degreesOfFreedom, prior.variance, sizeFactors);
    }

    private static DifferentialFit moderated(String method, FeatureFit[] fits,
            double[] raw, Prior prior, double[] sizeFactors, boolean countScale) {
        List<DifferentialResult> results = new ArrayList<>();
        for (int feature = 0; feature < fits.length; feature++) {
            FeatureFit fit = fits[feature];
            double variance = (prior.degreesOfFreedom * prior.variance
                + fit.degreesOfFreedom * Math.max(raw[feature], MINIMUM_VARIANCE))
                / (prior.degreesOfFreedom + fit.degreesOfFreedom);
            double standardError = Math.sqrt(Math.max(MINIMUM_VARIANCE,
                fit.unscaledContrastVariance * variance));
            double statistic = fit.effect / standardError;
            double degrees = prior.degreesOfFreedom + fit.degreesOfFreedom;
            double pValue = Math.min(1.0, 2.0 * T.cumulative(
                Math.abs(statistic), degrees, false, false));
            results.add(new DifferentialResult(feature, fit.effect,
                countScale ? fit.effect / LOG_TWO : fit.effect,
                standardError, statistic, degrees, pValue, fit.baseMean,
                raw[feature], variance, Double.NaN, fit.meanWeight, true, 1));
        }
        return new DifferentialFit(method, results, prior.degreesOfFreedom,
            prior.variance, sizeFactors);
    }

    private static FeatureFit gaussian(double[] response, double[][] design,
            double[] contrast, double[] weights) {
        OlsResult fit = Ols.fit(response, design, weights, null,
            OlsOptions.defaults(), BackendPolicy.CPU);
        double effect = dot(contrast, fit.coefficients());
        double[] covariance = fit.covariance();
        double variance = Math.max(fit.residualVariance(), MINIMUM_VARIANCE);
        double unscaled = quadratic(contrast, covariance) / variance;
        if (fit.residualVariance() < MINIMUM_VARIANCE) {
            int columns = design[0].length;
            double[] weighted = new double[design.length * columns];
            for (int row = 0; row < design.length; row++) for (int j = 0; j < columns; j++)
                weighted[row * columns + j] = design[row][j] * Math.sqrt(weights == null ? 1 : weights[row]);
            try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
                unscaled = quadratic(contrast, LeastSquaresSolver.solve(weighted,
                    new double[design.length], design.length, columns, false,
                    context.backend()).unscaledCovariance());
            }
        }
        double mean = Arrays.stream(response).average().orElseThrow();
        double meanWeight = weights == null ? 1.0
            : Arrays.stream(weights).average().orElseThrow();
        return new FeatureFit(effect, unscaled, fit.residualVariance(),
            fit.residualDegreesOfFreedom(), mean, meanWeight,
            fit.fittedValues());
    }

    private static NegativeBinomialFit negativeBinomial(double[] response,
            double[][] design, double[] contrast, double[] sizeFactors,
            double dispersion) {
        int observations = response.length;
        double[] offset = new double[observations];
        double[] transformed = new double[observations];
        for (int sample = 0; sample < observations; sample++) {
            offset[sample] = Math.log(sizeFactors[sample]);
            transformed[sample] = Math.log(response[sample] + 0.5);
        }
        OlsResult starting = Ols.fit(transformed, design, null, offset,
            OlsOptions.defaults(), BackendPolicy.CPU);
        double[] beta = starting.coefficients();
        OlsResult last = starting;
        boolean converged = false;
        int iteration = 0;
        for (; iteration < 100; iteration++) {
            double[] working = new double[observations];
            double[] weights = new double[observations];
            for (int sample = 0; sample < observations; sample++) {
                double eta = offset[sample] + dot(design[sample], beta);
                double mean = Math.exp(Math.max(-30.0, Math.min(30.0, eta)));
                working[sample] = eta + (response[sample] - mean) / mean;
                weights[sample] = Math.max(1e-10,
                    mean / (1.0 + dispersion * mean));
            }
            last = Ols.fit(working, design, weights, offset,
                OlsOptions.defaults(), BackendPolicy.CPU);
            double[] updated = last.coefficients();
            double change = 0.0;
            for (int column = 0; column < beta.length; column++)
                change = Math.max(change, Math.abs(updated[column] - beta[column]));
            beta = updated;
            if (change < 1e-8) { converged = true; iteration++; break; }
        }
        double effect = dot(contrast, beta);
        double standardError = Double.NaN;
        if (converged) {
            int columns = beta.length;
            double[] weightedDesign = new double[observations * columns];
            for (int row = 0; row < observations; row++) {
                double eta = offset[row] + dot(design[row], beta);
                if (!Double.isFinite(eta) || Math.abs(eta) >= 30) { converged = false; break; }
                double mean = Math.exp(eta);
                double rootWeight = Math.sqrt(mean / (1 + dispersion * mean));
                for (int column = 0; column < columns; column++)
                    weightedDesign[row * columns + column] = rootWeight * design[row][column];
            }
            if (converged) try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
                double[] inverseInformation = LeastSquaresSolver.solve(weightedDesign,
                    new double[observations], observations, columns, false,
                    context.backend()).unscaledCovariance();
                standardError = Math.sqrt(quadratic(contrast, inverseInformation));
                if (!(standardError > 0) || !Double.isFinite(standardError)) converged = false;
            }
        }
        if (!converged) standardError = Double.NaN;
        return new NegativeBinomialFit(effect, standardError,
            converged, iteration);
    }

    private static Prior variancePrior(double[] variance, int residualDf) {
        double[] log = Arrays.stream(variance)
            .map(value -> Math.log(Math.max(value, MINIMUM_VARIANCE))).toArray();
        double mean = Arrays.stream(log).average().orElseThrow();
        double observed = sampleVariance(log, mean);
        double residualTrigamma = trigamma(residualDf / 2.0);
        double priorDf = observed <= residualTrigamma + 1e-10
            ? 1e6 : 2.0 * inverseTrigamma(observed - residualTrigamma);
        double adjustedMean = mean - digamma(residualDf / 2.0)
            + Math.log(residualDf / 2.0);
        double priorVariance = Math.exp(adjustedMean
            + digamma(priorDf / 2.0) - Math.log(priorDf / 2.0));
        return new Prior(priorDf, Math.max(MINIMUM_VARIANCE, priorVariance));
    }

    private static Prior logPrior(double[] values, int residualDf) {
        double[] log = Arrays.stream(values).map(Math::log).toArray();
        double mean = Arrays.stream(log).average().orElseThrow();
        double observed = sampleVariance(log, mean);
        double sampling = trigamma(Math.max(1.0, residualDf / 2.0));
        double priorDf = observed <= sampling + 1e-10
            ? 1e6 : 2.0 * inverseTrigamma(observed - sampling);
        return new Prior(priorDf, Math.exp(mean));
    }

    private static double[] medianRatioSizeFactors(double[][] counts) {
        int samples = counts[0].length;
        List<List<Double>> ratios = new ArrayList<>();
        for (int sample = 0; sample < samples; sample++) ratios.add(new ArrayList<>());
        for (double[] feature : counts) {
            boolean positive = true;
            double logSum = 0.0;
            for (double value : feature) {
                if (!(value > 0.0)) { positive = false; break; }
                logSum += Math.log(value);
            }
            if (!positive) continue;
            double geometric = Math.exp(logSum / samples);
            for (int sample = 0; sample < samples; sample++)
                ratios.get(sample).add(feature[sample] / geometric);
        }
        double[] factors = new double[samples];
        if (ratios.get(0).isEmpty()) {
            for (int sample = 0; sample < samples; sample++) {
                for (double[] feature : counts) factors[sample] += feature[sample];
                factors[sample] = Math.max(1.0, factors[sample]);
            }
        } else {
            for (int sample = 0; sample < samples; sample++)
                factors[sample] = median(ratios.get(sample));
        }
        double geometric = geometricMean(factors);
        for (int sample = 0; sample < samples; sample++) factors[sample] /= geometric;
        return factors;
    }

    private static Dimensions validate(double[][] features, double[][] design,
            double[] contrast, boolean counts) {
        if (features == null || features.length < 2 || features[0] == null
                || features[0].length < 3)
            throw new IllegalArgumentException(
                "at least two features and three observations are required");
        int observations = features[0].length;
        for (double[] feature : features) {
            if (feature == null || feature.length != observations)
                throw new IllegalArgumentException("feature rows must be rectangular");
            for (double value : feature) {
                if (!Double.isFinite(value) || (counts && (value < 0.0
                        || Math.abs(value - Math.rint(value)) > 1e-8)))
                    throw new IllegalArgumentException(counts
                        ? "counts must be finite nonnegative integers"
                        : "features must be finite");
            }
        }
        if (design == null || design.length != observations || design[0] == null
                || design[0].length == 0)
            throw new IllegalArgumentException("design rows must align with observations");
        int parameters = design[0].length;
        if (observations <= parameters)
            throw new IllegalArgumentException("design needs residual degrees of freedom");
        for (double[] row : design) {
            if (row == null || row.length != parameters)
                throw new IllegalArgumentException("design must be rectangular");
            for (double value : row) if (!Double.isFinite(value))
                throw new IllegalArgumentException("design must be finite");
        }
        if (contrast == null || contrast.length != parameters)
            throw new IllegalArgumentException("contrast length must match design columns");
        for (double value : contrast) if (!Double.isFinite(value))
            throw new IllegalArgumentException("contrast must be finite");
        return new Dimensions(observations, parameters);
    }

    private static double quadratic(double[] vector, double[] matrix) {
        int size = vector.length;
        if (matrix.length != size * size)
            throw new IllegalArgumentException("covariance dimension differs from contrast");
        double result = 0.0;
        for (int row = 0; row < size; row++)
            for (int column = 0; column < size; column++)
                result = Math.fma(vector[row] * matrix[row * size + column],
                    vector[column], result);
        return result;
    }

    private static double dot(double[] left, double[] right) {
        if (left.length != right.length)
            throw new IllegalArgumentException("vector dimensions differ");
        double result = 0.0;
        for (int index = 0; index < left.length; index++)
            result = Math.fma(left[index], right[index], result);
        return result;
    }

    private static double sampleVariance(double[] values, double mean) {
        if (values.length < 2) return 0.0;
        double result = 0.0;
        for (double value : values) result += (value - mean) * (value - mean);
        return result / (values.length - 1.0);
    }

    private static double sum(double[] values) {
        double result = 0.0;
        for (double value : values) result += value;
        return result;
    }

    private static double geometricMean(double[] values) {
        double log = 0.0;
        for (double value : values) log += Math.log(Math.max(value, 1e-12));
        return Math.exp(log / values.length);
    }

    private static double median(List<Double> values) {
        double[] copy = values.stream().mapToDouble(Double::doubleValue).sorted().toArray();
        int middle = copy.length / 2;
        return copy.length % 2 == 0 ? (copy[middle - 1] + copy[middle]) / 2.0
            : copy[middle];
    }

    private static double[] ones(int size) {
        double[] result = new double[size];
        Arrays.fill(result, 1.0);
        return result;
    }

    private static double inverseTrigamma(double target) {
        double lower = 1e-6;
        double upper = 1e6;
        for (int iteration = 0; iteration < 100; iteration++) {
            double middle = Math.sqrt(lower * upper);
            if (trigamma(middle) > target) lower = middle;
            else upper = middle;
        }
        return Math.sqrt(lower * upper);
    }

    private static double digamma(double value) {
        double result = 0.0;
        double x = value;
        while (x < 8.0) { result -= 1.0 / x; x += 1.0; }
        double inverse = 1.0 / x;
        double squared = inverse * inverse;
        return result + Math.log(x) - 0.5 * inverse
            - squared * (1.0 / 12.0 - squared * (1.0 / 120.0
                - squared * (1.0 / 252.0)));
    }

    private static double trigamma(double value) {
        double result = 0.0;
        double x = value;
        while (x < 8.0) { result += 1.0 / (x * x); x += 1.0; }
        double inverse = 1.0 / x;
        double squared = inverse * inverse;
        return result + inverse + 0.5 * squared + inverse * squared / 6.0
            - inverse * squared * squared / 30.0
            + inverse * squared * squared * squared / 42.0;
    }

    private record Dimensions(int observations, int parameters) { }
    private record Prior(double degreesOfFreedom, double variance) { }
    private record FeatureFit(double effect, double unscaledContrastVariance,
        double variance, int degreesOfFreedom, double baseMean,
        double meanWeight, double[] fitted) { }
    private record NegativeBinomialFit(double effect, double standardError,
        boolean converged, int iterations) { }
}
