/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.Arrays;
import java.util.List;

/** Built-in location-scale and vector-response distributional families. */
public final class DistributionalFamilies {
    private static final DistributionalFamily GAUSSIAN_LOCATION_SCALE =
        new GaussianLocationScale();

    private DistributionalFamilies() { }

    /** Gaussian response with modeled mean and log standard deviation. */
    public static DistributionalFamily gaussianLocationScale() {
        return GAUSSIAN_LOCATION_SCALE;
    }

    /** Baseline-category multinomial logit with {@code categories - 1} predictors. */
    public static DistributionalFamily multinomial(int categories) {
        return new Multinomial(categories);
    }

    /** Gamma response with modeled log mean and log shape. */
    public static DistributionalFamily gammaMeanShape() {
        return new GammaMeanShapeFamily();
    }

    /** Beta response with modeled logit mean and log precision. */
    public static DistributionalFamily betaMeanPrecision() {
        return new BetaMeanPrecisionFamily();
    }

    /** NB2 counts with modeled log mean and log size. */
    public static DistributionalFamily negativeBinomialMeanDispersion() {
        return new NegativeBinomialMeanDispersionFamily();
    }

    /** Poisson counts with a modeled structural-zero probability. */
    public static DistributionalFamily zeroInflatedPoisson() {
        return new ZeroInflatedPoissonFamily();
    }

    /** Zero-hurdle probability plus a zero-truncated Poisson count model. */
    public static DistributionalFamily hurdlePoisson() {
        return new HurdlePoissonFamily();
    }

    /** Adjacent-category ordinal logits, compatible with VGAM's acat family. */
    public static DistributionalFamily ordinalAdjacentCategories(int categories) {
        return new AdjacentCategoryOrdinalFamily(categories);
    }

    private static final class GaussianLocationScale
            implements DistributionalFamily {
        private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
        private static final double MINIMUM_SCALE = 1e-10;

        @Override public String name() { return "gaussian-location-scale"; }
        @Override public int parameterCount() { return 2; }
        @Override public List<String> parameterNames() {
            return List.of("mu", "sigma");
        }
        @Override public void validateResponse(double response) {
            if (!Double.isFinite(response)) {
                throw new IllegalArgumentException(
                    "Gaussian responses must be finite");
            }
        }
        @Override public double[] initialParameters(double[] response) {
            double mean = Arrays.stream(response).average().orElseThrow();
            double sumSquares = 0.0;
            for (double value : response) {
                double difference = value - mean;
                sumSquares += difference * difference;
            }
            double scale = Math.sqrt(sumSquares / Math.max(1, response.length - 1));
            return new double[] {mean, Math.max(MINIMUM_SCALE, scale)};
        }
        @Override public double link(int parameter, double value) {
            checkParameter(parameter);
            if (parameter == 0) return value;
            if (!(value > 0.0)) {
                throw new IllegalArgumentException("sigma must be positive");
            }
            return Math.log(value);
        }
        @Override public double inverseLink(int parameter, double predictor) {
            checkParameter(parameter);
            return parameter == 0 ? predictor
                : Math.max(MINIMUM_SCALE, Math.exp(Math.min(350.0, predictor)));
        }
        @Override public double logLikelihood(
                double response, double[] parameters) {
            double residual = response - parameters[0];
            double variance = parameters[1] * parameters[1];
            return -0.5 * (LOG_TWO_PI + Math.log(variance)
                + residual * residual / variance);
        }
        @Override public void derivatives(
                double response,
                double[] parameters,
                double[] score,
                double[] information) {
            double residual = response - parameters[0];
            double variance = parameters[1] * parameters[1];
            score[0] = residual / variance;
            score[1] = -1.0 + residual * residual / variance;
            Arrays.fill(information, 0.0);
            information[0] = 1.0 / variance;
            information[3] = 2.0;
        }
        private static void checkParameter(int parameter) {
            if (parameter < 0 || parameter >= 2) {
                throw new IllegalArgumentException("unknown Gaussian parameter");
            }
        }
    }

    private static final class Multinomial implements DistributionalFamily {
        private final int categories;
        private final int predictors;
        private final List<String> names;

        private Multinomial(int categories) {
            if (categories < 2) {
                throw new IllegalArgumentException(
                    "multinomial family requires at least two categories");
            }
            this.categories = categories;
            this.predictors = categories - 1;
            java.util.ArrayList<String> values = new java.util.ArrayList<>();
            for (int category = 0; category < predictors; category++) {
                values.add("logit" + category);
            }
            this.names = List.copyOf(values);
        }

        @Override public String name() {
            return "multinomial-logit(" + categories + ")";
        }
        @Override public int parameterCount() { return predictors; }
        @Override public List<String> parameterNames() { return names; }
        @Override public void validateResponse(double response) {
            if (response != Math.rint(response)
                    || response < 0.0 || response >= categories) {
                throw new IllegalArgumentException(
                    "multinomial responses must be integer category indices");
            }
        }
        @Override public double[] initialParameters(double[] response) {
            double[] counts = new double[categories];
            Arrays.fill(counts, 0.5);
            for (double value : response) counts[(int) value] += 1.0;
            double[] result = new double[predictors];
            for (int category = 0; category < predictors; category++) {
                result[category] = Math.log(
                    counts[category] / counts[categories - 1]);
            }
            return result;
        }
        @Override public double link(int parameter, double value) {
            checkParameter(parameter);
            return value;
        }
        @Override public double inverseLink(int parameter, double predictor) {
            checkParameter(parameter);
            return predictor;
        }
        @Override public double logLikelihood(
                double response, double[] parameters) {
            Probabilities probabilities = probabilities(parameters);
            int category = (int) response;
            return category == categories - 1
                ? -probabilities.logDenominator()
                : parameters[category] - probabilities.logDenominator();
        }
        @Override public void derivatives(
                double response,
                double[] parameters,
                double[] score,
                double[] information) {
            Probabilities values = probabilities(parameters);
            int observed = (int) response;
            for (int row = 0; row < predictors; row++) {
                double probability = values.probabilities()[row];
                score[row] = (observed == row ? 1.0 : 0.0) - probability;
                for (int column = 0; column < predictors; column++) {
                    information[row * predictors + column] = probability
                        * ((row == column ? 1.0 : 0.0)
                            - values.probabilities()[column]);
                }
            }
        }
        @Override public double[] categoryProbabilities(double[] parameters) {
            if (parameters == null || parameters.length != predictors) {
                throw new IllegalArgumentException(
                    "wrong number of multinomial predictors");
            }
            Probabilities values = probabilities(parameters);
            double[] result = new double[categories];
            System.arraycopy(values.probabilities(), 0,
                result, 0, predictors);
            result[categories - 1] = Math.exp(-values.logDenominator());
            return result;
        }
        private Probabilities probabilities(double[] predictorsValues) {
            double maximum = 0.0;
            for (double value : predictorsValues) maximum = Math.max(maximum, value);
            double denominator = Math.exp(-maximum);
            double[] probabilities = new double[predictors];
            for (int category = 0; category < predictors; category++) {
                probabilities[category] =
                    Math.exp(predictorsValues[category] - maximum);
                denominator += probabilities[category];
            }
            for (int category = 0; category < predictors; category++) {
                probabilities[category] /= denominator;
            }
            return new Probabilities(
                probabilities, maximum + Math.log(denominator));
        }
        private void checkParameter(int parameter) {
            if (parameter < 0 || parameter >= predictors) {
                throw new IllegalArgumentException(
                    "unknown multinomial predictor");
            }
        }
        private record Probabilities(
                double[] probabilities, double logDenominator) { }
    }
}
