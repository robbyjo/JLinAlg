/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.Arrays;
import java.util.Comparator;
import jdistlib.Normal;
import jdistlib.T;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Working-independence marginal multinomial-logit GEE. */
public final class NominalGee {
    private NominalGee() { }

    /**
     * Fits category-specific marginal logits against reference category zero.
     * The cluster sandwich accounts for arbitrary within-cluster dependence.
     */
    public static NominalGeeResult fit(
            int[] response,
            double[][] design,
            int[] cluster,
            int categories,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || design == null || cluster == null
                || options == null || backendPolicy == null
                || response.length == 0 || response.length != design.length
                || response.length != cluster.length || design[0] == null
                || design[0].length == 0 || categories < 3) {
            throw new IllegalArgumentException("nominal GEE inputs are invalid");
        }
        if (options.correlation() != GeeCorrelation.INDEPENDENCE
                || options.method() != GeeMethod.ORDINARY) {
            throw new IllegalArgumentException(
                "nominal GEE currently uses ordinary working independence");
        }
        if (options.covariance() != GeeCovariance.NAIVE
                && options.covariance() != GeeCovariance.ROBUST
                && options.covariance() != GeeCovariance.DF_ADJUSTED) {
            throw new IllegalArgumentException(
                "nominal GEE supports naive, robust, and DF-adjusted covariance");
        }
        int n = response.length;
        int predictors = design[0].length;
        Integer[] order = new Integer[n];
        for (int row = 0; row < n; row++) {
            if (design[row] == null || design[row].length != predictors
                    || response[row] < 0 || response[row] >= categories) {
                throw new IllegalArgumentException(
                    "nominal responses or design rows are invalid");
            }
            order[row] = row;
        }
        Arrays.sort(order, Comparator.comparingInt(row -> cluster[row]));
        int[] sortedResponse = new int[n];
        int[] sortedCluster = new int[n];
        double[] sortedDesign = new double[n * predictors];
        int[] original = new int[n];
        for (int destination = 0; destination < n; destination++) {
            int source = order[destination];
            sortedResponse[destination] = response[source];
            sortedCluster[destination] = cluster[source];
            original[destination] = source;
            System.arraycopy(design[source], 0, sortedDesign,
                destination * predictors, predictors);
        }
        int[] starts = starts(sortedCluster);
        int parameters = predictors * (categories - 1);
        if (starts.length - 1 <= 1 || n <= parameters) {
            throw new IllegalArgumentException(
                "nominal GEE requires at least two clusters and n > parameters");
        }
        if ((options.covariance() == GeeCovariance.DF_ADJUSTED
                || options.inference() == GeeInference.CLUSTER_T)
                && starts.length - 1 <= parameters) {
            throw new IllegalArgumentException(
                "small-sample nominal inference requires clusters > parameters");
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] beta = new double[parameters];
            double[] probabilities = probabilities(sortedDesign, n, predictors,
                categories, beta);
            boolean converged = false;
            int iterations = 0;
            for (int iteration = 1;
                    iteration <= options.maximumIterations(); iteration++) {
                iterations = iteration;
                Components values = components(sortedResponse, sortedDesign,
                    starts, n, predictors, categories, probabilities);
                double[] step = solve(values.bread(), parameters,
                    values.score(), backend);
                double baseline = logLikelihood(sortedResponse, probabilities,
                    categories);
                double multiplier = 1.0;
                double[] candidate = null;
                double[] candidateProbabilities = null;
                while (multiplier >= 0x1.0p-20) {
                    double[] trial = beta.clone();
                    for (int index = 0; index < parameters; index++) {
                        trial[index] += multiplier * step[index];
                    }
                    double[] trialProbabilities = probabilities(sortedDesign,
                        n, predictors, categories, trial);
                    if (logLikelihood(sortedResponse, trialProbabilities, categories)
                            >= baseline - 1e-10 * (1.0 + Math.abs(baseline))) {
                        candidate = trial;
                        candidateProbabilities = trialProbabilities;
                        break;
                    }
                    multiplier *= 0.5;
                }
                if (candidate == null) break;
                double change = relativeChange(beta, candidate);
                beta = candidate;
                probabilities = candidateProbabilities;
                if (change <= options.relativeTolerance()) {
                    converged = true;
                    break;
                }
            }
            Components finalValues = components(sortedResponse, sortedDesign,
                starts, n, predictors, categories, probabilities);
            double[] inverseBread = inverse(finalValues.bread(), parameters, backend);
            double[] covariance = sandwich(inverseBread, finalValues.meat(), parameters);
            if (options.covariance() == GeeCovariance.NAIVE) {
                covariance = inverseBread;
            } else if (options.covariance() == GeeCovariance.DF_ADJUSTED) {
                double factor = (double) (starts.length - 1)
                    / (starts.length - 1 - parameters);
                for (int index = 0; index < covariance.length; index++) {
                    covariance[index] *= factor;
                }
            }
            double degrees = options.inference() == GeeInference.CLUSTER_T
                ? starts.length - 1 - parameters : Double.POSITIVE_INFINITY;
            if (options.inference() == GeeInference.CLUSTER_T && !(degrees > 0.0)) {
                throw new IllegalArgumentException(
                    "cluster-t inference requires clusters > parameters");
            }
            double[] se = new double[parameters];
            double[] statistic = new double[parameters];
            double[] pValue = new double[parameters];
            for (int index = 0; index < parameters; index++) {
                se[index] = Math.sqrt(Math.max(0.0,
                    covariance[index * parameters + index]));
                statistic[index] = beta[index] / se[index];
                pValue[index] = options.inference() == GeeInference.CLUSTER_T
                    ? 2.0 * T.cumulative(Math.abs(statistic[index]),
                        degrees, false, false)
                    : 2.0 * Normal.cumulative(Math.abs(statistic[index]),
                        0.0, 1.0, false, false);
            }
            double[] outputProbabilities = new double[n * categories];
            for (int sorted = 0; sorted < n; sorted++) {
                System.arraycopy(probabilities, sorted * categories,
                    outputProbabilities, original[sorted] * categories, categories);
            }
            return new NominalGeeResult(categories, predictors, beta, covariance,
                se, statistic, pValue, outputProbabilities, starts.length - 1,
                iterations, converged, degrees, context.provenance());
        }
    }

    private static double[] probabilities(
            double[] design, int rows, int predictors,
            int categories, double[] beta) {
        double[] result = new double[rows * categories];
        for (int row = 0; row < rows; row++) {
            double maximum = 0.0;
            for (int category = 1; category < categories; category++) {
                double eta = 0.0;
                int coefficient = (category - 1) * predictors;
                for (int column = 0; column < predictors; column++) {
                    eta += design[row * predictors + column]
                        * beta[coefficient + column];
                }
                result[row * categories + category] = eta;
                maximum = Math.max(maximum, eta);
            }
            double denominator = Math.exp(-maximum);
            for (int category = 1; category < categories; category++) {
                double value = Math.exp(result[row * categories + category] - maximum);
                result[row * categories + category] = value;
                denominator += value;
            }
            result[row * categories] = Math.exp(-maximum) / denominator;
            for (int category = 1; category < categories; category++) {
                result[row * categories + category] /= denominator;
            }
        }
        return result;
    }

    private static Components components(
            int[] response, double[] design, int[] starts,
            int rows, int predictors, int categories, double[] probabilities) {
        int parameters = predictors * (categories - 1);
        double[] bread = new double[parameters * parameters];
        double[] meat = new double[parameters * parameters];
        double[] score = new double[parameters];
        for (int cluster = 0; cluster < starts.length - 1; cluster++) {
            double[] clusterScore = new double[parameters];
            for (int row = starts[cluster]; row < starts[cluster + 1]; row++) {
                for (int firstCategory = 1;
                        firstCategory < categories; firstCategory++) {
                    double residual = (response[row] == firstCategory ? 1.0 : 0.0)
                        - probabilities[row * categories + firstCategory];
                    for (int firstColumn = 0;
                            firstColumn < predictors; firstColumn++) {
                        int first = (firstCategory - 1) * predictors + firstColumn;
                        clusterScore[first] += design[row * predictors + firstColumn]
                            * residual;
                        for (int secondCategory = 1;
                                secondCategory < categories; secondCategory++) {
                            double weight = probabilities[row * categories + firstCategory]
                                * ((firstCategory == secondCategory ? 1.0 : 0.0)
                                    - probabilities[row * categories + secondCategory]);
                            for (int secondColumn = 0;
                                    secondColumn < predictors; secondColumn++) {
                                int second = (secondCategory - 1) * predictors
                                    + secondColumn;
                                bread[first * parameters + second] += weight
                                    * design[row * predictors + firstColumn]
                                    * design[row * predictors + secondColumn];
                            }
                        }
                    }
                }
            }
            for (int first = 0; first < parameters; first++) {
                score[first] += clusterScore[first];
                for (int second = 0; second < parameters; second++) {
                    meat[first * parameters + second] += clusterScore[first]
                        * clusterScore[second];
                }
            }
        }
        return new Components(bread, meat, score);
    }

    private static double logLikelihood(
            int[] response, double[] probabilities, int categories) {
        double result = 0.0;
        for (int row = 0; row < response.length; row++) {
            result += Math.log(Math.max(1e-300,
                probabilities[row * categories + response[row]]));
        }
        return result;
    }

    private static int[] starts(int[] cluster) {
        int count = 1;
        for (int row = 1; row < cluster.length; row++) {
            if (cluster[row] != cluster[row - 1]) count++;
        }
        int[] result = new int[count + 1];
        int destination = 1;
        for (int row = 1; row < cluster.length; row++) {
            if (cluster[row] != cluster[row - 1]) result[destination++] = row;
        }
        result[count] = cluster.length;
        return result;
    }

    private static double[] solve(
            double[] matrix, int dimension,
            double[] right, ComputeBackend backend) {
        try {
            return backend.dpotrf(matrix, dimension).solve(right);
        } catch (IllegalArgumentException exception) {
            return backend.dsytrf(matrix, dimension).solve(right);
        }
    }

    private static double[] inverse(
            double[] matrix, int dimension, ComputeBackend backend) {
        double[] result = new double[dimension * dimension];
        CholeskyFactor factor = backend.dpotrf(matrix, dimension);
        for (int column = 0; column < dimension; column++) {
            double[] unit = new double[dimension];
            unit[column] = 1.0;
            double[] solved = factor.solve(unit);
            for (int row = 0; row < dimension; row++) {
                result[row * dimension + column] = solved[row];
            }
        }
        return result;
    }

    private static double[] sandwich(
            double[] inverseBread, double[] meat, int dimension) {
        double[] first = new double[dimension * dimension];
        double[] result = new double[dimension * dimension];
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                for (int inner = 0; inner < dimension; inner++) {
                    first[row * dimension + column] += inverseBread[row * dimension + inner]
                        * meat[inner * dimension + column];
                }
            }
        }
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                for (int inner = 0; inner < dimension; inner++) {
                    result[row * dimension + column] += first[row * dimension + inner]
                        * inverseBread[column * dimension + inner];
                }
            }
        }
        return result;
    }

    private static double relativeChange(double[] current, double[] candidate) {
        double maximum = 0.0;
        for (int index = 0; index < current.length; index++) {
            maximum = Math.max(maximum, Math.abs(candidate[index] - current[index])
                / (1.0 + Math.abs(current[index])));
        }
        return maximum;
    }

    private record Components(double[] bread, double[] meat, double[] score) { }
}
