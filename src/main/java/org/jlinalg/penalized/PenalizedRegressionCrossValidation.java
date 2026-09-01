/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import java.util.Random;
import org.jlinalg.internal.MatrixOps;

/** K-fold cross-validation utilities for Gaussian penalized regression. */
public final class PenalizedRegressionCrossValidation {
    private PenalizedRegressionCrossValidation() { }

    /** Fits user-supplied descending lambdas with deterministic shuffled folds. */
    public static PenalizedCrossValidationResult fit(
            double[] response,
            double[][] predictors,
            double[] lambdas,
            int folds,
            long randomSeed,
            ElasticNetOptions options) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] design = MatrixOps.rowMajor(predictors, response.length);
        int columns = predictors[0].length;
        MatrixOps.validateModelData(
            response, design, response.length, columns);
        validateFolds(response.length, folds);
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }

        PenalizedRegressionPath fullPath = PenalizedRegression.path(
            response, design, response.length, columns, lambdas, options);
        double[] validatedLambdas = fullPath.lambdas();
        int[] foldByObservation = folds(response.length, folds, randomSeed);
        double[][] foldErrors = new double[folds][validatedLambdas.length];
        double[] originalWeights = options.observationWeights();

        for (int fold = 0; fold < folds; fold++) {
            int trainingRows = 0;
            for (int value : foldByObservation) {
                if (value != fold) {
                    trainingRows++;
                }
            }
            double[] trainingResponse = new double[trainingRows];
            double[] trainingPredictors = new double[trainingRows * columns];
            double[] trainingWeights = originalWeights == null
                ? null : new double[trainingRows];
            int training = 0;
            for (int row = 0; row < response.length; row++) {
                if (foldByObservation[row] == fold) {
                    continue;
                }
                trainingResponse[training] = response[row];
                System.arraycopy(design, row * columns,
                    trainingPredictors, training * columns, columns);
                if (trainingWeights != null) {
                    trainingWeights[training] = originalWeights[row];
                }
                training++;
            }
            ElasticNetOptions trainingOptions = copyOptions(
                options, trainingWeights);
            PenalizedRegressionPath trainingPath = PenalizedRegression.path(
                trainingResponse, trainingPredictors,
                trainingRows, columns, validatedLambdas, trainingOptions);

            for (int lambdaIndex = 0;
                    lambdaIndex < validatedLambdas.length; lambdaIndex++) {
                PenalizedRegressionResult fit = trainingPath.fit(lambdaIndex);
                double[] coefficients = fit.coefficients();
                double sumSquared = 0.0;
                double sumWeight = 0.0;
                for (int row = 0; row < response.length; row++) {
                    if (foldByObservation[row] != fold) {
                        continue;
                    }
                    double prediction = fit.intercept();
                    for (int column = 0; column < columns; column++) {
                        prediction += design[row * columns + column]
                            * coefficients[column];
                    }
                    double weight = originalWeights == null
                        ? 1.0 : originalWeights[row];
                    double error = response[row] - prediction;
                    sumSquared += weight * error * error;
                    sumWeight += weight;
                }
                foldErrors[fold][lambdaIndex] = sumSquared / sumWeight;
            }
        }

        double[] means = new double[validatedLambdas.length];
        double[] standardErrors = new double[validatedLambdas.length];
        for (int lambdaIndex = 0;
                lambdaIndex < validatedLambdas.length; lambdaIndex++) {
            for (int fold = 0; fold < folds; fold++) {
                means[lambdaIndex] += foldErrors[fold][lambdaIndex];
            }
            means[lambdaIndex] /= folds;
            double sumSquares = 0.0;
            for (int fold = 0; fold < folds; fold++) {
                double centered = foldErrors[fold][lambdaIndex]
                    - means[lambdaIndex];
                sumSquares += centered * centered;
            }
            standardErrors[lambdaIndex] = Math.sqrt(
                sumSquares / (folds - 1) / folds);
        }
        int minimum = minimumIndex(means);
        double threshold = means[minimum] + standardErrors[minimum];
        int oneStandardError = minimum;
        for (int index = 0; index <= minimum; index++) {
            if (means[index] <= threshold) {
                oneStandardError = index;
                break;
            }
        }
        return new PenalizedCrossValidationResult(
            fullPath, means, standardErrors, minimum, oneStandardError,
            folds, randomSeed);
    }

    /** Generates the full-data automatic path, then cross-validates that path. */
    public static PenalizedCrossValidationResult automatic(
            double[] response,
            double[][] predictors,
            int lambdaCount,
            double minimumRatio,
            int folds,
            long randomSeed,
            ElasticNetOptions options) {
        PenalizedRegressionPath generated = PenalizedRegression.automaticPath(
            response, predictors, lambdaCount, minimumRatio, options);
        return fit(response, predictors, generated.lambdas(),
            folds, randomSeed, options);
    }

    private static ElasticNetOptions copyOptions(
            ElasticNetOptions source, double[] weights) {
        ElasticNetOptions.Builder builder = ElasticNetOptions.builder()
            .alpha(source.alpha())
            .fitIntercept(source.fitIntercept())
            .standardize(source.standardize())
            .maximumIterations(source.maximumIterations())
            .relativeTolerance(source.relativeTolerance());
        if (weights != null) {
            builder.observationWeights(weights);
        }
        if (source.penaltyFactors() != null) {
            builder.penaltyFactors(source.penaltyFactors());
        }
        return builder.build();
    }

    private static int[] folds(int rows, int foldCount, long seed) {
        int[] order = new int[rows];
        for (int index = 0; index < rows; index++) {
            order[index] = index;
        }
        Random random = new Random(seed);
        for (int index = rows - 1; index > 0; index--) {
            int other = random.nextInt(index + 1);
            int temporary = order[index];
            order[index] = order[other];
            order[other] = temporary;
        }
        int[] result = new int[rows];
        for (int position = 0; position < rows; position++) {
            result[order[position]] = position % foldCount;
        }
        return result;
    }

    private static void validateFolds(int rows, int folds) {
        if (folds < 2 || folds > rows) {
            throw new IllegalArgumentException(
                "folds must lie between two and the observation count");
        }
    }

    private static int minimumIndex(double[] values) {
        int result = 0;
        for (int index = 1; index < values.length; index++) {
            if (values[index] < values[result]) {
                result = index;
            }
        }
        return result;
    }
}
