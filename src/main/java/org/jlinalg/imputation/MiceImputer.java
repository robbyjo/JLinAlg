/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SplittableRandom;

/** Bootstrap chained-equations imputation for mixed numeric data.
 * Each update resamples observed rows before fitting its conditional model;
 * PMM also samples donors from that bootstrap, preserving marginal uncertainty.
 * Binary and categorical responses use a joint multinomial logistic model.
 * Rubin inference remains conditional on suitable MAR and model assumptions. */
public final class MiceImputer {
    private MiceImputer() { }

    /** NaN denotes missingness; infinities are rejected. */
    public static MiceResult impute(double[][] data, VariableType[] types,
            MiceOptions options) {
        Dimensions dimensions = validate(data, types, options);
        boolean[][] missing = new boolean[dimensions.rows][dimensions.columns];
        for (int row = 0; row < dimensions.rows; row++)
            for (int column = 0; column < dimensions.columns; column++)
                missing[row][column] = Double.isNaN(data[row][column]);
        List<double[][]> datasets = new ArrayList<>();
        SplittableRandom master = new SplittableRandom(options.seed());
        for (int imputation = 0; imputation < options.imputations(); imputation++) {
            SplittableRandom random = master.split();
            double[][] current = initialize(data, missing, random);
            for (int iteration = 0; iteration < options.iterations(); iteration++) {
                for (int column = 0; column < dimensions.columns; column++) {
                    if (!anyMissing(missing, column)) continue;
                    switch (types[column]) {
                        case CONTINUOUS -> continuous(current, data, missing,
                            column, types, options, random);
                        case BINARY, CATEGORICAL -> categorical(current, data, missing,
                            column, types, options.ridge(), random);
                    }
                }
            }
            datasets.add(copy(current));
        }
        return new MiceResult(datasets,
            diagnostics(datasets, missing, types), options);
    }

    private static void continuous(double[][] current, double[][] original,
            boolean[][] missing, int target, VariableType[] types, MiceOptions options,
            SplittableRandom random) {
        double[][] design = design(current, target, types);
        int[] observed = bootstrap(rows(missing, target, false), random);
        int[] absent = rows(missing, target, true);
        double[] response = new double[observed.length];
        double[][] x = new double[observed.length][];
        for (int index = 0; index < observed.length; index++) {
            response[index] = original[observed[index]][target];
            x[index] = design[observed[index]];
        }
        double[] coefficients = ridge(x, response, options.ridge());
        double[] predictedObserved = new double[observed.length];
        for (int index = 0; index < observed.length; index++)
            predictedObserved[index] = dot(x[index], coefficients);
        int donors = Math.min(options.predictiveMeanDonors(), observed.length);
        for (int row : absent) {
            double predicted = dot(design[row], coefficients);
            Integer[] order = new Integer[observed.length];
            for (int index = 0; index < order.length; index++) order[index] = index;
            // A stable distance sort after shuffling randomizes all ties,
            // including the boundary of the k-nearest donor pool.
            for (int i = order.length - 1; i > 0; i--) {
                int j = random.nextInt(i + 1);
                Integer swap = order[i]; order[i] = order[j]; order[j] = swap;
            }
            Arrays.sort(order, (left, right) -> Double.compare(
                Math.abs(predictedObserved[left] - predicted),
                Math.abs(predictedObserved[right] - predicted)));
            int donor = observed[order[random.nextInt(donors)]];
            current[row][target] = original[donor][target];
        }
    }

    private static void categorical(double[][] current, double[][] original,
            boolean[][] missing, int target, VariableType[] types, double ridge,
            SplittableRandom random) {
        double[] categories = categories(original, missing, target);
        double[][] design = design(current, target, types);
        int[] observed = bootstrap(rows(missing, target, false), random);
        double[][] x = new double[observed.length][];
        int[] response = new int[observed.length];
        for (int i = 0; i < observed.length; i++) {
            x[i] = design[observed[i]];
            response[i] = Arrays.binarySearch(categories, original[observed[i]][target]);
        }
        double[] coefficients = multinomial(x, response, categories.length, ridge);
        for (int row : rows(missing, target, true)) {
            double[] probability = probabilities(design[row], coefficients, categories.length);
            double draw = random.nextDouble();
            int selected = probability.length - 1;
            for (int category = 0; category < probability.length; category++) {
                draw -= probability[category];
                if (draw <= 0) { selected = category; break; }
            }
            current[row][target] = categories[selected];
        }
    }

    private static int[] bootstrap(int[] observed, SplittableRandom random) {
        int[] result = new int[observed.length];
        for (int i = 0; i < result.length; i++) result[i] = observed[random.nextInt(observed.length)];
        return result;
    }

    private static double[][] initialize(double[][] data, boolean[][] missing,
            SplittableRandom random) {
        double[][] result = copy(data);
        for (int column = 0; column < data[0].length; column++) {
            int[] observed = rows(missing, column, false);
            for (int row = 0; row < data.length; row++)
                if (missing[row][column])
                    result[row][column] = data[observed[random.nextInt(observed.length)]][column];
        }
        return result;
    }

    private static double[][] design(double[][] data, int target, VariableType[] types) {
        int rows = data.length;
        List<double[]> columns = new ArrayList<>();
        double[] intercept = new double[rows]; Arrays.fill(intercept, 1); columns.add(intercept);
        for (int column = 0; column < data[0].length; column++) {
            if (column == target) continue;
            if (types[column] == VariableType.CATEGORICAL) {
                final int index = column;
                double[] levels = Arrays.stream(data).mapToDouble(r -> r[index]).distinct().sorted().toArray();
                for (int k = 1; k < levels.length; k++) {
                    double[] dummy = new double[rows];
                    for (int row = 0; row < rows; row++) dummy[row] = data[row][column] == levels[k] ? 1 : 0;
                    columns.add(dummy);
                }
            } else {
                boolean constant = true;
                for (double[] row : data) constant &= row[column] == data[0][column];
                if (constant) continue;
                double mean = 0;
                for (double[] row : data) mean += row[column] / rows;
                double variance = 0;
                for (double[] row : data) variance += Math.pow(row[column] - mean, 2);
                if (variance == 0) continue;
                double scale = Math.sqrt(variance / (rows - 1));
                double[] values = new double[rows];
                for (int row = 0; row < rows; row++) values[row] = (data[row][column] - mean) / scale;
                columns.add(values);
            }
        }
        double[][] result = new double[rows][columns.size()];
        for (int row = 0; row < rows; row++) for (int c = 0; c < columns.size(); c++)
            result[row][c] = columns.get(c)[row];
        return result;
    }

    private static double[] ridge(double[][] design, double[] response,
            double ridge) {
        int columns = design[0].length;
        double[][] cross = new double[columns][columns];
        double[] target = new double[columns];
        for (int row = 0; row < design.length; row++) {
            for (int left = 0; left < columns; left++) {
                target[left] = Math.fma(design[row][left], response[row], target[left]);
                for (int right = 0; right < columns; right++)
                    cross[left][right] = Math.fma(design[row][left],
                        design[row][right], cross[left][right]);
            }
        }
        for (int column = 1; column < columns; column++) cross[column][column] += ridge;
        cross[0][0] += ridge * 1e-3;
        return solve(cross, target);
    }

    // Baseline-category multinomial likelihood, with one joint information
    // matrix. A weak ridge keeps bootstrap samples with separation finite.
    private static double[] multinomial(double[][] x, int[] y, int categories, double ridge) {
        int p = x[0].length, size = p * (categories - 1);
        double[] beta = new double[size];
        for (int iteration = 0; iteration < 200; iteration++) {
            double[][] information = new double[size][size];
            double[] score = new double[size];
            for (int row = 0; row < x.length; row++) {
                double[] probability = probabilities(x[row], beta, categories);
                for (int a = 0; a < categories - 1; a++) for (int j = 0; j < p; j++) {
                    int left = a * p + j;
                    score[left] += x[row][j] * ((y[row] == a ? 1 : 0) - probability[a]);
                    for (int b = 0; b < categories - 1; b++) for (int k = 0; k < p; k++)
                        information[left][b * p + k] += x[row][j] * x[row][k]
                            * probability[a] * ((a == b ? 1 : 0) - probability[b]);
                }
            }
            double maxScore = 0;
            for (int j = 0; j < size; j++) {
                double penalty = ridge * (j % p == 0 ? 1e-3 : 1);
                information[j][j] += penalty; score[j] -= penalty * beta[j];
                maxScore = Math.max(maxScore, Math.abs(score[j]));
            }
            if (maxScore < 1e-8) return beta;
            double[] step = solve(information, score);
            double old = objective(x, y, beta, categories, ridge);
            double fraction = 1; boolean accepted = false;
            for (int search = 0; search < 40; search++, fraction *= .5) {
                double[] candidate = beta.clone();
                for (int j = 0; j < size; j++) candidate[j] += fraction * step[j];
                if (objective(x, y, candidate, categories, ridge) <= old + 1e-12) {
                    beta = candidate; accepted = true; break;
                }
            }
            if (!accepted) throw new IllegalArgumentException("imputation multinomial line search failed");
        }
        throw new IllegalArgumentException("imputation multinomial model did not converge");
    }

    private static double[] probabilities(double[] x, double[] beta, int categories) {
        double[] logits = new double[categories]; double maximum = 0;
        for (int k = 0; k < categories - 1; k++) {
            for (int j = 0; j < x.length; j++) logits[k] += x[j] * beta[k * x.length + j];
            maximum = Math.max(maximum, logits[k]);
        }
        double sum = 0;
        for (int k = 0; k < categories; k++) { logits[k] = Math.exp(logits[k] - maximum); sum += logits[k]; }
        for (int k = 0; k < categories; k++) logits[k] /= sum;
        return logits;
    }

    private static double objective(double[][] x, int[] y, double[] beta, int categories, double ridge) {
        double value = 0;
        for (int i = 0; i < x.length; i++) {
            double[] logits = new double[categories]; double max = 0;
            for (int k = 0; k < categories - 1; k++) {
                for (int j = 0; j < x[i].length; j++) logits[k] += x[i][j] * beta[k * x[i].length + j];
                max = Math.max(max, logits[k]);
            }
            double sum = 0;
            for (double logit : logits) sum += Math.exp(logit - max);
            value += max - logits[y[i]] + Math.log(sum);
        }
        for (int j = 0; j < beta.length; j++)
            value += .5 * ridge * (j % x[0].length == 0 ? 1e-3 : 1) * beta[j] * beta[j];
        return value;
    }

    private static double[] solve(double[][] matrix, double[] right) {
        int size = right.length;
        double[][] augmented = new double[size][size + 1];
        for (int row = 0; row < size; row++) {
            System.arraycopy(matrix[row], 0, augmented[row], 0, size);
            augmented[row][size] = right[row];
        }
        for (int pivot = 0; pivot < size; pivot++) {
            int selected = pivot;
            for (int row = pivot + 1; row < size; row++)
                if (Math.abs(augmented[row][pivot]) > Math.abs(augmented[selected][pivot]))
                    selected = row;
            double[] swap = augmented[pivot]; augmented[pivot] = augmented[selected];
            augmented[selected] = swap;
            double diagonal = augmented[pivot][pivot];
            if (Math.abs(diagonal) < 1e-12)
                throw new IllegalArgumentException("imputation regression is numerically singular");
            for (int column = pivot; column <= size; column++)
                augmented[pivot][column] /= diagonal;
            for (int row = 0; row < size; row++) {
                if (row == pivot) continue;
                double factor = augmented[row][pivot];
                for (int column = pivot; column <= size; column++)
                    augmented[row][column] -= factor * augmented[pivot][column];
            }
        }
        double[] result = new double[size];
        for (int row = 0; row < size; row++) result[row] = augmented[row][size];
        return result;
    }

    private static List<ImputationDiagnostic> diagnostics(
            List<double[][]> datasets, boolean[][] missing, VariableType[] types) {
        List<ImputationDiagnostic> result = new ArrayList<>();
        for (int column = 0; column < types.length; column++) {
            int missingCount = rows(missing, column, true).length;
            double[] means = new double[datasets.size()];
            double[] variances = new double[datasets.size()];
            if (missingCount > 0) {
                int[] absent = rows(missing, column, true);
                for (int chain = 0; chain < datasets.size(); chain++) {
                    for (int row : absent) means[chain] += datasets.get(chain)[row][column]
                        / missingCount;
                    if (missingCount > 1) for (int row : absent)
                        variances[chain] += Math.pow(datasets.get(chain)[row][column]
                            - means[chain], 2.0) / (missingCount - 1.0);
                }
            }
            double meanOfMeans = Arrays.stream(means).average().orElse(0.0);
            double between = sampleVariance(means, meanOfMeans);
            double within = Arrays.stream(variances).average().orElse(0.0);
            result.add(new ImputationDiagnostic(column, types[column],
                missing.length - missingCount, missingCount, means, variances,
                between, within));
        }
        return List.copyOf(result);
    }

    private static double[] categories(double[][] data, boolean[][] missing,
            int column) {
        return java.util.stream.IntStream.range(0, data.length)
            .filter(row -> !missing[row][column]).mapToDouble(row -> data[row][column])
            .distinct().sorted().toArray();
    }

    private static int[] rows(boolean[][] missing, int column, boolean value) {
        return java.util.stream.IntStream.range(0, missing.length)
            .filter(row -> missing[row][column] == value).toArray();
    }

    private static boolean anyMissing(boolean[][] missing, int column) {
        for (boolean[] row : missing) if (row[column]) return true;
        return false;
    }

    private static double dot(double[] left, double[] right) {
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

    private static double[][] copy(double[][] source) {
        double[][] result = new double[source.length][];
        for (int row = 0; row < source.length; row++) result[row] = source[row].clone();
        return result;
    }

    private static Dimensions validate(double[][] data, VariableType[] types,
            MiceOptions options) {
        if (data == null || data.length < 3 || data[0] == null || data[0].length < 2)
            throw new IllegalArgumentException("MICE needs at least three rows and two variables");
        if (types == null || types.length != data[0].length || options == null)
            throw new IllegalArgumentException("one variable type and options are required");
        int columns = data[0].length;
        for (double[] row : data)
            if (row == null || row.length != columns)
                throw new IllegalArgumentException("imputation data must be rectangular");
        for (int column = 0; column < columns; column++) {
            if (types[column] == null) throw new IllegalArgumentException("variable types cannot be null");
            int observed = 0;
            boolean zero = false;
            boolean one = false;
            for (double[] row : data) {
                double value = row[column];
                if (Double.isInfinite(value))
                    throw new IllegalArgumentException("imputation data cannot contain infinities");
                if (Double.isNaN(value)) continue;
                observed++;
                if (types[column] != VariableType.CONTINUOUS
                        && value != Math.rint(value))
                    throw new IllegalArgumentException(
                        "categorical and binary values must be integer-coded");
                zero |= value == 0.0;
                one |= value == 1.0;
                if (types[column] == VariableType.BINARY && value != 0.0 && value != 1.0)
                    throw new IllegalArgumentException("binary variables must use 0 and 1");
            }
            if (observed < 2) throw new IllegalArgumentException(
                "every variable needs at least two observed values");
            if (types[column] == VariableType.BINARY && (!zero || !one))
                throw new IllegalArgumentException(
                    "binary variables need observed values in both classes");
            if (types[column] == VariableType.CATEGORICAL
                    && categories(data, missingFor(data), column).length < 2)
                throw new IllegalArgumentException(
                    "categorical variables need at least two observed categories");
        }
        return new Dimensions(data.length, columns);
    }

    private static boolean[][] missingFor(double[][] data) {
        boolean[][] result = new boolean[data.length][data[0].length];
        for (int row = 0; row < data.length; row++)
            for (int column = 0; column < data[row].length; column++)
                result[row][column] = Double.isNaN(data[row][column]);
        return result;
    }

    private record Dimensions(int rows, int columns) { }
}
