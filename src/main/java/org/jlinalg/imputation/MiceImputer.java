/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SplittableRandom;

/** Stochastic chained-equations imputation for mixed numeric data. */
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
                            column, options, random);
                        case BINARY -> binary(current, data, missing,
                            column, options.ridge(), random);
                        case CATEGORICAL -> categorical(current, data, missing,
                            column, options.ridge(), random);
                    }
                }
            }
            datasets.add(copy(current));
        }
        return new MiceResult(datasets,
            diagnostics(datasets, missing, types), options);
    }

    private static void continuous(double[][] current, double[][] original,
            boolean[][] missing, int target, MiceOptions options,
            SplittableRandom random) {
        double[][] design = design(current, target);
        int[] observed = rows(missing, target, false);
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
            Arrays.sort(order, (left, right) -> {
                int comparison = Double.compare(
                    Math.abs(predictedObserved[left] - predicted),
                    Math.abs(predictedObserved[right] - predicted));
                return comparison != 0 ? comparison : Integer.compare(left, right);
            });
            int donor = observed[order[random.nextInt(donors)]];
            current[row][target] = original[donor][target];
        }
    }

    private static void binary(double[][] current, double[][] original,
            boolean[][] missing, int target, double ridge,
            SplittableRandom random) {
        double[][] design = design(current, target);
        int[] observed = rows(missing, target, false);
        int[] absent = rows(missing, target, true);
        double[] response = new double[observed.length];
        double[][] x = new double[observed.length][];
        for (int index = 0; index < observed.length; index++) {
            response[index] = original[observed[index]][target];
            x[index] = design[observed[index]];
        }
        double[] coefficients = logistic(x, response, ridge);
        for (int row : absent) {
            double probability = logistic(dot(design[row], coefficients));
            current[row][target] = random.nextDouble() < probability ? 1.0 : 0.0;
        }
    }

    private static void categorical(double[][] current, double[][] original,
            boolean[][] missing, int target, double ridge,
            SplittableRandom random) {
        double[] categories = categories(original, missing, target);
        if (categories.length == 2) {
            double[][] design = design(current, target);
            int[] observed = rows(missing, target, false);
            int[] absent = rows(missing, target, true);
            double[] response = new double[observed.length];
            double[][] x = new double[observed.length][];
            for (int index = 0; index < observed.length; index++) {
                response[index] = original[observed[index]][target]
                    == categories[1] ? 1.0 : 0.0;
                x[index] = design[observed[index]];
            }
            double[] coefficients = logistic(x, response, ridge);
            for (int row : absent) {
                double probability = logistic(dot(design[row], coefficients));
                current[row][target] = random.nextDouble() < probability
                    ? categories[1] : categories[0];
            }
            return;
        }
        double[][] design = design(current, target);
        int[] observed = rows(missing, target, false);
        int[] absent = rows(missing, target, true);
        double[][] coefficients = new double[categories.length][];
        for (int category = 0; category < categories.length; category++) {
            double[] response = new double[observed.length];
            double[][] x = new double[observed.length][];
            for (int index = 0; index < observed.length; index++) {
                response[index] = original[observed[index]][target]
                    == categories[category] ? 1.0 : 0.0;
                x[index] = design[observed[index]];
            }
            coefficients[category] = logistic(x, response, ridge);
        }
        for (int row : absent) {
            double[] probability = new double[categories.length];
            double maximum = Double.NEGATIVE_INFINITY;
            for (int category = 0; category < categories.length; category++) {
                probability[category] = dot(design[row], coefficients[category]);
                maximum = Math.max(maximum, probability[category]);
            }
            double sum = 0.0;
            for (int category = 0; category < probability.length; category++) {
                probability[category] = Math.exp(probability[category] - maximum);
                sum += probability[category];
            }
            double draw = random.nextDouble() * sum;
            int selected = probability.length - 1;
            for (int category = 0; category < probability.length; category++) {
                draw -= probability[category];
                if (draw <= 0.0) { selected = category; break; }
            }
            current[row][target] = categories[selected];
        }
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

    private static double[][] design(double[][] data, int target) {
        int rows = data.length;
        int columns = data[0].length;
        double[][] result = new double[rows][columns];
        for (int row = 0; row < rows; row++) result[row][0] = 1.0;
        int output = 1;
        for (int column = 0; column < columns; column++) {
            if (column == target) continue;
            double mean = 0.0;
            for (double[] row : data) mean += row[column] / rows;
            double variance = 0.0;
            for (double[] row : data) variance += (row[column] - mean) * (row[column] - mean);
            double scale = Math.sqrt(variance / Math.max(1, rows - 1));
            if (!(scale > 0.0)) scale = 1.0;
            for (int row = 0; row < rows; row++)
                result[row][output] = (data[row][column] - mean) / scale;
            output++;
        }
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

    private static double[] logistic(double[][] design, double[] response,
            double ridge) {
        double[] coefficients = new double[design[0].length];
        for (int iteration = 0; iteration < 50; iteration++) {
            double[][] cross = new double[coefficients.length][coefficients.length];
            double[] target = new double[coefficients.length];
            for (int row = 0; row < design.length; row++) {
                double eta = dot(design[row], coefficients);
                double mean = logistic(eta);
                double weight = Math.max(1e-6, mean * (1.0 - mean));
                double working = eta + (response[row] - mean) / weight;
                for (int left = 0; left < coefficients.length; left++) {
                    target[left] = Math.fma(weight * design[row][left], working,
                        target[left]);
                    for (int right = 0; right < coefficients.length; right++)
                        cross[left][right] = Math.fma(weight * design[row][left],
                            design[row][right], cross[left][right]);
                }
            }
            for (int column = 1; column < coefficients.length; column++)
                cross[column][column] += ridge;
            cross[0][0] += ridge * 1e-3;
            double[] updated = solve(cross, target);
            double change = 0.0;
            for (int column = 0; column < coefficients.length; column++)
                change = Math.max(change, Math.abs(updated[column] - coefficients[column]));
            coefficients = updated;
            if (change < 1e-8) break;
        }
        return coefficients;
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

    private static double logistic(double value) {
        if (value >= 0.0) { double exponential = Math.exp(-value); return 1.0 / (1.0 + exponential); }
        double exponential = Math.exp(value); return exponential / (1.0 + exponential);
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
