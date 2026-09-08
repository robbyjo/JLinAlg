/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import jdistlib.Normal;

/** Internal validation, ranking, and scalar numerical helpers. */
final class StatisticsSupport {
    private StatisticsSupport() { }

    static double[] sample(double[] values, int minimum, String name) {
        if (values == null || values.length < minimum) {
            throw new IllegalArgumentException(
                name + " must contain at least " + minimum + " observations");
        }
        double[] copy = values.clone();
        for (double value : copy) {
            if (!Double.isFinite(value)) {
                throw new IllegalArgumentException(name + " must contain only finite values");
            }
        }
        return copy;
    }

    static void paired(double[] first, double[] second, int minimum) {
        sample(first, minimum, "first sample");
        sample(second, minimum, "second sample");
        if (first.length != second.length) {
            throw new IllegalArgumentException("paired samples must have equal lengths");
        }
    }

    static int[] groups(int[] groups, int length) {
        if (groups == null || groups.length != length) {
            throw new IllegalArgumentException(
                "group labels must have the same length as the observations");
        }
        return groups.clone();
    }

    static double mean(double[] values) {
        double scale = 0.0;
        for (double value : values) scale = Math.max(scale, Math.abs(value));
        if (scale == 0.0) return 0.0;
        double sum = 0.0, correction = 0.0;
        for (double value : values) {
            double adjusted = value / scale - correction;
            double next = sum + adjusted;
            correction = (next - sum) - adjusted;
            sum = next;
        }
        return (sum / values.length) * scale;
    }

    static double variance(double[] values) {
        double mean = mean(values);
        double sum = 0.0;
        for (double value : values) {
            double difference = value - mean;
            sum += difference * difference;
        }
        return sum / (values.length - 1.0);
    }

    static double[] differences(double[] first, double[] second) {
        paired(first, second, 2);
        double[] result = new double[first.length];
        for (int index = 0; index < result.length; index++) {
            result[index] = first[index] - second[index];
        }
        return result;
    }

    /** Applies a common finite unit scale in place to private sample copies. */
    static double rescale(double reference, double[]... samples) {
        double scale = Math.abs(reference);
        for (double[] values : samples) for (double value : values) scale = Math.max(scale, Math.abs(value));
        if (scale == 0.0) return 1.0;
        for (double[] values : samples) for (int i = 0; i < values.length; i++) values[i] /= scale;
        return scale;
    }

    static double[] ranks(double[] values) {
        Integer[] order = new Integer[values.length];
        for (int index = 0; index < order.length; index++) order[index] = index;
        Arrays.sort(order, (left, right) -> Double.compare(values[left], values[right]));
        double[] ranks = new double[values.length];
        int start = 0;
        while (start < order.length) {
            int end = start + 1;
            while (end < order.length
                    && values[order[start]] == values[order[end]]) {
                end++;
            }
            double rank = (start + 1.0 + end) / 2.0;
            for (int index = start; index < end; index++) ranks[order[index]] = rank;
            start = end;
        }
        return ranks;
    }

    static Map<Double, Integer> tieCounts(double[] values) {
        Map<Double, Integer> counts = new LinkedHashMap<>();
        for (double value : values) counts.merge(value == 0.0 ? 0.0 : value, 1, Integer::sum);
        counts.values().removeIf(count -> count == 1);
        return counts;
    }

    static double correlation(double[] first, double[] second) {
        double[] xValues = centeredUnitScale(first);
        double[] yValues = centeredUnitScale(second);
        double cross = 0.0;
        double firstSquare = 0.0;
        double secondSquare = 0.0;
        for (int index = 0; index < first.length; index++) {
            double x = xValues[index];
            double y = yValues[index];
            cross += x * y;
            firstSquare += x * x;
            secondSquare += y * y;
        }
        return Math.max(-1.0, Math.min(1.0,
            (cross / Math.sqrt(firstSquare)) / Math.sqrt(secondSquare)));
    }

    /** Center before scaling so a large shared offset does not lose the variation. */
    static double[] centeredUnitScale(double[] values) {
        double minimum = values[0], maximum = values[0];
        for (double value : values) {
            minimum = Math.min(minimum, value);
            maximum = Math.max(maximum, value);
        }
        double range = maximum - minimum;
        double center = Double.isFinite(range) ? minimum + range / 2.0
            : minimum / 2.0 + maximum / 2.0;
        double scale = Math.max(Math.abs(minimum - center), Math.abs(maximum - center));
        double[] result = new double[values.length];
        if (scale == 0.0) return result;
        double sum = 0.0, correction = 0.0;
        for (int i = 0; i < result.length; i++) {
            result[i] = (values[i] - center) / scale;
            double adjusted = result[i] - correction;
            double next = sum + adjusted;
            correction = (next - sum) - adjusted;
            sum = next;
        }
        double mean = sum / result.length;
        for (int i = 0; i < result.length; i++) result[i] -= mean;
        return result;
    }

    /** Counts concordant/discordant pairs in O(n log n), delaying updates within x ties. */
    static long[] kendallPairs(double[] first, double[] second) {
        Integer[] order = new Integer[first.length];
        for (int i = 0; i < order.length; i++) order[i] = i;
        Arrays.sort(order, (a, b) -> Double.compare(first[a], first[b]));
        double[] ranks = ranks(second);
        long[] tree = new long[2 * first.length + 2];
        long concordant = 0, discordant = 0, processed = 0;
        for (int start = 0; start < order.length;) {
            int end = start + 1;
            while (end < order.length && first[order[start]] == first[order[end]]) end++;
            for (int i = start; i < end; i++) {
                int rank = (int) (2 * ranks[order[i]]);
                long less = 0, lessOrEqual = 0;
                for (int j = rank - 1; j > 0; j -= j & -j) less += tree[j];
                for (int j = rank; j > 0; j -= j & -j) lessOrEqual += tree[j];
                concordant += less;
                discordant += processed - lessOrEqual;
            }
            for (int i = start; i < end; i++) {
                for (int j = (int) (2 * ranks[order[i]]); j < tree.length; j += j & -j) tree[j]++;
            }
            processed += end - start;
            start = end;
        }
        return new long[]{concordant, discordant};
    }

    static double normalPValue(double statistic, Alternative alternative) {
        return switch (alternative) {
            case LESS -> Normal.cumulative(statistic, 0.0, 1.0, true, false);
            case GREATER -> Normal.cumulative(statistic, 0.0, 1.0, false, false);
            case TWO_SIDED -> Math.min(1.0, 2.0 * Normal.cumulative(
                Math.abs(statistic), 0.0, 1.0, false, false));
        };
    }

    static double[] validateMatrix(double[][] matrix, int minimumRows,
            int minimumColumns, String name) {
        if (matrix == null || matrix.length < minimumRows || matrix[0] == null
                || matrix[0].length < minimumColumns) {
            throw new IllegalArgumentException(name + " has too few rows or columns");
        }
        int columns = matrix[0].length;
        double[] flat = new double[matrix.length * columns];
        for (int row = 0; row < matrix.length; row++) {
            if (matrix[row] == null || matrix[row].length != columns) {
                throw new IllegalArgumentException(name + " must be rectangular");
            }
            for (int column = 0; column < columns; column++) {
                double value = matrix[row][column];
                if (!Double.isFinite(value)) {
                    throw new IllegalArgumentException(name + " must contain finite values");
                }
                flat[row * columns + column] = value;
            }
        }
        return flat;
    }

    static long total(long[][] table, int minimumRows, int minimumColumns) {
        if (table == null || table.length < minimumRows || table[0] == null
                || table[0].length < minimumColumns) {
            throw new IllegalArgumentException("table has too few rows or columns");
        }
        int columns = table[0].length;
        long total = 0L;
        for (long[] row : table) {
            if (row == null || row.length != columns) {
                throw new IllegalArgumentException("table must be rectangular");
            }
            for (long count : row) {
                if (count < 0L || count > Long.MAX_VALUE - total) {
                    throw new IllegalArgumentException(
                        "table counts must be nonnegative and have a representable total");
                }
                total += count;
            }
        }
        return total;
    }

    static void probability(double probability, String name) {
        if (!(probability > 0.0 && probability < 1.0)
                || !Double.isFinite(probability)) {
            throw new IllegalArgumentException(name + " must be in (0, 1)");
        }
    }

    static double bisect(DoubleFunction function, double lower, double upper) {
        double lowerValue = function.value(lower);
        double upperValue = function.value(upper);
        if (Double.isNaN(lowerValue) || Double.isNaN(upperValue)
                || Math.signum(lowerValue) == Math.signum(upperValue)) {
            throw new IllegalArgumentException("unable to bracket numerical root");
        }
        double midpoint = Double.NaN;
        for (int iteration = 0; iteration < 200; iteration++) {
            midpoint = lower + (upper - lower) / 2.0;
            double value = function.value(midpoint);
            if (value == 0.0 || Math.abs(upper - lower)
                    <= 1e-12 * Math.max(1.0, Math.abs(midpoint))) {
                return midpoint;
            }
            if (Math.signum(value) == Math.signum(lowerValue)) {
                lower = midpoint;
                lowerValue = value;
            } else {
                upper = midpoint;
            }
        }
        return midpoint;
    }

    @FunctionalInterface
    interface DoubleFunction {
        double value(double value);
    }
}
