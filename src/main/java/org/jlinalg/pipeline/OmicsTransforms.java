/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import jdistlib.Normal;

/** Common deterministic row-wise preprocessing for TWAS, EWAS, and PWAS. */
public final class OmicsTransforms {
    private OmicsTransforms() { }

    public static OmicsTransform identity() {
        return values -> required(values).clone();
    }

    /** Clamps finite values to empirical linear-interpolation quantiles. */
    public static OmicsTransform winsorize(
            double lowerQuantile, double upperQuantile) {
        if (!Double.isFinite(lowerQuantile) || !Double.isFinite(upperQuantile)
                || lowerQuantile < 0 || upperQuantile > 1
                || lowerQuantile > upperQuantile)
            throw new IllegalArgumentException(
                "winsorization quantiles must satisfy 0 <= lower <= upper <= 1");
        return values -> {
            double[] result = required(values).clone();
            double[] finite = Arrays.stream(result)
                .filter(Double::isFinite).sorted().toArray();
            if (finite.length == 0) return result;
            double lower = quantile(finite, lowerQuantile);
            double upper = quantile(finite, upperQuantile);
            for (int index = 0; index < result.length; index++) {
                if (Double.isFinite(result[index]))
                    result[index] = Math.max(lower, Math.min(upper, result[index]));
            }
            return result;
        };
    }

    /** Applies log(1+x) to finite values; all must be greater than -1. */
    public static OmicsTransform log1p() {
        return values -> {
            double[] result = required(values).clone();
            for (int index = 0; index < result.length; index++) {
                if (!Double.isFinite(result[index])) continue;
                if (result[index] <= -1)
                    throw new IllegalArgumentException(
                        "log1p transform requires finite values greater than -1");
                result[index] = Math.log1p(result[index]);
            }
            return result;
        };
    }

    /** Applies log(x+shift); shifted finite values must be positive. */
    public static OmicsTransform shiftedLog(double shift) {
        if (!Double.isFinite(shift))
            throw new IllegalArgumentException("log shift must be finite");
        return values -> {
            double[] result = required(values).clone();
            for (int index = 0; index < result.length; index++) {
                if (!Double.isFinite(result[index])) continue;
                double shifted = result[index] + shift;
                if (!(shifted > 0))
                    throw new IllegalArgumentException(
                        "shifted log requires positive shifted values");
                result[index] = Math.log(shifted);
            }
            return result;
        };
    }

    /** Centers and scales finite values using the sample standard deviation. */
    public static OmicsTransform zScore() {
        return values -> {
            double[] result = required(values).clone();
            Summary summary = summary(result);
            if (summary.count() < 2 || !(summary.sampleVariance() > 0))
                throw new IllegalArgumentException(
                    "z-score transform requires at least two varying finite values");
            double standardDeviation = Math.sqrt(summary.sampleVariance());
            for (int index = 0; index < result.length; index++)
                if (Double.isFinite(result[index]))
                    result[index] = (result[index] - summary.mean())
                        / standardDeviation;
            return result;
        };
    }

    /** Converts methylation beta values to base-2 M values after clipping. */
    public static OmicsTransform mValue(double epsilon) {
        if (!(epsilon > 0.0 && epsilon < 0.5) || !Double.isFinite(epsilon))
            throw new IllegalArgumentException(
                "M-value epsilon must be finite and in (0, 0.5)");
        return values -> {
            double[] result = required(values).clone();
            for (int index = 0; index < result.length; index++) {
                if (!Double.isFinite(result[index])) continue;
                if (result[index] < 0.0 || result[index] > 1.0)
                    throw new IllegalArgumentException(
                        "M-value transform requires beta values in [0,1]");
                double beta = Math.max(epsilon,
                    Math.min(1.0 - epsilon, result[index]));
                result[index] = Math.log(beta / (1.0 - beta)) / Math.log(2.0);
            }
            return result;
        };
    }

    /**
     * Blom rank-based inverse-normal transform with average ranks for ties.
     * Missing values remain missing and do not enter the ranks.
     */
    public static OmicsTransform rankInverseNormal() {
        return values -> {
            double[] source = required(values);
            double[] result = new double[source.length];
            Arrays.fill(result, Double.NaN);
            List<Integer> order = new ArrayList<>();
            for (int index = 0; index < source.length; index++)
                if (Double.isFinite(source[index])) order.add(index);
            if (order.isEmpty()) return result;
            order.sort(Comparator.comparingDouble(index -> source[index]));
            int first = 0;
            while (first < order.size()) {
                int last = first + 1;
                double value = source[order.get(first)];
                while (last < order.size()
                        && Double.compare(source[order.get(last)], value) == 0)
                    last++;
                double averageRank = (first + 1 + last) * 0.5;
                double probability = (averageRank - 0.375)
                    / (order.size() + 0.25);
                double transformed = Normal.quantile(
                    probability, 0, 1, true, false);
                for (int index = first; index < last; index++)
                    result[order.get(index)] = transformed;
                first = last;
            }
            return result;
        };
    }

    public static OmicsTransform compose(OmicsTransform... transforms) {
        if (transforms == null)
            throw new IllegalArgumentException("transforms are required");
        OmicsTransform result = identity();
        for (OmicsTransform transform : transforms)
            result = result.andThen(transform);
        return result;
    }

    /** Applies one transform independently to every feature row. */
    public static double[][] rows(
            double[][] values, OmicsTransform transform) {
        if (values == null || transform == null)
            throw new IllegalArgumentException("values and transform are required");
        double[][] result = new double[values.length][];
        for (int row = 0; row < values.length; row++)
            result[row] = transform.apply(values[row]);
        return result;
    }

    private static double quantile(double[] sorted, double probability) {
        if (sorted.length == 1) return sorted[0];
        double index = probability * (sorted.length - 1);
        int lower = (int) Math.floor(index);
        int upper = (int) Math.ceil(index);
        double fraction = index - lower;
        return sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
    }

    private static Summary summary(double[] values) {
        int count = 0;
        double mean = 0;
        double sumSquares = 0;
        for (double value : values) {
            if (!Double.isFinite(value)) continue;
            count++;
            double difference = value - mean;
            mean += difference / count;
            sumSquares += difference * (value - mean);
        }
        return new Summary(count, mean,
            count > 1 ? Math.max(0, sumSquares / (count - 1)) : Double.NaN);
    }

    private static double[] required(double[] values) {
        if (values == null || values.length == 0)
            throw new IllegalArgumentException("omics row is required");
        return values;
    }

    private record Summary(int count, double mean, double sampleVariance) { }
}
