/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

/** Numerically stable aggregated Cauchy association test (ACAT). */
public final class Acat {
    private static final double LOG_10 = Math.log(10);

    private Acat() { }

    /** Combines p-values with equal weights. */
    public static Result combine(double[] pValues) {
        if (pValues == null)
            throw new IllegalArgumentException("p-values are required");
        double[] weights = new double[pValues.length];
        java.util.Arrays.fill(weights, 1);
        return combine(pValues, weights);
    }

    /**
     * Combines p-values with nonnegative weights. Weights are normalized
     * internally and need not sum to one.
     */
    public static Result combine(double[] pValues, double[] weights) {
        double[] normalized = normalizedWeights(pValues, weights);
        boolean zero = false, one = false;
        for (int index = 0; index < pValues.length; index++) {
            double p = pValues[index];
            if (!Double.isFinite(p) || p < 0 || p > 1)
                throw new IllegalArgumentException(
                    "ACAT p-values must be finite and within [0,1]");
            if (normalized[index] > 0) {
                zero |= p == 0;
                one |= p == 1;
            }
        }
        if (zero && one)
            throw new IllegalArgumentException(
                "ACAT cannot combine positive-weight p-values equal to both zero and one");
        if (zero)
            return new Result(Double.POSITIVE_INFINITY, 0,
                Double.NEGATIVE_INFINITY, positive(normalized));
        if (one)
            return new Result(Double.NEGATIVE_INFINITY, 1, 0,
                positive(normalized));

        // Scale extreme tails before division, retaining a representable p even
        // when the corresponding Cauchy statistic exceeds Double.MAX_VALUE.
        double scale = 1;
        for (int i = 0; i < pValues.length; i++) if (normalized[i] > 0
                && normalized[i] / Math.min(pValues[i], 1 - pValues[i]) > 1e150)
            scale = 1e-300;
        // Neumaier summation matters when moderately extreme tails oppose.
        double sum = 0, correction = 0;
        for (int index = 0; index < pValues.length; index++) {
            if (normalized[index] == 0) continue;
            double term = transform(pValues[index], normalized[index], scale);
            double next = sum + term;
            correction += Math.abs(sum) >= Math.abs(term)
                ? (sum - next) + term : (term - next) + sum;
            sum = next;
        }
        double scaled = sum + correction;
        double statistic = scaled / scale;
        if (Double.isNaN(statistic))
            throw new IllegalArgumentException(
                "ACAT statistic is unresolved because opposing tails overflow");
        double pValue = scaled > 0 ? Math.atan(scale / scaled) / Math.PI
            : scaled < 0 ? 1 - Math.atan(-scale / scaled) / Math.PI : 0.5;
        double log10 = pValue > 0 ? Math.log(pValue) / LOG_10
            : Double.NEGATIVE_INFINITY;
        return new Result(statistic, pValue, log10, positive(normalized));
    }

    static double[] normalizedWeights(double[] pValues, double[] weights) {
        if (pValues == null || weights == null || pValues.length == 0
                || weights.length != pValues.length)
            throw new IllegalArgumentException(
                "ACAT requires equally sized, nonempty p-value and weight arrays");
        double maximum = 0;
        for (double weight : weights) {
            if (!Double.isFinite(weight) || weight < 0)
                throw new IllegalArgumentException(
                    "ACAT weights must be finite and nonnegative");
            maximum = Math.max(maximum, weight);
        }
        if (!(maximum > 0))
            throw new IllegalArgumentException(
                "at least one ACAT weight must be positive");
        double scaledSum = 0;
        for (double weight : weights) scaledSum += weight / maximum;
        double[] normalized = new double[weights.length];
        for (int index = 0; index < weights.length; index++)
            normalized[index] = (weights[index] / maximum) / scaledSum;
        return normalized;
    }

    private static int positive(double[] values) {
        int result = 0;
        for (double value : values) if (value > 0) result++;
        return result;
    }

    private static double transform(double pValue, double weight, double scale) {
        if (pValue < 1e-8) return scaledReciprocal(weight, pValue, scale) / Math.PI;
        double upper = 1 - pValue;
        if (upper < 1e-8) return -scaledReciprocal(weight, upper, scale) / Math.PI;
        if (pValue < 0.25) return weight * (scale / Math.tan(Math.PI * pValue));
        if (pValue > 0.75) return -weight * (scale / Math.tan(Math.PI * upper));
        return weight * scale * Math.tan(Math.PI * (0.5 - pValue));
    }

    private static double scaledReciprocal(double weight, double probability, double scale) {
        double ratio = weight / probability;
        // Taking weight/probability first also preserves subnormal weights
        // attached to subnormal p-values when their ratio is ordinary sized.
        return Double.isInfinite(ratio) ? (scale / probability) * weight : ratio * scale;
    }

    /** ACAT statistic and its standard-Cauchy upper-tail probability. */
    public record Result(
            double statistic, double pValue, double log10PValue,
            int components) {
        public double negativeLog10PValue() { return -log10PValue; }
    }
}
