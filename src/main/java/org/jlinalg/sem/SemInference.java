/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import jdistlib.Normal;

/** Deterministic SEM inference utilities that operate on exported scores. */
public final class SemInference {
    private SemInference() { }

    /** Sobel/product-of-coefficients delta-method inference for an indirect effect. */
    public static IndirectEffect indirect(double a, double b, double varianceA,
                                          double varianceB, double covarianceAB,
                                          double confidenceLevel) {
        validateVariance(varianceA); validateVariance(varianceB);
        if (!Double.isFinite(a) || !Double.isFinite(b) || !Double.isFinite(covarianceAB)
                || !(confidenceLevel > 0.0) || !(confidenceLevel < 1.0))
            throw new IllegalArgumentException("indirect-effect inputs are invalid");
        double estimate = a * b;
        double variance = b * b * varianceA + a * a * varianceB
            + 2.0 * a * b * covarianceAB;
        double standardError = Math.sqrt(Math.max(0.0, variance));
        double critical = Normal.quantile(0.5 + confidenceLevel / 2.0, 0.0, 1.0, true, false);
        double z = standardError == 0.0 ? 0.0 : estimate / standardError;
        return new IndirectEffect(estimate, standardError, z,
            estimate - critical * standardError, estimate + critical * standardError,
            2.0 * Normal.cumulative(-Math.abs(z), 0.0, 1.0, true, false));
    }

    /** Sandwich covariance from case or cluster score rows and an inverse bread. */
    public static double[] sandwich(double[][] scores, int[] clusters, double[] bread) {
        if (scores == null || scores.length == 0 || clusters == null
                || clusters.length != scores.length || bread == null
                || bread.length == 0) {
            throw new IllegalArgumentException("robust SEM dimensions are invalid");
        }
        int parameters = bread.length;
        int dimension = (int) Math.sqrt(parameters);
        if (dimension * dimension != parameters) throw new IllegalArgumentException("bread must be square");
        double[] meat = new double[parameters];
        java.util.Map<Integer, double[]> grouped = new java.util.LinkedHashMap<>();
        for (int row = 0; row < scores.length; row++) {
            if (scores[row] == null || scores[row].length != dimension)
                throw new IllegalArgumentException("score rows must match bread");
            double[] sum = grouped.computeIfAbsent(clusters[row], ignored -> new double[dimension]);
            for (int column = 0; column < dimension; column++) sum[column] += scores[row][column];
        }
        for (double[] sum : grouped.values()) for (int row = 0; row < dimension; row++)
            for (int column = 0; column < dimension; column++) meat[row * dimension + column] += sum[row] * sum[column];
        double[] left = multiply(bread, meat, dimension), result = multiply(left, bread, dimension);
        symmetrize(result, dimension); return result;
    }

    /** One-score-per-parameter modification index approximation. */
    public static double[] modificationIndices(double[] score, double[] information) {
        if (score == null || information == null || score.length != information.length)
            throw new IllegalArgumentException("score and information dimensions must match");
        double[] result = new double[score.length];
        for (int i = 0; i < result.length; i++) {
            if (!(information[i] > 0.0) || !Double.isFinite(information[i])) result[i] = Double.NaN;
            else result[i] = score[i] * score[i] / information[i];
        }
        return result;
    }

    /** Estimates normal latent-response thresholds from ordinal category counts. */
    public static double[] ordinalThresholds(int[] categories, int categoryCount) {
        if (categories == null || categories.length < 2 || categoryCount < 2)
            throw new IllegalArgumentException("ordinal categories are invalid");
        int[] counts = new int[categoryCount];
        for (int category : categories) { if (category < 0 || category >= categoryCount) throw new IllegalArgumentException("ordinal category is outside its range"); counts[category]++; }
        double[] result = new double[categoryCount - 1]; int cumulative = 0;
        for (int index = 0; index < result.length; index++) {
            cumulative += counts[index]; double probability = (double) cumulative / categories.length;
            result[index] = Normal.quantile(Math.max(1e-12, Math.min(1.0 - 1e-12, probability)), 0.0, 1.0, true, false);
        }
        return result;
    }

    private static void validateVariance(double value) { if (!(value >= 0.0) || !Double.isFinite(value)) throw new IllegalArgumentException("variance must be finite and nonnegative"); }
    private static double[] multiply(double[] left, double[] right, int dimension) { double[] result = new double[left.length]; for (int row = 0; row < dimension; row++) for (int column = 0; column < dimension; column++) for (int k = 0; k < dimension; k++) result[row * dimension + column] += left[row * dimension + k] * right[k * dimension + column]; return result; }
    private static void symmetrize(double[] matrix, int dimension) { for (int row = 0; row < dimension; row++) for (int column = row + 1; column < dimension; column++) { double value = 0.5 * (matrix[row * dimension + column] + matrix[column * dimension + row]); matrix[row * dimension + column] = value; matrix[column * dimension + row] = value; } }

    public record IndirectEffect(double estimate, double standardError, double statistic,
                                 double confidenceLower, double confidenceUpper, double pValue) { }
}
