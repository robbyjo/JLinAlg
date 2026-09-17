/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.multipletesting;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/** Cross-fitted independent-hypothesis weighting followed by weighted BH. */
public final class IndependentHypothesisWeighting {
    private IndependentHypothesisWeighting() { }

    /** Cross-weighting controls and the threshold used to estimate bin null fractions. */
    public record Options(int bins, double lambda, double maximumWeight) {
        public Options {
            if (bins < 2) throw new IllegalArgumentException("IHW needs at least two bins");
            if (!(lambda > 0.0 && lambda < 1.0))
                throw new IllegalArgumentException("lambda must be in (0,1)");
            if (!(maximumWeight >= 1.0) || !Double.isFinite(maximumWeight))
                throw new IllegalArgumentException("maximumWeight must be finite and at least one");
        }
        public static Options defaults() { return new Options(5, 0.5, 5.0); }
    }

    /** Immutable per-hypothesis weights, covariate bins, and adjusted p-values. */
    public record Result(double[] weights, int[] bins, double[] adjustedPValues) {
        public Result {
            weights = weights.clone();
            bins = bins.clone();
            adjustedPValues = adjustedPValues.clone();
        }
        @Override public double[] weights() { return weights.clone(); }
        @Override public int[] bins() { return bins.clone(); }
        @Override public double[] adjustedPValues() { return adjustedPValues.clone(); }
    }

    /**
     * Learns bin weights only from other folds. The caller must ensure that the
     * covariate is independent of null p-values and that folds are independent.
     */
    public static Result adjust(double[] pValues, double[] covariate,
            int[] folds, Options options) {
        if (options == null) throw new IllegalArgumentException("options are required");
        if (pValues == null || covariate == null || folds == null
                || pValues.length == 0 || covariate.length != pValues.length
                || folds.length != pValues.length) {
            throw new IllegalArgumentException(
                "p-values, covariates, and folds must be nonempty and aligned");
        }
        Set<Integer> distinctFolds = new HashSet<>();
        for (int index = 0; index < pValues.length; index++) {
            if (!Double.isFinite(pValues[index]) || pValues[index] < 0.0
                    || pValues[index] > 1.0 || !Double.isFinite(covariate[index])) {
                throw new IllegalArgumentException(
                    "p-values must be in [0,1] and covariates must be finite");
            }
            distinctFolds.add(folds[index]);
        }
        if (distinctFolds.size() < 2)
            throw new IllegalArgumentException("IHW cross-fitting needs at least two folds");

        int[] bins = bins(covariate, options.bins());
        double[] weights = new double[pValues.length];
        for (int heldOut : distinctFolds) {
            double[] learned = learn(pValues, bins, folds, heldOut, options);
            double sum = 0.0;
            int count = 0;
            for (int index = 0; index < weights.length; index++) {
                if (folds[index] == heldOut) {
                    weights[index] = learned[bins[index]];
                    sum += weights[index];
                    count++;
                }
            }
            if (count == 0) continue;
            double scale = count / sum;
            for (int index = 0; index < weights.length; index++)
                if (folds[index] == heldOut) weights[index] *= scale;
        }
        return new Result(weights, bins,
            WeightedBenjaminiHochberg.adjust(pValues, weights));
    }

    private static double[] learn(double[] pValues, int[] bins, int[] folds,
            int heldOut, Options options) {
        int[] count = new int[options.bins()];
        int[] above = new int[options.bins()];
        for (int index = 0; index < pValues.length; index++) {
            if (folds[index] == heldOut) continue;
            count[bins[index]]++;
            if (pValues[index] > options.lambda()) above[bins[index]]++;
        }
        double[] weights = new double[options.bins()];
        double minimum = 1.0 / options.maximumWeight();
        for (int bin = 0; bin < weights.length; bin++) {
            if (count[bin] == 0) {
                weights[bin] = 1.0;
            } else {
                double pi0 = (above[bin] + 1.0)
                    / ((count[bin] + 1.0) * (1.0 - options.lambda()));
                pi0 = Math.max(minimum, Math.min(1.0, pi0));
                weights[bin] = Math.max(minimum,
                    Math.min(options.maximumWeight(), 1.0 / pi0));
            }
        }
        return weights;
    }

    private static int[] bins(double[] covariate, int binCount) {
        double[] sorted = covariate.clone();
        Arrays.sort(sorted);
        double[] cuts = new double[binCount - 1];
        for (int index = 1; index < binCount; index++)
            cuts[index - 1] = sorted[(int) ((long) index * sorted.length / binCount)];
        int[] result = new int[covariate.length];
        for (int index = 0; index < result.length; index++) {
            int bin = 0;
            while (bin < cuts.length && covariate[index] >= cuts[bin]) bin++;
            result[index] = bin;
        }
        return result;
    }
}
