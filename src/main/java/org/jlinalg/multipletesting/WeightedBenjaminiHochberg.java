/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.multipletesting;

import java.util.Arrays;

/** Weighted Benjamini-Hochberg adjustment with positive mean-one weights. */
public final class WeightedBenjaminiHochberg {
    private WeightedBenjaminiHochberg() { }

    /**
     * Adjusts a complete prespecified family. Weights are normalized to have
     * mean one before adjustment; larger weights give a hypothesis more power.
     */
    public static double[] adjust(double[] pValues, double[] weights) {
        validate(pValues, weights);
        int count = pValues.length;
        double mean = Arrays.stream(weights).sum() / count;
        Integer[] order = new Integer[count];
        double[] scaled = new double[count];
        for (int index = 0; index < count; index++) {
            order[index] = index;
            scaled[index] = Math.min(1.0, pValues[index] * mean / weights[index]);
        }
        Arrays.sort(order, (left, right) -> {
            int comparison = Double.compare(scaled[left], scaled[right]);
            return comparison != 0 ? comparison : Integer.compare(left, right);
        });
        double[] adjusted = new double[count];
        double next = 1.0;
        for (int rank = count; rank >= 1; rank--) {
            int index = order[rank - 1];
            next = Math.min(next, scaled[index] * count / rank);
            adjusted[index] = Math.min(1.0, next);
        }
        return adjusted;
    }

    /** Ordinary Benjamini-Hochberg adjustment. */
    public static double[] adjust(double[] pValues) {
        double[] weights = new double[pValues.length];
        Arrays.fill(weights, 1.0);
        return adjust(pValues, weights);
    }

    private static void validate(double[] pValues, double[] weights) {
        if (pValues == null || weights == null || pValues.length == 0
                || pValues.length != weights.length) {
            throw new IllegalArgumentException(
                "p-values and weights must be nonempty and have equal length");
        }
        for (int index = 0; index < pValues.length; index++) {
            if (!Double.isFinite(pValues[index]) || pValues[index] < 0.0
                    || pValues[index] > 1.0) {
                throw new IllegalArgumentException(
                    "p-values must be finite and in [0,1]");
            }
            if (!Double.isFinite(weights[index]) || !(weights[index] > 0.0)) {
                throw new IllegalArgumentException(
                    "multiple-testing weights must be finite and positive");
            }
        }
    }
}
