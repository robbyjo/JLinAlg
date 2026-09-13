/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.enrichment;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/** Gene-level inputs for missMethyl gsameth's equivalent-CpG / fractional-count method.
 * Annotation is supplied explicitly as unique feature-to-gene edges. Rank ties use
 * lexicographic gene ID order, equivalent to the R reference with LC_COLLATE=C.
 */
public final class Gsameth {
    private Gsameth() { }

    /** Effective coverage, fractional selected contribution, and fitted selection probability. */
    public record GeneWeight(double equivalentProbes, double selectedWeight, double probability) { }

    /** Fits the reference span-0.5 rank tricube probability weighting function. */
    public static Map<String, GeneWeight> weights(Map<String, Set<String>> eligibleMapping,
            Set<String> selectedFeatures) {
        if (!eligibleMapping.keySet().containsAll(selectedFeatures))
            throw new IllegalArgumentException("selected features must belong to the universe");
        Map<String, double[]> counts = new TreeMap<>();
        // Sort features too, making floating sums independent of input row order.
        for (String feature : new java.util.TreeSet<>(eligibleMapping.keySet())) {
            Set<String> genes = eligibleMapping.get(feature);
            if (genes.isEmpty()) continue;
            double weight = 1.0 / genes.size();
            for (String gene : genes) {
                double[] c = counts.computeIfAbsent(gene, ignored -> new double[4]);
                add(c, 0, weight);
                if (selectedFeatures.contains(feature)) add(c, 1, weight);
            }
        }
        List<String> ordered = new ArrayList<>(counts.keySet());
        ordered.sort(Comparator.comparingDouble(g -> counts.get(g)[0]));
        double[] selected = new double[ordered.size()];
        for (int i = 0; i < selected.length; i++) selected[i] = counts.get(ordered.get(i))[1] > 0 ? 1 : 0;
        double[] probability = tricube(selected);
        Map<String, GeneWeight> result = new TreeMap<>();
        for (int i = 0; i < ordered.size(); i++) {
            double[] c = counts.get(ordered.get(i));
            result.put(ordered.get(i), new GeneWeight(c[0], Math.min(1, c[1]), probability[i]));
        }
        return java.util.Collections.unmodifiableMap(result);
    }

    // Avoid moving an integer fractional-overlap threshold down by one because
    // repeated fractions (for example ten contributions of 0.1) lost low bits.
    private static void add(double[] counts, int index, double value) {
        double adjusted = value - counts[index + 2];
        double sum = counts[index] + adjusted;
        counts[index + 2] = (sum - counts[index]) - adjusted;
        counts[index] = sum;
    }

    // limma::tricubeMovingAverage: nearest odd width, positive kernel endpoints,
    // boundary renormalization. This is a moving average over ranks, not LOESS.
    static double[] tricube(double[] values) {
        int n = values.length, half = n / 4, width = 2 * half + 1;
        if (width > n) { width -= 2; half--; }
        if (half <= 0) return values.clone();
        double[] kernel = new double[half + 1];
        for (int d = 0; d <= half; d++) {
            double u = (double) d / half * width / (width + 1.0);
            kernel[d] = Math.pow(1 - u * u * u, 3);
        }
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0, total = 0;
            for (int j = Math.max(0, i - half); j <= Math.min(n - 1, i + half); j++) {
                double w = kernel[Math.abs(j - i)];
                sum += w * values[j]; total += w;
            }
            result[i] = sum / total;
        }
        return result;
    }
}
