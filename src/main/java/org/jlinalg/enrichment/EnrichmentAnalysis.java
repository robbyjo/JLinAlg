/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.enrichment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Database-independent one-sided gene-set enrichment with explicit universes. */
public final class EnrichmentAnalysis {
    private EnrichmentAnalysis() { }
    /** Correction family is every eligible term supplied to a run. */
    public enum Fdr { BH, BY, NONE }
    /** Ordinary ORA defaults, with terms filtered only by background membership. */
    public static List<EnrichmentResult> ora(Set<String> universe, Set<String> selected,
            List<GeneSet> sets, int minimumSize, int maximumSize, Fdr correction) {
        return run(universe, selected, sets, minimumSize, maximumSize, correction, null, 100_000_000L);
    }
    /** EWAS ORA using equivalent-probe bias and fractional gene counts. */
    public static List<EnrichmentResult> gsameth(Map<String, Set<String>> eligibleMapping,
            Set<String> selectedFeatures, List<GeneSet> sets, int minimumSize,
            int maximumSize, Fdr correction, long maximumUpdates) {
        Map<String, Gsameth.GeneWeight> weights = Gsameth.weights(eligibleMapping, selectedFeatures);
        Set<String> selected = new HashSet<>();
        weights.forEach((g, w) -> { if (w.selectedWeight() > 0) selected.add(g); });
        return run(weights.keySet(), selected, sets, minimumSize, maximumSize,
            correction, weights, maximumUpdates);
    }
    private static List<EnrichmentResult> run(Set<String> universe, Set<String> selected,
            List<GeneSet> sets, int minimumSize, int maximumSize, Fdr correction,
            Map<String, Gsameth.GeneWeight> weights, long maximumUpdates) {
        if (universe.isEmpty() || !universe.containsAll(selected)
                || minimumSize < 1 || maximumSize < minimumSize || correction == null)
            throw new IllegalArgumentException("invalid universe, selected subset, or size filters");
        int N = universe.size(), n = selected.size();
        List<EnrichmentResult> results = new ArrayList<>();
        Set<String> ids = new HashSet<>();
        for (GeneSet set : sets) {
            if (!ids.add(set.id())) throw new IllegalArgumentException("duplicate term ID: " + set.id());
            Set<String> members = new HashSet<>(set.genes()); members.retainAll(universe);
            int K = members.size();
            if (K < minimumSize || K > maximumSize) continue;
            List<String> hits = members.stream().filter(selected::contains).sorted().toList();
            int k = hits.size();
            double weighted = k, selectionOdds = 1;
            double logP;
            if (weights == null) logP = jdistlib.HyperGeometric.cumulative(k - 1.0, K, N - K, n, false, true);
            else {
                weighted = hits.stream().mapToDouble(g -> weights.get(g).selectedWeight()).sum();
                double in = 0, out = 0;
                for (var entry : weights.entrySet()) {
                    if (members.contains(entry.getKey())) in += entry.getValue().probability();
                    else out += entry.getValue().probability();
                }
                selectionOdds = K == N || n == 0 ? 1 : (in / K) / (out / (N - K));
                if (Double.isNaN(selectionOdds)) throw new ArithmeticException("undefined gsameth selection odds");
                // BiasedUrn's R interface truncates fractional x to integer before both calls.
                int threshold = (int) weighted;
                logP = Math.log(Wallenius.upperTail(K, N - K, n, threshold, selectionOdds, maximumUpdates));
            }
            if (Double.isNaN(logP) || logP > 1e-12)
                throw new ArithmeticException("invalid enrichment probability for " + set.id());
            double expected = (double) n * K / N;
            double denominator = (double) (n - k) * (K - k);
            double numerator = (double) k * (N - K - n + k);
            double oddsRatio = denominator == 0 ? numerator == 0 ? Double.NaN : Double.POSITIVE_INFINITY
                : numerator / denominator;
            results.add(new EnrichmentResult(set.id(), set.name(), N, n, K, k, weighted,
                expected, expected == 0 ? Double.NaN : k / expected, oddsRatio, selectionOdds,
                Math.exp(logP), logP, 0, 0, hits));
        }
        double[] adjusted = adjustLog(results.stream().mapToDouble(EnrichmentResult::logPValue).toArray(), correction);
        List<EnrichmentResult> output = new ArrayList<>();
        for (int i = 0; i < results.size(); i++) {
            var r = results.get(i);
            output.add(new EnrichmentResult(r.id(), r.name(), r.universeSize(), r.selectedSize(),
                r.setSize(), r.overlap(), r.weightedOverlap(), r.expectedOverlap(), r.foldEnrichment(),
                r.oddsRatio(), r.selectionOdds(), r.pValue(), r.logPValue(), Math.exp(adjusted[i]), adjusted[i], r.hitGenes()));
        }
        output.sort(Comparator.comparingDouble(EnrichmentResult::logPValue).thenComparing(EnrichmentResult::id));
        return List.copyOf(output);
    }
    /** Monotone BH/BY adjustment in log space, preserving input order. */
    public static double[] adjustLog(double[] logP, Fdr method) {
        if (method == null) throw new IllegalArgumentException("FDR method is required");
        double[] result = logP.clone();
        for (double p : logP) if (Double.isNaN(p) || p > 0)
            throw new IllegalArgumentException("log probabilities must be <= 0");
        if (method == Fdr.NONE || logP.length == 0) return result;
        Integer[] order = new Integer[logP.length];
        for (int i = 0; i < order.length; i++) order[i] = i;
        Arrays.sort(order, Comparator.comparingDouble(i -> logP[i]));
        double harmonic = 1;
        if (method == Fdr.BY) { harmonic = 0; for (int i = 1; i <= order.length; i++) harmonic += 1.0 / i; }
        double previous = 0, factor = Math.log(order.length) + Math.log(harmonic);
        for (int i = order.length - 1; i >= 0; i--) {
            previous = Math.min(previous, logP[order[i]] + factor - Math.log(i + 1.0));
            result[order[i]] = previous;
        }
        return result;
    }
}
