/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.ewas;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import jdistlib.Normal;
import org.jlinalg.multipletesting.WeightedBenjaminiHochberg;

/** Distance-aware signed Stouffer aggregation of neighboring EWAS probes. */
public final class RegionLevelEwas {
    private RegionLevelEwas() { }

    /** Region definition and exponential spatial-correlation scale in bases. */
    public record Options(long maximumGap, int minimumProbes,
            double correlationLength) {
        public Options {
            if (maximumGap < 1) throw new IllegalArgumentException("maximumGap must be positive");
            if (minimumProbes < 2) throw new IllegalArgumentException("minimumProbes must be at least two");
            if (!(correlationLength > 0.0) || !Double.isFinite(correlationLength))
                throw new IllegalArgumentException("correlationLength must be finite and positive");
        }
        public static Options defaults() { return new Options(500L, 3, 200.0); }
    }

    /**
     * Sorts probes by coordinate, forms nonoverlapping same-chromosome runs,
     * and adjusts every retained region as one complete BH family.
     */
    public static List<EwasRegion> scan(List<EwasProbe> probes,
            String genomeBuild, Options options) {
        if (probes == null || probes.isEmpty())
            throw new IllegalArgumentException("EWAS probes are required");
        if (genomeBuild == null || genomeBuild.isBlank())
            throw new IllegalArgumentException("genomeBuild is required");
        if (options == null) throw new IllegalArgumentException("options are required");
        List<EwasProbe> sorted = new ArrayList<>(probes);
        Set<String> ids = new HashSet<>();
        Set<String> coordinates = new HashSet<>();
        for (EwasProbe probe : sorted) {
            if (probe == null || probe.id() == null || probe.id().isBlank()
                    || probe.chromosome() == null || probe.chromosome().isBlank())
                throw new IllegalArgumentException("probe IDs and chromosomes must be nonblank");
            if (!ids.add(probe.id())) throw new IllegalArgumentException(
                "duplicate probe ID: " + probe.id());
            if (probe.position() < 1) throw new IllegalArgumentException(
                "probe positions must be positive");
            if (!Double.isFinite(probe.effect()) || !Double.isFinite(probe.pValue())
                    || probe.pValue() < 0.0 || probe.pValue() > 1.0)
                throw new IllegalArgumentException(
                    "probe effects must be finite and p-values must be in [0,1]");
            String coordinate = canonicalChromosome(probe.chromosome()) + ":" + probe.position();
            if (!coordinates.add(coordinate)) throw new IllegalArgumentException(
                "duplicate probe coordinate: " + coordinate);
        }
        sorted.sort(Comparator.comparing((EwasProbe probe) ->
            canonicalChromosome(probe.chromosome()), RegionLevelEwas::compareChromosomes)
            .thenComparingLong(EwasProbe::position).thenComparing(EwasProbe::id));
        List<List<EwasProbe>> candidates = new ArrayList<>();
        List<EwasProbe> current = new ArrayList<>();
        for (EwasProbe probe : sorted) {
            if (!current.isEmpty()) {
                EwasProbe previous = current.get(current.size() - 1);
                if (!canonicalChromosome(previous.chromosome()).equals(
                        canonicalChromosome(probe.chromosome()))
                        || probe.position() - previous.position() > options.maximumGap()) {
                    if (current.size() >= options.minimumProbes())
                        candidates.add(List.copyOf(current));
                    current.clear();
                }
            }
            current.add(probe);
        }
        if (current.size() >= options.minimumProbes()) candidates.add(List.copyOf(current));
        if (candidates.isEmpty()) throw new IllegalArgumentException(
            "no region met the minimum probe count and maximum-gap definition");

        double[] pValues = new double[candidates.size()];
        List<Unadjusted> unadjusted = new ArrayList<>();
        for (int index = 0; index < candidates.size(); index++) {
            Unadjusted region = aggregate(candidates.get(index), options.correlationLength());
            unadjusted.add(region);
            pValues[index] = region.pValue;
        }
        double[] adjusted = WeightedBenjaminiHochberg.adjust(pValues);
        List<EwasRegion> result = new ArrayList<>();
        for (int index = 0; index < unadjusted.size(); index++) {
            Unadjusted region = unadjusted.get(index);
            long span = Math.max(1L, region.end - region.start + 1L);
            result.add(new EwasRegion(genomeBuild,
                canonicalChromosome(region.probes.get(0).chromosome()),
                region.start, region.end,
                region.probes.stream().map(EwasProbe::id).toList(),
                region.meanEffect, region.statistic, region.pValue,
                adjusted[index], region.probes.size() * 1000.0 / span,
                region.maximumGap, options.correlationLength()));
        }
        return List.copyOf(result);
    }

    private static Unadjusted aggregate(List<EwasProbe> probes,
            double correlationLength) {
        int count = probes.size();
        double sumZ = 0.0;
        double effect = 0.0;
        double variance = 0.0;
        long maximumGap = 0L;
        for (int left = 0; left < count; left++) {
            EwasProbe probe = probes.get(left);
            double probability = Math.max(1e-300,
                Math.min(1.0 - 1e-16, probe.pValue()));
            double magnitude = Normal.quantile(1.0 - probability / 2.0,
                0.0, 1.0, true, false);
            sumZ += Math.signum(probe.effect()) * magnitude;
            effect += probe.effect() / count;
            if (left > 0) maximumGap = Math.max(maximumGap,
                probe.position() - probes.get(left - 1).position());
            for (int right = 0; right < count; right++) {
                variance += Math.exp(-Math.abs(probe.position()
                    - probes.get(right).position()) / correlationLength);
            }
        }
        double statistic = sumZ / Math.sqrt(variance);
        double pValue = Math.min(1.0, 2.0 * Normal.cumulative(
            Math.abs(statistic), 0.0, 1.0, false, false));
        return new Unadjusted(probes, probes.get(0).position(),
            probes.get(count - 1).position(), effect, statistic,
            pValue, maximumGap);
    }

    private static String canonicalChromosome(String chromosome) {
        String result = chromosome.trim();
        if (result.toLowerCase(Locale.ROOT).startsWith("chr")) result = result.substring(3);
        if (result.isBlank()) throw new IllegalArgumentException("chromosome is blank");
        return "chr" + result.toUpperCase(Locale.ROOT);
    }

    private static int compareChromosomes(String left, String right) {
        String a = left.substring(3);
        String b = right.substring(3);
        int ai = chromosomeNumber(a);
        int bi = chromosomeNumber(b);
        int comparison = Integer.compare(ai, bi);
        return comparison != 0 ? comparison : a.compareTo(b);
    }

    private static int chromosomeNumber(String value) {
        try { return Integer.parseInt(value); }
        catch (NumberFormatException ignored) {
            return switch (value) { case "X" -> 23; case "Y" -> 24;
                case "M", "MT" -> 25; default -> 1000; };
        }
    }

    private record Unadjusted(List<EwasProbe> probes, long start, long end,
        double meanEffect, double statistic, double pValue,
        long maximumGap) { }
}
