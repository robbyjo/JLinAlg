/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/** Deterministic conditional/secondary-signal clumping from an LD matrix. */
public final class SecondarySignalClumper {
    private SecondarySignalClumper() { }

    /**
     * Retains the strongest independent signal in each group, then continues
     * selecting secondary signals until the group limit is reached.
     */
    public static LdClumpResult clump(List<LdClumpCandidate> candidates,
            double[][] correlation, double rSquaredThreshold, int maximumSignalsPerGroup) {
        if (candidates == null || candidates.isEmpty() || correlation == null
                || correlation.length != candidates.size()
                || !(rSquaredThreshold >= 0.0 && rSquaredThreshold <= 1.0)
                || maximumSignalsPerGroup < 1)
            throw new IllegalArgumentException("secondary-signal clumping inputs are invalid");
        int n = candidates.size();
        for (int i = 0; i < n; i++) {
            if (correlation[i] == null || correlation[i].length != n
                    || Math.abs(correlation[i][i] - 1.0) > 1e-10)
                throw new IllegalArgumentException("LD correlation must be square with unit diagonal");
            for (int j = 0; j < i; j++) if (!Double.isFinite(correlation[i][j])
                    || Math.abs(correlation[i][j] - correlation[j][i]) > 1e-10
                    || Math.abs(correlation[i][j]) > 1.0 + 1e-12)
                throw new IllegalArgumentException("LD correlation must be symmetric and bounded");
        }
        List<Integer> order = new ArrayList<>(); for (int i = 0; i < n; i++) order.add(i);
        order.sort(Comparator.comparingDouble(i -> candidates.get(i).pValue()));
        List<LdClumpCandidate> retained = new ArrayList<>();
        List<LdClumpExclusion> exclusions = new ArrayList<>();
        Set<Integer> removed = new HashSet<>();
        java.util.Map<String, Integer> counts = new java.util.HashMap<>();
        for (int lead : order) {
            if (removed.contains(lead)) continue;
            LdClumpCandidate candidate = candidates.get(lead);
            int count = counts.getOrDefault(candidate.group(), 0);
            if (count >= maximumSignalsPerGroup) {
                exclusions.add(new LdClumpExclusion(candidate,
                    LdClumpExclusionReason.IN_LINKAGE_DISEQUILIBRIUM, null, Double.NaN));
                removed.add(lead);
                continue;
            }
            retained.add(candidate); counts.put(candidate.group(), count + 1);
            for (int other : order) {
                if (other == lead || removed.contains(other)
                        || !candidate.group().equals(candidates.get(other).group())) continue;
                double r2 = correlation[lead][other] * correlation[lead][other];
                if (r2 > rSquaredThreshold) {
                    removed.add(other);
                    exclusions.add(new LdClumpExclusion(candidates.get(other),
                        LdClumpExclusionReason.IN_LINKAGE_DISEQUILIBRIUM,
                        candidate.variantId(), r2));
                }
            }
        }
        return new LdClumpResult(retained, exclusions);
    }
}
