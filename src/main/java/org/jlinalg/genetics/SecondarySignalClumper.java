/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/** Marginal LD pruning and summary-regression conditional signal selection. */
public final class SecondarySignalClumper {
    private SecondarySignalClumper() { }

    /**
     * Forward selection with backward joint-p-value elimination. Recomputes all
     * conditional tests after every addition/removal. Groups limit lead counts;
     * all selected SNPs enter the conditioning model. Run different phenotypes
     * separately. P-values are nominal, not adjusted for adaptive selection.
     * Cycles/exhaustion return converged=false with the final recomputed tests.
     */
    public static ConditionalSignalSelection select(List<ConditionalAssociation> associations,
            double[][] correlation, double pValueThreshold, int maximumSignalsPerGroup,
            double confidenceLevel) {
        if (!(pValueThreshold > 0.0 && pValueThreshold < 1.0) || maximumSignalsPerGroup < 1)
            throw new IllegalArgumentException("selection threshold must be in (0,1) and signal limit positive");
        ConditionalAssociationModel model = new ConditionalAssociationModel(associations, correlation);
        List<Integer> selected = new ArrayList<>();
        List<String> changes = new ArrayList<>();
        Set<List<Integer>> seen = new HashSet<>();
        boolean converged = false;
        for (int iteration = 0; iteration < Math.max(100, 10 * associations.size()); iteration++) {
            List<Integer> state = selected.stream().sorted().toList();
            if (!seen.add(state)) break;
            List<ConditionalAssociationResult> tests = model.condition(
                selected.stream().mapToInt(Integer::intValue).toArray(), confidenceLevel);
            int worst = -1;
            for (int index : selected) if (tests.get(index).pValue() > pValueThreshold
                    && (worst < 0 || tests.get(index).pValue() > tests.get(worst).pValue())) worst = index;
            if (worst >= 0) {
                selected.remove(Integer.valueOf(worst));
                changes.add("remove " + associations.get(worst).variantId());
                continue;
            }
            java.util.Map<String, Integer> counts = new java.util.HashMap<>();
            for (int index : selected) counts.merge(associations.get(index).group(), 1, Integer::sum);
            int best = -1;
            for (int index = 0; index < associations.size(); index++) {
                if (selected.contains(index) || counts.getOrDefault(associations.get(index).group(), 0)
                        >= maximumSignalsPerGroup || tests.get(index).pValue() > pValueThreshold) continue;
                // Compare z, avoiding ties caused by underflow of tiny p-values.
                if (best < 0 || Math.abs(tests.get(index).statistic()) > Math.abs(tests.get(best).statistic())) best = index;
            }
            if (best < 0) { converged = true; break; }
            selected.add(best);
            changes.add("add " + associations.get(best).variantId());
        }
        return new ConditionalSignalSelection(selected.stream().map(i -> associations.get(i).variantId()).toList(),
            model.condition(selected.stream().mapToInt(Integer::intValue).toArray(), confidenceLevel), changes, converged);
    }

    /**
     * Marginal-p-value LD pruning only. For conditional inference use select;
     * p-values alone cannot reconstruct conditional effects or standard errors.
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
                    || !Double.isFinite(correlation[i][i])
                    || Math.abs(correlation[i][i] - 1.0) > 1e-10)
                throw new IllegalArgumentException("LD correlation must be square with unit diagonal");
            for (int j = 0; j < i; j++) if (!Double.isFinite(correlation[i][j])
                    || !Double.isFinite(correlation[j][i])
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
            removed.add(lead);
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
