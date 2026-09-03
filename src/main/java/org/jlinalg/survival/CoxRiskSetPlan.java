/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/** Immutable ordering and event groups reused across Cox evaluations. */
final class CoxRiskSetPlan {
    private final List<Stratum> strata;

    private CoxRiskSetPlan(List<Stratum> strata) {
        this.strata = List.copyOf(strata);
    }

    static CoxRiskSetPlan prepare(CoxSurvivalData survival) {
        for (double value : survival.startView())
            if (value != 0) return null;
        Set<Integer> labels = new TreeSet<>();
        for (int value : survival.strataView()) labels.add(value);
        List<Stratum> result = new ArrayList<>(labels.size());
        for (int label : labels) result.add(prepare(survival, label));
        return new CoxRiskSetPlan(result);
    }

    List<Stratum> strata() {
        return strata;
    }

    private static Stratum prepare(CoxSurvivalData survival, int label) {
        List<Integer> ordered = new ArrayList<>();
        List<Integer> events = new ArrayList<>();
        for (int row = 0; row < survival.observations(); row++) {
            if (survival.strataView()[row] != label) continue;
            ordered.add(row);
            if (survival.eventView()[row]) events.add(row);
        }
        Comparator<Integer> descending = Comparator.comparingDouble(
            (Integer row) -> survival.stopView()[row]).reversed();
        ordered.sort(descending);
        events.sort(descending);
        List<Double> times = new ArrayList<>();
        List<int[]> deaths = new ArrayList<>();
        for (int start = 0; start < events.size();) {
            double time = survival.stopView()[events.get(start)];
            int end = start + 1;
            while (end < events.size()
                    && survival.stopView()[events.get(end)] == time) end++;
            int[] group = new int[end - start];
            for (int index = start; index < end; index++)
                group[index - start] = events.get(index);
            times.add(time);
            deaths.add(group);
            start = end;
        }
        int[] orderedRows = ordered.stream().mapToInt(Integer::intValue).toArray();
        double[] eventTimes = times.stream().mapToDouble(Double::doubleValue).toArray();
        return new Stratum(label, orderedRows, eventTimes,
            deaths.toArray(int[][]::new));
    }

    record Stratum(
        int label, int[] orderedRows, double[] eventTimes, int[][] deaths) { }
}
