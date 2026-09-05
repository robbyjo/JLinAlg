/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/** Immutable entry/exit sweeps for start-stop Cox risk sets. */
final class CoxCountingProcessPlan {
    private final List<Stratum> strata;

    private CoxCountingProcessPlan(List<Stratum> strata) {
        this.strata = List.copyOf(strata);
    }

    static CoxCountingProcessPlan prepare(CoxSurvivalData survival) {
        Set<Integer> labels = new TreeSet<>();
        for (int label : survival.strataView()) labels.add(label);
        List<Stratum> result = new ArrayList<>(labels.size());
        for (int label : labels) result.add(prepare(survival, label));
        return new CoxCountingProcessPlan(result);
    }

    List<Stratum> strata() { return strata; }

    private static Stratum prepare(CoxSurvivalData survival, int label) {
        List<Integer> rows = new ArrayList<>();
        List<Integer> events = new ArrayList<>();
        for (int row = 0; row < survival.observations(); row++) {
            if (survival.strataView()[row] != label) continue;
            rows.add(row);
            if (survival.eventView()[row]) events.add(row);
        }
        Comparator<Integer> byStart = Comparator
            .comparingDouble((Integer row) -> survival.startView()[row])
            .thenComparingInt(Integer::intValue);
        Comparator<Integer> byStop = Comparator
            .comparingDouble((Integer row) -> survival.stopView()[row])
            .thenComparingInt(Integer::intValue);
        rows.sort(byStart);
        int[] starts = rows.stream().mapToInt(Integer::intValue).toArray();
        rows.sort(byStop);
        int[] stops = rows.stream().mapToInt(Integer::intValue).toArray();
        events.sort(byStop);
        List<Double> times = new ArrayList<>();
        List<int[]> deaths = new ArrayList<>();
        for (int first = 0; first < events.size();) {
            double time = survival.stopView()[events.get(first)];
            int end = first + 1;
            while (end < events.size()
                    && survival.stopView()[events.get(end)] == time) end++;
            int[] group = new int[end - first];
            for (int index = first; index < end; index++)
                group[index - first] = events.get(index);
            times.add(time);
            deaths.add(group);
            first = end;
        }
        return new Stratum(label, starts, stops,
            times.stream().mapToDouble(Double::doubleValue).toArray(),
            deaths.toArray(int[][]::new));
    }

    record Stratum(int label, int[] rowsByStart, int[] rowsByStop,
        double[] eventTimes, int[][] deaths) { }
}
