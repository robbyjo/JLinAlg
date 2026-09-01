/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.internal.MatrixOps;

/** Grouped random coefficients sharing one unstructured covariance matrix. */
public final class CorrelatedRandomEffectBlock {
    private final String name;
    private final int observations;
    private final int effectCount;
    private final int[] groupIndices;
    private final List<String> groupNames;
    private final List<String> effectNames;
    private final double[] effectDesign;

    private CorrelatedRandomEffectBlock(
            String name, int observations, int effectCount,
            int[] groupIndices, List<String> groupNames,
            List<String> effectNames, double[] effectDesign) {
        this.name = name;
        this.observations = observations;
        this.effectCount = effectCount;
        this.groupIndices = groupIndices;
        this.groupNames = groupNames;
        this.effectNames = effectNames;
        this.effectDesign = effectDesign;
    }

    public static CorrelatedRandomEffectBlock of(
            String name,
            List<String> groups,
            List<String> effectNames,
            double[][] observationEffectDesign) {
        if (name == null || name.isBlank() || groups == null
                || groups.isEmpty() || effectNames == null
                || effectNames.isEmpty() || observationEffectDesign == null
                || observationEffectDesign.length != groups.size())
            throw new IllegalArgumentException(
                "correlated random block dimensions are invalid");
        double[] design = MatrixOps.rowMajor(
            observationEffectDesign, groups.size());
        int effects = observationEffectDesign[0].length;
        if (effectNames.size() != effects
                || effectNames.stream().anyMatch(
                    value -> value == null || value.isBlank()))
            throw new IllegalArgumentException(
                "one nonblank name is required per random effect");
        Map<String, Integer> groupMap = new LinkedHashMap<>();
        int[] indices = new int[groups.size()];
        for (int row = 0; row < groups.size(); row++) {
            String group = groups.get(row);
            if (group == null || group.isBlank())
                throw new IllegalArgumentException(
                    "random-effect groups must not be blank");
            indices[row] = groupMap.computeIfAbsent(
                group, ignored -> groupMap.size());
        }
        return new CorrelatedRandomEffectBlock(name, groups.size(), effects,
            indices, List.copyOf(groupMap.keySet()), List.copyOf(effectNames),
            design);
    }

    public String name() { return name; }
    public int observations() { return observations; }
    public int effectCount() { return effectCount; }
    public int groupCount() { return groupNames.size(); }
    public List<String> groupNames() { return groupNames; }
    public List<String> effectNames() { return effectNames; }
    public int[] groupIndices() { return groupIndices.clone(); }
    public double[] effectDesign() { return effectDesign.clone(); }
}
