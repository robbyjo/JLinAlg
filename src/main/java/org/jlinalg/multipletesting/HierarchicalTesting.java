/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.multipletesting;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Hierarchy-respecting weighted Bonferroni testing under arbitrary dependence. */
public final class HierarchicalTesting {
    private HierarchicalTesting() { }

    /** Per-node normalized weights, gated adjusted p-values, and decisions. */
    public record Result(double[] weights, double[] adjustedPValues,
            boolean[] rejected) {
        public Result {
            weights = weights.clone();
            adjustedPValues = adjustedPValues.clone();
            rejected = rejected.clone();
        }
        @Override public double[] weights() { return weights.clone(); }
        @Override public double[] adjustedPValues() { return adjustedPValues.clone(); }
        @Override public boolean[] rejected() { return rejected.clone(); }
    }

    /**
     * Tests every named node in a prespecified forest. Empty parent IDs denote
     * roots. Leaf-count weights are normalized across all nodes, and a child
     * cannot be rejected unless every tested ancestor is rejected.
     */
    public static Result adjust(String[] ids, String[] parentIds,
            double[] pValues, double alpha) {
        if (ids == null || parentIds == null || pValues == null || ids.length == 0
                || parentIds.length != ids.length || pValues.length != ids.length) {
            throw new IllegalArgumentException(
                "hierarchy IDs, parents, and p-values must be nonempty and aligned");
        }
        if (!(alpha > 0.0 && alpha <= 1.0) || !Double.isFinite(alpha))
            throw new IllegalArgumentException("alpha must be in (0,1]");
        Map<String, Integer> index = new HashMap<>();
        for (int node = 0; node < ids.length; node++) {
            if (ids[node] == null || ids[node].isBlank()
                    || index.put(ids[node], node) != null) {
                throw new IllegalArgumentException("hierarchy IDs must be unique and nonblank");
            }
            if (!Double.isFinite(pValues[node]) || pValues[node] < 0.0
                    || pValues[node] > 1.0)
                throw new IllegalArgumentException("p-values must be finite and in [0,1]");
        }
        int[] parent = new int[ids.length];
        Arrays.fill(parent, -1);
        List<List<Integer>> children = new ArrayList<>();
        for (int node = 0; node < ids.length; node++) children.add(new ArrayList<>());
        for (int node = 0; node < ids.length; node++) {
            String namedParent = parentIds[node];
            if (namedParent == null || namedParent.isBlank()) continue;
            Integer found = index.get(namedParent);
            if (found == null) throw new IllegalArgumentException(
                "parent is absent from the complete family: " + namedParent);
            if (found == node) throw new IllegalArgumentException("a hierarchy node cannot parent itself");
            parent[node] = found;
            children.get(found).add(node);
        }
        int[] leaves = new int[ids.length];
        int[] depth = new int[ids.length];
        for (int node = 0; node < ids.length; node++) {
            depth[node] = depth(node, parent, new HashSet<>());
            leaves[node] = leaves(node, children, new HashSet<>());
        }
        double denominator = Arrays.stream(leaves).sum();
        double[] weights = new double[ids.length];
        for (int node = 0; node < ids.length; node++) weights[node] = leaves[node] / denominator;
        Integer[] order = new Integer[ids.length];
        for (int node = 0; node < ids.length; node++) order[node] = node;
        Arrays.sort(order, (left, right) -> {
            int comparison = Integer.compare(depth[left], depth[right]);
            return comparison != 0 ? comparison : Integer.compare(left, right);
        });
        double[] adjusted = new double[ids.length];
        boolean[] rejected = new boolean[ids.length];
        for (int node : order) {
            double value = Math.min(1.0, pValues[node] / weights[node]);
            if (parent[node] >= 0) value = Math.max(value, adjusted[parent[node]]);
            adjusted[node] = value;
            rejected[node] = value <= alpha;
        }
        return new Result(weights, adjusted, rejected);
    }

    private static int depth(int node, int[] parent, Set<Integer> path) {
        if (!path.add(node)) throw new IllegalArgumentException("hierarchy contains a cycle");
        int result = parent[node] < 0 ? 0 : 1 + depth(parent[node], parent, path);
        path.remove(node);
        return result;
    }

    private static int leaves(int node, List<List<Integer>> children,
            Set<Integer> path) {
        if (!path.add(node)) throw new IllegalArgumentException("hierarchy contains a cycle");
        if (children.get(node).isEmpty()) {
            path.remove(node);
            return 1;
        }
        int result = 0;
        for (int child : children.get(node)) result += leaves(child, children, path);
        path.remove(node);
        return result;
    }
}
