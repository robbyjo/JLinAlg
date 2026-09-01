/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/** Deterministic sample-ID alignment between data sources. */
public record SampleAlignment(
        List<String> sampleIds, int[] leftOrder, int[] rightOrder) {
    public SampleAlignment {
        sampleIds = List.copyOf(sampleIds);
        leftOrder = leftOrder.clone();
        rightOrder = rightOrder.clone();
        if (sampleIds.size() != leftOrder.length || sampleIds.size() != rightOrder.length)
            throw new IllegalArgumentException("alignment lengths must match");
    }

    @Override public int[] leftOrder() { return leftOrder.clone(); }
    @Override public int[] rightOrder() { return rightOrder.clone(); }

    /** Aligns in left/source order, retaining only IDs present in both lists. */
    public static SampleAlignment intersect(
            List<String> left, List<String> right) {
        Map<String, Integer> rightIndex = indexes(right, "right");
        indexes(left, "left");
        List<String> ids = new ArrayList<>();
        List<Integer> leftRows = new ArrayList<>();
        List<Integer> rightRows = new ArrayList<>();
        for (int index = 0; index < left.size(); index++) {
            Integer match = rightIndex.get(left.get(index));
            if (match != null) {
                ids.add(left.get(index));
                leftRows.add(index);
                rightRows.add(match);
            }
        }
        if (ids.isEmpty())
            throw new IllegalArgumentException("data sources have no sample IDs in common");
        return new SampleAlignment(ids, integers(leftRows), integers(rightRows));
    }

    /** Requires every requested analysis ID and returns its source index. */
    public static int[] requireOrder(
            List<String> source, List<String> requested) {
        Map<String, Integer> sourceIndex = indexes(source, "source");
        indexes(requested, "requested");
        int[] order = new int[requested.size()];
        for (int index = 0; index < order.length; index++) {
            Integer match = sourceIndex.get(requested.get(index));
            if (match == null)
                throw new IllegalArgumentException(
                    "sample is absent from source: " + requested.get(index));
            order[index] = match;
        }
        return order;
    }

    private static Map<String, Integer> indexes(
            List<String> values, String name) {
        if (values == null || values.isEmpty())
            throw new IllegalArgumentException(name + " sample IDs are required");
        Map<String, Integer> result = new HashMap<>();
        HashSet<String> seen = new HashSet<>();
        for (int index = 0; index < values.size(); index++) {
            String value = values.get(index);
            if (value == null || value.isBlank() || !seen.add(value))
                throw new IllegalArgumentException(
                    name + " sample IDs must be unique and nonblank");
            result.put(value, index);
        }
        return result;
    }

    private static int[] integers(List<Integer> values) {
        int[] result = new int[values.size()];
        for (int index = 0; index < result.length; index++)
            result[index] = values.get(index);
        return result;
    }
}
