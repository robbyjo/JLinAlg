/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;

import java.util.Arrays;

/** Span-0.5 local-constant tricube trend, interpolated on at most 2049 query
 * quantiles. Small inputs are evaluated exactly. Sorting is performed once. */
final class VoomTrend {
    private final double[] location, value, knots, fitted;
    VoomTrend(double[] location, double[] value, double[] queries) {
        Integer[] order = new Integer[location.length];
        for (int i = 0; i < order.length; i++) order[i] = i;
        Arrays.sort(order, java.util.Comparator.comparingDouble(i -> location[i]));
        this.location = new double[order.length]; this.value = new double[order.length];
        for (int i = 0; i < order.length; i++) {
            this.location[i] = location[order[i]]; this.value[i] = value[order[i]];
        }
        if (location.length <= 128) { knots = fitted = new double[0]; return; }
        double[] sorted = queries.clone(); Arrays.sort(sorted);
        int count = Math.min(2049, sorted.length);
        double[] candidate = new double[count];
        for (int i = 0; i < count; i++)
            candidate[i] = sorted[(int) ((long) i * (sorted.length - 1) / (count - 1))];
        knots = Arrays.stream(candidate).distinct().toArray();
        fitted = new double[knots.length];
        for (int i = 0; i < knots.length; i++) fitted[i] = exact(knots[i]);
    }
    double at(double x) {
        if (knots.length == 0) return exact(x);
        int index = Arrays.binarySearch(knots, x);
        if (index >= 0) return fitted[index];
        int upper = -index - 1;
        if (upper == 0 || upper == knots.length) return exact(x);
        double fraction = (x - knots[upper - 1]) / (knots[upper] - knots[upper - 1]);
        return fitted[upper - 1] + fraction * (fitted[upper] - fitted[upper - 1]);
    }
    private double exact(double x) {
        int n = location.length, neighbors = Math.min(n, Math.max(3, (n + 1) / 2));
        int low = 0, high = n - neighbors;
        while (low < high) {
            int middle = (low + high) >>> 1;
            if (x - location[middle] > location[middle + neighbors] - x) low = middle + 1;
            else high = middle;
        }
        double radius = Math.max(1e-12, Math.max(Math.abs(x - location[low]),
            Math.abs(x - location[low + neighbors - 1])));
        double sum = 0, weightSum = 0;
        for (int i = 0; i < n; i++) {
            double ratio = Math.abs(location[i] - x) / radius;
            if (ratio >= 1) continue;
            double kernel = 1 - ratio * ratio * ratio;
            double weight = kernel * kernel * kernel;
            sum += weight * value[i]; weightSum += weight;
        }
        return Math.max(1e-6, weightSum > 0 ? sum / weightSum
            : Arrays.stream(value).average().orElseThrow());
    }
}
