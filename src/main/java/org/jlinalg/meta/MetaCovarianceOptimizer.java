/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;
import java.util.Comparator;
import java.util.function.ToDoubleFunction;

/** Deterministic Nelder-Mead in scaled Cholesky coordinates, with a stationarity check. */
final class MetaCovarianceOptimizer {
    private MetaCovarianceOptimizer() { }
    static Result minimize(ToDoubleFunction<double[]> f, double[] start, int maxIterations, double tolerance) {
        int d = start.length;
        double[][] points = new double[d + 1][]; double[] values = new double[d + 1];
        for (int i = 0; i <= d; i++) {
            points[i] = start.clone();
            if (i > 0) points[i][i - 1] += 0.2;
            values[i] = f.applyAsDouble(points[i]);
        }
        for (int iteration = 0; iteration < maxIterations; iteration++) {
            Integer[] order = new Integer[d + 1]; for (int i = 0; i <= d; i++) order[i] = i;
            final double[] oldValues = values;
            Arrays.sort(order, Comparator.comparingDouble(i -> oldValues[i]));
            double[][] sorted = new double[d + 1][]; double[] sortedValues = new double[d + 1];
            for (int i = 0; i <= d; i++) { sorted[i] = points[order[i]]; sortedValues[i] = values[order[i]]; }
            points = sorted; values = sortedValues;
            double spread = 0;
            for (int i = 1; i <= d; i++) for (int j = 0; j < d; j++)
                spread = Math.max(spread, Math.abs(points[i][j] - points[0][j]) / (1 + Math.abs(points[0][j])));
            if (spread < Math.sqrt(tolerance) && values[d] - values[0] < tolerance * (1 + Math.abs(values[0]))) {
                // A collapsed simplex alone is insufficient: verify finite-difference stationarity.
                double score = 0;
                for (int j = 0; j < d; j++) {
                    double[] a = points[0].clone(), b = a.clone();
                    double h = 1e-5 * (1 + Math.abs(a[j])); a[j] += h; b[j] -= h;
                    score = Math.max(score, Math.abs((f.applyAsDouble(a) - f.applyAsDouble(b)) / (2 * h)));
                }
                return new Result(points[0], values[0], score < 10 * Math.sqrt(tolerance) * (1 + Math.abs(values[0])), iteration + 1);
            }
            double[] center = new double[d];
            for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) center[j] += points[i][j] / d;
            double[] reflected = move(center, points[d], -1); double fr = f.applyAsDouble(reflected);
            if (fr < values[0]) {
                double[] expanded = move(center, points[d], -2); double fe = f.applyAsDouble(expanded);
                points[d] = fe < fr ? expanded : reflected; values[d] = Math.min(fe, fr);
            } else if (fr < values[d - 1]) { points[d] = reflected; values[d] = fr; }
            else {
                double[] contracted = move(center, points[d], fr < values[d] ? -0.5 : 0.5);
                double fc = f.applyAsDouble(contracted);
                if (fc < Math.min(fr, values[d])) { points[d] = contracted; values[d] = fc; }
                else for (int i = 1; i <= d; i++) { points[i] = move(points[0], points[i], 0.5); values[i] = f.applyAsDouble(points[i]); }
            }
        }
        int best = 0; for (int i = 1; i <= d; i++) if (values[i] < values[best]) best = i;
        return new Result(points[best], values[best], false, maxIterations);
    }
    private static double[] move(double[] center, double[] point, double multiplier) {
        double[] out = new double[point.length];
        for (int i = 0; i < out.length; i++) out[i] = center[i] + multiplier * (point[i] - center[i]);
        return out;
    }
    record Result(double[] point, double objective, boolean converged, int iterations) { }
}
