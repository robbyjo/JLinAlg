/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.glm.GlmFamily;

/** Fixed Gauss-Hermite quadrature alternative for one-dimensional random intercepts. */
public final class GlmmQuadrature {
    private static final double[] NODES = {-3.436159118837738, -2.532731674232790, -1.756683649299882, -1.036610829789514, -0.342901327223705, 0.342901327223705, 1.036610829789514, 1.756683649299882, 2.532731674232790, 3.436159118837738};
    private static final double[] WEIGHTS = {7.640432855232620e-6, 0.001343645746781233, 0.03387439445548106, 0.2401386110823147, 0.6108626337353258, 0.6108626337353258, 0.2401386110823147, 0.03387439445548106, 0.001343645746781233, 7.640432855232620e-6};
    private GlmmQuadrature() { }

    public static GlmmQuadratureResult fit(double[] response, double[][] fixedEffects,
                                           List<String> groups, GlmFamily family) {
        if (response == null || fixedEffects == null || groups == null || family == null
                || response.length != fixedEffects.length || groups.size() != response.length
                || response.length < 4 || fixedEffects[0] == null)
            throw new IllegalArgumentException("quadrature GLMM inputs are invalid");
        int rows = response.length, columns = fixedEffects[0].length; double[] x = new double[rows * columns];
        for (int row = 0; row < rows; row++) { if (fixedEffects[row] == null || fixedEffects[row].length != columns || groups.get(row) == null || groups.get(row).isBlank()) throw new IllegalArgumentException("quadrature rows are invalid"); System.arraycopy(fixedEffects[row], 0, x, row * columns, columns); family.validateResponse(response[row], 1.0); }
        List<int[]> groupRows = groups(groups); double[] beta = new double[columns]; double logSd = -0.5, current = logLikelihood(response, x, groupRows, family, beta, Math.exp(logSd)); double step = 0.2; int evaluations = 1; boolean converged = false;
        while (evaluations < 500 && step > 1e-7) { boolean improved = false; for (int parameter = 0; parameter <= columns; parameter++) for (int direction : new int[] {-1, 1}) { double[] candidate = beta.clone(); double candidateLogSd = logSd; if (parameter == columns) candidateLogSd += direction * step; else candidate[parameter] += direction * step; candidateLogSd = Math.max(-5.0, Math.min(3.0, candidateLogSd)); double value = logLikelihood(response, x, groupRows, family, candidate, Math.exp(candidateLogSd)); evaluations++; if (value > current) { beta = candidate; logSd = candidateLogSd; current = value; improved = true; } } if (!improved) step *= 0.5; if (step <= 1e-7) converged = true; }
        double sd = Math.exp(logSd), fitted = 0.0; double[] means = new double[rows]; for (int row = 0; row < rows; row++) { for (int node = 0; node < NODES.length; node++) means[row] += WEIGHTS[node] * family.inverseLink(dot(x, row, beta, columns) + Math.sqrt(2.0) * sd * NODES[node]) / Math.sqrt(Math.PI); }
        return new GlmmQuadratureResult(family.name(), beta, sd, means, current, NODES.length, evaluations, converged);
    }

    private static double logLikelihood(double[] y, double[] x, List<int[]> groupRows, GlmFamily family, double[] beta, double sd) { double result = 0.0; for (int[] rows : groupRows) { double[] terms = new double[NODES.length]; for (int node = 0; node < NODES.length; node++) { double random = Math.sqrt(2.0) * sd * NODES[node], value = Math.log(WEIGHTS[node] / Math.sqrt(Math.PI)); for (int row : rows) { double eta = dot(x, row, beta, beta.length) + random; value += family.logLikelihood(y[row], family.inverseLink(eta), 1.0, 1.0); } terms[node] = value; } result += logSumExp(terms); } return result; }
    private static List<int[]> groups(List<String> labels) { java.util.LinkedHashMap<String, List<Integer>> map = new java.util.LinkedHashMap<>(); for (int row = 0; row < labels.size(); row++) map.computeIfAbsent(labels.get(row), ignored -> new ArrayList<>()).add(row); List<int[]> result = new ArrayList<>(); for (List<Integer> values : map.values()) { int[] rows = new int[values.size()]; for (int index = 0; index < rows.length; index++) rows[index] = values.get(index); result.add(rows); } return result; }
    private static double dot(double[] x, int row, double[] beta, int columns) { double result = 0.0; for (int column = 0; column < columns; column++) result += x[row * columns + column] * beta[column]; return result; }
    private static double logSumExp(double[] values) { double maximum = values[0]; for (double value : values) maximum = Math.max(maximum, value); double sum = 0.0; for (double value : values) sum += Math.exp(value - maximum); return maximum + Math.log(sum); }
}
