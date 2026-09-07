/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** A grouped sparse random block with a caller-supplied within-group covariance shape. */
public final class SparseCorrelatedRandomEffectBlock {
    private final String name;
    private final RandomEffectTerm term;
    private final SparsePrecisionMatrix precision;

    private SparseCorrelatedRandomEffectBlock(String name, RandomEffectTerm term,
            SparsePrecisionMatrix precision) {
        this.name = name;
        this.term = term;
        this.precision = precision;
    }

    /**
     * Creates a block whose covariance is {@code variance * covarianceShape}
     * independently for each group. The sparse fitter estimates one scale per block.
     */
    public static SparseCorrelatedRandomEffectBlock of(String name,
            List<String> groups, List<String> effectNames,
            double[][] observationEffectDesign, double[][] covarianceShape) {
        if (name == null || name.isBlank() || groups == null || groups.isEmpty()
                || effectNames == null || effectNames.isEmpty()
                || observationEffectDesign == null
                || observationEffectDesign.length != groups.size())
            throw new IllegalArgumentException("sparse correlated block dimensions are invalid");
        int effects = effectNames.size();
        if (covarianceShape == null || covarianceShape.length != effects)
            throw new IllegalArgumentException("covariance shape must match effect count");
        Map<String, Integer> groupIndex = new LinkedHashMap<>();
        for (String group : groups) {
            if (group == null || group.isBlank()) throw new IllegalArgumentException("group labels must not be blank");
            groupIndex.computeIfAbsent(group, ignored -> groupIndex.size());
        }
        int coefficients = groupIndex.size() * effects;
        int[] rowStarts = new int[groups.size() + 1];
        int[] columns = new int[groups.size() * effects];
        double[] values = new double[columns.length];
        for (int row = 0; row < groups.size(); row++) {
            if (observationEffectDesign[row] == null
                    || observationEffectDesign[row].length != effects)
                throw new IllegalArgumentException("observation effect designs must have equal width");
            rowStarts[row] = row * effects;
            int offset = groupIndex.get(groups.get(row)) * effects;
            for (int effect = 0; effect < effects; effect++) {
                columns[row * effects + effect] = offset + effect;
                values[row * effects + effect] = observationEffectDesign[row][effect];
            }
        }
        rowStarts[groups.size()] = columns.length;
        List<String> coefficientNames = new ArrayList<>(coefficients);
        for (String group : groupIndex.keySet()) for (String effect : effectNames)
            coefficientNames.add(group + ":" + effect);
        validateCovariance(covarianceShape, effects);
        double[][] inverse = inverse(covarianceShape);
        int dimension = coefficients;
        int entries = groupIndex.size() * effects * (effects + 1) / 2;
        int[] precisionRows = new int[dimension + 1];
        int[] precisionColumns = new int[entries];
        double[] precisionValues = new double[entries];
        int position = 0;
        for (int group = 0; group < groupIndex.size(); group++) {
            for (int row = 0; row < effects; row++) {
                int globalRow = group * effects + row;
                precisionRows[globalRow] = position;
                for (int column = 0; column <= row; column++) {
                    precisionColumns[position] = group * effects + column;
                    precisionValues[position++] = inverse[row][column];
                }
            }
        }
        precisionRows[dimension] = position;
        return new SparseCorrelatedRandomEffectBlock(name,
            RandomEffectTerm.ofSparseCsr(name, groups.size(), coefficients,
                rowStarts, columns, values, coefficientNames),
            new SparsePrecisionMatrix(dimension, precisionRows,
                precisionColumns, precisionValues));
    }

    public String name() { return name; }
    public RandomEffectTerm term() { return term; }
    public SparsePrecisionMatrix precision() { return precision; }

    private static void validateCovariance(double[][] matrix, int n) {
        for (int i = 0; i < n; i++) {
            if (matrix[i] == null || matrix[i].length != n)
                throw new IllegalArgumentException("covariance shape must be square");
            for (int j = 0; j < n; j++) if (!Double.isFinite(matrix[i][j])
                    || Math.abs(matrix[i][j] - matrix[j][i]) > 1e-10)
                throw new IllegalArgumentException("covariance shape must be finite and symmetric");
        }
    }

    private static double[][] inverse(double[][] matrix) {
        int n = matrix.length;
        double[][] a = new double[n][2 * n];
        for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) a[i][j] = matrix[i][j]; a[i][n + i] = 1.0; }
        for (int column = 0; column < n; column++) {
            int pivot = column;
            for (int row = column + 1; row < n; row++) if (Math.abs(a[row][column]) > Math.abs(a[pivot][column])) pivot = row;
            if (Math.abs(a[pivot][column]) < 1e-12) throw new IllegalArgumentException("covariance shape must be positive definite");
            double[] swap = a[column]; a[column] = a[pivot]; a[pivot] = swap;
            double scale = a[column][column]; for (int j = 0; j < 2 * n; j++) a[column][j] /= scale;
            for (int row = 0; row < n; row++) if (row != column) { double factor = a[row][column]; for (int j = 0; j < 2 * n; j++) a[row][j] -= factor * a[column][j]; }
        }
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) System.arraycopy(a[i], n, result[i], 0, n);
        return result;
    }
}
