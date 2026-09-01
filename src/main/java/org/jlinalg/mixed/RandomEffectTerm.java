/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.jlinalg.internal.MatrixOps;

/** A named random-effect design matrix whose coefficients share one variance. */
public final class RandomEffectTerm {
    private final String name;
    private final int observations;
    private final int coefficients;
    private final double[] denseDesign;
    private final int[] rowPointers;
    private final int[] columnIndices;
    private final double[] sparseValues;
    private final List<String> coefficientNames;

    private RandomEffectTerm(
            String name, int observations, int coefficients,
            double[] design, List<String> coefficientNames) {
        this(name, observations, coefficients, design,
            null, null, null, coefficientNames);
    }

    private RandomEffectTerm(
            String name, int observations, int coefficients,
            double[] denseDesign,
            int[] rowPointers,
            int[] columnIndices,
            double[] sparseValues,
            List<String> coefficientNames) {
        this.name = Objects.requireNonNull(name, "name");
        if (name.isBlank()) {
            throw new IllegalArgumentException("random-effect term name must not be blank");
        }
        boolean dense = denseDesign != null;
        if (observations < 1 || coefficients < 1
                || (dense && denseDesign.length != observations * coefficients)
                || (!dense && !validSparse(observations, coefficients,
                    rowPointers, columnIndices, sparseValues))
                || coefficientNames.size() != coefficients) {
            throw new IllegalArgumentException("random-effect term dimensions are invalid");
        }
        this.observations = observations;
        this.coefficients = coefficients;
        this.denseDesign = dense
            ? MatrixOps.finiteCopy(denseDesign, "randomEffectDesign") : null;
        this.rowPointers = rowPointers == null ? null : rowPointers.clone();
        this.columnIndices = columnIndices == null ? null : columnIndices.clone();
        this.sparseValues = sparseValues == null ? null
            : MatrixOps.finiteCopy(sparseValues, "sparseValues");
        this.coefficientNames = List.copyOf(coefficientNames);
        for (String coefficientName : this.coefficientNames) {
            if (coefficientName == null || coefficientName.isBlank()) {
                throw new IllegalArgumentException(
                    "random-effect coefficient names must not be blank");
            }
        }
    }

    /** Creates a term from a conventional row-by-coefficient design matrix. */
    public static RandomEffectTerm of(
            String name,
            double[][] design,
            List<String> coefficientNames) {
        if (design == null || design.length == 0 || design[0] == null) {
            throw new IllegalArgumentException("random-effect design is required");
        }
        return new RandomEffectTerm(name, design.length, design[0].length,
            MatrixOps.rowMajor(design, design.length), coefficientNames);
    }

    /** Creates a random-effect term from a validated CSR design matrix. */
    public static RandomEffectTerm ofSparseCsr(
            String name,
            int observations,
            int coefficients,
            int[] rowPointers,
            int[] columnIndices,
            double[] values,
            List<String> coefficientNames) {
        return new RandomEffectTerm(name, observations, coefficients,
            null, rowPointers, columnIndices, values, coefficientNames);
    }

    /** Creates an lme4-style {@code (1 | group)} random-intercept term. */
    public static RandomEffectTerm randomIntercept(
            String name, List<String> groups) {
        return grouped(name, groups, null);
    }

    /** Creates an independent {@code (0 + covariate | group)} random-slope term. */
    public static RandomEffectTerm randomSlope(
            String name, List<String> groups, double[] covariate) {
        if (covariate == null || groups == null
                || covariate.length != groups.size()) {
            throw new IllegalArgumentException(
                "one random-slope covariate value is required per observation");
        }
        return grouped(name, groups,
            MatrixOps.finiteCopy(covariate, "covariate"));
    }

    private static RandomEffectTerm grouped(
            String name, List<String> groups, double[] covariate) {
        if (groups == null || groups.isEmpty()) {
            throw new IllegalArgumentException("group labels are required");
        }
        Map<String, Integer> indexByGroup = new LinkedHashMap<>();
        for (String group : groups) {
            if (group == null || group.isBlank()) {
                throw new IllegalArgumentException("group labels must not be blank");
            }
            indexByGroup.computeIfAbsent(group, ignored -> indexByGroup.size());
        }
        int coefficients = indexByGroup.size();
        int[] rowPointers = new int[groups.size() + 1];
        int[] columnIndices = new int[groups.size()];
        double[] values = new double[groups.size()];
        for (int observation = 0; observation < groups.size(); observation++) {
            int coefficient = indexByGroup.get(groups.get(observation));
            rowPointers[observation] = observation;
            columnIndices[observation] = coefficient;
            values[observation] = covariate == null ? 1.0 : covariate[observation];
        }
        rowPointers[groups.size()] = groups.size();
        return new RandomEffectTerm(name, groups.size(), coefficients,
            null, rowPointers, columnIndices, values,
            new ArrayList<>(indexByGroup.keySet()));
    }

    public String name() { return name; }
    public int observations() { return observations; }
    public int coefficients() { return coefficients; }
    /** Materializes and returns the row-major design matrix. */
    public double[] design() {
        if (denseDesign != null) {
            return denseDesign.clone();
        }
        double[] result = new double[observations * coefficients];
        for (int row = 0; row < observations; row++) {
            for (int index = rowPointers[row]; index < rowPointers[row + 1]; index++) {
                result[row * coefficients + columnIndices[index]] = sparseValues[index];
            }
        }
        return result;
    }
    public List<String> coefficientNames() { return coefficientNames; }

    public boolean sparse() { return denseDesign == null; }
    public int nonzeroCount() {
        return sparse() ? sparseValues.length : denseDesign.length;
    }
    public int[] rowPointers() {
        return rowPointers == null ? null : rowPointers.clone();
    }
    public int[] columnIndices() {
        return columnIndices == null ? null : columnIndices.clone();
    }
    public double[] sparseValues() {
        return sparseValues == null ? null : sparseValues.clone();
    }

    double[] denseDesignView() { return denseDesign; }
    int[] rowPointersView() { return rowPointers; }
    int[] columnIndicesView() { return columnIndices; }
    double[] sparseValuesView() { return sparseValues; }

    private static boolean validSparse(
            int observations,
            int coefficients,
            int[] rowPointers,
            int[] columnIndices,
            double[] values) {
        if (rowPointers == null || columnIndices == null || values == null
                || rowPointers.length != observations + 1
                || columnIndices.length != values.length
                || rowPointers[0] != 0
                || rowPointers[observations] != values.length) {
            return false;
        }
        for (int row = 0; row < observations; row++) {
            if (rowPointers[row] > rowPointers[row + 1]) {
                return false;
            }
            int previous = -1;
            for (int index = rowPointers[row]; index < rowPointers[row + 1]; index++) {
                int column = columnIndices[index];
                if (column < 0 || column >= coefficients || column <= previous) {
                    return false;
                }
                previous = column;
            }
        }
        return true;
    }
}
