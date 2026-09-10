/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import org.jlinalg.model.MissingDataPolicy;

/** Immutable named numeric and categorical columns used by formula compilation. */
public final class ModelTable {
    private final int rows;
    private final Map<String, double[]> numeric;
    private final Map<String, String[]> categorical;

    private ModelTable(Builder builder) {
        this.rows = builder.rows;
        this.numeric = copyNumeric(builder.numeric);
        this.categorical = copyCategorical(builder.categorical);
    }

    public static Builder builder(int rows) { return new Builder(rows); }
    public int rows() { return rows; }

    /** Returns a defensive copy of a named numeric column. */
    public double[] numericColumn(String name) { return numeric(name).clone(); }

    double[] numeric(String name) {
        double[] value = numeric.get(name);
        if (value == null) {
            throw new IllegalArgumentException("numeric column is absent: " + name);
        }
        return value;
    }

    String[] categorical(String name) {
        String[] value = categorical.get(name);
        if (value == null) {
            throw new IllegalArgumentException("categorical column is absent: " + name);
        }
        return value;
    }

    boolean isNumeric(String name) { return numeric.containsKey(name); }
    boolean isCategorical(String name) { return categorical.containsKey(name); }

    int[] completeRows(Set<String> columns, MissingDataPolicy policy) {
        if (columns == null || columns.isEmpty() || policy == null) {
            throw new IllegalArgumentException(
                "formula columns and missing-data policy are required");
        }
        for (String column : columns) {
            if (!isNumeric(column) && !isCategorical(column)) {
                throw new IllegalArgumentException("column is absent: " + column);
            }
        }
        int[] retained = new int[rows];
        int count = 0;
        for (int row = 0; row < rows; row++) {
            String missing = null;
            for (String column : columns) {
                if (isNumeric(column) && !Double.isFinite(numeric.get(column)[row])
                        || isCategorical(column) && categorical.get(column)[row] == null) {
                    missing = column;
                    break;
                }
            }
            if (missing == null) {
                retained[count++] = row;
            } else if (policy == MissingDataPolicy.ERROR) {
                throw new IllegalArgumentException(
                    "model table contains a missing value in column "
                        + missing + " at row " + row);
            }
        }
        if (count == 0) {
            throw new IllegalArgumentException("no complete observations remain");
        }
        return java.util.Arrays.copyOf(retained, count);
    }

    ModelTable retainRows(int[] retained) {
        if (retained == null || retained.length == 0) {
            throw new IllegalArgumentException("retained rows are required");
        }
        if (retained.length == rows) {
            boolean identity = true;
            for (int row = 0; row < rows; row++) identity &= retained[row] == row;
            if (identity) return this;
        }
        Builder builder = builder(retained.length);
        for (Map.Entry<String, double[]> entry : numeric.entrySet()) {
            double[] values = new double[retained.length];
            for (int row = 0; row < retained.length; row++) {
                values[row] = entry.getValue()[checkedRow(retained, row)];
            }
            builder.numeric(entry.getKey(), values);
        }
        for (Map.Entry<String, String[]> entry : categorical.entrySet()) {
            String[] values = new String[retained.length];
            for (int row = 0; row < retained.length; row++) {
                values[row] = entry.getValue()[checkedRow(retained, row)];
            }
            builder.categorical(entry.getKey(), values);
        }
        return builder.build();
    }

    private int checkedRow(int[] retained, int position) {
        int row = retained[position];
        if (row < 0 || row >= rows || position > 0 && row <= retained[position - 1]) {
            throw new IllegalArgumentException(
                "retained rows must be strictly increasing and in range");
        }
        return row;
    }

    private static Map<String, double[]> copyNumeric(Map<String, double[]> source) {
        Map<String, double[]> result = new LinkedHashMap<>();
        source.forEach((name, value) -> result.put(name, value.clone()));
        return Map.copyOf(result);
    }

    private static Map<String, String[]> copyCategorical(Map<String, String[]> source) {
        Map<String, String[]> result = new LinkedHashMap<>();
        source.forEach((name, value) -> result.put(name, value.clone()));
        return Map.copyOf(result);
    }

    /** Builder retaining caller-specified column order. */
    public static final class Builder {
        private final int rows;
        private final Map<String, double[]> numeric = new LinkedHashMap<>();
        private final Map<String, String[]> categorical = new LinkedHashMap<>();

        private Builder(int rows) {
            if (rows < 1) throw new IllegalArgumentException("rows must be positive");
            this.rows = rows;
        }

        public Builder numeric(String name, double... values) {
            validate(name, values == null ? -1 : values.length);
            numeric.put(name, values.clone());
            categorical.remove(name);
            return this;
        }

        public Builder categorical(String name, String... values) {
            validate(name, values == null ? -1 : values.length);
            categorical.put(name, values.clone());
            numeric.remove(name);
            return this;
        }

        private void validate(String name, int length) {
            if (name == null || name.isBlank()) {
                throw new IllegalArgumentException("column name must not be blank");
            }
            if (length != rows) {
                throw new IllegalArgumentException(
                    "column length must equal table rows: " + name);
            }
        }

        public ModelTable build() { return new ModelTable(this); }
    }
}
