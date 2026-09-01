/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/** Compiler for an intentionally small, fast subset of R fixed-effect formulas. */
public final class Formula {
    private Formula() { }

    /** Compiles treatment-coded fixed effects, interactions, and offsets. */
    public static CompiledFormula compile(String formula, ModelTable table) {
        return compile(formula, table, FormulaOptions.defaults());
    }

    /**
     * Supports {@code y ~ x + factor + x:factor + x*factor + offset(exposure)},
     * with {@code 0} or {@code -1} suppressing the intercept.
     */
    public static CompiledFormula compile(
            String formula, ModelTable table, FormulaOptions options) {
        if (formula == null || table == null || options == null) {
            throw new IllegalArgumentException("formula, table, and options are required");
        }
        String[] sides = formula.replaceAll("\\s+", "").split("~", -1);
        if (sides.length != 2 || sides[0].isEmpty() || sides[1].isEmpty()) {
            throw new IllegalArgumentException("formula must have the form response ~ terms");
        }
        double[] response = table.numeric(sides[0]).clone();
        boolean intercept = true;
        String offsetName = null;
        Set<String> terms = new LinkedHashSet<>();
        String rightHandSide = sides[1].replace("-1", "+-1");
        for (String token : rightHandSide.split("\\+")) {
            if (token.isEmpty()) continue;
            if (token.equals("1")) {
                intercept = true;
            } else if (token.equals("0") || token.equals("-1")) {
                intercept = false;
            } else if (token.startsWith("offset(") && token.endsWith(")")) {
                offsetName = token.substring(7, token.length() - 1);
            } else if (token.contains("*")) {
                String[] factors = token.split("\\*", -1);
                if (factors.length != 2 || factors[0].isEmpty() || factors[1].isEmpty()) {
                    throw new IllegalArgumentException("invalid interaction: " + token);
                }
                terms.add(factors[0]);
                terms.add(factors[1]);
                terms.add(factors[0] + ":" + factors[1]);
            } else {
                terms.add(token);
            }
        }

        List<Encoded> encoded = new ArrayList<>();
        if (intercept) {
            double[] values = new double[table.rows()];
            java.util.Arrays.fill(values, 1.0);
            encoded.add(new Encoded(List.of("(Intercept)"), values, 1));
        }
        for (String term : terms) {
            String[] factors = term.split(":", -1);
            Encoded value = encode(factors[0], table, options.contrastCoding());
            for (int index = 1; index < factors.length; index++) {
                value = interact(value,
                    encode(factors[index], table, options.contrastCoding()),
                    table.rows());
            }
            encoded.add(value);
        }

        int columns = encoded.stream().mapToInt(Encoded::columns).sum();
        if (columns == 0) {
            throw new IllegalArgumentException("formula produces no fixed-effect columns");
        }
        double[] design = new double[table.rows() * columns];
        List<String> names = new ArrayList<>(columns);
        int destinationColumn = 0;
        for (Encoded value : encoded) {
            names.addAll(value.names());
            for (int row = 0; row < table.rows(); row++) {
                System.arraycopy(value.values(), row * value.columns(),
                    design, row * columns + destinationColumn, value.columns());
            }
            destinationColumn += value.columns();
        }
        double[] weights = options.weightColumn() == null ? null
            : table.numeric(options.weightColumn()).clone();
        double[] offset = offsetName == null ? null
            : table.numeric(offsetName).clone();
        return new CompiledFormula(response, design, table.rows(), columns,
            names, weights, offset);
    }

    private static Encoded encode(
            String name, ModelTable table, ContrastCoding coding) {
        if (table.isNumeric(name)) {
            return new Encoded(List.of(name), table.numeric(name).clone(), 1);
        }
        if (!table.isCategorical(name)) {
            throw new IllegalArgumentException("column is absent: " + name);
        }
        String[] values = table.categorical(name);
        List<String> levels = new ArrayList<>();
        for (String value : values) {
            if (value != null && !levels.contains(value)) levels.add(value);
        }
        if (levels.size() < 2) {
            throw new IllegalArgumentException(
                "categorical term requires at least two observed levels: " + name);
        }
        int columns = levels.size() - 1;
        double[] matrix = new double[values.length * columns];
        List<String> names = new ArrayList<>(columns);
        for (int column = 0; column < columns; column++) {
            String level = coding == ContrastCoding.TREATMENT
                ? levels.get(column + 1) : levels.get(column);
            names.add(name + level);
        }
        for (int row = 0; row < values.length; row++) {
            String value = values[row];
            if (value == null) {
                java.util.Arrays.fill(matrix, row * columns,
                    (row + 1) * columns, Double.NaN);
                continue;
            }
            if (coding == ContrastCoding.TREATMENT) {
                int level = levels.indexOf(value) - 1;
                if (level >= 0) matrix[row * columns + level] = 1.0;
            } else {
                int level = levels.indexOf(value);
                if (level == levels.size() - 1) {
                    java.util.Arrays.fill(matrix, row * columns,
                        (row + 1) * columns, -1.0);
                } else {
                    matrix[row * columns + level] = 1.0;
                }
            }
        }
        return new Encoded(names, matrix, columns);
    }

    private static Encoded interact(Encoded left, Encoded right, int rows) {
        int columns = left.columns() * right.columns();
        double[] values = new double[rows * columns];
        List<String> names = new ArrayList<>(columns);
        for (String leftName : left.names()) {
            for (String rightName : right.names()) {
                names.add(leftName + ":" + rightName);
            }
        }
        for (int row = 0; row < rows; row++) {
            int output = row * columns;
            for (int leftColumn = 0; leftColumn < left.columns(); leftColumn++) {
                for (int rightColumn = 0;
                        rightColumn < right.columns(); rightColumn++) {
                    values[output++] = left.values()[row * left.columns() + leftColumn]
                        * right.values()[row * right.columns() + rightColumn];
                }
            }
        }
        return new Encoded(names, values, columns);
    }

    private record Encoded(List<String> names, double[] values, int columns) { }
}
