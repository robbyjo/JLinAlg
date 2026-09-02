/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.CompiledFormula;
import org.jlinalg.formula.Formula;
import org.jlinalg.formula.ModelTable;
import org.jlinalg.glm.GlmFamily;

/** Formula adapter for fitting GEE from a reusable compiled model matrix. */
public final class GeeFormula {
    private GeeFormula() { }

    /**
     * Fits a formula containing {@code cluster(id)} and optional
     * {@code wave(visit)} declarations on its right-hand side.
     */
    public static GeeResult fit(
            String formula,
            ModelTable table,
            GlmFamily family,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (formula == null || table == null) {
            throw new IllegalArgumentException("formula and table are required");
        }
        String compact = formula.replaceAll("\\s+", "");
        String clusterColumn = specialColumn(compact, "cluster", true);
        String waveColumn = specialColumn(compact, "wave", false);
        String meanFormula = removeSpecial(compact, "cluster", clusterColumn);
        if (waveColumn != null) {
            meanFormula = removeSpecial(meanFormula, "wave", waveColumn);
        }
        CompiledFormula compiled = Formula.compile(meanFormula, table);
        return fit(compiled, integerColumn(table, clusterColumn),
            waveColumn == null ? null : integerColumn(table, waveColumn),
            family, options, backendPolicy);
    }

    /** Fits a GEE without reparsing the supplied compiled formula. */
    public static GeeResult fit(
            CompiledFormula formula,
            int[] cluster,
            int[] repeated,
            GlmFamily family,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (formula == null) {
            throw new IllegalArgumentException("formula is required");
        }
        return Gee.fit(formula.response(), formula.design(),
            formula.rows(), formula.columns(), cluster, repeated,
            family, formula.weights(), formula.offset(), options, backendPolicy);
    }

    private static String specialColumn(
            String formula, String function, boolean required) {
        String marker = function + "(";
        int start = formula.indexOf(marker);
        if (start < 0) {
            if (required) {
                throw new IllegalArgumentException(
                    "GEE formula requires " + function + "(column)");
            }
            return null;
        }
        if (formula.indexOf(marker, start + marker.length()) >= 0) {
            throw new IllegalArgumentException(
                "GEE formula may contain only one " + function + " declaration");
        }
        int end = formula.indexOf(')', start + marker.length());
        if (end < 0) {
            throw new IllegalArgumentException("unterminated " + function + " declaration");
        }
        String column = formula.substring(start + marker.length(), end);
        if (column.isEmpty() || column.indexOf('(') >= 0) {
            throw new IllegalArgumentException("invalid " + function + " column");
        }
        return column;
    }

    private static String removeSpecial(
            String formula, String function, String column) {
        String result = formula.replace(function + "(" + column + ")", "");
        while (result.contains("++")) result = result.replace("++", "+");
        result = result.replace("~+", "~");
        if (result.endsWith("+")) result = result.substring(0, result.length() - 1);
        return result;
    }

    private static int[] integerColumn(ModelTable table, String name) {
        double[] values = table.numericColumn(name);
        int[] result = new int[values.length];
        for (int row = 0; row < values.length; row++) {
            if (!Double.isFinite(values[row]) || values[row] != Math.rint(values[row])
                    || values[row] < Integer.MIN_VALUE
                    || values[row] > Integer.MAX_VALUE) {
                throw new IllegalArgumentException(
                    "cluster and wave columns must contain finite integers: " + name);
            }
            result[row] = (int) values[row];
        }
        return result;
    }
}
