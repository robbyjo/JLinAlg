/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.coloc;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.jlinalg.susie.CredibleSet;
import org.jlinalg.susie.SusieResult;

/**
 * Per-signal, per-variant log Bayes factors used by SuSiE colocalization.
 * Values are stored row-major with one row per credible signal.
 */
public final class ColocSusieInput {
    private final List<String> variantNames;
    private final double[] logBayesFactors;
    private final int[] effectIndices;

    /** Creates an input whose source effect indices are {@code 0..L-1}. */
    public ColocSusieInput(
            List<String> variantNames, double[][] logBayesFactors) {
        this(variantNames, flatten(logBayesFactors),
            sequence(logBayesFactors == null ? 0 : logBayesFactors.length));
    }

    /**
     * Creates an input from row-major log Bayes factors and their source
     * zero-based SuSiE effect indices.
     */
    public ColocSusieInput(
            List<String> variantNames, double[] logBayesFactors,
            int[] effectIndices) {
        if (variantNames == null || variantNames.isEmpty()) {
            throw new IllegalArgumentException("variant names are required");
        }
        List<String> names = List.copyOf(variantNames);
        Set<String> unique = new HashSet<>(names.size());
        for (String name : names) {
            if (name == null || name.isBlank() || !unique.add(name)) {
                throw new IllegalArgumentException(
                    "variant names must be nonblank and unique");
            }
        }
        if (logBayesFactors == null || effectIndices == null
                || logBayesFactors.length != effectIndices.length * names.size()) {
            throw new IllegalArgumentException(
                "log Bayes factor dimensions are invalid");
        }
        for (int effect : effectIndices) {
            if (effect < 0) {
                throw new IllegalArgumentException(
                    "effect indices must be nonnegative");
            }
        }
        for (int signal = 0; signal < effectIndices.length; signal++) {
            boolean finite = false;
            int offset = signal * names.size();
            for (int variant = 0; variant < names.size(); variant++) {
                double value = logBayesFactors[offset + variant];
                if (Double.isNaN(value) || value == Double.POSITIVE_INFINITY) {
                    throw new IllegalArgumentException(
                        "log Bayes factors must be finite or negative infinity");
                }
                finite |= Double.isFinite(value);
            }
            if (!finite) {
                throw new IllegalArgumentException(
                    "each signal needs at least one finite log Bayes factor");
            }
        }
        this.variantNames = names;
        this.logBayesFactors = logBayesFactors.clone();
        this.effectIndices = effectIndices.clone();
    }

    /** Extracts the credible-effect rows required by {@code coloc.susie}. */
    public static ColocSusieInput from(SusieResult result) {
        if (result == null) {
            throw new IllegalArgumentException("SuSiE result is required");
        }
        List<CredibleSet> sets = result.credibleSets();
        List<String> names = result.variableNames();
        double[] all = result.logBayesFactors();
        double[] selected = new double[sets.size() * names.size()];
        int[] indices = new int[sets.size()];
        for (int signal = 0; signal < sets.size(); signal++) {
            int effect = sets.get(signal).effectIndex();
            indices[signal] = effect;
            System.arraycopy(all, effect * names.size(), selected,
                signal * names.size(), names.size());
        }
        return new ColocSusieInput(names, selected, indices);
    }

    public List<String> variantNames() { return variantNames; }
    public double[] logBayesFactors() { return logBayesFactors.clone(); }
    public int[] effectIndices() { return effectIndices.clone(); }
    public int signals() { return effectIndices.length; }
    public int variants() { return variantNames.size(); }

    double[] rawLogBayesFactors() { return logBayesFactors; }
    int[] rawEffectIndices() { return effectIndices; }

    private static double[] flatten(double[][] matrix) {
        if (matrix == null) {
            throw new IllegalArgumentException("log Bayes factors are required");
        }
        if (matrix.length == 0) return new double[0];
        if (matrix[0] == null || matrix[0].length == 0) {
            throw new IllegalArgumentException(
                "log Bayes factor matrix must have variants");
        }
        int columns = matrix[0].length;
        double[] result = new double[matrix.length * columns];
        for (int row = 0; row < matrix.length; row++) {
            if (matrix[row] == null || matrix[row].length != columns) {
                throw new IllegalArgumentException(
                    "log Bayes factor matrix must be rectangular");
            }
            System.arraycopy(matrix[row], 0, result, row * columns, columns);
        }
        return result;
    }

    private static int[] sequence(int length) {
        int[] result = new int[length];
        for (int index = 0; index < length; index++) result[index] = index;
        return result;
    }
}
