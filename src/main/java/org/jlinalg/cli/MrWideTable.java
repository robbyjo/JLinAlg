/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Pattern;
import org.jlinalg.mr.MultivariateInstrument;

/** Strict wide-table reader shared by the explicitly distinct MVMR commands. */
final class MrWideTable {
    private MrWideTable() { }

    static List<MultivariateInstrument> read(Path path,
            List<String> exposures, List<String> outcomes) throws IOException {
        List<String> lines = Files.readAllLines(path);
        if (lines.size() < 2)
            throw new IOException("MR input needs a header and at least one row");
        String separator = path.toString().toLowerCase(Locale.ROOT).endsWith(".csv")
            ? "," : "\t";
        String[] header = lines.get(0).split(Pattern.quote(separator), -1);
        Map<String, Integer> columns = new HashMap<>();
        for (int index = 0; index < header.length; index++) {
            String key = normalize(header[index]);
            if (columns.putIfAbsent(key, index) != null)
                throw new IllegalArgumentException(
                    "duplicate normalized MR column: " + header[index]);
        }
        int variant = column(columns, "variant_id", "snp", "rsid");
        int[] bx = new int[exposures.size()], sx = new int[exposures.size()];
        for (int index = 0; index < exposures.size(); index++) {
            String name = normalize(exposures.get(index));
            bx[index] = named(columns, "beta_exposure", name,
                exposures.size() == 1, "exposure_effect");
            sx[index] = named(columns, "se_exposure", name,
                exposures.size() == 1, "exposure_se");
        }
        int[] by = new int[outcomes.size()], sy = new int[outcomes.size()];
        for (int index = 0; index < outcomes.size(); index++) {
            String name = normalize(outcomes.get(index));
            by[index] = named(columns, "beta_outcome", name,
                outcomes.size() == 1, "outcome_effect");
            sy[index] = named(columns, "se_outcome", name,
                outcomes.size() == 1, "outcome_se");
        }
        List<MultivariateInstrument> result = new ArrayList<>();
        for (int line = 1; line < lines.size(); line++) {
            if (lines.get(line).isBlank()) continue;
            String[] values = lines.get(line).split(Pattern.quote(separator), -1);
            if (values.length != header.length)
                throw new IOException("inconsistent MR input row " + (line + 1));
            double[] exposureEffects = values(values, bx);
            double[] exposureErrors = values(values, sx);
            double[] outcomeEffects = values(values, by);
            double[] outcomeErrors = values(values, sy);
            result.add(new MultivariateInstrument(values[variant].trim(),
                exposureEffects, exposureErrors, outcomeEffects, outcomeErrors));
        }
        return List.copyOf(result);
    }

    static double[][] readMatrix(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path);
        lines = lines.stream().filter(line -> !line.isBlank()).toList();
        if (lines.isEmpty()) throw new IOException("matrix is empty: " + path);
        String separator = path.toString().toLowerCase(Locale.ROOT).endsWith(".csv")
            ? "," : "\t";
        int size = lines.size();
        double[][] result = new double[size][size];
        for (int row = 0; row < size; row++) {
            String[] values = lines.get(row).split(Pattern.quote(separator), -1);
            if (values.length != size)
                throw new IOException("matrix must be square: " + path);
            for (int column = 0; column < size; column++)
                result[row][column] = number(values[column]);
        }
        return result;
    }

    static List<String> names(String value, String option) {
        if (value == null || value.isBlank())
            throw new IllegalArgumentException(option + " is required");
        List<String> result = Arrays.stream(value.split(",", -1))
            .map(String::trim).toList();
        if (result.stream().anyMatch(String::isBlank)
                || new java.util.HashSet<>(result).size() != result.size())
            throw new IllegalArgumentException(
                option + " must contain unique nonblank names");
        return result;
    }

    private static int named(Map<String, Integer> columns, String prefix,
            String name, boolean allowCanonical, String... canonicalAliases) {
        List<String> candidates = new ArrayList<>();
        candidates.add(prefix + "_" + name);
        candidates.add(prefix + "." + name);
        if (allowCanonical) {
            candidates.add(prefix);
            candidates.addAll(List.of(canonicalAliases));
        }
        return column(columns, candidates.toArray(String[]::new));
    }

    private static int column(Map<String, Integer> columns, String... names) {
        for (String name : names) {
            Integer index = columns.get(normalize(name));
            if (index != null) return index;
        }
        throw new IllegalArgumentException(
            "MR input is missing one of " + Arrays.toString(names));
    }

    private static double[] values(String[] row, int[] columns) {
        double[] result = new double[columns.length];
        for (int index = 0; index < columns.length; index++)
            result[index] = number(row[columns[index]]);
        return result;
    }

    private static double number(String value) {
        try {
            double result = Double.parseDouble(value.trim());
            if (!Double.isFinite(result))
                throw new NumberFormatException("non-finite");
            return result;
        } catch (NumberFormatException failure) {
            throw new IllegalArgumentException(
                "invalid numeric MR value: " + value, failure);
        }
    }

    private static String normalize(String value) {
        return value.trim().toLowerCase(Locale.ROOT)
            .replace('-', '_').replace('.', '_');
    }
}
