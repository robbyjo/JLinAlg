/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.formula.ModelTable;

/** In-memory phenotype table; large omics matrices remain block streamed. */
final class PhenotypeData {
    private static final Set<String> MISSING = Set.of(
        "", ".", "na", "n/a", "null", "nan");
    private final DelimitedData table;
    private final String idColumn;
    private final int idIndex;
    private final Map<String, Integer> rowById;

    private PhenotypeData(DelimitedData table, String idColumn) {
        this.table = table;
        this.idColumn = idColumn;
        idIndex = table.column(idColumn);
        rowById = new HashMap<>();
        for (int row = 0; row < table.rows().size(); row++) {
            String id = table.rows().get(row)[idIndex].trim();
            if (id.isEmpty() || rowById.put(id, row) != null)
                throw new IllegalArgumentException(
                    "phenotype IDs must be unique and nonblank: " + idColumn);
        }
    }

    static PhenotypeData read(Path path, String idColumn) throws IOException {
        return new PhenotypeData(DelimitedData.read(path), idColumn);
    }

    Prepared prepare(
            List<String> requestedIds, String response,
            boolean encodeBinomial, String requestedCase,
            String requestedControl) {
        List<String> ids = requestedIds == null
            ? originalIds() : List.copyOf(requestedIds);
        if (requestedIds != null) {
            if (ids.size() != rowById.size())
                throw new IllegalArgumentException(
                    "omics and phenotype ID sets differ; exact alignment is required");
            for (String id : ids)
                if (!rowById.containsKey(id))
                    throw new IllegalArgumentException(
                        "omics sample is absent from phenotype table: " + id);
        }
        int[] order = new int[ids.size()];
        for (int index = 0; index < ids.size(); index++)
            order[index] = rowById.get(ids.get(index));

        BinaryMapping binary = binaryMapping(
            response, order, requestedCase, requestedControl);
        ModelTable.Builder builder = ModelTable.builder(order.length);
        for (int column = 0; column < table.header().size(); column++) {
            String name = table.header().get(column);
            if (column == idIndex) {
                String[] values = strings(column, order);
                builder.categorical(name, values);
                continue;
            }
            if (encodeBinomial && name.equals(response)) {
                if (binary == null)
                    throw new IllegalArgumentException(
                        "binomial response must contain two recognized levels; "
                            + "use --case-value and --control-value");
                double[] values = new double[order.length];
                for (int index = 0; index < order.length; index++) {
                    String raw = normalized(table.rows().get(order[index])[column]);
                    values[index] = missing(raw) ? Double.NaN
                        : raw.equals(binary.caseValue()) ? 1.0 : 0.0;
                }
                builder.numeric(name, values);
            } else if (numeric(column, order)) {
                double[] values = new double[order.length];
                for (int index = 0; index < order.length; index++) {
                    String raw = normalized(table.rows().get(order[index])[column]);
                    values[index] = missing(raw) ? Double.NaN
                        : Double.parseDouble(raw);
                }
                builder.numeric(name, values);
            } else {
                builder.categorical(name, strings(column, order));
            }
        }
        int[] groups = binary == null ? null : groups(response, order, binary);
        return new Prepared(ids, builder.build(), binary, groups);
    }

    List<String> originalIds() {
        List<String> result = new ArrayList<>(table.rows().size());
        for (String[] row : table.rows()) result.add(row[idIndex].trim());
        return result;
    }

    List<String> alignedValues(List<String> ids, String columnName) {
        int column = table.column(columnName);
        List<String> result = new ArrayList<>(ids.size());
        for (String id : ids) {
            Integer row = rowById.get(id);
            if (row == null)
                throw new IllegalArgumentException(
                    "phenotype ID is absent while aligning " + columnName
                        + ": " + id);
            String value = table.rows().get(row)[column].trim();
            if (missing(normalized(value)))
                throw new IllegalArgumentException(
                    "GRM matching column contains a missing value: "
                        + columnName + " for " + id);
            result.add(value);
        }
        return List.copyOf(result);
    }

    private BinaryMapping binaryMapping(
            String response, int[] order, String requestedCase,
            String requestedControl) {
        if (response == null || !table.header().contains(response)) return null;
        int column = table.column(response);
        LinkedHashSet<String> levels = new LinkedHashSet<>();
        for (int row : order) {
            String value = normalized(table.rows().get(row)[column]);
            if (!missing(value)) levels.add(value);
        }
        if (levels.size() != 2) return null;
        if (requestedCase != null) {
            String positive = normalized(requestedCase);
            String negative = requestedControl == null
                ? levels.stream().filter(value -> !value.equals(positive))
                    .findFirst().orElseThrow()
                : normalized(requestedControl);
            if (!levels.contains(positive) || !levels.contains(negative)
                    || positive.equals(negative))
                throw new IllegalArgumentException(
                    "case/control values do not match the response levels");
            return new BinaryMapping(positive, negative);
        }
        Map<String, String> pairs = Map.ofEntries(
            Map.entry("1", "0"), Map.entry("yes", "no"),
            Map.entry("y", "n"), Map.entry("true", "false"),
            Map.entry("case", "control"),
            Map.entry("affected", "unaffected"));
        for (Map.Entry<String, String> pair : pairs.entrySet())
            if (levels.contains(pair.getKey()) && levels.contains(pair.getValue()))
                return new BinaryMapping(pair.getKey(), pair.getValue());
        return null;
    }

    private int[] groups(String response, int[] order, BinaryMapping mapping) {
        int column = table.column(response);
        int[] groups = new int[order.length];
        java.util.Arrays.fill(groups, -1);
        for (int index = 0; index < order.length; index++) {
            String value = normalized(table.rows().get(order[index])[column]);
            if (value.equals(mapping.caseValue())) groups[index] = 1;
            else if (value.equals(mapping.controlValue())) groups[index] = 0;
        }
        return groups;
    }

    private boolean numeric(int column, int[] order) {
        boolean observed = false;
        for (int row : order) {
            String value = normalized(table.rows().get(row)[column]);
            if (missing(value)) continue;
            observed = true;
            try {
                double parsed = Double.parseDouble(value);
                if (!Double.isFinite(parsed)) return false;
            } catch (NumberFormatException exception) {
                return false;
            }
        }
        return observed;
    }

    private String[] strings(int column, int[] order) {
        String[] result = new String[order.length];
        for (int index = 0; index < order.length; index++) {
            String value = table.rows().get(order[index])[column].trim();
            result[index] = missing(normalized(value)) ? null : value;
        }
        return result;
    }

    private static boolean missing(String value) {
        return MISSING.contains(value.toLowerCase(Locale.ROOT));
    }
    private static String normalized(String value) {
        return value.trim().toLowerCase(Locale.ROOT);
    }

    record BinaryMapping(String caseValue, String controlValue) { }
    record Prepared(
        List<String> ids, ModelTable modelTable,
        BinaryMapping binaryMapping, int[] caseControlGroups) {
        Prepared {
            ids = List.copyOf(ids);
            caseControlGroups = caseControlGroups == null ? null
                : caseControlGroups.clone();
        }
        @Override public int[] caseControlGroups() {
            return caseControlGroups == null ? null : caseControlGroups.clone();
        }
    }
}
