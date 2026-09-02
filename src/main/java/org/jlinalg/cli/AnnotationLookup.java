/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Selected annotation columns keyed without expanding result row counts. */
final class AnnotationLookup {
    private final List<String> columns;
    private final Map<String, String[]> values;

    private AnnotationLookup(
            List<String> columns, Map<String, String[]> values) {
        this.columns = List.copyOf(columns);
        this.values = Map.copyOf(values);
    }

    static AnnotationLookup empty() {
        return new AnnotationLookup(List.of(), Map.of());
    }

    static AnnotationLookup read(
            Path path, String requestedId, List<String> requestedColumns)
            throws IOException {
        DelimitedData table = DelimitedData.read(path);
        String id = requestedId == null ? table.header().get(0) : requestedId;
        int idIndex = table.column(id);
        List<String> columns;
        if (requestedColumns.size() == 1
                && requestedColumns.get(0).equalsIgnoreCase("all")) {
            columns = table.header().stream()
                .filter(column -> !column.equals(id)).toList();
        } else {
            columns = requestedColumns.isEmpty()
                ? List.of() : List.copyOf(requestedColumns);
        }
        int[] indices = new int[columns.size()];
        for (int index = 0; index < columns.size(); index++)
            indices[index] = table.column(columns.get(index));
        Map<String, String[]> values = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            String key = row[idIndex].trim();
            if (key.isEmpty() || values.containsKey(key))
                throw new IllegalArgumentException(
                    "annotation IDs must be unique and nonblank: " + id);
            String[] selected = new String[indices.length];
            for (int index = 0; index < indices.length; index++)
                selected[index] = row[indices[index]];
            values.put(key, selected);
        }
        return new AnnotationLookup(columns, values);
    }

    List<String> columns() { return columns; }
    String[] values(String id) {
        String[] result = values.get(id);
        return result == null ? new String[columns.size()] : result.clone();
    }
    long size() { return values.size(); }
}
