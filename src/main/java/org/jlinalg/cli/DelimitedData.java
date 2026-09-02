/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Small delimited-table reader used for phenotype and annotation metadata. */
final class DelimitedData {
    private final List<String> header;
    private final List<String[]> rows;

    private DelimitedData(List<String> header, List<String[]> rows) {
        this.header = List.copyOf(header);
        this.rows = List.copyOf(rows);
    }

    static DelimitedData read(Path path) throws IOException {
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        char delimiter = name.endsWith(".csv") ? ',' : '\t';
        try (BufferedReader input = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            String first = input.readLine();
            if (first == null) throw new IOException("table is empty: " + path);
            List<String> header = parse(first, delimiter, 1, path);
            Map<String, Boolean> unique = new LinkedHashMap<>();
            for (int index = 0; index < header.size(); index++) {
                String value = header.get(index).trim();
                if (value.isEmpty() || unique.put(value, Boolean.TRUE) != null)
                    throw new IOException(
                        "table column names must be unique and nonblank: " + path);
                header.set(index, value);
            }
            List<String[]> rows = new ArrayList<>();
            long lineNumber = 1;
            for (String line; (line = input.readLine()) != null;) {
                lineNumber++;
                if (line.isBlank()) continue;
                List<String> fields = parse(line, delimiter, lineNumber, path);
                if (fields.size() != header.size())
                    throw new IOException("expected " + header.size()
                        + " fields but found " + fields.size() + " at line "
                        + lineNumber + " in " + path);
                rows.add(fields.toArray(String[]::new));
            }
            if (rows.isEmpty()) throw new IOException(
                "table contains no data rows: " + path);
            return new DelimitedData(header, rows);
        }
    }

    List<String> header() { return header; }
    List<String[]> rows() { return rows; }

    int column(String name) {
        int index = header.indexOf(name);
        if (index < 0)
            throw new IllegalArgumentException("column is absent: " + name);
        return index;
    }

    private static List<String> parse(
            String line, char delimiter, long lineNumber, Path path)
            throws IOException {
        List<String> result = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean quoted = false;
        for (int index = 0; index < line.length(); index++) {
            char value = line.charAt(index);
            if (value == '"') {
                if (quoted && index + 1 < line.length()
                        && line.charAt(index + 1) == '"') {
                    current.append('"');
                    index++;
                } else {
                    quoted = !quoted;
                }
            } else if (value == delimiter && !quoted) {
                result.add(current.toString());
                current.setLength(0);
            } else {
                current.append(value);
            }
        }
        if (quoted)
            throw new IOException("unterminated quote at line " + lineNumber
                + " in " + path);
        result.add(current.toString());
        return result;
    }
}
