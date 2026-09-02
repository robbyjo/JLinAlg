/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.OffsetDateTime;
import java.util.LinkedHashMap;
import java.util.Map;

/** Minimal stable JSON provenance manifest without an additional dependency. */
final class ManifestWriter {
    private final Map<String, String> values = new LinkedHashMap<>();

    ManifestWriter put(String name, Object value) {
        values.put(name, value == null ? null : value.toString());
        return this;
    }

    void write(Path path) throws IOException {
        Path absolute = path.toAbsolutePath().normalize();
        Path parent = absolute.getParent();
        if (parent != null) Files.createDirectories(parent);
        StringBuilder json = new StringBuilder("{\n");
        json.append("  \"schema_version\": 1,\n");
        json.append("  \"created\": \"")
            .append(escape(OffsetDateTime.now().toString())).append("\"");
        for (Map.Entry<String, String> entry : values.entrySet()) {
            json.append(",\n  \"").append(escape(entry.getKey()))
                .append("\": ");
            if (entry.getValue() == null) json.append("null");
            else json.append('\"').append(escape(entry.getValue())).append('\"');
        }
        json.append("\n}\n");
        Files.writeString(absolute, json, StandardCharsets.UTF_8);
    }

    private static String escape(String value) {
        return value.replace("\\", "\\\\").replace("\"", "\\\"")
            .replace("\r", "\\r").replace("\n", "\\n");
    }
}
