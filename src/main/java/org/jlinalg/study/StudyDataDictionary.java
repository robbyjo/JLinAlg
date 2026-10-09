/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

/** Local, explicitly typed dictionary mappings; prose never changes an analysis. */
final class StudyDataDictionary {
    private StudyDataDictionary() { }

    static void apply(Path root, Map<String, Object> dataset) {
        Map<String, Object> files = StudyIO.optionalMap(dataset.get("files"), "files");
        if (!files.containsKey("data_dictionary")) return;
        Path path = root.resolve(String.valueOf(files.get("data_dictionary"))).normalize();
        try {
            Map<String, Object> dictionary = StudyIO.load(path);
            Map<String, Object> columns = StudyIO.optionalMap(dictionary.get("columns"), "dictionary.columns");
            Map<String, Object> explicit = StudyIO.optionalMap(dataset.get("columns"), "columns");
            Map<String, Object> previous = StudyIO.optionalMap(StudyIO.optionalMap(dataset.get("metadata"), "metadata")
                .get("resolved_data_dictionary"), "resolved_data_dictionary");
            if (!previous.isEmpty()) {
                Map<String, Object> oldExplicit = StudyIO.optionalMap(previous.get("explicit_columns"), "explicit_columns");
                Map<String, Object> oldColumns = StudyIO.optionalMap(StudyIO.optionalMap(previous.get("content"), "content")
                    .get("columns"), "columns");
                oldColumns.forEach((key, value) -> {
                    if (!oldExplicit.containsKey(key) && value.equals(explicit.get(key))) explicit.remove(key);
                });
            }
            Map<String, Object> overrides = new LinkedHashMap<>();
            explicit.forEach((key, value) -> {
                if (columns.containsKey(key) && !columns.get(key).equals(value))
                    overrides.put(key, Map.of("dictionary", columns.get(key), "study", value));
            });
            columns.putAll(explicit);
            dataset.put("columns", columns);
            Map<String, Object> metadata = StudyIO.optionalMap(dataset.get("metadata"), "metadata");
            metadata.put("resolved_data_dictionary", Map.of(
                "path", files.get("data_dictionary"), "sha256", StudyIO.hash(path),
                "content", dictionary, "explicit_overrides", overrides, "explicit_columns", explicit));
            dataset.put("metadata", metadata);
        } catch (IOException failure) {
            throw new IllegalArgumentException("Cannot read data dictionary: " + path, failure);
        }
    }

    static void verify(Path root, Map<String, Object> study) throws IOException {
        for (Object item : StudyIO.list(study.get("datasets"), "datasets")) {
            Map<String, Object> dataset = StudyIO.map(item, "dataset");
            Map<String, Object> metadata = StudyIO.optionalMap(dataset.get("metadata"), "metadata");
            Map<String, Object> dictionary = StudyIO.optionalMap(
                metadata.get("resolved_data_dictionary"), "resolved_data_dictionary");
            if (!dictionary.isEmpty() && !dictionary.get("sha256").equals(
                    StudyIO.hash(root.resolve(String.valueOf(dictionary.get("path"))))))
                throw new IOException("Data dictionary changed after plan approval for dataset "
                    + dataset.get("id") + "; use a versioned dictionary and an approved amendment");
        }
    }
}
