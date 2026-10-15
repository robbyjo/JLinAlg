/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import org.yaml.snakeyaml.LoaderOptions;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.SafeConstructor;

/** Bounded safe YAML reader for study-pipeline artifacts. */
final class StudyYaml {
    static final int CODE_POINT_LIMIT = 16_000_000;

    private StudyYaml() { }

    static Object load(Reader input) {
        return reader().load(input);
    }

    static Object load(String input) {
        return reader().load(input);
    }

    static Map<String, Object> loadMap(Path path) throws IOException {
        return loadMap(path, CODE_POINT_LIMIT);
    }

    static Map<String, Object> loadArtifact(Path path) throws IOException {
        if (Files.size(path) > 128L * 1024 * 1024)
            throw new IOException("Generated study artifact exceeds 128 MiB: " + path);
        return loadMap(path, 128 * 1024 * 1024);
    }

    private static Map<String, Object> loadMap(Path path, int codePointLimit) throws IOException {
        Object value;
        try (Reader input = Files.newBufferedReader(path,
                StandardCharsets.UTF_8)) {
            value = reader(codePointLimit).load(input);
        }
        if (!(value instanceof Map<?, ?> raw))
            throw new IOException("Expected YAML mapping: " + path);
        Map<String, Object> result = new LinkedHashMap<>();
        raw.forEach((key, item) -> result.put(String.valueOf(key), item));
        return result;
    }

    private static Yaml reader() {
        return reader(CODE_POINT_LIMIT);
    }

    private static Yaml reader(int codePointLimit) {
        LoaderOptions options = new LoaderOptions();
        options.setAllowDuplicateKeys(false);
        options.setMaxAliasesForCollections(0);
        options.setNestingDepthLimit(80);
        options.setCodePointLimit(codePointLimit);
        return new Yaml(new SafeConstructor(options));
    }
}
