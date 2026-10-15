/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.yaml.snakeyaml.DumperOptions;
import org.yaml.snakeyaml.LoaderOptions;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.SafeConstructor;

/** Deterministic, safe serialization primitives for study documents. */
final class StudyIO {
    private StudyIO() { }

    static Yaml reader() {
        return reader(16_000_000);
    }

    private static Yaml reader(int codePointLimit) {
        LoaderOptions options = new LoaderOptions();
        options.setAllowDuplicateKeys(false);
        options.setMaxAliasesForCollections(0);
        options.setNestingDepthLimit(80);
        options.setCodePointLimit(codePointLimit);
        return new Yaml(new SafeConstructor(options));
    }

    static Yaml writer() {
        DumperOptions options = new DumperOptions();
        options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK);
        options.setPrettyFlow(true);
        options.setIndent(2);
        options.setIndicatorIndent(0);
        options.setWidth(100);
        options.setSplitLines(true);
        return new Yaml(options);
    }

    static Map<String, Object> load(Path path) throws IOException {
        return load(path, 16_000_000);
    }

    /** Generated, ledger/hash-checked evidence can exceed a study contract's budget.
     * Retains the same alias, duplicate-key, depth and safe-constructor protections. */
    static Map<String, Object> loadArtifact(Path path) throws IOException {
        if (Files.size(path) > 128L * 1024 * 1024)
            throw new IOException("Generated study artifact exceeds 128 MiB: " + path);
        return load(path, 128 * 1024 * 1024);
    }

    private static Map<String, Object> load(Path path, int codePointLimit) throws IOException {
        if (!Files.isRegularFile(path))
            throw new IOException("Study document is absent: " + path);
        Object value;
        try (Reader input = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            value = reader(codePointLimit).load(input);
        } catch (RuntimeException failure) {
            throw new IllegalArgumentException(
                "Invalid YAML in " + path + ": " + failure.getMessage(), failure);
        }
        return map(value, path.toString());
    }

    static Map<String, Object> map(Object value, String name) {
        if (!(value instanceof Map<?, ?> raw))
            throw new IllegalArgumentException(name + " must be a mapping");
        Map<String, Object> result = new LinkedHashMap<>();
        for (Map.Entry<?, ?> entry : raw.entrySet()) {
            if (!(entry.getKey() instanceof String key) || key.isBlank())
                throw new IllegalArgumentException(name + " requires nonblank string keys");
            result.put(key, normalize(entry.getValue(), name + "." + key));
        }
        return result;
    }

    static Map<String, Object> optionalMap(Object value, String name) {
        return value == null ? new LinkedHashMap<>() : map(value, name);
    }

    static List<Object> list(Object value, String name) {
        if (!(value instanceof Collection<?> raw))
            throw new IllegalArgumentException(name + " must be a list");
        List<Object> result = new ArrayList<>(raw.size());
        int index = 0;
        for (Object item : raw) result.add(normalize(item, name + "[" + index++ + "]"));
        return result;
    }

    static List<Object> optionalList(Object value, String name) {
        return value == null ? new ArrayList<>() : list(value, name);
    }

    private static Object normalize(Object value, String name) {
        if (value == null || value instanceof String || value instanceof Boolean
                || value instanceof Number) return value;
        if (value instanceof Map<?, ?>) return map(value, name);
        if (value instanceof Collection<?>) return list(value, name);
        // SnakeYAML may materialize timestamp-looking scalars as dates. Treat
        // them as their stable textual form rather than allowing Java tags.
        if (value instanceof java.util.Date date) return date.toInstant().toString();
        throw new IllegalArgumentException(name + " has unsupported value type "
            + value.getClass().getSimpleName());
    }

    @SuppressWarnings("unchecked")
    static Map<String, Object> copyMap(Map<String, Object> source) {
        return (Map<String, Object>) copy(source);
    }

    private static Object copy(Object value) {
        if (value instanceof Map<?, ?> raw) {
            Map<String, Object> result = new LinkedHashMap<>();
            raw.forEach((key, item) -> result.put(String.valueOf(key), copy(item)));
            return result;
        }
        if (value instanceof Collection<?> raw) {
            List<Object> result = new ArrayList<>(raw.size());
            raw.forEach(item -> result.add(copy(item)));
            return result;
        }
        return value;
    }

    static String canonicalJson(Object value) {
        StringBuilder output = new StringBuilder();
        json(value, output);
        return output.toString();
    }

    private static void json(Object value, StringBuilder output) {
        if (value == null) {
            output.append("null");
        } else if (value instanceof Boolean || value instanceof Number) {
            if (value instanceof Double && !Double.isFinite((Double) value)
                    || value instanceof Float && !Float.isFinite((Float) value))
                throw new IllegalArgumentException("Canonical documents cannot contain nonfinite numbers");
            output.append(value);
        } else if (value instanceof Map<?, ?> raw) {
            output.append('{');
            boolean first = true;
            TreeMap<String, Object> sorted = new TreeMap<>();
            raw.forEach((key, item) -> sorted.put(String.valueOf(key), item));
            for (Map.Entry<String, Object> entry : sorted.entrySet()) {
                if (!first) output.append(',');
                first = false;
                string(entry.getKey(), output);
                output.append(':');
                json(entry.getValue(), output);
            }
            output.append('}');
        } else if (value instanceof Collection<?> items) {
            output.append('[');
            boolean first = true;
            for (Object item : items) {
                if (!first) output.append(',');
                first = false;
                json(item, output);
            }
            output.append(']');
        } else {
            string(String.valueOf(value), output);
        }
    }

    private static void string(String value, StringBuilder output) {
        output.append('"');
        for (int i = 0; i < value.length(); i++) {
            char character = value.charAt(i);
            switch (character) {
                case '"' -> output.append('\\').append('"');
                case '\\' -> output.append('\\').append('\\');
                case '\b' -> output.append("\\b");
                case '\f' -> output.append("\\f");
                case '\n' -> output.append("\\n");
                case '\r' -> output.append("\\r");
                case '\t' -> output.append("\\t");
                default -> {
                    if (character < 0x20)
                        output.append(String.format("\\u%04x", (int) character));
                    else output.append(character);
                }
            }
        }
        output.append('"');
    }

    static String hash(Object value) {
        return hashBytes(canonicalJson(value).getBytes(StandardCharsets.UTF_8));
    }

    static String hash(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (var input = Files.newInputStream(path)) {
                byte[] buffer = new byte[1024 * 1024];
                int count;
                while ((count = input.read(buffer)) >= 0)
                    digest.update(buffer, 0, count);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException impossible) {
            throw new AssertionError(impossible);
        }
    }

    private static String hashBytes(byte[] bytes) {
        try {
            return HexFormat.of().formatHex(
                MessageDigest.getInstance("SHA-256").digest(bytes));
        } catch (NoSuchAlgorithmException impossible) {
            throw new AssertionError(impossible);
        }
    }

    static String yaml(Map<String, Object> value) {
        return writer().dump(sort(value));
    }

    private static Object sort(Object value) {
        if (value instanceof Map<?, ?> raw) {
            Map<String, Object> result = new LinkedHashMap<>();
            raw.entrySet().stream()
                .sorted(Comparator.comparing(entry -> String.valueOf(entry.getKey())))
                .forEach(entry -> result.put(String.valueOf(entry.getKey()),
                    sort(entry.getValue())));
            return result;
        }
        if (value instanceof Collection<?> raw) {
            List<Object> result = new ArrayList<>(raw.size());
            raw.forEach(item -> result.add(sort(item)));
            return result;
        }
        return value;
    }

    static void writeYaml(Path path, Map<String, Object> value) throws IOException {
        atomicWrite(path, yaml(value));
    }

    static void atomicWrite(Path path, String content) throws IOException {
        Path absolute = path.toAbsolutePath().normalize();
        Path parent = absolute.getParent();
        if (parent != null) Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, "." + absolute.getFileName(), ".tmp");
        try {
            Files.writeString(temporary, content, StandardCharsets.UTF_8);
            try {
                Files.move(temporary, absolute, StandardCopyOption.ATOMIC_MOVE,
                    StandardCopyOption.REPLACE_EXISTING);
            } catch (AtomicMoveNotSupportedException unsupported) {
                Files.move(temporary, absolute, StandardCopyOption.REPLACE_EXISTING);
            }
        } finally {
            Files.deleteIfExists(temporary);
        }
    }

    static String requireString(Map<String, Object> map, String key, String path) {
        Object value = map.get(key);
        if (!(value instanceof String text) || text.isBlank())
            throw new IllegalArgumentException(path + "." + key
                + " requires a nonblank string");
        return text;
    }

    static String string(Map<String, Object> map, String key, String fallback) {
        Object value = map.get(key);
        if (value == null) return fallback;
        if (!(value instanceof String text) || text.isBlank())
            throw new IllegalArgumentException(key + " requires a nonblank string");
        return text;
    }
}
