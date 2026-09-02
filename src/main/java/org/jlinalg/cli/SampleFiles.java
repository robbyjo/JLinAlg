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
import java.util.List;
import java.util.Locale;

/** BGEN companion and one-column sample-ID file parsing. */
final class SampleFiles {
    private SampleFiles() { }

    static List<String> read(Path path) throws IOException {
        List<String> result = new ArrayList<>();
        try (BufferedReader input = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            String first = input.readLine();
            if (first == null) throw new IOException(
                "sample file is empty: " + path);
            boolean oxford = first.trim().toLowerCase(Locale.ROOT)
                .startsWith("id_1");
            if (!oxford) add(result, first);
            if (oxford) input.readLine();
            for (String line; (line = input.readLine()) != null;) add(result, line);
        }
        if (result.isEmpty())
            throw new IOException("sample file contains no IDs: " + path);
        return List.copyOf(result);
    }

    private static void add(List<String> destination, String line) {
        if (line.isBlank()) return;
        String[] fields = line.trim().split("[\\t ,]+", -1);
        destination.add(fields.length > 1 ? fields[1] : fields[0]);
    }
}
