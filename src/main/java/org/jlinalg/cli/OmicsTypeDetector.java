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
import java.util.Locale;
import java.util.regex.Pattern;
import org.jlinalg.pipeline.DataFormat;

/** Conservative format-then-identifier omics schema inference. */
final class OmicsTypeDetector {
    private static final Pattern GWAS = Pattern.compile(
        "(?i)(rs\\d+|(?:chr)?[0-9XYM]+:\\d+(?::[^:]+:[^:]+)?)");
    private static final Pattern EWAS = Pattern.compile("(?i)cg\\d+");
    private static final Pattern EXPRESSION = Pattern.compile(
        "(?i)ENS(?:G|T|E)\\d+(?:\\.\\d+)?");

    private OmicsTypeDetector() { }

    static Detection detect(Path path, String requested) throws IOException {
        if (!requested.equals("auto"))
            return new Detection(requested, "explicit", 1.0);
        DataFormat format = DataFormat.infer(path);
        if (format == DataFormat.VCF || format == DataFormat.VCF_GZ
                || format == DataFormat.BCF || format == DataFormat.BGEN)
            return new Detection("gwas", "file-format", 1.0);
        char delimiter = format == DataFormat.CSV ? ',' : '\t';
        int gwas = 0;
        int ewas = 0;
        int expression = 0;
        int examined = 0;
        try (BufferedReader input = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            input.readLine();
            for (String line; examined < 1000
                    && (line = input.readLine()) != null;) {
                if (line.isBlank()) continue;
                String id = firstField(line, delimiter).trim();
                if (GWAS.matcher(id).matches()) gwas++;
                if (EWAS.matcher(id).matches()) ewas++;
                if (EXPRESSION.matcher(id).matches()) expression++;
                examined++;
            }
        }
        if (examined == 0)
            throw new IOException("omics input contains no data rows: " + path);
        int maximum = Math.max(gwas, Math.max(ewas, expression));
        double confidence = maximum / (double) examined;
        if (confidence < 0.8)
            return new Detection("generic", "identifier-pattern", confidence);
        String type = maximum == gwas ? "gwas"
            : maximum == ewas ? "ewas" : "expression";
        return new Detection(type, "identifier-pattern", confidence);
    }

    private static String firstField(String line, char delimiter) {
        boolean quoted = false;
        StringBuilder result = new StringBuilder();
        for (int index = 0; index < line.length(); index++) {
            char value = line.charAt(index);
            if (value == '"') {
                quoted = !quoted;
            } else if (value == delimiter && !quoted) {
                break;
            } else {
                result.append(value);
            }
        }
        return result.toString();
    }

    record Detection(String type, String source, double confidence) {
        Detection {
            type = type.toLowerCase(Locale.ROOT);
        }
    }
}
