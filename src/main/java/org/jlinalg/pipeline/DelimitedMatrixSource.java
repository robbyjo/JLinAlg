/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;

/** Streaming {@code row_id,sample1,...} CSV/TSV numeric matrix. */
public final class DelimitedMatrixSource implements NumericMatrixSource {
    private final Path path;
    private final char delimiter;
    private final NumericMatrixMetadata metadata;

    public DelimitedMatrixSource(Path path, DataFormat format) throws IOException {
        if (format != DataFormat.CSV && format != DataFormat.TSV)
            throw new IllegalArgumentException("matrix source requires CSV or TSV");
        this.path = path.toAbsolutePath().normalize();
        delimiter = format == DataFormat.CSV ? ',' : '\t';
        Scan scan = scan();
        metadata = new NumericMatrixMetadata(
            this.path, scan.rows(), scan.samples(), format);
    }

    public static DelimitedMatrixSource open(Path path) throws IOException {
        return new DelimitedMatrixSource(path, DataFormat.infer(path));
    }

    @Override public NumericMatrixMetadata metadata() { return metadata; }

    @Override
    public NumericBlockReader open(int[] requestedOrder) throws IOException {
        int[] order = DelimitedVariantSource.normalizeOrder(requestedOrder,
            metadata.sampleIds().size());
        BufferedReader input = Files.newBufferedReader(path, StandardCharsets.UTF_8);
        if (input.readLine() == null) {
            input.close();
            throw new IOException("matrix is empty: " + path);
        }
        return new Reader(input, order);
    }

    private Scan scan() throws IOException {
        try (BufferedReader input =
                 Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String first = input.readLine();
            if (first == null) throw new IOException("matrix is empty: " + path);
            List<String> header = parse(first, 1);
            if (header.size() < 2)
                throw new IOException(
                    "matrix header requires a row ID and samples: " + path);
            List<String> samples = new ArrayList<>();
            HashSet<String> seenSamples = new HashSet<>();
            for (int index = 1; index < header.size(); index++) {
                String sample = header.get(index).trim();
                if (sample.isEmpty() || !seenSamples.add(sample))
                    throw new IOException(
                        "sample IDs must be unique and nonblank: " + path);
                samples.add(sample);
            }
            long rows = 0;
            for (String line; (line = input.readLine()) != null;) {
                if (line.isBlank()) continue;
                long lineNumber = rows + 2;
                List<String> fields = parse(line, lineNumber);
                if (fields.size() != header.size())
                    throw new IOException("expected " + header.size()
                        + " fields but found " + fields.size()
                        + " at line " + lineNumber + " in " + path);
                String id = fields.get(0).trim();
                if (id.isEmpty())
                    throw new IOException(
                        "row IDs must be nonblank at line "
                        + lineNumber + " in " + path);
                rows++;
            }
            if (rows == 0) throw new IOException("matrix has no rows: " + path);
            return new Scan(rows, List.copyOf(samples));
        }
    }

    private final class Reader implements NumericBlockReader {
        private final BufferedReader input;
        private final int[] order;
        private long sourceIndex;
        private long lineNumber = 1;
        private boolean closed;

        private Reader(BufferedReader input, int[] order) {
            this.input = input;
            this.order = order;
        }

        @Override
        public NumericBlock read(int maximumRows) throws IOException {
            if (maximumRows < 1)
                throw new IllegalArgumentException("maximum rows must be positive");
            if (closed) throw new IOException("matrix reader is closed");
            List<NumericRow> rows = new ArrayList<>(maximumRows);
            long first = sourceIndex;
            for (String line; rows.size() < maximumRows
                    && (line = input.readLine()) != null;) {
                lineNumber++;
                if (line.isBlank()) continue;
                List<String> fields = parse(line, lineNumber);
                int expected = metadata.sampleIds().size() + 1;
                if (fields.size() != expected)
                    throw new IOException("expected " + expected
                        + " fields but found " + fields.size()
                        + " at line " + lineNumber + " in " + path);
                double[] values = new double[order.length];
                for (int index = 0; index < order.length; index++) {
                    String token = fields.get(order[index] + 1).trim();
                    values[index] = missing(token) ? Double.NaN
                        : number(token, lineNumber, order[index]);
                }
                rows.add(new NumericRow(fields.get(0), values));
                sourceIndex++;
            }
            return rows.isEmpty() ? null : new NumericBlock(first, rows);
        }

        @Override public void close() throws IOException {
            if (!closed) {
                closed = true;
                input.close();
            }
        }
    }

    private List<String> parse(String line, long lineNumber) throws IOException {
        List<String> fields = new ArrayList<>();
        StringBuilder value = new StringBuilder();
        boolean quoted = false;
        for (int index = 0; index < line.length(); index++) {
            char current = line.charAt(index);
            if (current == '"') {
                if (quoted && index + 1 < line.length()
                        && line.charAt(index + 1) == '"') {
                    value.append('"');
                    index++;
                } else {
                    quoted = !quoted;
                }
            } else if (current == delimiter && !quoted) {
                fields.add(value.toString());
                value.setLength(0);
            } else {
                value.append(current);
            }
        }
        if (quoted)
            throw new IOException("unterminated quoted field at line "
                + lineNumber + " in " + path);
        fields.add(value.toString());
        return fields;
    }

    private double number(String token, long line, int sample) throws IOException {
        try {
            double value = Double.parseDouble(token);
            if (!Double.isFinite(value)) throw new NumberFormatException("not finite");
            return value;
        } catch (NumberFormatException exception) {
            throw new IOException("invalid value for sample column "
                + (sample + 1) + " at line " + line + " in " + path,
                exception);
        }
    }

    private static boolean missing(String token) {
        String value = token.toLowerCase(Locale.ROOT);
        return value.isEmpty() || value.equals(".") || value.equals("na")
            || value.equals("n/a") || value.equals("null") || value.equals("nan");
    }

    private record Scan(long rows, List<String> samples) { }
}
