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

/**
 * Streaming CSV/TSV additive-dosage matrix.
 *
 * <p>The compact layout is {@code id,sample1,...}. A genomic layout may use
 * {@code id,chromosome,position,ref,alt,sample1,...}. Quoted fields and doubled
 * quote escapes are supported. Missing dosages are blank, {@code .},
 * {@code NA}, {@code N/A}, {@code null}, or {@code NaN}.</p>
 */
public final class DelimitedVariantSource implements VariantSource {
    private final Path path;
    private final char delimiter;
    private final VariantSourceMetadata metadata;
    private final boolean genomicLayout;
    private final int sampleStart;

    public DelimitedVariantSource(Path path, DataFormat format) throws IOException {
        if (format != DataFormat.CSV && format != DataFormat.TSV)
            throw new IllegalArgumentException("delimited source requires CSV or TSV");
        this.path = path.toAbsolutePath().normalize();
        delimiter = format == DataFormat.CSV ? ',' : '\t';
        Scan scan = scan();
        genomicLayout = scan.genomicLayout();
        sampleStart = scan.sampleStart();
        metadata = new VariantSourceMetadata(this.path, scan.rows(),
            scan.samples(), format);
    }

    public static DelimitedVariantSource open(Path path) throws IOException {
        return new DelimitedVariantSource(path, DataFormat.infer(path));
    }

    @Override public VariantSourceMetadata metadata() { return metadata; }

    @Override
    public VariantBlockReader open(int[] requestedOrder) throws IOException {
        int[] order = normalizeOrder(requestedOrder, metadata.sampleIds().size());
        BufferedReader input = Files.newBufferedReader(path, StandardCharsets.UTF_8);
        String header = input.readLine();
        if (header == null) {
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
            boolean genomic = isGenomicHeader(header);
            int firstSample = genomic ? 5 : 1;
            if (header.size() <= firstSample)
                throw new IOException(
                    "matrix header must contain at least one sample: " + path);
            List<String> samples = new ArrayList<>();
            HashSet<String> seenSamples = new HashSet<>();
            for (int index = firstSample; index < header.size(); index++) {
                String sample = header.get(index).trim();
                if (sample.isEmpty() || !seenSamples.add(sample))
                    throw new IOException(
                        "sample IDs must be unique and nonblank: " + path);
                samples.add(sample);
            }
            HashSet<String> seenVariants = new HashSet<>();
            long rows = 0;
            for (String line; (line = input.readLine()) != null;) {
                long lineNumber = rows + 2;
                if (line.isBlank()) continue;
                List<String> fields = parse(line, lineNumber);
                if (fields.size() != header.size())
                    throw new IOException("expected " + header.size()
                        + " fields but found " + fields.size()
                        + " at line " + lineNumber + " in " + path);
                String id = fields.get(0).trim();
                if (id.isEmpty() || !seenVariants.add(id))
                    throw new IOException(
                        "variant IDs must be unique and nonblank at line "
                        + lineNumber + " in " + path);
                if (genomic) parsePosition(fields.get(2), lineNumber);
                rows++;
            }
            if (rows == 0) throw new IOException("matrix has no variants: " + path);
            return new Scan(rows, List.copyOf(samples), genomic, firstSample);
        }
    }

    private final class Reader implements VariantBlockReader {
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
        public VariantBlock read(int maximumVariants) throws IOException {
            if (maximumVariants < 1)
                throw new IllegalArgumentException("maximum variants must be positive");
            if (closed) throw new IOException("variant reader is closed");
            List<VariantRecord> result = new ArrayList<>(maximumVariants);
            long first = sourceIndex;
            for (String line; result.size() < maximumVariants
                    && (line = input.readLine()) != null;) {
                lineNumber++;
                if (line.isBlank()) continue;
                List<String> fields = parse(line, lineNumber);
                int expected = sampleStart + metadata.sampleIds().size();
                if (fields.size() != expected)
                    throw new IOException("expected " + expected
                        + " fields but found " + fields.size()
                        + " at line " + lineNumber + " in " + path);
                double[] dosages = new double[order.length];
                for (int index = 0; index < order.length; index++) {
                    String token = fields.get(sampleStart + order[index]).trim();
                    dosages[index] = missing(token)
                        ? Double.NaN : number(token, lineNumber, order[index]);
                }
                String id = fields.get(0).trim();
                Coordinates coordinates = genomicLayout
                    ? new Coordinates(fields.get(1).trim(),
                        parsePosition(fields.get(2), lineNumber),
                        fields.get(3).trim(), fields.get(4).trim())
                    : coordinates(id);
                result.add(new VariantRecord(id, coordinates.chromosome(),
                    coordinates.position(), coordinates.reference(),
                    coordinates.alternate(), dosages, Double.NaN));
                sourceIndex++;
            }
            return result.isEmpty() ? null : new VariantBlock(first, result);
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

    private static boolean isGenomicHeader(List<String> header) {
        if (header.size() < 6) return false;
        String chromosome = normalized(header.get(1));
        String position = normalized(header.get(2));
        return (chromosome.equals("chrom") || chromosome.equals("chromosome")
            || chromosome.equals("chr"))
            && (position.equals("position") || position.equals("pos"))
            && normalized(header.get(3)).equals("ref")
            && normalized(header.get(4)).equals("alt");
    }

    private static String normalized(String value) {
        return value.trim().toLowerCase(Locale.ROOT)
            .replace("_", "").replace("-", "");
    }

    private static Coordinates coordinates(String id) {
        String[] pieces = id.split(":", -1);
        if (pieces.length < 4) return new Coordinates("", 0, "", "");
        try {
            return new Coordinates(pieces[0], Long.parseLong(pieces[1]),
                pieces[2], pieces[3]);
        } catch (NumberFormatException exception) {
            return new Coordinates("", 0, "", "");
        }
    }

    private static long parsePosition(String token, long line) throws IOException {
        try {
            long value = Long.parseLong(token.trim());
            if (value < 1) throw new NumberFormatException("not positive");
            return value;
        } catch (NumberFormatException exception) {
            throw new IOException("invalid genomic position at line " + line,
                exception);
        }
    }

    private double number(String token, long line, int sample) throws IOException {
        try {
            double value = Double.parseDouble(token);
            if (!Double.isFinite(value)) throw new NumberFormatException("not finite");
            return value;
        } catch (NumberFormatException exception) {
            throw new IOException("invalid dosage for sample column "
                + (sample + 1) + " at line " + line + " in " + path,
                exception);
        }
    }

    private static boolean missing(String token) {
        String value = token.toLowerCase(Locale.ROOT);
        return value.isEmpty() || value.equals(".") || value.equals("na")
            || value.equals("n/a") || value.equals("null") || value.equals("nan");
    }

    static int[] normalizeOrder(int[] requested, int samples) {
        if (requested == null) {
            int[] identity = new int[samples];
            for (int index = 0; index < samples; index++) identity[index] = index;
            return identity;
        }
        if (requested.length == 0)
            throw new IllegalArgumentException("sample order cannot be empty");
        boolean[] seen = new boolean[samples];
        for (int value : requested) {
            if (value < 0 || value >= samples || seen[value])
                throw new IllegalArgumentException(
                    "sample order contains an invalid or duplicate index");
            seen[value] = true;
        }
        return requested.clone();
    }

    private record Scan(long rows, List<String> samples,
        boolean genomicLayout, int sampleStart) { }
    private record Coordinates(String chromosome, long position,
        String reference, String alternate) { }
}
