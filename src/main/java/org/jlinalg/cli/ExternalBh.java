/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

/** Disk-backed Benjamini-Hochberg adjustment and final result assembly. */
final class ExternalBh implements AutoCloseable {
    private final Path output;
    private final Path raw;
    private final Path work;
    private final Path qValues;
    private final BufferedWriter rawWriter;
    private final DataOutputStream initialQ;
    private final List<PValue> buffered = new ArrayList<>();
    private final List<Path> chunks = new ArrayList<>();
    private final int chunkCapacity;
    private long rows;
    private long tests;
    private boolean finished;
    private boolean closed;
    private boolean headerWritten;
    private final boolean overwrite;

    ExternalBh(Path output, boolean overwrite) throws IOException {
        this(output, overwrite, (int) Math.max(10_000,
            Math.min(500_000, Runtime.getRuntime().maxMemory() / 512)));
    }

    /** Explicit run size for deterministic bounded-memory tests/benchmarks. */
    ExternalBh(Path output, boolean overwrite, int chunkCapacity) throws IOException {
        if (chunkCapacity < 1) throw new IllegalArgumentException("BH run size must be positive");
        this.overwrite = overwrite;
        this.output = output.toAbsolutePath().normalize();
        if (Files.exists(this.output) && !overwrite)
            throw new IOException("output exists; use --overwrite: " + output);
        Path parent = this.output.getParent();
        if (parent != null) Files.createDirectories(parent);
        raw = Path.of(this.output + ".partial");
        if (Files.exists(raw) && !overwrite)
            throw new IOException(
                "partial output exists; use --resume or --overwrite: " + raw);
        work = Files.createTempDirectory(parent, ".jlinalg-bh-");
        qValues = work.resolve("q.bin");
        rawWriter = Files.newBufferedWriter(raw, StandardCharsets.UTF_8,
            overwrite ? java.nio.file.StandardOpenOption.CREATE
                : java.nio.file.StandardOpenOption.CREATE_NEW,
            java.nio.file.StandardOpenOption.TRUNCATE_EXISTING);
        initialQ = new DataOutputStream(new BufferedOutputStream(
            Files.newOutputStream(qValues)));
        this.chunkCapacity = chunkCapacity;
    }

    void writeHeader(List<String> fields) throws IOException {
        if (closed || headerWritten) throw new IOException("BH header already written or output closed");
        rawWriter.write(join(fields));
        rawWriter.newLine();
        headerWritten = true;
    }

    void write(List<String> fields, double pValue) throws IOException {
        if (closed || !headerWritten) throw new IOException("BH output requires an open header");
        rawWriter.write(join(fields));
        rawWriter.newLine();
        initialQ.writeDouble(Double.NaN);
        if (Double.isFinite(pValue) && pValue >= 0.0 && pValue <= 1.0) {
            buffered.add(new PValue(pValue, rows));
            tests++;
            if (buffered.size() >= chunkCapacity) flushChunk();
        }
        rows++;
    }

    long tests() { return tests; }

    void finish() throws IOException {
        if (finished) return;
        if (closed || !headerWritten) throw new IOException("BH output is closed or lacks a header");
        closed = true;
        rawWriter.close();
        initialQ.close();
        flushChunk();
        Path sorted = work.resolve("sorted.bin");
        merge(sorted);
        assign(sorted);
        Path complete = work.resolve("complete.tsv");
        try (BufferedReader input = Files.newBufferedReader(
                    raw, StandardCharsets.UTF_8);
             BufferedWriter writer = Files.newBufferedWriter(
                    complete, StandardCharsets.UTF_8);
             DataInputStream qInput = new DataInputStream(
                    new BufferedInputStream(Files.newInputStream(qValues)))) {
            String header = readRecord(input);
            writer.write(header);
            writer.write("\tfdr_bh");
            writer.newLine();
            for (String line; (line = readRecord(input)) != null;) {
                writer.write(line);
                writer.write('\t');
                double q = qInput.readDouble();
                if (Double.isFinite(q)) writer.write(Double.toString(q));
                writer.newLine();
            }
        }
        // Without overwrite, the move must itself reject a concurrent output;
        // ATOMIC_MOVE may replace it even without REPLACE_EXISTING.
        if (!overwrite) Files.move(complete, output);
        else try {
            Files.move(complete, output, StandardCopyOption.REPLACE_EXISTING,
                StandardCopyOption.ATOMIC_MOVE);
        } catch (AtomicMoveNotSupportedException exception) {
            Files.move(complete, output, StandardCopyOption.REPLACE_EXISTING);
        }
        Files.deleteIfExists(raw);
        cleanup();
        finished = true;
    }

    @Override
    public void close() throws IOException {
        try {
            if (!closed) {
                closed = true;
                try { rawWriter.close(); } finally { initialQ.close(); }
            }
        } finally { if (!finished) cleanup(); }
    }

    private void flushChunk() throws IOException {
        if (buffered.isEmpty()) return;
        buffered.sort(Comparator.comparingDouble(PValue::p)
            .thenComparingLong(PValue::row));
        Path chunk = work.resolve("chunk-" + chunks.size() + ".bin");
        try (DataOutputStream output = new DataOutputStream(
                new BufferedOutputStream(Files.newOutputStream(chunk)))) {
            output.writeLong(buffered.size());
            for (PValue value : buffered) {
                output.writeDouble(value.p());
                output.writeLong(value.row());
            }
        }
        chunks.add(chunk);
        buffered.clear();
    }

    private void merge(Path sorted) throws IOException {
        List<Path> runs = new ArrayList<>(chunks);
        // Limit open files independently of the number of input records.
        while (runs.size() > 64) {
            List<Path> next = new ArrayList<>();
            for (int first = 0; first < runs.size(); first += 64) {
                Path combined = work.resolve("chunk-" + chunks.size() + ".bin");
                chunks.add(combined);
                mergeRuns(runs.subList(first, Math.min(first + 64, runs.size())), combined, true);
                next.add(combined);
            }
            for (Path run : runs) Files.deleteIfExists(run);
            runs = next;
        }
        mergeRuns(runs, sorted, false);
    }

    private void mergeRuns(List<Path> runs, Path sorted, boolean header) throws IOException {
        PriorityQueue<Cursor> queue = new PriorityQueue<>(
            Comparator.comparingDouble((Cursor value) -> value.current.p())
                .thenComparingLong(value -> value.current.row()));
        List<Cursor> cursors = new ArrayList<>();
        try {
            long count = 0;
            for (Path chunk : runs) {
                Cursor cursor = new Cursor(chunk);
                cursors.add(cursor);
                if (cursor.current != null) queue.add(cursor);
                count = Math.addExact(count, cursor.remaining + (cursor.current == null ? 0 : 1));
            }
            try (DataOutputStream output = new DataOutputStream(
                    new BufferedOutputStream(Files.newOutputStream(sorted)))) {
                if (header) output.writeLong(count);
                while (!queue.isEmpty()) {
                    Cursor cursor = queue.remove();
                    output.writeDouble(cursor.current.p());
                    output.writeLong(cursor.current.row());
                    cursor.advance();
                    if (cursor.current != null) queue.add(cursor);
                }
            }
        } finally {
            for (Cursor cursor : cursors) cursor.close();
        }
    }

    private void assign(Path sorted) throws IOException {
        if (tests == 0) return;
        try (RandomAccessFile values = new RandomAccessFile(sorted.toFile(), "r");
             RandomAccessFile q = new RandomAccessFile(qValues.toFile(), "rw")) {
            double running = 1.0;
            byte[] record = new byte[16];
            ByteBuffer decoded = ByteBuffer.wrap(record);
            byte[] encodedQ = new byte[8];
            ByteBuffer qBuffer = ByteBuffer.wrap(encodedQ);
            for (long rank = tests; rank >= 1; rank--) {
                values.seek((rank - 1) * 16);
                values.readFully(record);
                double p = decoded.getDouble(0);
                long row = decoded.getLong(8);
                running = Math.min(running, p * tests / rank);
                q.seek(row * 8);
                qBuffer.putDouble(0, Math.min(1.0, running));
                q.write(encodedQ);
            }
        }
    }

    private void cleanup() throws IOException {
        for (Path chunk : chunks) Files.deleteIfExists(chunk);
        Files.deleteIfExists(work.resolve("sorted.bin"));
        Files.deleteIfExists(qValues);
        Files.deleteIfExists(work.resolve("complete.tsv"));
        Files.deleteIfExists(work);
    }

    /** Reads one quoted TSV record, retaining embedded CR/LF verbatim. */
    private static String readRecord(BufferedReader input) throws IOException {
        StringBuilder record = new StringBuilder();
        boolean quoted = false;
        for (int value; (value = input.read()) != -1;) {
            char c = (char) value;
            if (c == '"') quoted = !quoted; // doubled quotes toggle twice
            if (!quoted && (c == '\n' || c == '\r')) {
                if (c == '\r') {
                    input.mark(1);
                    if (input.read() != '\n') input.reset();
                }
                return record.toString();
            }
            record.append(c);
        }
        if (quoted) throw new IOException("unterminated quoted BH record");
        return record.isEmpty() ? null : record.toString();
    }

    private static String join(List<String> fields) {
        StringBuilder result = new StringBuilder();
        for (int index = 0; index < fields.size(); index++) {
            if (index > 0) result.append('\t');
            result.append(escape(fields.get(index)));
        }
        return result.toString();
    }

    private static String escape(String value) {
        String text = value == null ? "" : value;
        if (text.indexOf('\t') < 0 && text.indexOf('\n') < 0
                && text.indexOf('\r') < 0 && text.indexOf('"') < 0)
            return text;
        return "\"" + text.replace("\"", "\"\"") + "\"";
    }

    private record PValue(double p, long row) { }

    private static final class Cursor implements AutoCloseable {
        private final DataInputStream input;
        private long remaining;
        private PValue current;
        private Cursor(Path path) throws IOException {
            input = new DataInputStream(new BufferedInputStream(
                Files.newInputStream(path)));
            remaining = input.readLong();
            advance();
        }
        private void advance() throws IOException {
            if (remaining == 0) {
                current = null;
                return;
            }
            try {
                current = new PValue(input.readDouble(), input.readLong());
                remaining--;
            } catch (EOFException exception) {
                throw new IOException("truncated BH sort chunk", exception);
            }
        }
        @Override public void close() throws IOException { input.close(); }
    }
}
