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

    ExternalBh(Path output, boolean overwrite) throws IOException {
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
        rawWriter = Files.newBufferedWriter(raw, StandardCharsets.UTF_8);
        initialQ = new DataOutputStream(new BufferedOutputStream(
            Files.newOutputStream(qValues)));
        long memory = Runtime.getRuntime().maxMemory();
        chunkCapacity = (int) Math.max(10_000,
            Math.min(500_000, memory / 512));
    }

    void writeHeader(List<String> fields) throws IOException {
        rawWriter.write(String.join("\t", fields));
        rawWriter.newLine();
    }

    void write(List<String> fields, double pValue) throws IOException {
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
        finished = true;
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
            String header = input.readLine();
            writer.write(header);
            writer.write("\tfdr_bh");
            writer.newLine();
            for (String line; (line = input.readLine()) != null;) {
                writer.write(line);
                writer.write('\t');
                double q = qInput.readDouble();
                if (Double.isFinite(q)) writer.write(Double.toString(q));
                writer.newLine();
            }
        }
        Files.move(complete, output, StandardCopyOption.REPLACE_EXISTING,
            StandardCopyOption.ATOMIC_MOVE);
        Files.deleteIfExists(raw);
        cleanup();
    }

    @Override
    public void close() throws IOException {
        if (!finished) {
            rawWriter.close();
            initialQ.close();
        }
    }

    private void flushChunk() throws IOException {
        if (buffered.isEmpty()) return;
        buffered.sort(Comparator.comparingDouble(PValue::p)
            .thenComparingLong(PValue::row));
        Path chunk = work.resolve("chunk-" + chunks.size() + ".bin");
        try (DataOutputStream output = new DataOutputStream(
                new BufferedOutputStream(Files.newOutputStream(chunk)))) {
            output.writeInt(buffered.size());
            for (PValue value : buffered) {
                output.writeDouble(value.p());
                output.writeLong(value.row());
            }
        }
        chunks.add(chunk);
        buffered.clear();
    }

    private void merge(Path sorted) throws IOException {
        PriorityQueue<Cursor> queue = new PriorityQueue<>(
            Comparator.comparingDouble((Cursor value) -> value.current.p())
                .thenComparingLong(value -> value.current.row()));
        List<Cursor> cursors = new ArrayList<>();
        try {
            for (Path chunk : chunks) {
                Cursor cursor = new Cursor(chunk);
                cursors.add(cursor);
                if (cursor.current != null) queue.add(cursor);
            }
            try (DataOutputStream output = new DataOutputStream(
                    new BufferedOutputStream(Files.newOutputStream(sorted)))) {
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
            for (long rank = tests; rank >= 1; rank--) {
                values.seek((rank - 1) * 16);
                double p = values.readDouble();
                long row = values.readLong();
                running = Math.min(running, p * tests / rank);
                q.seek(row * 8);
                q.writeDouble(Math.min(1.0, running));
            }
        }
    }

    private void cleanup() throws IOException {
        for (Path chunk : chunks) Files.deleteIfExists(chunk);
        Files.deleteIfExists(work.resolve("sorted.bin"));
        Files.deleteIfExists(qValues);
        Files.deleteIfExists(work);
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
        private int remaining;
        private PValue current;
        private Cursor(Path path) throws IOException {
            input = new DataInputStream(new BufferedInputStream(
                Files.newInputStream(path)));
            remaining = input.readInt();
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
