/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.UUID;

/** Human-readable, timestamped lab-notebook log. */
final class RunLog implements AutoCloseable {
    private final BufferedWriter writer;
    private final RunTiming timing = new RunTiming();
    private final String runId;
    private final boolean plain;
    private boolean completed;

    private RunLog(BufferedWriter writer, boolean plain) throws IOException {
        this.writer = writer;
        this.plain = plain;
        runId = UUID.randomUUID().toString();
        line("INFO", "run_id=" + runId);
        line("INFO", "started=" + timing.started());
    }

    static RunLog open(Path path, boolean append) throws IOException {
        return open(path, append, false);
    }

    static RunLog openPlain(Path path, boolean append) throws IOException {
        return open(path, append, true);
    }

    private static RunLog open(Path path, boolean append, boolean plain) throws IOException {
        Path absolute = path.toAbsolutePath().normalize();
        Path parent = absolute.getParent();
        if (parent != null) Files.createDirectories(parent);
        BufferedWriter writer = Files.newBufferedWriter(absolute,
            StandardCharsets.UTF_8, StandardOpenOption.CREATE,
            append ? StandardOpenOption.APPEND
                : StandardOpenOption.TRUNCATE_EXISTING);
        try { return new RunLog(writer, plain); }
        catch (IOException failure) { writer.close(); throw failure; }
    }

    void info(String message) throws IOException { line("INFO", message); }
    void warning(String message) throws IOException { line("WARN", message); }
    void error(String message) throws IOException { line("ERROR", message); }
    String runId() { return runId; }

    void metadata(String text) throws IOException {
        writer.write(text);
        if (!text.endsWith("\n")) writer.newLine();
        writer.flush();
    }

    void complete(String status) throws IOException {
        if (completed) return;
        RunTiming.End end = timing.finish();
        line("INFO", "finished=" + end.finished());
        line("INFO", "elapsed_ms=" + end.elapsedMillis());
        line("INFO", "elapsed=" + end.elapsed());
        line("INFO", "status=" + status);
        completed = true;
    }

    @Override public void close() throws IOException {
        try { if (!completed) complete("failed"); }
        finally { writer.close(); }
    }

    private void line(String level, String message) throws IOException {
        if (plain) {
            writer.write(level.equals("INFO") ? message : level.toLowerCase(java.util.Locale.ROOT) + "=" + message);
            writer.newLine();
            writer.flush();
            return;
        }
        writer.write(timestamp(OffsetDateTime.now()));
        writer.write(" [");
        writer.write(level);
        writer.write("] ");
        writer.write(message);
        writer.newLine();
        writer.flush();
    }

    private static String timestamp(OffsetDateTime value) {
        return DateTimeFormatter.ISO_OFFSET_DATE_TIME.format(
            value.withOffsetSameInstant(ZoneOffset.UTC));
    }
}
