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
import java.time.Duration;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.UUID;

/** Human-readable, timestamped lab-notebook log. */
final class RunLog implements AutoCloseable {
    private final BufferedWriter writer;
    private final OffsetDateTime started;
    private final String runId;

    private RunLog(BufferedWriter writer) throws IOException {
        this.writer = writer;
        started = OffsetDateTime.now();
        runId = UUID.randomUUID().toString();
        line("INFO", "run_id=" + runId);
        line("INFO", "started=" + timestamp(started));
    }

    static RunLog open(Path path, boolean append) throws IOException {
        Path absolute = path.toAbsolutePath().normalize();
        Path parent = absolute.getParent();
        if (parent != null) Files.createDirectories(parent);
        BufferedWriter writer = Files.newBufferedWriter(absolute,
            StandardCharsets.UTF_8, StandardOpenOption.CREATE,
            append ? StandardOpenOption.APPEND
                : StandardOpenOption.TRUNCATE_EXISTING);
        return new RunLog(writer);
    }

    void info(String message) throws IOException { line("INFO", message); }
    void warning(String message) throws IOException { line("WARN", message); }
    void error(String message) throws IOException { line("ERROR", message); }
    String runId() { return runId; }

    void complete(String status) throws IOException {
        OffsetDateTime ended = OffsetDateTime.now();
        line("INFO", "finished=" + timestamp(ended));
        line("INFO", "elapsed_ms="
            + Duration.between(started, ended).toMillis());
        line("INFO", "status=" + status);
    }

    @Override public void close() throws IOException { writer.close(); }

    private void line(String level, String message) throws IOException {
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
