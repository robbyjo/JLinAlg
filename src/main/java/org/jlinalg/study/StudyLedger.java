/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Append-only, hash-chained project event ledger. */
final class StudyLedger {
    private StudyLedger() { }

    static List<Map<String, Object>> read(Path project) throws IOException {
        Path ledger = project.resolve("ledger/events.jsonl");
        if (!Files.exists(ledger)) return List.of();
        List<Map<String, Object>> result = new ArrayList<>();
        int lineNumber = 0;
        for (String line : Files.readAllLines(ledger, StandardCharsets.UTF_8)) {
            lineNumber++;
            if (line.isBlank()) continue;
            Object parsed;
            try {
                parsed = StudyIO.reader().load(line);
            } catch (RuntimeException failure) {
                throw new IOException("Invalid event ledger JSON at line "
                    + lineNumber + ": " + failure.getMessage(), failure);
            }
            result.add(StudyIO.map(parsed, "event line " + lineNumber));
        }
        verify(result);
        return List.copyOf(result);
    }

    static Map<String, Object> append(Path project, String type,
            StudyProject.Actor actor, Map<String, Object> payload)
            throws IOException {
        Path directory = project.resolve("ledger");
        Files.createDirectories(directory);
        Path lockPath = directory.resolve(".events.lock");
        try (FileChannel channel = FileChannel.open(lockPath,
                StandardOpenOption.CREATE, StandardOpenOption.WRITE);
             FileLock lock = channel.lock()) {
            if (!lock.isValid()) throw new IOException("Could not lock the event ledger");
            List<Map<String, Object>> existing = new ArrayList<>(read(project));
            String previous = existing.isEmpty() ? null
                : String.valueOf(existing.get(existing.size() - 1).get("event_hash"));
            Map<String, Object> event = new LinkedHashMap<>();
            event.put("schema_version", 1);
            event.put("event_id", String.format("event-%06d", existing.size() + 1));
            event.put("timestamp", Instant.now().toString());
            event.put("event_type", type);
            event.put("actor", actor.asMap());
            event.put("payload", StudyIO.copyMap(payload));
            event.put("previous_event_hash", previous);
            event.put("event_hash", StudyIO.hash(event));
            existing.add(event);
            StringBuilder text = new StringBuilder();
            for (Map<String, Object> item : existing)
                text.append(StudyIO.canonicalJson(item)).append('\n');
            StudyIO.atomicWrite(directory.resolve("events.jsonl"), text.toString());
            return event;
        }
    }

    static void verify(List<Map<String, Object>> events) throws IOException {
        String previous = null;
        for (int index = 0; index < events.size(); index++) {
            Map<String, Object> event = events.get(index);
            String expectedId = String.format("event-%06d", index + 1);
            if (!expectedId.equals(event.get("event_id")))
                throw new IOException("Event ledger sequence is invalid at "
                    + expectedId);
            Object linked = event.get("previous_event_hash");
            if (!(previous == null ? linked == null : previous.equals(linked)))
                throw new IOException("Event ledger chain is broken at "
                    + expectedId);
            Object recorded = event.get("event_hash");
            Map<String, Object> unhashed = StudyIO.copyMap(event);
            unhashed.remove("event_hash");
            String calculated = StudyIO.hash(unhashed);
            if (!calculated.equals(recorded))
                throw new IOException("Event ledger hash mismatch at "
                    + expectedId);
            previous = calculated;
        }
    }
}
