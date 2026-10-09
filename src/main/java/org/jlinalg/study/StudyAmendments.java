/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/** Validation, impact preview, and immutable application of amendments. */
public final class StudyAmendments {
    private StudyAmendments() { }

    private static final Pattern ID =
        Pattern.compile("[A-Za-z][A-Za-z0-9._-]{0,127}");
    private static final Set<String> COLLECTIONS = Set.of(
        "cohorts", "cohort_relationships", "datasets", "hypotheses",
        "analysis_requests");

    public static Map<String, Object> preview(Path amendment, Path project)
            throws IOException {
        Loaded loaded = load(amendment, project);
        Map<String, Object> impact = impact(loaded.changes);
        impact.put("amendment_id", loaded.id);
        impact.put("base_study_hash", loaded.baseHash);
        impact.put("resulting_resolved_study_hash",
            StudyProject.resolvedHash(loaded.resolved.document()));
        return impact;
    }

    public static Map<String, Object> apply(Path amendment, Path project,
            StudyProject.Actor approver) throws IOException {
        return StudyProjectLock.with(project, () -> applyLocked(amendment, project, approver));
    }

    private static Map<String, Object> applyLocked(Path amendment, Path project,
            StudyProject.Actor approver) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        Map<String, Object> raw = StudyIO.load(amendment);
        Map<String, Object> header = StudyIO.map(raw.get("amendment"),
            "amendment");
        String proposedId = StudyIO.requireString(header, "id", "amendment");
        if (applied(root).contains(proposedId))
            throw new IllegalArgumentException(
                "Amendment ID has already been applied: " + proposedId);
        Loaded loaded = load(amendment, root);
        if ("scientist_answers".equals(StudyIO.optionalMap(
                StudyIO.map(loaded.amendment.get("amendment"), "amendment").get("source"), "source").get("type"))
                && (!"human".equals(approver.type()) || "unspecified".equals(approver.id())))
            throw new IllegalArgumentException("Answer-derived amendments require an explicit human approver");

        Path destination = root.resolve("amendments/amendment-" + loaded.id + ".yaml");
        if (java.nio.file.Files.exists(destination))
            throw new IOException(
                "Accepted amendment path already exists but is not in the ledger: "
                    + destination);
        StudyIO.writeYaml(destination, loaded.amendment);
        StudyIO.writeYaml(root.resolve("resolved-study.yaml"),
            loaded.resolved.document());
        Map<String, Object> payload = impact(loaded.changes);
        payload.put("amendment_id", loaded.id);
        payload.put("amendment_sha256", StudyIO.hash(destination));
        payload.put("base_study_hash", loaded.baseHash);
        payload.put("resolved_study_hash",
            StudyProject.resolvedHash(loaded.resolved.document()));
        payload.put("state", "ready_for_planning");
        payload.put("phase", "planning");
        payload.put("summary", "Approved and applied amendment " + loaded.id);
        List<Object> questionUpdates = StudyQuestions.amendmentUpdates(loaded.amendment);
        if (!questionUpdates.isEmpty()) payload.put("question_updates", questionUpdates);
        StudyProject.append(root, "amendment_applied", approver, payload);
        StudyEvidenceIndex.rebuild(root);
        return payload;
    }

    private static Loaded load(Path path, Path project) throws IOException {
        StudyProject.verify(project);
        Map<String, Object> amendment = StudyIO.load(path);
        StudyQuestions.validateAmendment(project.toAbsolutePath().normalize(), amendment);
        for (String key : amendment.keySet())
            if (!Set.of("schema_version", "amendment", "changes").contains(key))
                throw new IllegalArgumentException(
                    "Unknown amendment field: " + key);
        if (!Integer.valueOf(1).equals(amendment.get("schema_version")))
            throw new IllegalArgumentException(
                "Amendment requires schema_version: 1");
        Map<String, Object> header = StudyIO.map(
            amendment.get("amendment"), "amendment");
        for (String key : header.keySet())
            if (!Set.of("id", "title", "base_study_hash", "proposed_at",
                    "proposed_by", "source", "rationale", "revision").contains(key))
                throw new IllegalArgumentException(
                    "Unknown amendment metadata field: " + key);
        String id = StudyIO.requireString(header, "id", "amendment");
        if (!ID.matcher(id).matches())
            throw new IllegalArgumentException("Invalid amendment ID: " + id);
        String file = path.getFileName().toString();
        if (!file.equals("amendment-" + id + ".yaml")
                && !file.equals("amendment-" + id + ".yml"))
            throw new IllegalArgumentException(
                "Amendment filename must be amendment-" + id + ".yaml");
        String baseHash = StudyIO.requireString(
            header, "base_study_hash", "amendment");
        StudyContract current = StudyProject.resolved(project);
        String currentHash = StudyProject.resolvedHash(current.document());
        if (!baseHash.equals(currentHash))
            throw new IllegalArgumentException(
                "Amendment base_study_hash does not match the current resolved study; expected "
                    + currentHash);
        Map<String, Object> changes = StudyIO.map(
            amendment.get("changes"), "changes");
        if (changes.isEmpty())
            throw new IllegalArgumentException("Amendment changes cannot be empty");
        for (String key : changes.keySet())
            if (!Set.of("add", "update", "retire", "configure").contains(key))
                throw new IllegalArgumentException(
                    "Unknown amendment change section: " + key);

        Map<String, Object> effective = current.document();
        effective.remove("assumptions");
        effective.remove("_resolution");
        applyChanges(effective, changes);
        StudyValidation validation = StudyValidator.validateDocument(
            effective, project.toAbsolutePath().normalize(), false);
        if (!validation.valid())
            throw new IllegalArgumentException(
                "Amended study is invalid:\n" + validation.report());
        List<String> amendments = new ArrayList<>(applied(project));
        amendments.add(id);
        StudyContract source = new StudyContract(
            project.resolve("resolved-study.yaml"), effective);
        StudyContract resolved = StudyResolver.resolve(source, amendments);
        return new Loaded(id, baseHash, amendment, changes, resolved);
    }

    private static void applyChanges(Map<String, Object> target,
            Map<String, Object> changes) {
        Map<String, Object> configuration = StudyIO.optionalMap(changes.get("configure"), "changes.configure");
        for (Map.Entry<String, Object> entry : configuration.entrySet()) {
            if (!Set.of("assistance", "resources", "project").contains(entry.getKey()))
                throw new IllegalArgumentException("Unsupported configuration amendment: " + entry.getKey());
            Map<String, Object> existing = StudyIO.optionalMap(target.get(entry.getKey()), entry.getKey());
            deepMerge(existing, StudyIO.map(entry.getValue(), entry.getKey()));
            target.put(entry.getKey(), existing);
        }
        Map<String, Object> add = StudyIO.optionalMap(changes.get("add"), "changes.add");
        Map<String, Object> update = StudyIO.optionalMap(
            changes.get("update"), "changes.update");
        Map<String, Object> retire = StudyIO.optionalMap(
            changes.get("retire"), "changes.retire");
        knownCollections(add, "changes.add");
        knownCollections(update, "changes.update");
        knownCollections(retire, "changes.retire");

        for (Map.Entry<String, Object> entry : add.entrySet()) {
            List<Object> current = new ArrayList<>(StudyIO.optionalList(
                target.get(entry.getKey()), entry.getKey()));
            Set<String> ids = ids(current, entry.getKey());
            for (Object value : StudyIO.list(entry.getValue(),
                    "changes.add." + entry.getKey())) {
                Map<String, Object> item = StudyIO.map(value, entry.getKey());
                String id = objectId(item, entry.getKey());
                if (!ids.add(id))
                    throw new IllegalArgumentException(
                        "Cannot add duplicate " + entry.getKey() + " ID: " + id);
                current.add(item);
            }
            target.put(entry.getKey(), current);
        }

        for (Map.Entry<String, Object> entry : update.entrySet()) {
            List<Object> current = new ArrayList<>(StudyIO.optionalList(
                target.get(entry.getKey()), entry.getKey()));
            for (Object value : StudyIO.list(entry.getValue(),
                    "changes.update." + entry.getKey())) {
                Map<String, Object> patch = StudyIO.map(value, entry.getKey());
                String id = StudyIO.requireString(patch, "id", entry.getKey());
                int index = index(current, id, entry.getKey());
                Map<String, Object> item = StudyIO.map(current.get(index), entry.getKey());
                deepMerge(item, patch);
                current.set(index, item);
            }
            target.put(entry.getKey(), current);
        }

        for (Map.Entry<String, Object> entry : retire.entrySet()) {
            List<Object> current = new ArrayList<>(StudyIO.optionalList(
                target.get(entry.getKey()), entry.getKey()));
            for (Object value : StudyIO.list(entry.getValue(),
                    "changes.retire." + entry.getKey())) {
                Map<String, Object> retirement = StudyIO.map(value, entry.getKey());
                String id = StudyIO.requireString(retirement, "id", entry.getKey());
                String reason = StudyIO.requireString(
                    retirement, "reason", entry.getKey());
                int index = index(current, id, entry.getKey());
                Map<String, Object> item = StudyIO.map(current.get(index), entry.getKey());
                item.put("status", "retired");
                item.put("retired_reason", reason);
                current.set(index, item);
            }
            target.put(entry.getKey(), current);
        }
    }

    private static void deepMerge(Map<String, Object> target,
            Map<String, Object> patch) {
        for (Map.Entry<String, Object> entry : patch.entrySet()) {
            if (entry.getKey().equals("id")) continue;
            Object value = entry.getValue();
            if (value instanceof Map<?, ?> raw) {
                Map<String, Object> map = StudyIO.map(raw, entry.getKey());
                if (map.size() == 1 && map.containsKey("append")) {
                    List<Object> result = new ArrayList<>(StudyIO.optionalList(
                        target.get(entry.getKey()), entry.getKey()));
                    for (Object appended : StudyIO.list(
                            map.get("append"), entry.getKey() + ".append"))
                        if (!result.contains(appended)) result.add(appended);
                    target.put(entry.getKey(), result);
                } else if (target.get(entry.getKey()) instanceof Map<?, ?> existing) {
                    Map<String, Object> child = StudyIO.map(existing, entry.getKey());
                    deepMerge(child, map);
                    target.put(entry.getKey(), child);
                } else target.put(entry.getKey(), map);
            } else target.put(entry.getKey(), value);
        }
    }

    private static int index(List<Object> items, String id, String collection) {
        for (int index = 0; index < items.size(); index++) {
            Map<String, Object> item = StudyIO.map(items.get(index), collection);
            if (id.equals(item.get("id"))) return index;
        }
        throw new IllegalArgumentException(
            "Cannot update or retire absent " + collection + " ID: " + id);
    }

    private static Set<String> ids(List<Object> items, String collection) {
        Set<String> result = new HashSet<>();
        for (Object value : items) {
            Map<String, Object> item = StudyIO.map(value, collection);
            result.add(objectId(item, collection));
        }
        return result;
    }

    private static String objectId(Map<String, Object> item, String collection) {
        if (item.get("id") instanceof String id && !id.isBlank()) return id;
        if (collection.equals("cohort_relationships")) {
            List<Object> between = StudyIO.list(item.get("between"),
                "cohort_relationships.between");
            return String.valueOf(item.get("type")) + ":" + between;
        }
        return StudyIO.requireString(item, "id", collection);
    }

    private static void knownCollections(Map<String, Object> section,
            String path) {
        for (String key : section.keySet())
            if (!COLLECTIONS.contains(key))
                throw new IllegalArgumentException(
                    path + " has unknown collection: " + key);
    }

    private static Map<String, Object> impact(Map<String, Object> changes) {
        Map<String, Object> result = new LinkedHashMap<>();
        if (changes.containsKey("configure")) result.put("configure", new ArrayList<>(
            StudyIO.map(changes.get("configure"), "configure").keySet()));
        for (String action : List.of("add", "update", "retire")) {
            Map<String, Object> section = StudyIO.optionalMap(
                changes.get(action), "changes." + action);
            Map<String, Object> counts = new LinkedHashMap<>();
            for (Map.Entry<String, Object> entry : section.entrySet())
                counts.put(entry.getKey(),
                    StudyIO.list(entry.getValue(), entry.getKey()).size());
            result.put(action, counts);
        }
        return result;
    }

    private static Set<String> applied(Path project) throws IOException {
        Set<String> result = new HashSet<>();
        for (Map<String, Object> event : StudyLedger.read(project))
            if ("amendment_applied".equals(event.get("event_type"))) {
                Map<String, Object> payload = StudyIO.map(
                    event.get("payload"), "event.payload");
                result.add(String.valueOf(payload.get("amendment_id")));
            }
        return result;
    }

    private record Loaded(String id, String baseHash,
            Map<String, Object> amendment, Map<String, Object> changes,
            StudyContract resolved) { }
}
