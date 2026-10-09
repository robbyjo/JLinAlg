/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/** Immutable unusual-result observations and optional amendment proposals. */
public final class StudyObservations {
    private StudyObservations() { }

    private static final Pattern ID =
        Pattern.compile("[A-Za-z][A-Za-z0-9._-]{0,127}");

    public static Map<String, Object> record(Path project, String id,
            String type, Path artifact, List<String> hypotheses,
            String description, List<String> explanations,
            String proposedDirection, String amendmentId,
            StudyProject.Actor actor) throws IOException {
        return StudyProjectLock.with(project, () -> recordLocked(project, id, type, artifact,
            hypotheses, description, explanations, proposedDirection, amendmentId, actor));
    }

    private static Map<String, Object> recordLocked(Path project, String id,
            String type, Path artifact, List<String> hypotheses,
            String description, List<String> explanations,
            String proposedDirection, String amendmentId,
            StudyProject.Actor actor) throws IOException {
        if (!ID.matcher(id).matches())
            throw new IllegalArgumentException("Invalid observation ID: " + id);
        Path root = project.toAbsolutePath().normalize();
        Path source = artifact.isAbsolute() ? artifact.normalize()
            : root.resolve(artifact).normalize();
        if (!Files.isRegularFile(source))
            throw new IOException("Observed artifact is absent: " + source);
        if (amendmentId != null && !ID.matcher(amendmentId).matches())
            throw new IllegalArgumentException(
                "Invalid proposed amendment ID: " + amendmentId);
        Path proposedAmendment = amendmentId == null ? null
            : root.resolve("amendments/proposed/amendment-"
                + amendmentId + ".yaml");
        if (proposedAmendment != null && Files.exists(proposedAmendment))
            throw new IOException(
                "Proposed amendment already exists: " + amendmentId);
        Map<String, Object> resolved = StudyProject.resolved(root).document();
        java.util.Set<String> knownHypotheses = new java.util.HashSet<>();
        for (Object value : StudyIO.list(resolved.get("hypotheses"), "hypotheses"))
            knownHypotheses.add(String.valueOf(
                StudyIO.map(value, "hypothesis").get("id")));
        for (String hypothesis : hypotheses)
            if (!knownHypotheses.contains(hypothesis))
                throw new IllegalArgumentException(
                    "Observation references an unknown hypothesis: " + hypothesis);
        Path output = root.resolve("observations/observation-" + id + ".yaml");
        if (Files.exists(output))
            throw new IOException("Observation ID already exists: " + id);

        Map<String, Object> observation = new LinkedHashMap<>();
        observation.put("schema_version", 1);
        observation.put("id", id);
        observation.put("type", required(type, "observation type"));
        observation.put("description", required(description, "description"));
        observation.put("artifact", portable(root, source));
        observation.put("artifact_sha256", StudyIO.hash(source));
        observation.put("hypotheses", List.copyOf(hypotheses));
        observation.put("competing_explanations", List.copyOf(explanations));
        observation.put("proposed_direction",
            required(proposedDirection, "proposed direction"));
        observation.put("status", "proposed");
        observation.put("recorded_at", Instant.now().toString());
        observation.put("actor", actor.asMap());
        StudyIO.writeYaml(output, observation);

        Path proposal = null;
        if (amendmentId != null) {
            proposal = proposal(root, amendmentId, observation);
        }
        Map<String, Object> payload = new LinkedHashMap<>();
        payload.put("observation_id", id);
        payload.put("observation_type", type);
        payload.put("observation_sha256", StudyIO.hash(output));
        payload.put("artifact_sha256", observation.get("artifact_sha256"));
        if (proposal != null) {
            payload.put("proposed_amendment",
                portable(root, proposal));
            payload.put("proposed_amendment_sha256", StudyIO.hash(proposal));
        }
        payload.put("state", "waiting_for_input");
        payload.put("phase", "review");
        payload.put("summary", "Recorded unusual-result observation " + id);
        StudyProject.append(root, "observation_recorded", actor, payload);
        return observation;
    }

    private static Path proposal(Path root, String amendmentId,
            Map<String, Object> observation) throws IOException {
        Path path = root.resolve("amendments/proposed/amendment-"
            + amendmentId + ".yaml");
        if (Files.exists(path))
            throw new IOException("Proposed amendment already exists: " + amendmentId);
        Map<String, Object> header = new LinkedHashMap<>();
        header.put("id", amendmentId);
        header.put("title", "Follow up observation " + observation.get("id"));
        header.put("base_study_hash",
            StudyProject.resolvedHash(StudyProject.resolved(root).document()));
        header.put("proposed_at", Instant.now().toString());
        header.put("proposed_by", observation.get("actor"));
        header.put("source", Map.of("type", "diagnostic_observation",
            "id", observation.get("id")));
        header.put("rationale", observation.get("description"));
        Map<String, Object> request = new LinkedHashMap<>();
        request.put("id", "request-" + observation.get("id"));
        request.put("type", observation.get("proposed_direction"));
        request.put("for_hypotheses", observation.get("hypotheses"));
        request.put("description", "Investigate " + observation.get("description"));
        request.put("origin_observation", observation.get("id"));
        Map<String, Object> changes = Map.of("add",
            Map.of("analysis_requests", List.of(request)));
        Map<String, Object> document = new LinkedHashMap<>();
        document.put("schema_version", 1);
        document.put("amendment", header);
        document.put("changes", changes);
        StudyIO.writeYaml(path, document);
        return path;
    }

    private static String portable(Path root, Path path) {
        try {
            return root.relativize(path.toAbsolutePath().normalize())
                .toString().replace('\\', '/');
        } catch (IllegalArgumentException differentRoot) {
            return path.toAbsolutePath().normalize().toString();
        }
    }

    private static String required(String value, String name) {
        if (value == null || value.isBlank())
            throw new IllegalArgumentException(name + " is required");
        return value;
    }
}
