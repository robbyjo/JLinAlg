/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Project state, ledger, and derived human/machine audit views. */
public final class StudyProject {
    private StudyProject() { }

    public record Actor(String type, String id, String model) {
        public Actor {
            type = required(type, "actor type").toLowerCase(Locale.ROOT);
            id = required(id, "actor id");
            if (!List.of("human", "jlinalg", "llm", "external").contains(type))
                throw new IllegalArgumentException(
                    "actor type must be human, jlinalg, llm, or external");
            if (type.equals("llm") && (model == null || model.isBlank()))
                throw new IllegalArgumentException("LLM actors require a model identifier");
        }

        public static Actor human(String id) {
            return new Actor("human", id, null);
        }

        public Map<String, Object> asMap() {
            Map<String, Object> result = new LinkedHashMap<>();
            result.put("type", type);
            result.put("id", id);
            if (model != null) result.put("model", model);
            result.put("software", "JLinAlg");
            result.put("software_version", version());
            return result;
        }

        private static String required(String value, String name) {
            if (value == null || value.isBlank())
                throw new IllegalArgumentException(name + " is required");
            return value;
        }
    }

    public static Map<String, Object> resolve(Path source, Path project,
            Actor actor) throws IOException {
        return StudyProjectLock.with(project, () -> resolveLocked(source, project, actor));
    }

    private static Map<String, Object> resolveLocked(Path source, Path project,
            Actor actor) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        Files.createDirectories(root);
        StudyContract human = StudyValidator.requireValid(source, false);
        Map<String, Object> current = status(root);
        if (current.get("study_hash") != null) {
            if (!human.hash().equals(current.get("study_hash")))
                throw new IllegalArgumentException(
                    "A resolved project cannot replace its base study directly; submit an amendment");
            StudyContract stored = resolved(root);
            String storedHash = resolvedHash(stored.document());
            if (!storedHash.equals(current.get("resolved_study_hash")))
                throw new IllegalArgumentException(
                    "The stored resolved study does not match the project ledger");
            return current;
        }
        if (!"beginning".equals(current.get("state"))) {
            throw new IllegalArgumentException(
                "The project ledger has no valid base study; verification is required");
        }
        StudyContract resolved = StudyResolver.resolve(human, List.of());
        Map<String, Object> document = resolved.document();
        String resolvedHash = resolvedHash(document);

        if (resolvedHash.equals(current.get("resolved_study_hash"))
                && Files.isRegularFile(root.resolve("resolved-study.yaml")))
            return current;

        StudyIO.writeYaml(root.resolve("resolved-study.yaml"), document);
        if (!source.toAbsolutePath().normalize().equals(root.resolve("study.yaml"))
                && !Files.exists(root.resolve("study.yaml")))
            StudyIO.atomicWrite(root.resolve("study.yaml"),
                Files.readString(source));
        Map<String, Object> payload = new LinkedHashMap<>();
        payload.put("project_id", resolved.projectId());
        payload.put("study_hash", human.hash());
        payload.put("resolved_study_hash", resolvedHash);
        payload.put("state", "ready_for_planning");
        payload.put("phase", "planning");
        payload.put("summary", "Validated and resolved the study contract");
        StudyLedger.append(root, "study_resolved", actor, payload);
        regenerate(root);
        return status(root);
    }

    public static Map<String, Object> status(Path project) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        List<Map<String, Object>> events = StudyLedger.read(root);
        if (events.isEmpty()) return beginning(root);
        return deriveStatus(root, events);
    }

    public static void regenerate(Path project) throws IOException {
        StudyProjectLock.with(project, () -> { regenerateLocked(project); return null; });
    }

    private static void regenerateLocked(Path project) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        List<Map<String, Object>> events = StudyLedger.read(root);
        if (events.isEmpty()) return;
        Map<String, Object> status = deriveStatus(root, events);
        StudyIO.writeYaml(root.resolve("state/project-status.yaml"), status);
        StudyIO.atomicWrite(root.resolve("SESSION_HISTORY.md"),
            sessionHistory(events));
        StudyQuestions.materialize(root, events);
        if (Files.isRegularFile(root.resolve("resolved-study.yaml"))) {
            Map<String, Object> resolved = StudyIO.load(
                root.resolve("resolved-study.yaml"));
            StudyIO.atomicWrite(root.resolve("AGENTS.md"),
                agents(resolved, status));
        }
    }

    public static void verify(Path project) throws IOException {
        verify(project, true);
    }

    /** Repair derived notebooks only after all authoritative records pass verification. */
    public static void repairViews(Path project) throws IOException {
        StudyProjectLock.with(project, () -> {
            verify(project, false);
            regenerate(project);
            StudyEvidenceIndex.rebuild(project);
            return null;
        });
    }

    private static void verify(Path project, boolean checkViews) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        List<Map<String, Object>> events = StudyLedger.read(root);
        if (events.isEmpty()) return;
        Map<String, Object> status = deriveStatus(root, events);
        Path baseStudy = root.resolve("study.yaml");
        if (!Files.isRegularFile(baseStudy))
            throw new IOException("study.yaml is missing");
        String baseHash = StudyIO.hash(StudyIO.load(baseStudy));
        if (!baseHash.equals(status.get("study_hash")))
            throw new IOException("study.yaml hash does not match the project ledger");
        if (!Files.isRegularFile(root.resolve("resolved-study.yaml")))
            throw new IOException("resolved-study.yaml is missing");
        Map<String, Object> resolved = StudyIO.load(root.resolve("resolved-study.yaml"));
        String expected = resolvedHash(resolved);
        Object recorded = status.get("resolved_study_hash");
        if (recorded != null && !expected.equals(recorded))
            throw new IOException("resolved-study.yaml hash does not match project status");
        StudyQuestions.verifyArtifacts(root, events);
        if (checkViews) {
            Path statusPath = root.resolve("state/project-status.yaml");
            if (!Files.isRegularFile(statusPath))
                throw new IOException("state/project-status.yaml is missing");
            if (!StudyIO.canonicalJson(StudyIO.load(statusPath))
                    .equals(StudyIO.canonicalJson(status)))
                throw new IOException("project-status.yaml is stale or modified");
            Path history = root.resolve("SESSION_HISTORY.md");
            if (!Files.isRegularFile(history)
                    || !Files.readString(history).equals(sessionHistory(events)))
                throw new IOException("SESSION_HISTORY.md is missing, stale, or modified");
            Path agents = root.resolve("AGENTS.md");
            if (!Files.isRegularFile(agents)
                    || !Files.readString(agents).equals(agents(resolved, status)))
                throw new IOException("AGENTS.md is missing, stale, or modified");
            if (StudyQuestions.hasQuestions(events)) {
                Map<String, Object> questions = StudyQuestions.view(events);
                Path statePath = root.resolve("state/questions.yaml");
                Path notebook = root.resolve("QUESTIONS.md");
                if (!Files.isRegularFile(statePath) || !StudyIO.canonicalJson(questions)
                        .equals(StudyIO.canonicalJson(StudyIO.load(statePath))))
                    throw new IOException("state/questions.yaml is missing, stale, or modified");
                if (!Files.isRegularFile(notebook) || !StudyQuestions.notebook(questions).equals(Files.readString(notebook)))
                    throw new IOException("QUESTIONS.md is missing, stale, or modified");
            }
        }

        Map<String, String> plans = new LinkedHashMap<>();
        Map<String, String> manifests = new LinkedHashMap<>();
        for (Map<String, Object> event : events) {
            Map<String, Object> payload = StudyIO.optionalMap(
                event.get("payload"), "event.payload");
            String type = String.valueOf(event.get("event_type"));
            if (type.equals("plan_created"))
                plans.put(String.valueOf(payload.get("plan_id")),
                    String.valueOf(payload.get("plan_hash")));
            if (type.startsWith("operation_"))
                manifests.put(payload.get("run_id") + "/"
                    + payload.get("operation_id"),
                    String.valueOf(payload.get("manifest_sha256")));
            if (type.equals("amendment_applied"))
                verifyAmendment(root, payload);
            if (type.equals("observation_recorded"))
                verifyObservation(root, payload);
        }
        for (Map.Entry<String, String> plan : plans.entrySet()) {
            Path path = root.resolve("plans").resolve(plan.getKey())
                .resolve("analysis-plan.yaml");
            if (!Files.isRegularFile(path))
                throw new IOException("Recorded plan is missing: " + plan.getKey());
            Map<String, Object> document = StudyIO.load(path);
            if (!plan.getKey().equals(document.get("plan_id"))
                    || !plan.getValue().equals(document.get("plan_hash"))
                    || !plan.getValue().equals(
                        StudyPlanner.hashWithoutPlanHash(document)))
                throw new IOException("Recorded plan is modified: " + plan.getKey());
        }
        Path runs = root.resolve("runs");
        if (Files.isDirectory(runs)) try (var paths = Files.walk(runs)) {
            for (Path manifestPath : paths.filter(path ->
                    path.getFileName().toString().equals("manifest.yaml")
                        && !path.getParent().getFileName().toString().startsWith(".operation-")).toList()) {
                Map<String, Object> manifest = StudyIO.load(manifestPath);
                String manifestKey = manifest.get("run_id") + "/"
                    + manifest.get("operation_id");
                String manifestHash = manifests.remove(manifestKey);
                if (manifestHash == null
                        || !manifestHash.equals(StudyIO.hash(manifestPath)))
                    throw new IOException(
                        "Operation manifest is not anchored in the ledger: "
                            + manifestPath);
                Map<String, Object> outputs = StudyIO.map(
                    manifest.get("output_sha256"), "output_sha256");
                Path directory = manifestPath.getParent();
                for (Map.Entry<String, Object> output : outputs.entrySet()) {
                    Path artifact = directory.resolve(output.getKey()).normalize();
                    if (!artifact.startsWith(directory)
                            || !Files.isRegularFile(artifact)
                            || !String.valueOf(output.getValue())
                                .equals(StudyIO.hash(artifact)))
                        throw new IOException("Artifact hash mismatch: " + artifact);
                }
            }
        }
        if (!manifests.isEmpty())
            throw new IOException(
                "Ledger references missing operation manifests: "
                    + manifests.keySet());
    }

    private static void verifyAmendment(Path root, Map<String, Object> payload)
            throws IOException {
        String id = String.valueOf(payload.get("amendment_id"));
        Path path = root.resolve("amendments/amendment-" + id + ".yaml");
        if (!Files.isRegularFile(path)
                || !String.valueOf(payload.get("amendment_sha256"))
                    .equals(StudyIO.hash(path)))
            throw new IOException("Accepted amendment is missing or modified: " + id);
    }

    private static void verifyObservation(Path root, Map<String, Object> payload)
            throws IOException {
        String id = String.valueOf(payload.get("observation_id"));
        Path path = root.resolve("observations/observation-" + id + ".yaml");
        if (!Files.isRegularFile(path)
                || !String.valueOf(payload.get("observation_sha256"))
                    .equals(StudyIO.hash(path)))
            throw new IOException("Observation is missing or modified: " + id);
        Map<String, Object> observation = StudyIO.load(path);
        Path artifact = Path.of(String.valueOf(observation.get("artifact")));
        if (!artifact.isAbsolute()) artifact = root.resolve(artifact).normalize();
        if (!Files.isRegularFile(artifact)
                || !String.valueOf(observation.get("artifact_sha256"))
                    .equals(StudyIO.hash(artifact)))
            throw new IOException("Observed artifact is missing or modified: " + artifact);
        if (payload.get("proposed_amendment") != null) {
            Path proposal = root.resolve(String.valueOf(
                payload.get("proposed_amendment"))).normalize();
            if (!proposal.startsWith(root) || !Files.isRegularFile(proposal)
                    || !String.valueOf(payload.get("proposed_amendment_sha256"))
                        .equals(StudyIO.hash(proposal)))
                throw new IOException(
                    "Proposed amendment is missing or modified: " + proposal);
        }
    }

    static Map<String, Object> append(Path project, String type, Actor actor,
            Map<String, Object> payload) throws IOException {
        return StudyProjectLock.with(project, () -> appendLocked(project, type, actor, payload));
    }

    private static Map<String, Object> appendLocked(Path project, String type, Actor actor,
            Map<String, Object> payload) throws IOException {
        Map<String, Object> event = StudyLedger.append(
            project.toAbsolutePath().normalize(), type, actor, payload);
        regenerate(project);
        return event;
    }

    static StudyContract resolved(Path project) throws IOException {
        Path path = project.toAbsolutePath().normalize().resolve("resolved-study.yaml");
        Map<String, Object> document = StudyIO.load(path);
        StudyValidation validation = StudyValidator.validateDocument(
            document, path.getParent(), false);
        if (!validation.valid())
            throw new IllegalArgumentException(
                "Resolved study validation failed:\n" + validation.report());
        return new StudyContract(path, document);
    }

    static String resolvedHash(Map<String, Object> document) {
        Map<String, Object> resolution = StudyIO.optionalMap(
            document.get("_resolution"), "_resolution");
        Object recorded = resolution.get("resolved_study_hash");
        String calculated = StudyResolver.hashWithoutResolvedHash(document);
        if (recorded != null && !calculated.equals(recorded))
            throw new IllegalArgumentException(
                "resolved-study.yaml has an invalid resolved_study_hash");
        return calculated;
    }

    private static Map<String, Object> beginning(Path root) {
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("schema_version", 1);
        result.put("project", root.toString());
        result.put("state", "beginning");
        result.put("phase", "initialization");
        result.put("message",
            "No event ledger exists; no compiled or executed project state is claimed.");
        return result;
    }

    private static Map<String, Object> deriveStatus(Path root,
            List<Map<String, Object>> events) throws IOException {
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("schema_version", 1);
        String state = "beginning", phase = "initialization";
        Map<String, Object> activePlan = null;
        List<Object> amendments = new ArrayList<>();
        List<Object> pending = new ArrayList<>();
        Map<String, Object> operationCounts = new LinkedHashMap<>();
        int completed = 0, reused = 0, waiting = 0, failed = 0, skipped = 0;
        java.util.Set<String> acceptedAssumptions = new java.util.HashSet<>();
        for (Map<String, Object> event : events) {
            Map<String, Object> payload = StudyIO.optionalMap(
                event.get("payload"), "event.payload");
            copy(payload, result, "project_id", "study_hash",
                "resolved_study_hash", "run_id");
            if (payload.get("state") instanceof String value) state = value;
            if (payload.get("phase") instanceof String value) phase = value;
            String type = String.valueOf(event.get("event_type"));
            if (type.equals("amendment_applied"))
                amendments.add(payload.get("amendment_id"));
            if (type.equals("amendment_applied") && activePlan != null) {
                activePlan.put("status", "superseded");
                pending.removeIf(item -> item instanceof Map<?, ?> map
                    && "plan_approval".equals(map.get("type")));
            }
            if (type.equals("plan_created")) {
                activePlan = new LinkedHashMap<>();
                copy(payload, activePlan, "plan_id", "plan_hash");
                activePlan.put("status", "proposed");
                pending.add(Map.of("id", "approve-" + payload.get("plan_id"),
                    "type", "plan_approval"));
            } else if (type.equals("plan_approved") && activePlan != null) {
                activePlan.put("status", "approved");
                pending.removeIf(item -> item instanceof Map<?, ?> map
                    && "plan_approval".equals(map.get("type")));
                for (Object id : StudyIO.optionalList(
                        payload.get("accepted_assumptions"), "accepted_assumptions"))
                    acceptedAssumptions.add(String.valueOf(id));
            } else if (type.equals("observation_recorded")) {
                pending.add(Map.of("id", payload.get("observation_id"),
                    "type", "observation_disposition"));
            } else if (type.equals("operation_completed")) completed++;
            else if (type.equals("operation_reused")) reused++;
            else if (type.equals("operation_waiting")) waiting++;
            else if (type.equals("operation_failed")) failed++;
            else if (type.equals("operation_skipped")) skipped++;
        }
        Map<String, Object> last = events.get(events.size() - 1);
        result.put("state", state);
        result.put("phase", phase);
        result.put("updated_at", last.get("timestamp"));
        result.put("last_event_id", last.get("event_id"));
        if (activePlan != null) result.put("active_plan", activePlan);
        result.put("applied_amendments", amendments);
        result.put("pending_decisions", pending);
        if (StudyQuestions.hasQuestions(events)) {
            Map<String, Object> questions = StudyQuestions.view(events);
            result.put("scientist_questions", Map.of("counts", questions.get("counts"),
                "notebook", "QUESTIONS.md", "state", "state/questions.yaml"));
        }
        operationCounts.put("completed", completed);
        operationCounts.put("reused", reused);
        operationCounts.put("waiting", waiting);
        operationCounts.put("failed", failed);
        if (skipped > 0) operationCounts.put("skipped", skipped);
        result.put("operations", operationCounts);
        if (Files.isRegularFile(root.resolve("resolved-study.yaml"))) {
            Map<String, Object> study = StudyIO.load(root.resolve("resolved-study.yaml"));
            Map<String, Object> hypotheses = new LinkedHashMap<>();
            for (Object value : StudyIO.list(study.get("hypotheses"), "hypotheses")) {
                Map<String, Object> hypothesis = StudyIO.map(value, "hypothesis");
                Map<String, Object> summary = new LinkedHashMap<>();
                summary.put("status", hypothesis.getOrDefault("status", "approved"));
                hypotheses.put(String.valueOf(hypothesis.get("id")), summary);
            }
            result.put("hypotheses", hypotheses);
            List<Object> assumptions = StudyIO.optionalList(
                study.get("assumptions"), "assumptions");
            long unresolved = assumptions.stream()
                .filter(item -> item instanceof Map<?, ?> map
                    && "pending".equals(map.get("review_status"))
                    && !acceptedAssumptions.contains(String.valueOf(map.get("id"))))
                .count();
            result.put("pending_assumptions", unresolved);
        }
        return result;
    }

    private static void copy(Map<String, Object> from, Map<String, Object> to,
            String... keys) {
        for (String key : keys)
            if (from.containsKey(key)) to.put(key, from.get(key));
    }

    private static String sessionHistory(List<Map<String, Object>> events) {
        StringBuilder result = new StringBuilder("# Session history\n\n");
        result.append("Generated from ledger/events.jsonl; do not edit this file to change project state.\n\n");
        for (Map<String, Object> event : events) {
            Map<String, Object> actor = StudyIO.map(event.get("actor"), "actor");
            Map<String, Object> payload = StudyIO.optionalMap(
                event.get("payload"), "payload");
            result.append("## ").append(event.get("event_id")).append(" — ")
                .append(event.get("event_type")).append("\n\n")
                .append("- Timestamp: ").append(event.get("timestamp")).append("\n")
                .append("- Actor: ").append(actor.get("type")).append(':')
                .append(actor.get("id"));
            if (actor.get("model") != null)
                result.append(" (").append(actor.get("model")).append(')');
            result.append("\n- Event hash: ").append(event.get("event_hash"))
                .append("\n");
            if (payload.get("summary") != null)
                result.append("- Summary: ").append(payload.get("summary")).append("\n");
            result.append('\n');
        }
        if (StudyQuestions.hasQuestions(events))
            result.append("## Documented scientist questions\n\nSee QUESTIONS.md for exact prompts, reasons, answer options, evidence, answers, and decisions.\n\n")
                .append(StudyQuestions.notebook(StudyQuestions.view(events)).replaceAll("(?m)^#", "##"));
        return result.toString();
    }

    private static String agents(Map<String, Object> study,
            Map<String, Object> status) {
        Map<String, Object> project = StudyIO.map(study.get("project"), "project");
        StringBuilder result = new StringBuilder("# Project agent contract\n\n");
        result.append("This file is generated from the approved resolved study and project ledger.\n")
            .append("It is context for an agent, not authority to alter scientific state.\n\n")
            .append("## Project\n\n")
            .append("- ID: ").append(project.get("id")).append("\n")
            .append("- Title: ").append(project.get("title")).append("\n")
            .append("- Objective: ").append(project.get("objective")).append("\n")
            .append("- State: ").append(status.get("state")).append("\n\n")
            .append("## Active hypotheses\n\n");
        for (Object value : StudyIO.list(study.get("hypotheses"), "hypotheses")) {
            Map<String, Object> hypothesis = StudyIO.map(value, "hypothesis");
            if (!"retired".equals(hypothesis.get("status")))
                result.append("- ").append(hypothesis.get("id")).append(": ")
                    .append(hypothesis.get("description")).append("\n");
        }
        result.append("\n## Rules\n\n")
            .append("- Do not modify resolved-study.yaml, the event ledger, status, or session history manually.\n")
            .append("- Propose scientific changes through an amendment-XXX.yaml file.\n")
            .append("- Do not execute an unapproved plan or send restricted data to a remote service.\n")
            .append("- Cite stable result IDs; never invent or transcribe numerical results without lineage.\n")
            .append("- Record manual, JLinAlg, external-software, and LLM actions with their actual actor type.\n");
        return result.toString();
    }

    private static String version() {
        Package pkg = StudyProject.class.getPackage();
        return pkg.getImplementationVersion() == null
            ? "development" : pkg.getImplementationVersion();
    }
}
