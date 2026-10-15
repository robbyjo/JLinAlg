/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.UUID;
import java.util.regex.Pattern;

/** Durable local scientist questions. Answers never authorize analysis or egress. */
public final class StudyQuestions {
    private StudyQuestions() { }
    private static final Pattern ID = Pattern.compile("[A-Za-z][A-Za-z0-9._-]{0,127}");
    private static final Set<String> TERMINAL = Set.of("resolved", "superseded", "withdrawn");
    private static final StudyProject.Actor SYSTEM = new StudyProject.Actor("jlinalg", "question-coordinator", null);

    /** Read-only: no provider calls, question generation, or view writes. */
    public static Map<String, Object> read(Path project) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        List<Map<String, Object>> events = StudyLedger.read(root);
        verifyArtifacts(root, events);
        return view(events);
    }

    static Map<String, Object> view(List<Map<String, Object>> events) {
        Map<String, Map<String, Object>> records = new LinkedHashMap<>();
        String projectId = "unresolved";
        for (Map<String, Object> event : events) {
            Map<String, Object> payload = map(event.get("payload"));
            if (payload.get("project_id") != null) projectId = String.valueOf(payload.get("project_id"));
            if ("questions_recorded".equals(event.get("event_type"))) {
                for (Object value : list(payload.get("questions"))) {
                    Map<String, Object> question = map(value);
                    question.put("state", "open");
                    question.put("documented_at", event.get("timestamp"));
                    question.put("documented_by", event.get("actor"));
                    question.put("event_id", event.get("event_id"));
                    question.put("answers", new ArrayList<>());
                    records.put(key(question), question);
                }
            }
            for (Object value : list(payload.get("question_updates"))) {
                Map<String, Object> update = map(value);
                Map<String, Object> question = records.get(key(update));
                if (question == null) throw new IllegalArgumentException("Question update has no originating record");
                question.put("state", update.get("state"));
                question.put("state_reason", update.get("reason"));
                question.put("updated_at", event.get("timestamp"));
                if (update.get("amendment_id") != null) question.put("applied_amendment", update.get("amendment_id"));
            }
            if ("scientist_answers_recorded".equals(event.get("event_type"))) {
                for (Object value : list(payload.get("responses"))) {
                    Map<String, Object> response = map(value);
                    Map<String, Object> question = records.get(key(response));
                    if (question == null) throw new IllegalArgumentException("Answer has no originating question");
                    response.put("answer_id", payload.get("answer_id"));
                    response.put("recorded_at", event.get("timestamp"));
                    response.put("actor", event.get("actor"));
                    response.put("artifact", payload.get("answer_file"));
                    if (payload.get("amendment") != null) response.put("amendment", payload.get("amendment"));
                    List<Object> answers = list(question.get("answers"));
                    answers.add(response);
                    question.put("answers", answers);
                    question.put("latest_answer", response);
                    question.put("state", response.get("resulting_state"));
                    question.put("state_reason", response.get("decision"));
                    question.put("updated_at", event.get("timestamp"));
                }
            }
        }
        Map<String, Object> counts = new TreeMap<>();
        for (Map<String, Object> question : records.values())
            counts.merge(String.valueOf(question.get("state")), 1, (a, b) -> ((Number) a).intValue() + 1);
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("schema_version", 1);
        result.put("project_id", projectId);
        result.put("questions", new ArrayList<>(records.values()));
        result.put("counts", counts);
        result.put("authority", "Local attribution only; answers do not grant execution or transmission permission.");
        return result;
    }

    /** Materialize typed planner/runtime needs, preserving every historical revision. */
    static void reconcile(Path project, StudyPlan plan, String runId) throws IOException {
        Map<String, Object> study = StudyProject.resolved(project).document();
        String studyHash = StudyProject.resolvedHash(study);
        List<Map<String, Object>> events = StudyLedger.read(project);
        verifyArtifacts(project, events);
        Map<String, Map<String, Object>> latest = latest(view(events));
        Map<String, Map<String, Object>> candidates = new TreeMap<>();
        Map<String, String> datasetHashes = new TreeMap<>();
        for (Object raw : list(study.get("datasets"))) {
            Map<String, Object> dataset = map(raw);
            for (Map.Entry<String, Object> file : map(dataset.get("files")).entrySet()) {
                Path path = project.resolve(String.valueOf(file.getValue())).normalize();
                if (Files.isRegularFile(path)) datasetHashes.put("dataset:" + dataset.get("id") + ":" + file.getKey(), StudyIO.hash(path));
            }
        }
        for (Object raw : list(plan.document().get("operations"))) {
            Map<String, Object> operation = map(raw);
            List<StudyNeed> needs = new ArrayList<>(StudyQuestionRules.prerequisites(study, operation));
            Map<String, Object> inputHashes = new TreeMap<>();
            for (String id : StudyQuestionRules.strings(operation.get("datasets")))
                datasetHashes.forEach((name, hash) -> { if (name.startsWith("dataset:" + id + ":")) inputHashes.put(name, hash); });
            if (runId != null) {
                Path manifest = project.resolve("runs").resolve(runId).resolve("operations")
                    .resolve(String.valueOf(operation.get("id"))).resolve("manifest.yaml");
                if (Files.isRegularFile(manifest)) {
                    Map<String, Object> record = StudyIO.loadArtifact(manifest);
                    for (Object value : list(map(record.get("details")).get("needs")))
                        needs.add(StudyNeed.from(map(value)));
                }
            }
            for (StudyNeed need : needs) {
                Map<String, Object> candidate = candidates.computeIfAbsent(need.key(), ignored -> {
                    Map<String, Object> value = new LinkedHashMap<>();
                    value.put("need", need.document());
                    value.put("hypotheses", new java.util.TreeSet<String>());
                    value.put("operations", new java.util.TreeSet<String>());
                    value.put("datasets", new java.util.TreeSet<String>());
                    value.put("input_sha256", new TreeMap<String, Object>());
                    return value;
                });
                if (!StudyIO.canonicalJson(candidate.get("need")).equals(StudyIO.canonicalJson(need.document())))
                    throw new IOException("Conflicting question contracts for " + need.key());
                add(candidate, "hypotheses", String.valueOf(operation.get("hypothesis")));
                add(candidate, "operations", String.valueOf(operation.get("id")));
                for (String id : StudyQuestionRules.strings(operation.get("datasets"))) add(candidate, "datasets", id);
                @SuppressWarnings("unchecked") Map<String, Object> hashes = (Map<String, Object>) candidate.get("input_sha256");
                // Dataset hashes are stable across run IDs. Conservative study hashing remains intentional.
                inputHashes.forEach((name, hash) -> { if (name.startsWith("dataset:")) hashes.put(name, hash); });
            }
        }
        List<Object> created = new ArrayList<>();
        List<Object> updates = new ArrayList<>();
        Map<String, String> files = new LinkedHashMap<>();
        for (Map.Entry<String, Map<String, Object>> entry : candidates.entrySet()) {
            Map<String, Object> candidate = StudyIO.copyMap(entry.getValue());
            candidate.put("resolved_study_hash", studyHash);
            String fingerprint = StudyIO.hash(candidate);
            Map<String, Object> prior = latest.get(entry.getKey());
            if (prior != null && fingerprint.equals(prior.get("premises_hash"))
                    && !TERMINAL.contains(prior.get("state"))) continue;
            int revision = prior == null ? 1 : ((Number) prior.get("revision")).intValue() + 1;
            if (prior != null && !TERMINAL.contains(prior.get("state")))
                updates.add(update(prior, "superseded", "The question's scientific premises or scope changed."));
            Map<String, Object> question = new LinkedHashMap<>();
            question.put("schema_version", 1);
            question.put("project_id", map(study.get("project")).get("id"));
            question.put("question_id", "Q-" + StudyIO.hash(entry.getKey()).substring(0, 20));
            question.put("revision", revision);
            question.put("semantic_key", entry.getKey());
            question.put("premises_hash", fingerprint);
            question.put("resolved_study_hash", studyHash);
            question.put("plan_id", plan.id());
            question.put("need", candidate.get("need"));
            question.put("scope", Map.of("hypotheses", candidate.get("hypotheses"),
                "operations", candidate.get("operations"), "datasets", candidate.get("datasets")));
            question.put("input_sha256", candidate.get("input_sha256"));
            List<Object> evidence = new ArrayList<>();
            if (runId != null) for (Object operation : list(candidate.get("operations"))) {
                String path = "runs/" + runId + "/operations/" + operation + "/manifest.yaml";
                if (Files.isRegularFile(project.resolve(path)))
                    evidence.add(Map.of("manifest", path, "sha256", StudyIO.hash(project.resolve(path))));
            }
            question.put("evidence", evidence);
            question.put("question_hash", StudyIO.hash(question));
            String path = "questions/" + question.get("question_id") + "/revision-" + String.format("%04d", revision) + ".yaml";
            files.put(path, StudyIO.yaml(question));
            Map<String, Object> recorded = StudyIO.copyMap(question);
            recorded.put("artifact", path);
            created.add(recorded);
        }
        for (Map.Entry<String, Map<String, Object>> entry : latest.entrySet()) {
            Map<String, Object> prior = entry.getValue();
            if (candidates.containsKey(entry.getKey()) || TERMINAL.contains(prior.get("state"))) continue;
            if (!studyHash.equals(prior.get("resolved_study_hash")))
                updates.add(update(prior, "superseded", "The approved study changed; reassess under the new plan."));
            else if (runId != null && scopeCompleted(project, runId, prior))
                updates.add(update(prior, "withdrawn", "The current run no longer reports this need; this is not evidence of a biological negative."));
        }
        if (created.isEmpty() && updates.isEmpty()) return;
        Map<String, Object> payload = new LinkedHashMap<>();
        payload.put("questions", created);
        payload.put("question_updates", updates);
        payload.put("summary", "Documented " + created.size() + " scientist questions and " + updates.size() + " lifecycle updates");
        commit(project, "questions_recorded", SYSTEM, payload, files);
    }

    private static boolean scopeCompleted(Path project, String runId, Map<String, Object> question) throws IOException {
        for (Object operation : list(map(question.get("scope")).get("operations"))) {
            Path manifest = project.resolve("runs").resolve(runId).resolve("operations").resolve(String.valueOf(operation)).resolve("manifest.yaml");
            if (!Files.isRegularFile(manifest) || !Set.of("complete", "not_applicable").contains(StudyIO.loadArtifact(manifest).get("state"))) return false;
        }
        return true;
    }

    @SuppressWarnings("unchecked") private static void add(Map<String, Object> map, String key, String value) {
        ((Set<String>) map.get(key)).add(value);
    }

    /** Export blank responses. Nothing is preaccepted, and existing files are protected. */
    public static void export(Path project, Path output, boolean includeAnswered) throws IOException {
        Map<String, Object> view = read(project);
        List<Object> responses = new ArrayList<>();
        for (Map<String, Object> question : latest(view).values()) {
            if (TERMINAL.contains(question.get("state")) || !includeAnswered && !"open".equals(question.get("state"))) continue;
            Map<String, Object> response = new LinkedHashMap<>();
            response.put("question_id", question.get("question_id"));
            response.put("question_revision", question.get("revision"));
            response.put("question_hash", question.get("question_hash"));
            response.put("disposition", null);
            response.put("value", null);
            response.put("rationale", null);
            response.put("evidence", List.of());
            if (question.containsKey("latest_answer")) response.put("supersedes", map(question.get("latest_answer")).get("answer_id"));
            responses.add(response);
        }
        Map<String, Object> template = new LinkedHashMap<>();
        template.put("schema_version", 1);
        template.put("id", "A-" + UUID.randomUUID());
        template.put("project_id", view.get("project_id"));
        template.put("amendment_id", null);
        template.put("responses", responses);
        Path destination = output.toAbsolutePath().normalize();
        // Export outside the project is allowed, but never over an existing file.
        if (destination.getParent() != null) Files.createDirectories(destination.getParent());
        Files.writeString(destination, StudyIO.yaml(template), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
    }

    /** Validate a complete answer batch without any write or external request. */
    public static Map<String, Object> preview(Path project, Path answers, StudyProject.Actor actor) throws IOException {
        return prepare(project.toAbsolutePath().normalize(), answers, actor).receipt();
    }

    /** Atomically accept a validated batch; scientific changes remain proposals. */
    public static Map<String, Object> record(Path project, Path answers, StudyProject.Actor actor) throws IOException {
        return StudyProjectLock.with(project, () -> {
            Prepared prepared = prepare(project.toAbsolutePath().normalize(), answers, actor);
            if (prepared.duplicate()) return prepared.receipt();
            commit(project, "scientist_answers_recorded", actor, prepared.receipt(), prepared.files());
            StudyEvidenceIndex.rebuild(project);
            return prepared.receipt();
        });
    }

    private record Prepared(Map<String, Object> receipt, Map<String, String> files, boolean duplicate) { }

    private static Prepared prepare(Path root, Path input, StudyProject.Actor actor) throws IOException {
        if (!"human".equals(actor.type()) || "unspecified".equals(actor.id()))
            throw new IllegalArgumentException("Scientist answers require an explicit local human actor; model output is not approval");
        if (Files.size(input) > 1_048_576) throw new IOException("Answer file exceeds 1 MiB");
        StudyProject.verify(root);
        String imported = Files.readString(input);
        Map<String, Object> raw = StudyIO.map(StudyIO.reader().load(imported), "answers");
        fields(raw, Set.of("schema_version", "id", "project_id", "amendment_id", "meeting_id", "responses"), "answer file");
        if (!Integer.valueOf(1).equals(raw.get("schema_version"))) throw new IllegalArgumentException("Answers require schema_version: 1");
        String id = identifier(raw, "id");
        String requestHash = StudyIO.hash(raw);
        List<Map<String, Object>> events = StudyLedger.read(root);
        for (Map<String, Object> event : events) {
            Map<String, Object> payload = map(event.get("payload"));
            if (!"scientist_answers_recorded".equals(event.get("event_type")) || !id.equals(payload.get("answer_id"))) continue;
            Map<String, Object> recordedActor = map(event.get("actor"));
            if (!requestHash.equals(payload.get("request_hash")) || !actor.id().equals(recordedActor.get("id"))
                    || !actor.type().equals(recordedActor.get("type")) || !java.util.Objects.equals(actor.model(), recordedActor.get("model")))
                throw new IllegalArgumentException("Answer ID already exists with different content or attribution: " + id);
            return new Prepared(payload, Map.of(), true);
        }
        Map<String, Object> study = StudyProject.resolved(root).document();
        String studyHash = StudyProject.resolvedHash(study);
        if (!map(study.get("project")).get("id").equals(raw.get("project_id")))
            throw new IllegalArgumentException("Answer project_id does not match this project");
        Map<String, Map<String, Object>> current = new LinkedHashMap<>();
        latest(view(events)).values().forEach(question -> current.put(String.valueOf(question.get("question_id")), question));
        List<Object> responses = list(raw.get("responses"));
        if (responses.isEmpty() || responses.size() > 1000) throw new IllegalArgumentException("Supply between 1 and 1000 responses");
        Set<String> seen = new java.util.HashSet<>();
        Map<String, String> files = new LinkedHashMap<>();
        List<Object> accepted = new ArrayList<>();
        List<Object> bindings = new ArrayList<>();
        Map<String, Object> amended = StudyIO.copyMap(study);
        Map<String, Map<String, Object>> datasetPatches = new LinkedHashMap<>();
        Map<String, Map<String, Object>> hypothesisPatches = new LinkedHashMap<>();
        List<Object> relationshipAdds = new ArrayList<>();
        for (Object value : responses) {
            Map<String, Object> response = map(value);
            fields(response, Set.of("question_id", "question_revision", "question_hash", "disposition", "value", "rationale", "evidence", "supersedes"), "response");
            String questionId = identifier(response, "question_id");
            if (!seen.add(questionId)) throw new IllegalArgumentException("Multiple answers to the same question in one batch: " + questionId);
            Map<String, Object> question = current.get(questionId);
            if (question == null || TERMINAL.contains(question.get("state"))
                    || !question.get("revision").equals(response.get("question_revision"))
                    || !question.get("question_hash").equals(response.get("question_hash"))
                    || !studyHash.equals(question.get("resolved_study_hash")))
                throw new IllegalArgumentException("Stale or unknown question: " + questionId + "; replan/run and export current questions");
            checkInputPremises(root, question);
            Map<String, Object> previous = map(question.get("latest_answer"));
            if (!previous.isEmpty() && !previous.get("answer_id").equals(response.get("supersedes")))
                throw new IllegalArgumentException("A recorded answer already exists; explicitly name supersedes: " + previous.get("answer_id"));
            if (previous.isEmpty() && response.get("supersedes") != null)
                throw new IllegalArgumentException("Cannot supersede an absent answer");
            String disposition = StudyIO.requireString(response, "disposition", "response");
            if (!Set.of("answer", "unknown", "defer", "decline").contains(disposition))
                throw new IllegalArgumentException("disposition must be answer, unknown, defer, or decline");
            String rationale = StudyIO.requireString(response, "rationale", "response");
            bounded(rationale);
            Map<String, Object> need = map(question.get("need"));
            if (disposition.equals("decline") && Boolean.TRUE.equals(need.get("blocking")))
                throw new IllegalArgumentException("Cannot decline a required scientific prerequisite; use unknown or defer");
            Map<String, Object> normalized = new LinkedHashMap<>();
            normalized.put("question_id", questionId);
            normalized.put("revision", question.get("revision"));
            normalized.put("question_hash", question.get("question_hash"));
            normalized.put("prompt", need.get("prompt"));
            normalized.put("disposition", disposition);
            normalized.put("rationale", rationale);
            if (!previous.isEmpty()) normalized.put("supersedes", response.get("supersedes"));
            List<Object> evidence = evidence(root, response, files);
            normalized.put("evidence", evidence);
            String state;
            if (disposition.equals("answer")) {
                String supplied = StudyIO.requireString(response, "value", "response");
                bounded(supplied);
                Map<String, Object> contract = map(need.get("response"));
                if ("single_choice".equals(contract.get("type")) && !list(contract.get("options")).contains(supplied))
                    throw new IllegalArgumentException("Answer is not one of the declared choices: " + questionId);
                if ("local_file".equals(contract.get("type"))) {
                    Path path = localFile(root, supplied);
                    supplied = root.relativize(path).toString().replace('\\', '/');
                    normalized.put("resource_sha256", StudyIO.hash(path));
                }
                normalized.put("value", supplied);
                if (!map(need.get("target")).isEmpty()) {
                    if (evidence.isEmpty()) throw new IllegalArgumentException("A factual answer needs supporting evidence or a source note: " + questionId);
                    applyAnswer(amended, question, supplied, datasetPatches, hypothesisPatches, relationshipAdds);
                    bindings.add(Map.of("question_id", questionId, "revision", question.get("revision"), "question_hash", question.get("question_hash")));
                    state = "awaiting_approval";
                    normalized.put("decision", "Propose amendment only; original inputs and execution permissions are unchanged.");
                } else {
                    state = "needs_clarification";
                    normalized.put("decision", "Scientist response recorded for review; supply a supported amendment or missing input. No arbitrary text is executed.");
                }
            } else {
                if (response.get("value") != null) throw new IllegalArgumentException("Non-answer dispositions must not include a value");
                state = disposition.equals("unknown") ? "needs_clarification" : disposition.equals("defer") ? "deferred" : "declined";
                normalized.put("decision", disposition.equals("decline") ? "Optional investigation not pursued; not a biological negative."
                    : "Prerequisite remains unresolved; no permission granted. Revisit explicitly or when premises change.");
            }
            normalized.put("resulting_state", state);
            accepted.add(normalized);
        }
        Map<String, Object> changes = new LinkedHashMap<>();
        Map<String, Object> updates = new LinkedHashMap<>();
        if (!datasetPatches.isEmpty()) updates.put("datasets", new ArrayList<>(datasetPatches.values()));
        if (!hypothesisPatches.isEmpty()) updates.put("hypotheses", new ArrayList<>(hypothesisPatches.values()));
        if (!updates.isEmpty()) changes.put("update", updates);
        if (!relationshipAdds.isEmpty()) changes.put("add", Map.of("cohort_relationships", relationshipAdds));
        Map<String, Object> receipt = new LinkedHashMap<>();
        receipt.put("answer_id", id);
        receipt.put("request_hash", requestHash);
        receipt.put("responses", accepted);
        receipt.put("summary", "Recorded scientist answer batch " + id + "; no analysis or remote transmission authorized");
        if (raw.get("meeting_id") != null) { bounded(String.valueOf(raw.get("meeting_id"))); receipt.put("meeting_id", raw.get("meeting_id")); }
        if (!changes.isEmpty()) {
            String amendmentId = identifier(raw, "amendment_id");
            StudyValidation validation = StudyValidator.validateDocument(amended, root, false);
            if (!validation.valid()) throw new IllegalArgumentException("Answer would produce an invalid study: " + validation.report());
            Map<String, Object> header = new LinkedHashMap<>();
            header.put("id", amendmentId);
            header.put("title", "Scientist answers " + id);
            header.put("base_study_hash", studyHash);
            header.put("proposed_by", actor.asMap());
            header.put("rationale", "Proposed from recorded scientist answers; review evidence, effect on hypotheses, and the resulting plan before approval.");
            header.put("source", Map.of("type", "scientist_answers", "answer_id", id, "questions", bindings));
            Map<String, Object> amendment = Map.of("schema_version", 1, "amendment", header, "changes", changes);
            String path = "amendments/proposed/amendment-" + amendmentId + ".yaml";
            if (Files.exists(root.resolve("amendments/amendment-" + amendmentId + ".yaml")))
                throw new IllegalArgumentException("Amendment ID is already applied: " + amendmentId);
            files.put(path, StudyIO.yaml(amendment));
            receipt.put("amendment", Map.of("id", amendmentId, "path", path, "document_hash", StudyIO.hash(amendment)));
        } else if (raw.get("amendment_id") != null) {
            throw new IllegalArgumentException("No supported contract changes in this batch; omit amendment_id");
        }
        String path = "answers/" + id + ".yaml";
        files.put(path, StudyIO.yaml(receipt));
        files.put("answers/imports/" + id + ".yaml", imported);
        receipt.put("answer_file", path);
        // Check all destination collisions before publishing anything, including orphaned retries.
        for (Map.Entry<String, String> file : files.entrySet()) checkImmutable(root, file.getKey(), file.getValue());
        return new Prepared(receipt, files, false);
    }

    private static void applyAnswer(Map<String, Object> study, Map<String, Object> question,
            String value, Map<String, Map<String, Object>> datasets,
            Map<String, Map<String, Object>> hypotheses, List<Object> relationships) {
        Map<String, Object> target = map(map(question.get("need")).get("target"));
        String collection = String.valueOf(target.get("collection"));
        if (collection.equals("cohort_relationships")) {
            Map<String, Object> relation = Map.of("id", "relationship-" + question.get("question_id"),
                "between", target.get("between"), "type", value);
            relationships.add(relation);
            List<Object> existing = list(study.get("cohort_relationships"));
            existing.add(relation);
            study.put("cohort_relationships", existing);
            return;
        }
        String id = String.valueOf(target.get("id"));
        String field = String.valueOf(target.get("field"));
        Set<String> allowed = collection.equals("datasets")
            ? Set.of("genome_build", "coordinate_system", "files.annotation", "files.phenotype", "units")
            : collection.equals("hypotheses") ? Set.of("outcome") : Set.of();
        if (!allowed.contains(field)) throw new IllegalArgumentException("Unsupported answer target: " + collection + "." + field);
        Map<String, Map<String, Object>> patches = collection.equals("datasets") ? datasets : hypotheses;
        Map<String, Object> patch = patches.computeIfAbsent(id, ignored -> new LinkedHashMap<>(Map.of("id", id)));
        setField(patch, field, value);
        List<Object> objects = list(study.get(collection));
        for (int i = 0; i < objects.size(); i++) {
            Map<String, Object> object = map(objects.get(i));
            if (id.equals(object.get("id"))) { setField(object, field, value); objects.set(i, object); }
        }
        study.put(collection, objects);
    }

    private static void setField(Map<String, Object> object, String field, String value) {
        String[] parts = field.split("\\.");
        if (parts.length == 1) object.put(field, value);
        else { Map<String, Object> nested = map(object.get(parts[0])); nested.put(parts[1], value); object.put(parts[0], nested); }
    }

    private static List<Object> evidence(Path root, Map<String, Object> response, Map<String, String> files) throws IOException {
        List<Object> result = new ArrayList<>();
        for (Object item : list(response.get("evidence"))) {
            Map<String, Object> evidence = map(item);
            fields(evidence, Set.of("file", "selector", "citation", "note"), "evidence");
            if (evidence.isEmpty()) throw new IllegalArgumentException("Empty evidence reference");
            for (Object value : evidence.values()) {
                if (!(value instanceof String text)) throw new IllegalArgumentException("Evidence fields must be nonblank strings");
                bounded(text);
            }
            if (evidence.get("file") != null) {
                Path source = localFile(root, String.valueOf(evidence.get("file")));
                if (Files.size(source) > 1_048_576) throw new IOException("Evidence snapshots must be UTF-8 metadata of at most 1 MiB");
                String hash = StudyIO.hash(source);
                String snapshot = "answers/evidence/" + hash + ".txt";
                files.put(snapshot, Files.readString(source));
                evidence.put("sha256", hash);
                evidence.put("snapshot", snapshot);
            } else if (evidence.get("citation") == null && evidence.get("note") == null) {
                throw new IllegalArgumentException("Evidence needs a local file, citation, or explicit source note");
            }
            evidence.put("verification", "scientist_supplied; citation claims are not automatically verified");
            result.add(evidence);
        }
        return result;
    }

    private static void checkInputPremises(Path root, Map<String, Object> question) throws IOException {
        Map<String, Object> study = StudyProject.resolved(root).document();
        for (Map.Entry<String, Object> entry : map(question.get("input_sha256")).entrySet()) {
            String[] parts = entry.getKey().split(":", 3);
            if (parts.length != 3 || !parts[0].equals("dataset")) throw new IOException("Invalid question input reference");
            Map<String, Object> dataset = StudyQuestionRules.object(study, "datasets", parts[1]);
            Object file = map(dataset.get("files")).get(parts[2]);
            if (file == null || !Files.isRegularFile(root.resolve(String.valueOf(file)))
                    || !entry.getValue().equals(StudyIO.hash(root.resolve(String.valueOf(file)))))
                throw new IllegalArgumentException("Question inputs changed; rerun before answering " + question.get("question_id"));
        }
    }

    /** Validate a proposed amendment's attribution and reject superseded answers. */
    static void validateAmendment(Path root, Map<String, Object> amendment) throws IOException {
        Map<String, Object> source = map(map(amendment.get("amendment")).get("source"));
        if (!"scientist_answers".equals(source.get("type"))) return;
        List<Map<String, Object>> events = StudyLedger.read(root);
        verifyArtifacts(root, events);
        Map<String, Map<String, Object>> current = new LinkedHashMap<>();
        latest(view(events)).values().forEach(question -> current.put(String.valueOf(question.get("question_id")), question));
        Map<String, Object> receipt = events.stream().filter(event -> "scientist_answers_recorded".equals(event.get("event_type")))
            .map(event -> map(event.get("payload"))).filter(payload -> source.get("answer_id").equals(payload.get("answer_id")))
            .findFirst().orElseThrow(() -> new IllegalArgumentException("Amendment references an unrecorded answer"));
        if (!StudyIO.hash(amendment).equals(map(receipt.get("amendment")).get("document_hash")))
            throw new IllegalArgumentException("Answer-derived amendment changed; record a new reviewed answer/proposal");
        for (Object value : list(source.get("questions"))) {
            Map<String, Object> binding = map(value);
            Map<String, Object> question = current.get(String.valueOf(binding.get("question_id")));
            if (question == null || !"awaiting_approval".equals(question.get("state"))
                    || !binding.get("question_hash").equals(question.get("question_hash"))
                    || !source.get("answer_id").equals(map(question.get("latest_answer")).get("answer_id")))
                throw new IllegalArgumentException("Amendment is based on a stale or superseded scientist answer");
            checkInputPremises(root, question);
            Map<String, Object> answer = map(question.get("latest_answer"));
            if (answer.get("resource_sha256") != null
                    && !answer.get("resource_sha256").equals(StudyIO.hash(localFile(root, String.valueOf(answer.get("value"))))))
                throw new IllegalArgumentException("Answered local resource changed before amendment approval");
        }
    }

    static List<Object> amendmentUpdates(Map<String, Object> amendment) {
        Map<String, Object> header = map(amendment.get("amendment"));
        Map<String, Object> source = map(header.get("source"));
        List<Object> result = new ArrayList<>();
        if ("scientist_answers".equals(source.get("type"))) for (Object value : list(source.get("questions"))) {
            Map<String, Object> item = map(value);
            item.put("state", "resolved");
            item.put("reason", "Scientist-approved amendment applied; downstream validation and plan approval remain required.");
            item.put("amendment_id", header.get("id"));
            result.add(item);
        }
        return result;
    }

    static boolean hasQuestions(List<Map<String, Object>> events) {
        return events.stream().anyMatch(event -> "questions_recorded".equals(event.get("event_type")));
    }

    static void materialize(Path root, List<Map<String, Object>> events) throws IOException {
        if (!hasQuestions(events)) return;
        Map<String, Object> view = view(events);
        StudyIO.writeYaml(root.resolve("state/questions.yaml"), view);
        StudyIO.atomicWrite(root.resolve("QUESTIONS.md"), notebook(view));
    }

    static String notebook(Map<String, Object> view) {
        StringBuilder text = new StringBuilder("# Scientist questions and answers\n\n")
            .append("Generated from ledger/events.jsonl. Questions are documented for review; this does not assert they were read.\n")
            .append("Answers do not grant execution or disclosure permission. Keep this notebook local unless sharing is authorized.\n\n");
        for (Object raw : list(view.get("questions"))) {
            Map<String, Object> question = map(raw);
            Map<String, Object> need = map(question.get("need"));
            text.append("## ").append(question.get("question_id")).append(" revision ").append(question.get("revision")).append("\n\n")
                .append("- State: ").append(question.get("state")).append("\n")
                .append("- Documented at: ").append(question.get("documented_at")).append("\n")
                .append("- Documented by: ").append(md(String.valueOf(question.get("documented_by")))).append("\n")
                .append("- Record: ").append(question.get("artifact")).append("\n")
                .append("- Question hash: ").append(question.get("question_hash")).append("\n")
                .append("- Scope: ").append(md(String.valueOf(question.get("scope")))).append("\n\n")
                .append(md(String.valueOf(need.get("prompt")))).append("\n\n")
                .append("Why: ").append(md(String.valueOf(need.get("why")))).append("\n\n")
                .append("Response contract: ").append(md(String.valueOf(need.get("response")))).append("\n\n")
                .append("Consequences: ").append(md(String.valueOf(need.get("consequences")))).append("\n\n")
                .append("Evidence input hashes: ").append(md(String.valueOf(question.get("input_sha256")))).append("\n\n");
            text.append("Evidence manifests: ").append(md(String.valueOf(question.get("evidence")))).append("\n\n");
            if (question.containsKey("state_reason")) text.append("Disposition: ").append(md(String.valueOf(question.get("state_reason")))).append("\n\n");
            for (Object rawAnswer : list(question.get("answers"))) {
                Map<String, Object> answer = map(rawAnswer);
                text.append("### Answer ").append(answer.get("answer_id")).append("\n\n");
                for (Map.Entry<String, Object> field : answer.entrySet())
                    text.append("- ").append(field.getKey()).append(": ").append(md(String.valueOf(field.getValue()))).append('\n');
                text.append('\n');
            }
        }
        return text.toString();
    }

    static void verifyArtifacts(Path root, List<Map<String, Object>> events) throws IOException {
        for (Map<String, Object> event : events) {
            if (!Set.of("questions_recorded", "scientist_answers_recorded").contains(event.get("event_type"))) continue;
            Map<String, Object> payload = map(event.get("payload"));
            for (Object raw : list(payload.get("question_artifacts"))) {
                Map<String, Object> artifact = map(raw);
                Path path = artifactPath(root, String.valueOf(artifact.get("path")));
                if (!Files.isRegularFile(path) || !artifact.get("sha256").equals(StudyIO.hash(path)))
                    throw new IOException("Question/answer artifact is missing or modified: " + path);
            }
        }
    }

    private static void commit(Path project, String type, StudyProject.Actor actor,
            Map<String, Object> payload, Map<String, String> files) throws IOException {
        List<Object> artifacts = new ArrayList<>();
        for (Map.Entry<String, String> file : files.entrySet()) checkImmutable(project, file.getKey(), file.getValue());
        for (Map.Entry<String, String> file : files.entrySet()) {
            Path path = artifactPath(project, file.getKey());
            if (!Files.exists(path)) StudyIO.atomicWrite(path, file.getValue());
            artifacts.add(Map.of("path", file.getKey(), "sha256", StudyIO.hash(path)));
        }
        payload.put("question_artifacts", artifacts);
        StudyProject.append(project, type, actor, payload);
    }

    private static void checkImmutable(Path root, String relative, String text) throws IOException {
        Path path = artifactPath(root, relative);
        if (Files.exists(path) && (!Files.isRegularFile(path) || !Files.readString(path).equals(text)))
            throw new IOException("Refusing to replace an immutable or conflicting question artifact: " + relative);
    }

    private static Path artifactPath(Path root, String relative) throws IOException {
        Path absoluteRoot = root.toAbsolutePath().normalize();
        Path path = absoluteRoot.resolve(relative).normalize();
        if (Path.of(relative).isAbsolute() || !path.startsWith(absoluteRoot) || path.equals(absoluteRoot))
            throw new IOException("Unsafe question artifact path");
        for (Path parent = path; parent != null && parent.startsWith(absoluteRoot); parent = parent.getParent())
            if (Files.exists(parent, LinkOption.NOFOLLOW_LINKS) && !parent.toRealPath().startsWith(absoluteRoot.toRealPath()))
                throw new IOException("Question artifact escapes the project through a link");
        return path;
    }

    private static Path localFile(Path root, String relative) throws IOException {
        Path path = artifactPath(root, relative);
        if (!Files.isRegularFile(path)) throw new IOException("Local resource is absent: " + relative);
        return path;
    }

    private static Map<String, Map<String, Object>> latest(Map<String, Object> view) {
        Map<String, Map<String, Object>> result = new LinkedHashMap<>();
        for (Object raw : list(view.get("questions"))) {
            Map<String, Object> question = map(raw);
            result.put(String.valueOf(question.get("semantic_key")), question);
        }
        return result;
    }

    private static Map<String, Object> update(Map<String, Object> question, String state, String reason) {
        return Map.of("question_id", question.get("question_id"), "revision", question.get("revision"), "state", state, "reason", reason);
    }
    private static String key(Map<String, Object> value) { return value.get("question_id") + ":" + value.get("revision"); }
    private static String identifier(Map<String, Object> value, String key) {
        String result = StudyIO.requireString(value, key, "answers");
        if (!ID.matcher(result).matches()) throw new IllegalArgumentException("Invalid identifier: " + key);
        return result;
    }
    private static void bounded(String value) {
        if (value.isBlank() || value.length() > 8000) throw new IllegalArgumentException("Answer text must be nonblank and at most 8000 characters");
    }
    private static void fields(Map<String, Object> value, Set<String> allowed, String name) {
        for (String key : value.keySet()) if (!allowed.contains(key)) throw new IllegalArgumentException("Unknown " + name + " field: " + key);
    }
    private static String md(String value) {
        return value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace("\\", "\\\\").replace("`", "\\`").replace("[", "\\[")
            .replace("]", "\\]").replace("*", "\\*").replace("_", "\\_")
            .replace("\r", " ").replace("\n", " ");
    }
    private static Map<String, Object> map(Object value) { return StudyIO.optionalMap(value, "question record"); }
    private static List<Object> list(Object value) { return StudyIO.optionalList(value, "question list"); }
}
