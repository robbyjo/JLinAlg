/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Approved-plan execution with verified resume, content caching, and manifests. */
public final class StudyExecutionEngine {
    private StudyExecutionEngine() { }

    public record RunResult(String runId, String state,
            int completed, int reused, int waiting, int failed, int skipped) { }

    public static RunResult run(Path project, StudyProject.Actor initiator,
            StudyOperationExecutor executor, boolean resume) throws IOException {
        return run(project, initiator, executor, resume, false);
    }

    public static RunResult run(Path project, StudyProject.Actor initiator,
            StudyOperationExecutor executor, boolean resume, boolean offline)
            throws IOException {
        Path root = project.toAbsolutePath().normalize();
        return StudyProjectLock.with(root, () -> runLocked(root, initiator, executor, resume, offline));
    }

    private static RunResult runLocked(Path root, StudyProject.Actor initiator,
            StudyOperationExecutor executor, boolean resume, boolean offline)
            throws IOException {
        StudyProject.verify(root);
        StudyPlan plan = StudyPlanner.loadApproved(root);
        List<Path> priorManifests = new ArrayList<>();
        if (Files.isDirectory(root.resolve("runs"))) try (var paths = Files.walk(root.resolve("runs"))) {
            priorManifests.addAll(paths.filter(path -> path.getFileName().toString().equals("manifest.yaml")
                    && !path.getParent().getFileName().toString().startsWith(".operation-"))
                .sorted((left, right) -> {
                    try { return Files.getLastModifiedTime(right).compareTo(Files.getLastModifiedTime(left)); }
                    catch (IOException failure) { throw new java.io.UncheckedIOException(failure); }
                }).toList());
        }
        Map<String, Object> planDocument = plan.document();
        Map<String, Object> study = StudyProject.resolved(root).document();
        Map<String, Object> software = softwareIdentity();
        StudyDataDictionary.verify(root, study);
        if (offline) {
            Map<String, Object> resources = StudyIO.optionalMap(study.get("resources"), "resources");
            resources.put("network", "offline");
            study.put("resources", resources);
        }
        String baseRunId = "run-" + plan.id();
        String runId = baseRunId;
        for (int attempt = 2; Files.exists(root.resolve("runs").resolve(runId)); attempt++)
            runId = baseRunId + "-attempt-" + String.format("%04d", attempt);
        Path run = root.resolve("runs").resolve(runId);
        Path operationsRoot = run.resolve("operations");
        Files.createDirectories(operationsRoot);

        Map<String, Object> started = new LinkedHashMap<>();
        started.put("run_id", runId);
        started.put("plan_id", plan.id());
        started.put("plan_hash", plan.hash());
        started.put("state", "running");
        started.put("execution_mode", resume ? "resume" : "rerun");
        started.put("offline", offline);
        started.put("phase", "execution");
        started.put("summary", "Started approved plan " + plan.id());
        StudyProject.append(root, "run_started", initiator, started);

        Map<String, String> states = new LinkedHashMap<>();
        Map<String, Path> directories = new LinkedHashMap<>();
        int completed = 0, reused = 0, waiting = 0, failed = 0, skipped = 0;
        List<Object> operations = StudyIO.list(
            planDocument.get("operations"), "plan.operations");
        for (Object value : operations) {
            Map<String, Object> operation = StudyIO.map(value, "operation");
            String id = StudyIO.requireString(operation, "id", "operation");
            String type = StudyIO.requireString(operation, "type", "operation");
            Path destination = operationsRoot.resolve(id);
            List<String> dependencies = strings(
                operation.get("depends_on"), "operation.depends_on");
            boolean toleratesIncomplete = Boolean.TRUE.equals(
                operation.get("allow_incomplete_dependencies"));
            boolean blocked = !toleratesIncomplete
                && dependencies.stream().anyMatch(dependency ->
                !"complete".equals(states.get(dependency))
                    && !"reused".equals(states.get(dependency)));
            Map<String, Path> dependencyDirectories = new LinkedHashMap<>();
            dependencies.forEach(dependency -> {
                Path directory = directories.get(dependency);
                if (directory != null) dependencyDirectories.put(dependency, directory);
            });
            Map<String, String> inputs = inputHashes(
                root, study, operation, dependencyDirectories);
            Map<String, Object> cacheIdentity = new LinkedHashMap<>();
            cacheIdentity.put("operation", operation);
            cacheIdentity.put("resolved_study_hash",
                planDocument.get("resolved_study_hash"));
            cacheIdentity.put("input_sha256", inputs);
            cacheIdentity.put("implementation", 3);
            cacheIdentity.put("software", software);
            String cacheKey = StudyIO.hash(cacheIdentity);
            List<StudyNeed> needs = StudyQuestionRules.prerequisites(study, operation);
            Path reusable = null;
            if (resume && !blocked && needs.isEmpty()) for (Path prior : priorManifests) {
                if (verified(prior.getParent(), cacheKey)) {
                    reusable = prior.getParent();
                    break;
                }
            }

            if (Files.exists(destination))
                throw new IOException("Unverified operation output already exists: "
                    + destination);
            if (reusable != null) {
                copyDirectory(reusable, destination);
                Map<String, Object> manifest = StudyIO.load(
                    destination.resolve("manifest.yaml"));
                manifest.put("run_id", runId);
                manifest.put("reused_from_cache", true);
                manifest.put("reused_from", root.relativize(reusable).toString().replace('\\', '/'));
                StudyIO.writeYaml(destination.resolve("manifest.yaml"), manifest);
                states.put(id, "reused");
                directories.put(id, destination);
                reused++;
                event(root, "operation_reused", runId, operation, cacheKey,
                    "Reused verified prior operation with identical inputs", destination);
                continue;
            }

            Files.createDirectories(operationsRoot);
            Path staging = Files.createTempDirectory(operationsRoot, ".operation-");
            StudyOperationExecutor.Result result;
            try {
                if (!needs.isEmpty()) {
                    result = StudyOperationExecutor.Result.waiting(
                        "Scientific prerequisites require scientist input; see QUESTIONS.md", Map.of(), needs);
                } else if (blocked) {
                    result = StudyOperationExecutor.Result.waiting(
                        "A prerequisite operation is incomplete",
                        Map.of("dependencies", dependencies));
                } else {
                    result = executor.execute(new StudyOperationExecutor.Context(
                        root, staging, operation, study, planDocument,
                        dependencyDirectories));
                }
                Map<String, Object> manifest = manifest(runId, operation,
                    cacheKey, inputs, result, staging);
                manifest.put("software", software);
                StudyIO.writeYaml(staging.resolve("manifest.yaml"), manifest);
                publish(staging, destination);
                directories.put(id, destination);
                if (result.state() == StudyOperationExecutor.State.COMPLETE) {
                    states.put(id, "complete");
                    completed++;
                    event(root, "operation_completed", runId, operation,
                        cacheKey, result.summary(), destination);
                } else if (result.state() == StudyOperationExecutor.State.NOT_APPLICABLE) {
                    states.put(id, "not_applicable");
                    skipped++;
                    event(root, "operation_skipped", runId, operation,
                        cacheKey, result.summary(), destination);
                } else if (result.state() == StudyOperationExecutor.State.FAILED) {
                    states.put(id, "failed");
                    failed++;
                    event(root, "operation_failed", runId, operation,
                        cacheKey, result.summary(), destination);
                } else {
                    states.put(id, "waiting");
                    waiting++;
                    event(root, "operation_waiting", runId, operation,
                        cacheKey, result.summary(), destination);
                }
            } catch (Exception failure) {
                failed++;
                states.put(id, "failed");
                Map<String, Object> details = new LinkedHashMap<>();
                details.put("exception", failure.getClass().getName());
                details.put("message", String.valueOf(failure.getMessage()));
                result = new StudyOperationExecutor.Result(
                    StudyOperationExecutor.State.FAILED,
                    "Operation failed: " + failure.getMessage(), List.of(), details);
                Map<String, Object> manifest = manifest(runId, operation,
                    cacheKey, inputs, result, staging);
                manifest.put("software", software);
                StudyIO.writeYaml(staging.resolve("manifest.yaml"), manifest);
                if (Files.exists(destination)) throw failure instanceof IOException io
                    ? io : new IOException("Could not finalize published operation", failure);
                publish(staging, destination);
                directories.put(id, destination);
                event(root, "operation_failed", runId, operation,
                    cacheKey, result.summary(), destination);
            } finally {
                if (Files.exists(staging)) deleteTree(staging);
            }
        }

        String state = failed > 0 ? "failed"
            : waiting > 0 ? "waiting_for_input" : "completed";
        Map<String, Object> runManifest = new LinkedHashMap<>();
        runManifest.put("schema_version", 1);
        runManifest.put("run_id", runId);
        runManifest.put("plan_id", plan.id());
        runManifest.put("plan_hash", plan.hash());
        runManifest.put("state", state);
        runManifest.put("execution_mode", resume ? "resume" : "rerun");
        runManifest.put("offline", offline);
        runManifest.put("completed", completed);
        runManifest.put("reused", reused);
        runManifest.put("waiting", waiting);
        runManifest.put("failed", failed);
        runManifest.put("skipped", skipped);
        runManifest.put("operation_states", states);
        StudyIO.writeYaml(run.resolve("run-manifest.yaml"), runManifest);
        Map<String, Object> finished = new LinkedHashMap<>(runManifest);
        finished.put("phase", state.equals("completed") ? "reporting" : "execution");
        finished.put("summary", "Finished run " + runId + " with state " + state);
        StudyProject.append(root, state.equals("failed")
            ? "run_failed" : "run_completed", initiator, finished);
        StudyQuestions.reconcile(root, plan, runId);
        int indexed = StudyEvidenceIndex.rebuild(root);
        runManifest.put("indexed_results", indexed);
        StudyIO.writeYaml(run.resolve("run-manifest.yaml"), runManifest);
        return new RunResult(runId, state, completed, reused, waiting, failed, skipped);
    }

    private static Map<String, Object> softwareIdentity() throws IOException {
        try {
            Path source = Path.of(StudyExecutionEngine.class.getProtectionDomain().getCodeSource().getLocation().toURI());
            String hash;
            if (Files.isRegularFile(source)) hash = StudyIO.hash(source);
            else {
                Map<String, String> classes = new LinkedHashMap<>();
                try (var paths = Files.walk(source)) {
                    for (Path path : paths.filter(Files::isRegularFile).filter(path -> path.toString().endsWith(".class")).sorted().toList())
                        classes.put(source.relativize(path).toString().replace('\\', '/'), StudyIO.hash(path));
                }
                hash = StudyIO.hash(classes);
            }
            return Map.of("implementation_sha256", hash, "java_version", System.getProperty("java.version"),
                "identity_scope", Files.isRegularFile(source) ? "executable_jar" : "compiled_main_classes");
        } catch (java.net.URISyntaxException invalid) {
            throw new IOException("Cannot identify running JLinAlg implementation", invalid);
        }
    }

    private static Map<String, Object> manifest(String runId,
            Map<String, Object> operation, String cacheKey,
            Map<String, String> inputs, StudyOperationExecutor.Result result,
            Path staging) throws IOException {
        Map<String, Object> manifest = new LinkedHashMap<>();
        manifest.put("schema_version", 1);
        manifest.put("run_id", runId);
        manifest.put("operation_id", operation.get("id"));
        manifest.put("operation_type", operation.get("type"));
        manifest.put("hypothesis", operation.get("hypothesis"));
        manifest.put("state", result.state().value());
        manifest.put("summary", result.summary());
        manifest.put("cache_key", cacheKey);
        manifest.put("input_sha256", inputs);
        manifest.put("declared_outputs", result.outputs());
        manifest.put("details", result.details());
        Map<String, String> outputs = new LinkedHashMap<>();
        try (var paths = Files.walk(staging)) {
            for (Path path : paths.filter(Files::isRegularFile)
                    .sorted().toList()) {
                String relative = staging.relativize(path).toString()
                    .replace('\\', '/');
                if (!relative.equals("manifest.yaml"))
                    outputs.put(relative, StudyIO.hash(path));
            }
        }
        for (String output : result.outputs()) {
            Path expected = staging.resolve(output).normalize();
            if (!expected.startsWith(staging) || !Files.isRegularFile(expected))
                throw new IOException("Declared operation output is absent or unsafe: "
                    + output);
        }
        manifest.put("output_sha256", outputs);
        return manifest;
    }

    private static Map<String, String> inputHashes(Path project,
            Map<String, Object> study, Map<String, Object> operation,
            Map<String, Path> dependencies) throws IOException {
        Map<String, String> result = new LinkedHashMap<>();
        if (Set.of("advisor_review", "literature_support", "manuscript_draft").contains(operation.get("type"))) {
            Map<String, Object> assistance = StudyIO.optionalMap(study.get("assistance"), "assistance");
            for (Object raw : StudyIO.optionalList(assistance.get("response_files"), "assistance.response_files")) {
                Path receipt = project.resolve(String.valueOf(raw)).normalize();
                // Hash the receipt and all fixed-format companions before cache reuse, not only during execution.
                for (String name : List.of("receipt.yaml", "packet.yaml", "prompt.txt", "response.yaml", "index.html", "SESSION_HISTORY.md")) {
                    Path file = receipt.getParent().resolve(name);
                    result.put("browser_response:" + raw + ":" + name, Files.isRegularFile(file) ? StudyIO.hash(file) : "missing");
                }
            }
        }
        Map<String, Map<String, Object>> datasets = new LinkedHashMap<>();
        for (Object value : StudyIO.list(study.get("datasets"), "datasets")) {
            Map<String, Object> dataset = StudyIO.map(value, "dataset");
            datasets.put(String.valueOf(dataset.get("id")), dataset);
        }
        for (String id : strings(operation.get("datasets"), "operation.datasets")) {
            Map<String, Object> dataset = datasets.get(id);
            if (dataset == null) continue;
            Map<String, Object> files = StudyIO.map(dataset.get("files"), "dataset.files");
            for (Map.Entry<String, Object> file : files.entrySet()) {
                Path path = project.resolve(String.valueOf(file.getValue())).normalize();
                if (Files.isRegularFile(path))
                    result.put("dataset:" + id + ":" + file.getKey(),
                        StudyIO.hash(path));
            }
        }
        for (Map.Entry<String, Path> dependency : dependencies.entrySet()) {
            Path manifest = dependency.getValue().resolve("manifest.yaml");
            if (Files.isRegularFile(manifest)) {
                Map<String, Object> source = StudyIO.load(manifest);
                // Run IDs, reuse flags and timestamps do not change scientific inputs.
                Map<String, Object> identity = new LinkedHashMap<>();
                for (String key : List.of("state", "cache_key", "output_sha256", "details"))
                    identity.put(key, source.get(key));
                result.put("operation:" + dependency.getKey(),
                    StudyIO.hash(identity));
            }
        }
        if (Set.of("fine_mapping", "colocalization", "mendelian_randomization").contains(operation.get("type"))) {
            Map<String, Object> resources = StudyIO.optionalMap(study.get("resources"), "resources");
            Map<String, Object> ld = StudyIO.optionalMap(resources.get("ld"), "resources.ld");
            Path resourceRoot = project.resolve(String.valueOf(resources.getOrDefault("directory", "resources"))).normalize();
            Path reference = ld.containsKey("location") ? project.resolve(String.valueOf(ld.get("location"))).normalize()
                : resourceRoot.resolve("ld").resolve(String.valueOf(ld.getOrDefault("database", "1000g-phase3")));
            Set<String> populations = new java.util.LinkedHashSet<>();
            for (Object raw : StudyIO.list(study.get("hypotheses"), "hypotheses")) {
                Map<String, Object> hypothesis = StudyIO.map(raw, "hypothesis");
                if (!hypothesis.get("id").equals(operation.get("hypothesis"))) continue;
                Map<String, Object> bio = StudyIO.optionalMap(StudyIO.map(hypothesis.get("analysis"), "analysis").get("bioinformatics"), "bioinformatics");
                for (Object pair : StudyIO.optionalList(bio.get("genetic_pairs"), "genetic_pairs")) {
                    Map<String, Object> panel = StudyIO.optionalMap(StudyIO.map(pair, "pair").get("ld"), "ld");
                    if (!panel.containsKey("file_role") && !panel.containsKey("url")) populations.add(String.valueOf(panel.getOrDefault("population", ld.getOrDefault("population", "auto"))));
                }
            }
            for (String population : populations) if (population.matches("[A-Za-z0-9_-]+")) {
                Path prefix = reference.resolve("panels").resolve(population).resolve("genotypes");
                for (String suffix : List.of(".bed", ".bim", ".fam")) {
                    Path source = Path.of(prefix + suffix);
                    if (Files.isRegularFile(source)) result.put("ld:" + population + suffix, StudyIO.hash(source));
                }
            }
        }
        return result;
    }

    private static boolean verified(Path directory, String cacheKey)
            throws IOException {
        Path manifestPath = directory.resolve("manifest.yaml");
        if (!Files.isRegularFile(manifestPath)) return false;
        Map<String, Object> manifest = StudyIO.load(manifestPath);
        if (!"complete".equals(manifest.get("state"))
                || !cacheKey.equals(manifest.get("cache_key"))
                || Boolean.TRUE.equals(StudyIO.optionalMap(manifest.get("details"),
                    "details").get("retryable"))) return false;
        Map<String, Object> outputs = StudyIO.map(
            manifest.get("output_sha256"), "output_sha256");
        for (Map.Entry<String, Object> output : outputs.entrySet()) {
            Path path = directory.resolve(output.getKey()).normalize();
            if (!path.startsWith(directory) || !Files.isRegularFile(path)
                    || !String.valueOf(output.getValue()).equals(StudyIO.hash(path)))
                return false;
        }
        return true;
    }

    private static void event(Path project, String type, String runId,
            Map<String, Object> operation, String cacheKey, String summary,
            Path operationDirectory)
            throws IOException {
        Map<String, Object> payload = new LinkedHashMap<>();
        payload.put("run_id", runId);
        payload.put("operation_id", operation.get("id"));
        payload.put("operation_type", operation.get("type"));
        payload.put("hypothesis", operation.get("hypothesis"));
        payload.put("cache_key", cacheKey);
        payload.put("manifest_sha256", StudyIO.hash(
            operationDirectory.resolve("manifest.yaml")));
        payload.put("summary", summary);
        StudyProject.append(project, type,
            new StudyProject.Actor("jlinalg", "study-executor", null), payload);
        Map<String, Object> details = StudyIO.optionalMap(StudyIO.load(
            operationDirectory.resolve("manifest.yaml")).get("details"), "details");
        if (!type.equals("operation_reused") && "human_attested_not_observed".equals(details.get("browser_submission"))) {
            Map<String, Object> interaction = new LinkedHashMap<>(payload);
            interaction.put("summary", "Imported human-attested browser response from " + details.get("attested_provider")
                + " / " + details.get("llm_model") + "; JLinAlg did not observe external submission. " + summary);
            interaction.put("receipt", project.relativize(operationDirectory.resolve("browser-response/receipt.yaml")).toString().replace('\\', '/'));
            Map<String, Object> actor = StudyIO.map(details.get("browser_actor"), "browser_actor");
            StudyProject.append(project, "browser_handoff_response_recorded", StudyProject.Actor.human(String.valueOf(actor.get("id"))), interaction);
        }
        if (!type.equals("operation_reused") && Boolean.TRUE.equals(details.get("request_attempted"))) {
            Map<String, Object> interaction = new LinkedHashMap<>(payload);
            interaction.put("summary", "Recorded " + details.get("llm_provider") + " / "
                + details.get("llm_model") + " interaction and JLinAlg disposition: " + summary);
            interaction.put("interaction_directory", project.relativize(operationDirectory.resolve("llm/interaction-001")).toString().replace('\\', '/'));
            StudyProject.append(project, "llm_interaction_recorded", new StudyProject.Actor(
                "llm", String.valueOf(details.get("llm_provider")), String.valueOf(details.get("llm_model"))), interaction);
        }
    }

    private static List<String> strings(Object value, String name) {
        List<String> result = new ArrayList<>();
        for (Object item : StudyIO.optionalList(value, name))
            result.add(String.valueOf(item));
        return result;
    }

    private static void publish(Path staging, Path destination)
            throws IOException {
        try {
            Files.move(staging, destination, StandardCopyOption.ATOMIC_MOVE);
        } catch (AtomicMoveNotSupportedException unsupported) {
            Files.move(staging, destination);
        }
    }

    private static void copyDirectory(Path source, Path destination)
            throws IOException {
        if (Files.exists(destination))
            throw new IOException("Destination already exists: " + destination);
        try (var paths = Files.walk(source)) {
            for (Path path : paths.sorted().toList()) {
                Path target = destination.resolve(source.relativize(path));
                if (Files.isDirectory(path)) Files.createDirectories(target);
                else Files.copy(path, target, StandardCopyOption.COPY_ATTRIBUTES);
            }
        }
    }

    private static void deleteTree(Path directory) throws IOException {
        try (var paths = Files.walk(directory)) {
            for (Path path : paths.sorted(Comparator.reverseOrder()).toList())
                Files.deleteIfExists(path);
        }
    }
}
