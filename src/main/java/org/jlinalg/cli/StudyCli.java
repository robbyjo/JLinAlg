/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyAmendments;
import org.jlinalg.study.StudyEvidenceIndex;
import org.jlinalg.study.StudyExecutionEngine;
import org.jlinalg.study.StudyObservations;
import org.jlinalg.study.StudyOperationRegistry;
import org.jlinalg.study.StudyPlan;
import org.jlinalg.study.StudyPlanner;
import org.jlinalg.study.StudyProject;
import org.jlinalg.study.StudyQuestions;
import org.jlinalg.study.StudyValidation;
import org.jlinalg.study.StudyValidator;
import org.yaml.snakeyaml.Yaml;

/** Auditable deterministic study-contract and end-to-end pipeline CLI. */
final class StudyCli {
    private StudyCli() { }

    static int run(String[] arguments, PrintStream output, PrintStream error) {
        try {
            if (arguments.length == 0 || Set.of("help", "--help", "-h")
                    .contains(arguments[0])) {
                output.println(help());
                return 0;
            }
            String action = arguments[0];
            Options options = Options.parse(
                java.util.Arrays.copyOfRange(arguments, 1, arguments.length));
            return switch (action) {
                case "init" -> init(options, output);
                case "validate" -> validate(options, output);
                case "resolve" -> resolve(options, output);
                case "status" -> status(options, output);
                case "amend" -> amend(options, output);
                case "plan" -> plan(options, output);
                case "approve" -> approve(options, output);
                case "run" -> execute(options, output);
                case "auto", "execute" -> auto(options, output);
                case "observe" -> observe(options, output);
                case "index" -> index(options, output);
                case "verify" -> verify(options, output);
                case "questions" -> questions(options, output);
                case "answer" -> answer(options, output);
                case "operations" -> operations(options, output);
                case "export" -> export(options, output);
                case "manuscript" -> manuscript(options, output);
                case "references" -> {
                    Path result = StudyCitationSources.fetch(options.project(Path.of(".")), options.onePositional("references"),
                        Path.of(options.required("out")), options.actor(), options.flag("offline"), new StudyPublicEvidence());
                    output.println("Audited public reference bundle: " + result.resolve("metadata/manifest.yaml"));
                    output.println("Inspect unresolved retrievals and source support before citing. No analyses or manuscript text transmitted.");
                    yield 0;
                }
                case "manuscript-editorial" -> editorial(options, output);
                case "handoff-import" -> {
                    Path result = StudyBrowserHandoff.importReply(Path.of(options.required("handoff")),
                        options.onePositional("handoff-import"), Path.of(options.required("out")), options.actor());
                    output.println("Browser transcript recorded: " + result.resolve("index.html"));
                    output.println("Human-attested, not provider-authenticated. No edits or analyses applied. Use receipt.yaml in assistance.response_files to resume.");
                    yield 0;
                }
                case "handoff-verify" -> {
                    StudyBrowserHandoff.verifyReceipt(options.onePositional("handoff-verify"));
                    output.println("Browser receipt integrity verified. Provider authenticity and scientific correctness are not established.");
                    yield 0;
                }
                case "handoff-open" -> {
                    options.noPositionals("handoff-open");
                    output.println(StudyBrowserLauncher.open(Path.of(options.required("handoff")), options.one("browser", "auto")));
                    yield 0;
                }
                case "manuscript-verify" -> {
                    options.noPositionals("manuscript-verify");
                    StudyManuscriptRevision.verifyBundle(options.project(Path.of(".")));
                    output.println("Manuscript revision integrity verified; scientific review remains required.");
                    yield 0;
                }
                default -> throw new IllegalArgumentException(
                    "Unknown study action: " + action);
            };
        } catch (IOException | RuntimeException failure) {
            error.println("jlinalg: " + failure.getMessage());
            return 2;
        }
    }

    private static int init(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("init");
        String template = options.one("template", "bulk-transcriptomics");
        Path path = Path.of(options.one("out", "study.yaml"))
            .toAbsolutePath().normalize();
        if (Files.exists(path) && !options.flag("force"))
            throw new IOException("Study file exists; use --force to replace: " + path);
        Path parent = path.getParent();
        if (parent != null) Files.createDirectories(parent);
        Files.writeString(path, template(template), StandardCharsets.UTF_8,
            options.flag("force")
                ? new StandardOpenOption[] {StandardOpenOption.CREATE,
                    StandardOpenOption.TRUNCATE_EXISTING}
                : new StandardOpenOption[] {StandardOpenOption.CREATE_NEW});
        output.println("Created " + path);
        output.println("Project status remains beginning until study resolve succeeds.");
        return 0;
    }

    private static int export(Options options, PrintStream output) throws IOException {
        options.noPositionals("export");
        Path archive = Path.of(options.required("out"));
        String hash = org.jlinalg.study.StudyEvidenceArchive.export(options.project(Path.of(".")), archive,
            options.actor(), options.flag("acknowledge-sensitive"));
        output.println("Local evidence archive: " + archive.toAbsolutePath());
        output.println("SHA-256: " + hash);
        output.println("Contains potentially sensitive evidence. No external submission was made; inspect EVIDENCE_BUNDLE.json for exclusions and audit gaps.");
        return 0;
    }

    private static int manuscript(Options options, PrintStream output) throws IOException {
        Path revision = options.onePositional("manuscript");
        Path result = StudyManuscriptRevision.create(options.project(Path.of(".")), revision,
            Path.of(options.required("out")), options.actor());
        output.println("Local manuscript revision: " + result.resolve("artifacts/manuscript.html"));
        output.println("No analyses rerun, source project modified, or network requests made. Scientific review remains required.");
        return 0;
    }

    private static int editorial(Options options, PrintStream output) throws IOException {
        options.noPositionals("manuscript-editorial");
        String provider = options.one("provider", "browser");
        boolean local = StudyAssistance.isLocal(provider);
        if (!Set.of("browser", "openai", "gemini", "local", "llama_cpp", "openai_compatible").contains(provider)) throw new IllegalArgumentException("Unsupported editorial provider");
        if (provider.equals("browser")) for (String key : List.of("model", "api-key-env", "oauth-token-env", "auth", "endpoint"))
            if (options.one(key, null) != null) throw new IllegalArgumentException("Browser handoff does not use --" + key + "; choose an API provider explicitly");
        if (!Set.of("openai_compatible", "browser").contains(provider) && (options.one("auth", null) != null || options.one("oauth-token-env", null) != null))
            throw new IllegalArgumentException("--auth and --oauth-token-env require --provider openai_compatible");
        if (!local && !Set.of("openai_compatible", "browser").contains(provider) && options.one("endpoint", null) != null)
            throw new IllegalArgumentException("--endpoint requires --provider openai_compatible (or explicit local adapter)");
        if (local && options.one("api-key-env", null) != null) throw new IllegalArgumentException("Local adapter does not use API credentials");
        Map<String, Object> config = new java.util.LinkedHashMap<>();
        config.put("provider", provider);
        if (!provider.equals("browser")) config.put("model", options.required("model"));
        config.put("timeout_seconds", local ? 300 : 60);
        if (local) config.put("endpoint", StudyAssistance.localEndpoint(options.required("endpoint")).toString());
        else if (provider.equals("openai_compatible")) {
            config.put("endpoint", StudyAssistance.compatibleEndpoint(options.required("endpoint")).toString());
            String auth = options.one("auth", "api_key"); config.put("auth", auth);
            if (!Set.of("api_key", "oauth").contains(auth)) throw new IllegalArgumentException("--auth must be api_key or oauth");
            if (auth.equals("oauth") && options.one("api-key-env", null) != null || auth.equals("api_key") && options.one("oauth-token-env", null) != null)
                throw new IllegalArgumentException("Choose only the credential option belonging to --auth");
            if (auth.equals("oauth")) config.put("oauth_token_env", options.required("oauth-token-env"));
            else config.put("api_key_env", options.one("api-key-env", "OPENAI_API_KEY"));
        } else if (!provider.equals("browser")) config.put("api_key_env", options.one("api-key-env", provider.equals("openai") ? "OPENAI_API_KEY" : "GEMINI_API_KEY"));
        Path out = StudyEditorialReview.request(options.project(Path.of(".")), Path.of(options.required("out")),
            options.required("approved-payload-sha256"), config, options.actor(), options.flag("offline"), new StudyAssistance());
        output.println("Editorial attempt recorded: " + out.resolve("artifacts/index.html"));
        output.println("Inspect the recorded disposition: a complete audit bundle does not mean the provider returned a usable review. No manuscript edits were applied.");
        return 0;
    }

    private static int validate(Options options, PrintStream output)
            throws IOException {
        Path study = options.studyPath();
        StudyValidation validation = StudyValidator.validate(
            study, options.flag("check-files"));
        String report = validation.report();
        output.print(report);
        if (!report.endsWith("\n")) output.println();
        if (validation.valid())
            output.println("study_hash: "
                + StudyValidator.requireValid(study,
                    options.flag("check-files")).hash());
        return validation.valid() ? 0 : 2;
    }

    private static int resolve(Options options, PrintStream output)
            throws IOException {
        Path study = options.studyPath();
        Path project = options.project(study.getParent());
        Map<String, Object> status = StudyProject.resolve(
            study, project, options.actor());
        output.print(new Yaml().dump(status));
        return 0;
    }

    private static int status(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("status");
        output.print(new Yaml().dump(StudyProject.status(
            options.project(Path.of("")))));
        return 0;
    }

    private static int amend(Options options, PrintStream output)
            throws IOException {
        Path amendment = options.onePositional("amend");
        Path project = options.project(amendment.getParent());
        if (options.flag("check")) {
            output.print(new Yaml().dump(
                StudyAmendments.preview(amendment, project)));
            return 0;
        }
        String approvedBy = options.required("approved-by");
        output.print(new Yaml().dump(StudyAmendments.apply(
            amendment, project, StudyProject.Actor.human(approvedBy))));
        return 0;
    }

    private static int plan(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("plan");
        StudyPlan plan = StudyPlanner.create(
            options.project(Path.of("")), options.actor());
        output.print(new Yaml().dump(plan.document()));
        return 0;
    }

    private static int approve(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("approve");
        Path project = options.project(Path.of(""));
        String plan = options.values.get("plan") == null
            ? activePlan(project) : options.one("plan", null);
        StudyPlanner.approve(project, plan,
            StudyProject.Actor.human(options.required("approved-by")),
            options.flag("accept-assumptions"));
        output.println("Approved " + plan);
        return 0;
    }

    private static int execute(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("run");
        StudyExecutionEngine.RunResult result = StudyExecutionEngine.run(
            options.project(Path.of("")), options.actor(),
            new StudyCliExecutor(), !options.flag("rerun"), options.flag("offline"));
        output.print(new Yaml().dump(Map.of(
            "run_id", result.runId(), "state", result.state(),
            "completed", result.completed(), "reused", result.reused(),
            "waiting", result.waiting(), "failed", result.failed(), "skipped", result.skipped())));
        if (result.failed() > 0) return 2;
        if (result.waiting() > 0 && !options.flag("allow-partial")) return 3;
        return 0;
    }

    private static int auto(Options options, PrintStream output)
            throws IOException {
        Path study = options.studyPath();
        Path project = options.project(study.getParent());
        StudyValidation validation = StudyValidator.validate(study, true);
        if (!validation.valid())
            throw new IllegalArgumentException(
                "Study validation failed:\n" + validation.report());

        String actorId = options.one("actor-id", "command-line-user");
        StudyProject.Actor actor = new StudyProject.Actor(
            options.one("actor-type", "human"), actorId,
            options.one("model", null));
        StudyProject.resolve(study, project, actor);
        StudyPlan plan = StudyPlanner.create(project, actor);
        Map<String, Object> status = StudyProject.status(project);
        Map<String, Object> active = status.get("active_plan") instanceof Map<?, ?>
            ? map(status.get("active_plan")) : Map.of();
        if ("proposed".equals(active.get("status"))) {
            // Shared-parent relationships remain unresolved scientist questions.
            // Approving independent branches does not approve those assumptions.
            StudyPlanner.approve(project, plan.id(),
                StudyProject.Actor.human(options.one("approved-by", actorId)),
                !list(plan.document().get("pending_decisions")).isEmpty());
        } else if (!"approved".equals(active.get("status"))) {
            throw new IllegalArgumentException(
                "The active plan cannot be executed automatically: "
                    + active.get("status"));
        }

        StudyExecutionEngine.RunResult result = StudyExecutionEngine.run(
            project, actor, new StudyCliExecutor(), !options.flag("rerun"), options.flag("offline"));
        StudyProject.verify(project);
        int indexed = StudyEvidenceIndex.rebuild(project);
        Map<String, Object> summary = new LinkedHashMap<>();
        summary.put("project", project.toString());
        summary.put("run_id", result.runId());
        summary.put("state", result.state());
        summary.put("completed", result.completed());
        summary.put("reused", result.reused());
        summary.put("waiting", result.waiting());
        summary.put("failed", result.failed());
        summary.put("skipped", result.skipped());
        summary.put("indexed_results", indexed);
        summary.put("llm_context", project.resolve("index/LLM_CONTEXT.md")
            .toString());
        output.print(new Yaml().dump(summary));
        return result.failed() > 0 ? 2 : 0;
    }

    private static int observe(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("observe");
        Map<String, Object> observation = StudyObservations.record(
            options.project(Path.of("")), options.required("id"),
            options.required("type"), Path.of(options.required("artifact")),
            options.many("hypothesis"), options.required("description"),
            options.many("explanation"),
            options.required("proposed-direction"),
            options.values.get("amendment-id") == null ? null
                : options.one("amendment-id", null),
            options.actor());
        output.print(new Yaml().dump(observation));
        return 0;
    }

    private static int index(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("index");
        int count = StudyEvidenceIndex.rebuild(options.project(Path.of("")));
        output.println("Indexed " + count + " operation results");
        return 0;
    }

    private static int verify(Options options, PrintStream output)
            throws IOException {
        options.noPositionals("verify");
        Path project = options.project(Path.of(""));
        if (options.flag("repair-views")) StudyProject.repairViews(project);
        StudyProject.verify(project);
        int count = StudyEvidenceIndex.rebuild(project);
        output.println("Project evidence verified; indexed results: " + count);
        return 0;
    }

    private static int questions(Options options, PrintStream output) throws IOException {
        options.noPositionals("questions");
        Path project = options.project(Path.of(""));
        if (options.values.containsKey("out")) {
            Path destination = Path.of(options.required("out"));
            StudyQuestions.export(project, destination, options.flag("include-answered"));
            output.println("Exported blank answers to " + destination + "; no answers are preaccepted.");
        } else output.print(new Yaml().dump(StudyQuestions.read(project)));
        return 0;
    }

    private static int answer(Options options, PrintStream output) throws IOException {
        Path file = options.onePositional("answer");
        if (options.flag("check") == options.flag("record"))
            throw new IllegalArgumentException("Choose exactly one of --check or --record; neither executes analyses");
        options.required("actor-id");
        Path project = options.project(Path.of(""));
        output.print(new Yaml().dump(options.flag("check")
            ? StudyQuestions.preview(project, file, options.actor())
            : StudyQuestions.record(project, file, options.actor())));
        return 0;
    }

    private static int operations(Options options, PrintStream output) {
        options.noPositionals("operations");
        Map<String, Object> result = new LinkedHashMap<>();
        StudyOperationRegistry.descriptors().forEach((name, descriptor) -> {
            Map<String, Object> item = new LinkedHashMap<>();
            item.put("executor", descriptor.executor());
            item.put("accepted_stages", descriptor.acceptedStages());
            item.put("modalities", descriptor.modalities());
            item.put("outputs", descriptor.outputs());
            item.put("diagnostics", descriptor.diagnostics());
            result.put(name, item);
        });
        output.print(new Yaml().dump(result));
        return 0;
    }

    private static String activePlan(Path project) throws IOException {
        Map<String, Object> status = StudyProject.status(project);
        Object value = status.get("active_plan");
        if (!(value instanceof Map<?, ?> active)
                || active.get("plan_id") == null)
            throw new IllegalArgumentException(
                "No active plan; run study plan first");
        return String.valueOf(active.get("plan_id"));
    }

    private static String template(String name) {
        String type;
        String stage;
        String primary;
        String file;
        switch (name) {
            case "bulk-transcriptomics" -> {
                type = "transcriptomics"; stage = "processed_matrix";
                primary = "matrix"; file = "data/expression.tsv";
            }
            case "epigenomics" -> {
                type = "epigenomics"; stage = "processed_matrix";
                primary = "matrix"; file = "data/methylation.tsv";
            }
            case "gwas" -> {
                type = "gwas"; stage = "association_results";
                primary = "results"; file = "data/gwas-results.tsv";
            }
            case "single-cell" -> {
                type = "single_cell_rna"; stage = "raw_quantified";
                primary = "matrix"; file = "data/counts.tsv";
            }
            case "spatial" -> {
                type = "spatial_transcriptomics"; stage = "raw_quantified";
                primary = "matrix"; file = "data/counts.tsv";
            }
            case "microbiome" -> {
                type = "microbiome"; stage = "processed_matrix";
                primary = "matrix"; file = "data/taxa.tsv";
            }
            case "multimodal" -> {
                type = "multimodal"; stage = "processed_matrix";
                primary = "matrix"; file = "data/features.tsv";
            }
            case "association-results" -> {
                type = "transcriptomics"; stage = "association_results";
                primary = "results"; file = "data/association-results.tsv";
            }
            default -> throw new IllegalArgumentException(
                "Unknown template: " + name);
        }
        String files = primary.equals("matrix")
            ? "      matrix: " + file + "\n"
                + "      phenotype: data/phenotypes.tsv\n"
            : "      results: " + file + "\n";
        String columns = primary.equals("matrix")
            ? "      sample_id: sample_id\n"
                + "      feature_id: feature_id\n"
                + "      phenotypes: {trait: trait}\n"
            : "      feature_id: feature_id\n"
                + "      effect: beta\n"
                + "      standard_error: se\n"
                + "      p_value: p_value\n";
        return "schema_version: 1\n\n"
            + "project:\n"
            + "  id: example-study\n"
            + "  title: Example molecular study\n"
            + "  description: Test a prespecified molecular association and layer supporting evidence.\n"
            + "  objective: publication\n\n"
            + "datasets:\n"
            + "  - id: discovery-data\n"
            + "    cohort: discovery\n"
            + "    type: " + type + "\n"
            + "    tissue: specify_tissue\n"
            + "    stage: " + stage + "\n"
            + "    files:\n" + files
            + "    columns:\n" + columns + "\n"
            + "hypotheses:\n"
            + "  - id: H1\n"
            + "    description: Molecular features are associated with the prespecified trait.\n"
            + "    outcome: trait\n"
            + "    covariates: []\n"
            + "    inputs:\n"
            + "      discovery: [discovery-data]\n"
            + "    analysis:\n"
            + "      mode: auto\n";
    }

    private static String help() {
        return """
            Usage: java -jar jlinalg-<version>.jar study ACTION [options]

              study init [--template bulk-transcriptomics|epigenomics|gwas|single-cell|spatial|microbiome|multimodal|association-results] [--out study.yaml]
              study validate [study.yaml] [--check-files]
              study resolve [study.yaml] [--project DIRECTORY]
              study status [--project DIRECTORY]
              study amend amendment-ID.yaml [--project DIRECTORY] [--check | --approved-by PERSON]
              study plan [--project DIRECTORY]
              study approve [--project DIRECTORY] [--plan ID] --approved-by PERSON [--accept-assumptions]
              study run [--project DIRECTORY] [--rerun] [--offline] [--allow-partial]
              study auto|execute [study.yaml] [--project DIRECTORY]
                [--approved-by PERSON] [--accept-assumptions] [--rerun] [--offline]
              study observe --id ID --type TYPE --artifact FILE --hypothesis H1
                --description TEXT --proposed-direction TYPE
                [--explanation TEXT ...] [--amendment-id ID]
              study index|verify|operations [--project DIRECTORY]
              study verify [--project DIRECTORY] [--repair-views]
              study questions [--project DIRECTORY] [--out answers.yaml] [--include-answered]
              study answer answers.yaml [--project DIRECTORY] --actor-id PERSON --check|--record
              study export --project DIRECTORY --out evidence.zip --acknowledge-sensitive
              study manuscript revision.yaml --project DIRECTORY --out NEW_BUNDLE
              study references references.yaml --project DIRECTORY --out NEW_REFERENCE_BUNDLE [--offline]
              study manuscript-verify --project REVISION_BUNDLE
              study manuscript-editorial --project REVISION_BUNDLE --out NEW_REVIEW_BUNDLE
                --approved-payload-sha256 HASH [--provider browser|openai|gemini|openai_compatible|local|llama_cpp]
                [--model MODEL] [--api-key-env VARIABLE] [--endpoint URL] [--auth api_key|oauth]
                [--oauth-token-env VARIABLE] [--offline]
                Default provider browser prepares a local handoff; it never submits a prompt.
              study handoff-import RESPONSE.yaml --handoff PACKET_DIRECTORY --out NEW_RECEIPT_DIRECTORY --actor-id SCIENTIST
              study handoff-verify RECEIPT_DIRECTORY/receipt.yaml
              study handoff-open --handoff PACKET_DIRECTORY [--browser auto|chrome|chromium|firefox|edge|links|none]

            Actor options: [--actor-type human|jlinalg|llm|external]
              [--actor-id ID] [--actor-model ACTOR_MODEL_ID]
              --model configures an API model and remains the legacy actor-model fallback.

            validate, status, operations, amendment --check, and plan compilation
            do not execute analyses. Missing project-status.yaml means beginning.
            study auto validates files, resolves, plans, records command-line
            approval, executes every eligible branch, verifies, and indexes the
            project. Missing scientific prerequisites are indexed as blockers.
            Repeated runs resume by default, retrying unfinished operations.
            --rerun bypasses completed-result reuse and retains previous attempts.
            --offline prohibits built-in network requests; --resume remains an alias.
            Questions and answers are documented locally in QUESTIONS.md and SESSION_HISTORY.md.
            Answers create reviewable proposals, not execution or transmission permission.
            --repair-views regenerates derived views only after authoritative evidence verifies.
            """;
    }

    private static final class Options {
        private static final Set<String> FLAGS = Set.of(
            "check-files", "check", "force", "explain",
            "accept-assumptions", "resume", "allow-partial", "rerun", "offline",
            "record", "include-answered", "repair-views", "acknowledge-sensitive");
        private static final Set<String> VALUES = Set.of(
            "template", "out", "project", "actor-type", "actor-id", "actor-model", "model",
            "approved-by", "plan", "artifact", "type", "id", "description",
            "explanation", "hypothesis", "proposed-direction", "amendment-id",
            "provider", "api-key-env", "approved-payload-sha256", "endpoint", "auth", "oauth-token-env", "handoff", "browser");
        private final Map<String, List<String>> values = new LinkedHashMap<>();
        private final Set<String> flags = new java.util.HashSet<>();
        private final List<String> positional = new ArrayList<>();

        static Options parse(String[] arguments) {
            Options result = new Options();
            for (int index = 0; index < arguments.length; index++) {
                String argument = arguments[index];
                if (!argument.startsWith("--")) {
                    result.positional.add(argument);
                    continue;
                }
                String key = argument.substring(2);
                if (FLAGS.contains(key)) {
                    if (!result.flags.add(key))
                        throw new IllegalArgumentException(
                            "Duplicate option: " + argument);
                } else if (VALUES.contains(key)) {
                    if (++index == arguments.length
                            || arguments[index].startsWith("--"))
                        throw new IllegalArgumentException(
                            argument + " requires a value");
                    result.values.computeIfAbsent(key,
                        ignored -> new ArrayList<>()).add(arguments[index]);
                    if (!Set.of("explanation", "hypothesis").contains(key)
                            && result.values.get(key).size() > 1)
                        throw new IllegalArgumentException(
                            "Duplicate option: " + argument);
                } else throw new IllegalArgumentException(
                    "Unknown study option: " + argument);
            }
            return result;
        }

        boolean flag(String key) {
            return flags.contains(key);
        }

        String one(String key, String fallback) {
            List<String> found = values.get(key);
            return found == null ? fallback : found.get(0);
        }

        String required(String key) {
            String result = one(key, null);
            if (result == null || result.isBlank())
                throw new IllegalArgumentException("--" + key + " is required");
            return result;
        }

        List<String> many(String key) {
            return List.copyOf(values.getOrDefault(key, List.of()));
        }

        Path onePositional(String action) {
            if (positional.size() != 1)
                throw new IllegalArgumentException(
                    "study " + action + " requires exactly one file");
            return Path.of(positional.get(0)).toAbsolutePath().normalize();
        }

        Path studyPath() {
            if (positional.size() > 1)
                throw new IllegalArgumentException(
                    "Only one study YAML file may be supplied");
            return Path.of(positional.isEmpty() ? "study.yaml" : positional.get(0))
                .toAbsolutePath().normalize();
        }

        Path project(Path fallback) {
            return Path.of(one("project", fallback.toString()))
                .toAbsolutePath().normalize();
        }

        StudyProject.Actor actor() {
            return new StudyProject.Actor(one("actor-type", "human"),
                one("actor-id", "unspecified"), one("actor-model", one("model", null)));
        }

        void noPositionals(String action) {
            if (!positional.isEmpty())
                throw new IllegalArgumentException(
                    "study " + action + " does not accept positional arguments");
        }
    }

    private static Map<String, Object> map(Object value) {
        Map<String, Object> result = new LinkedHashMap<>();
        if (value instanceof Map<?, ?> raw)
            raw.forEach((key, item) -> result.put(String.valueOf(key), item));
        return result;
    }

    private static List<Object> list(Object value) {
        return value instanceof List<?> items
            ? new ArrayList<>(items) : List.of();
    }
}
