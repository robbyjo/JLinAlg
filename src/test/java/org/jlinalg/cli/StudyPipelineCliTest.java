/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyPlan;
import org.jlinalg.study.StudyPlanner;
import org.jlinalg.study.StudyProject;
import org.jlinalg.study.StudyResolver;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyPipelineCliTest {
    @TempDir Path directory;
    private String stdout;
    private String stderr;

    @Test void effectComparisonLabelsTwoPhenotypesFromOneCohortAndFormatsCorrelation()
            throws Exception {
        Path upstream = directory.resolve("same-cohort-input");
        Path standard = upstream.resolve("standard-results");
        Path output = directory.resolve("same-cohort-output");
        Files.createDirectories(standard);
        Files.createDirectories(output);
        Files.writeString(standard.resolve("bmi.tsv"), """
            feature_id\tbeta\tse\tp_value
            A\t1\t.1\t.01
            B\t2\t.1\t.01
            C\t3\t.1\t.01
            """);
        Files.writeString(standard.resolve("waist.tsv"), """
            feature_id\tbeta\tse\tp_value
            A\t2\t.1\t.01
            B\t4\t.1\t.01
            C\t6\t.1\t.01
            """);
        Map<String, Object> study = Map.of("datasets", List.of(
            Map.of("id", "bmi", "cohort", "FHS"),
            Map.of("id", "waist", "cohort", "FHS")));
        var context = new org.jlinalg.study.StudyOperationExecutor.Context(
            directory, output, Map.of("type", "effect_comparison"), study,
            Map.of(), Map.of("harmonized", upstream));
        var result = new StudyCliExecutor().execute(context);
        assertEquals(org.jlinalg.study.StudyOperationExecutor.State.COMPLETE,
            result.state());
        String comparison = Files.readString(output.resolve(
            "effect-comparisons.tsv"));
        assertTrue(comparison.contains("FHS\tFHS\twithin_cohort\t3\t1.0\t1.0000"),
            comparison);
        assertTrue(Files.isRegularFile(output.resolve(
            "effect-scatterplots/01-bmi-vs-waist.png")));
        assertTrue(Files.isRegularFile(output.resolve(
            "effect-scatterplots/01-bmi-vs-waist.pdf")));
    }

    @Test void candidatePrioritizationIgnoresStudyLevelEvidenceRowsWithoutFeatures()
            throws Exception {
        Path meta = directory.resolve("meta"), literature = directory.resolve("literature");
        Path output = directory.resolve("prioritization");
        Files.createDirectories(meta); Files.createDirectories(literature); Files.createDirectories(output);
        Files.writeString(meta.resolve("meta-analysis.tsv"),
            "feature_id\tbeta\tp_value\ti_squared\tstatus\nGENE1\t.2\t.001\t0\tok\n");
        Files.writeString(meta.resolve("manifest.yaml"), """
            operation_type: meta_analysis
            state: complete
            output_sha256: {meta-analysis.tsv: fixture}
            """);
        Files.writeString(literature.resolve("literature-search.yaml"), """
            evidence:
              - {kind: literature, feature_id: '', status: offline_policy}
              - {kind: literature, feature_id: GENE1, status: retrieved}
            """);
        Files.writeString(literature.resolve("manifest.yaml"), """
            operation_type: literature_search
            state: waiting_for_input
            output_sha256: {literature-search.yaml: fixture}
            """);
        var context = new org.jlinalg.study.StudyOperationExecutor.Context(
            directory, output, Map.of("type", "candidate_prioritization"),
            Map.of(), Map.of(), Map.of("meta", meta, "literature", literature));
        var result = new StudyCliExecutor().execute(context);
        assertEquals(org.jlinalg.study.StudyOperationExecutor.State.COMPLETE,
            result.state());
        String table = Files.readString(output.resolve("candidate-prioritization.tsv"));
        assertTrue(table.contains("GENE1"));
        assertTrue(table.contains("literature:retrieved"));
    }

    @Test void candidatePrioritizationReadsLargeGeneratedGeneticEvidence()
            throws Exception {
        Path meta = directory.resolve("large-meta"), genetics = directory.resolve("large-genetics");
        Path output = directory.resolve("large-prioritization");
        Files.createDirectories(meta); Files.createDirectories(genetics); Files.createDirectories(output);
        Files.writeString(meta.resolve("meta-analysis.tsv"),
            "feature_id\tbeta\tp_value\ti_squared\tstatus\nGENE1\t.2\t.001\t0\tok\n");
        Files.writeString(meta.resolve("manifest.yaml"), """
            operation_type: meta_analysis
            state: complete
            output_sha256: {meta-analysis.tsv: fixture}
            """);
        Files.writeString(genetics.resolve("genetic-evidence.yaml"),
            "pairs:\n  - {kind: mr, feature_id: GENE1, status: complete}\npadding: "
                + "x".repeat(3_300_000) + "\n");
        Files.writeString(genetics.resolve("manifest.yaml"), """
            operation_type: mendelian_randomization
            state: complete
            output_sha256: {genetic-evidence.yaml: fixture}
            """);
        var context = new org.jlinalg.study.StudyOperationExecutor.Context(
            directory, output, Map.of("type", "candidate_prioritization"),
            Map.of(), Map.of(), Map.of("meta", meta, "genetics", genetics));
        var result = new StudyCliExecutor().execute(context);
        assertEquals(org.jlinalg.study.StudyOperationExecutor.State.COMPLETE,
            result.state());
        assertTrue(Files.readString(output.resolve("candidate-prioritization.tsv"))
            .contains("mr:complete"));
    }

    @Test void preapprovedSensitivityRunsOfflineWithoutAnLlmAndKeepsPrimaryResults() throws Exception {
        Path project = directory.resolve("sensitivity");
        Files.createDirectories(project);
        for (int i = 1; i <= 4; i++) Files.writeString(project.resolve("c" + i + ".tsv"),
            "feature_id\tbeta\tse\tp_value\nGENE1\t" + (i * 10) + "\t1\t.001\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project: {id: sensitivity, title: Sensitivity fixture, description: Local exploratory follow-up.}
            resources: {downloads: deny}
            assistance:
              enabled: false
              allow_analyses: [fixed_random_comparison, leave_one_cohort_out]
            datasets:
              - {id: d1, cohort: c1, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: c1.tsv}}
              - {id: d2, cohort: c2, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: c2.tsv}}
              - {id: d3, cohort: c3, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: c3.tsv}}
              - {id: d4, cohort: c4, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: c4.tsv}}
            hypotheses:
              - id: H1
                description: Check cohort influence.
                analysis:
                  mode: manual
                  goals: [meta_analysis]
                  operations:
                    meta_analysis: {min_cohorts: 4}
                    sensitivity_analysis: {min_cohorts: 3}
            """);
        assertEquals(0, run("study", "auto", study.toString(), "--offline"), stderr);
        assertTrue(stdout.contains("state: completed"), stdout);
        Path decision = find(project.resolve("runs"), "decision.yaml");
        String disposition = Files.readString(decision);
        assertTrue(disposition.contains("leave_one_cohort_out"));
        assertTrue(disposition.contains("fixed_random_comparison"));
        try (var paths = Files.walk(decision.getParent())) {
            assertEquals(6, paths.filter(path -> path.getFileName().toString().equals("meta-analysis.tsv")).count());
        }
        String fixed = Files.readString(decision.getParent().resolve("fixed_random_comparison/analysis-1/meta-analysis.tsv"));
        String random = Files.readString(decision.getParent().resolve("fixed_random_comparison/analysis-2/meta-analysis.tsv"));
        assertNotEquals(fixed, random);
        assertEquals(0, run("study", "verify", "--project", project.toString()), stderr);
        assertEquals(0, run("study", "auto", study.toString(), "--offline"), stderr);
        assertTrue(stdout.contains("reused: 6"), stdout);
    }

    @Test void automaticLeaveOneCohortOutFlagsJointReplicationAndHeterogeneityImprovement()
            throws Exception {
        Path project = directory.resolve("automatic-odd-cohort");
        Files.createDirectories(project);
        Files.writeString(project.resolve("discovery.tsv"), """
            feature_id\tbeta\tse\tp_value
            GENE1\t1.0\t0.1\t0.001
            GENE2\t1.0\t0.1\t0.001
            """);
        Files.writeString(project.resolve("replication.tsv"), """
            feature_id\tbeta\tse\tp_value
            GENE1\t0.9\t0.1\t0.001
            GENE2\t1.1\t0.1\t0.001
            """);
        Files.writeString(project.resolve("odd.tsv"), """
            feature_id\tbeta\tse\tp_value
            GENE1\t-2.0\t0.1\t0.001
            GENE2\t-2.0\t0.1\t0.001
            """);
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project: {id: automatic-odd-cohort, title: Influence fixture, description: Detect an influential discordant cohort.}
            resources: {downloads: deny}
            cohorts:
              - {id: discovery}
              - {id: replication}
              - {id: odd}
            cohort_relationships:
              - {between: [discovery, replication], type: independent}
              - {between: [discovery, odd], type: independent}
              - {between: [replication, odd], type: independent}
            datasets:
              - {id: d1, cohort: discovery, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: discovery.tsv}}
              - {id: d2, cohort: replication, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: replication.tsv}}
              - {id: d3, cohort: odd, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: odd.tsv}}
            hypotheses:
              - id: H1
                description: Compare discovery with two replication cohorts.
                inputs: {discovery: [d1], replication: [d2, d3]}
                analysis: {mode: manual, goals: [discovery_replication, meta_analysis]}
            """);
        assertEquals(0, run("study", "auto", study.toString(), "--offline"), stderr);
        Path comparison = find(project.resolve("runs"), "sensitivity-comparison.tsv");
        String table = Files.readString(comparison);
        assertTrue(table.contains("leave_one_cohort_out\todd\t"), table);
        assertTrue(table.contains("\ttrue\tfeature_in_sensitivity_section_primary_unchanged"), table);
        String decision = Files.readString(comparison.getParent().resolve("decision.yaml"));
        assertTrue(decision.contains("omitted_cohort: odd"), decision);
        assertTrue(decision.contains("primary results unchanged"), decision);
        assertTrue(Files.isRegularFile(comparison.getParent()
            .resolve("leave_one_cohort_out/analysis-3/replication.tsv")));
        assertEquals(0, run("study", "verify", "--project", project.toString()), stderr);
    }

    @Test void automaticColocalizationPlansPhenotypeAdviceBeforeGwasDiscovery() throws Exception {
        Path project = directory.resolve("automatic-colocalization-plan");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/results.tsv"), "feature_id\tbeta\tse\tp_value\nprobe1\t.1\t.01\t1e-9\n");
        Files.writeString(project.resolve("data/irqtl.tsv"), "fixture\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project:
              id: coloc-plan
              title: Colocalization plan
              description: Plan related phenotype discovery before automatic colocalization.
            assistance:
              enabled: true
              provider: browser
              related_phenotypes: true
            datasets:
              - id: fhs-bmi
                cohort: fhs
                type: transcriptomics
                stage: association_results
                tissue: Whole Blood
                files: {results: data/results.tsv, irqtl: data/irqtl.tsv}
                columns: {feature_id: feature_id, effect: beta, standard_error: se, p_value: p_value}
            hypotheses:
              - id: H1
                description: Test shared genetic signals with BMI-related phenotypes.
                outcome: body mass index
                inputs: {discovery: [fhs-bmi]}
                analysis:
                  mode: manual
                  bioinformatics:
                    colocalization:
                      exposure_dataset: fhs-bmi
                      exposure_ancestry: European
                      ld_population: EUR
            """);
        StudyPlan plan = StudyPlanner.compile(StudyResolver.resolve(study));
        @SuppressWarnings("unchecked")
        List<Map<String, Object>> operations = (List<Map<String, Object>>) (List<?>) plan.document().get("operations");
        Map<String, Object> advice = operations.stream().filter(operation -> operation.get("type").equals("related_phenotype_advice")).findFirst().orElseThrow();
        Map<String, Object> discovery = operations.stream().filter(operation -> operation.get("type").equals("source_discovery")).findFirst().orElseThrow();
        Map<String, Object> coloc = operations.stream().filter(operation -> operation.get("type").equals("colocalization")).findFirst().orElseThrow();
        assertTrue(((List<?>) discovery.get("depends_on")).contains(advice.get("id")));
        assertTrue(((List<?>) coloc.get("depends_on")).contains(discovery.get("id")));
        assertTrue(((List<?>) coloc.get("depends_on")).stream().map(Object::toString).anyMatch(id -> id.endsWith("irqtl_candidate_review")));
    }

    @Test void annotatedDashboardShowsGenesRequirementsSensitivityAndInapplicableScoring() throws Exception {
        Path project = directory.resolve("annotation-dashboard");
        Files.createDirectories(project);
        Files.writeString(project.resolve("a.tsv"), "feature_id\tbeta\tse\tp_value\np1\t1\t.1\t.00001\np2\t2\t.1\t.000001\n");
        Files.writeString(project.resolve("b.tsv"), "feature_id\tbeta\tse\tp_value\np1\t3\t.1\t.00001\np2\t4\t.1\t.000001\n");
        Files.writeString(project.resolve("probes.tsv"), "feature_id\tgene_id\tgene_symbol\tchr\tposition\np1\tENSG1\tGENE1\t1\t100\np2\tENSG2\t<script>GENE2</script>\t2\t200\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project: {id: annotation-dashboard, title: Probe annotation, description: Test auditable feature annotation and recommendations.}
            resources: {downloads: deny}
            datasets:
              - id: a
                type: transcriptomics
                stage: association_results
                feature_id_type: probe
                genome_build: GRCh38
                coordinate_system: 1-based-closed
                files: {results: a.tsv, annotation: probes.tsv}
              - id: b
                type: transcriptomics
                stage: association_results
                feature_id_type: probe
                genome_build: GRCh38
                coordinate_system: 1-based-closed
                files: {results: b.tsv, annotation: probes.tsv}
            hypotheses:
              - id: H1
                description: Compare splicing effects.
                inputs: {discovery: [a], replication: [b]}
                analysis: {mode: manual, goals: [discovery_replication, meta_analysis, genetic_support]}
            """);
        assertEquals(0, run("study", "auto", study.toString()), stderr);
        assertTrue(stdout.contains("skipped: 1"), stdout);
        String html = Files.readString(project.resolve("index/index.html"));
        assertTrue(html.contains("Major findings and next steps"));
        assertTrue(html.contains("GENE1"));
        assertTrue(html.contains("&lt;script&gt;GENE2&lt;/script&gt;"));
        assertFalse(html.contains("<script>GENE2</script>"));
        assertTrue(html.contains("What is needed next?"));
        assertTrue(html.contains("leave-one-cohort-out"));
        assertTrue(html.contains("not_applicable"));
        assertTrue(html.contains("not_run"));
        assertTrue(html.contains("BH_FDR"));
        assertTrue(Files.readString(project.resolve("index/LLM_CONTEXT.md")).contains("GENE1"));
        assertEquals(0, run("study", "verify", "--project", project.toString()));
        assertEquals(0, run("study", "auto", study.toString()), stderr);
        assertTrue(stdout.contains("skipped: 1"), stdout);
    }

    @Test void missingStatusMeansBeginningAndInitDoesNotClaimProgress()
            throws Exception {
        Path project = directory.resolve("new-project");
        assertEquals(0, run("study", "status", "--project", project.toString()));
        assertTrue(stdout.contains("beginning"));
        assertFalse(Files.exists(project.resolve("state/project-status.yaml")));

        Path study = project.resolve("study.yaml");
        assertEquals(0, run("study", "init", "--template",
            "association-results", "--out", study.toString()));
        assertEquals(0, run("study", "validate", study.toString()));
        assertTrue(stdout.contains("valid: no validation findings"));
        List<String> validationLines = stdout.lines().toList();
        assertEquals("valid: no validation findings", validationLines.get(0));
        assertTrue(validationLines.get(1).startsWith("study_hash:"));
        assertFalse(Files.exists(project.resolve("ledger/events.jsonl")));
    }

    @Test void validationRejectsUnknownFieldsDuplicateIdsAndMissingFiles()
            throws Exception {
        Path study = directory.resolve("bad.yaml");
        Files.writeString(study, """
            schema_version: 1
            unknown: value
            project:
              id: bad-study
              title: Bad
              description: Invalid contract.
            datasets:
              - id: data
                type: transcriptomics
                stage: processed_matrix
                files: {matrix: absent.tsv}
              - id: data
                type: transcriptomics
                stage: association_results
                files: {results: absent.tsv}
            hypotheses:
              - id: H1
                description: Test.
                inputs: {discovery: [missing]}
                analysis:
                  goals: [meta_anlysis]
                  harmonization: {invalid_rows: ignore, common_universe: union}
                  replication: {adjustment_scope: all, require_same_sign: sometimes}
            analysis_requests:
              - id: request-1
                type: sensitivity_analysis
                description: Invalid hypothesis reference.
                for_hypothesis: UNKNOWN
            """);
        assertEquals(2, run("study", "validate", study.toString(),
            "--check-files"));
        assertTrue(stdout.contains("$.unknown"));
        assertTrue(stdout.contains("duplicate identifier"));
        assertTrue(stdout.contains("unknown dataset"));
        assertTrue(stdout.contains("unregistered analysis goal"));
        assertTrue(stdout.contains("harmonization.invalid_rows"));
        assertTrue(stdout.contains("harmonization.common_universe"));
        assertTrue(stdout.contains("replication.adjustment_scope"));
        assertTrue(stdout.contains("replication.require_same_sign"));
        assertTrue(stdout.contains("unknown hypothesis"));
        assertTrue(stdout.contains("file is absent"));
    }

    @Test void autoModeExtendsExplicitGoalsAndPlansIndependentEvidenceBranches()
            throws Exception {
        Path project = directory.resolve("auto-plan");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/a.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Files.writeString(project.resolve("data/b.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t.9\t.3\t.02\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project:
              id: auto-plan
              title: Automatic evidence graph
              description: Infer every eligible publication and drug-discovery branch.
              objective: [publication, drug_discovery_prioritization]
            cohorts: [{id: a}, {id: b}]
            cohort_relationships:
              - {between: [a, b], type: independent}
            datasets:
              - id: a-results
                cohort: a
                type: transcriptomics
                stage: association_results
                annotation_release: GENCODE_v48
                files: {results: data/a.tsv}
                columns: {feature_id: feature_id, effect: beta, standard_error: se, p_value: p_value}
              - id: b-results
                cohort: b
                type: transcriptomics
                stage: association_results
                annotation_release: GENCODE_v48
                files: {results: data/b.tsv}
                columns: {feature_id: feature_id, effect: beta, standard_error: se, p_value: p_value}
            hypotheses:
              - id: H1
                description: Test and layer evidence for the association.
                inputs: {discovery: [a-results], replication: [b-results]}
                analysis:
                  mode: auto
                  goals: [discovery_replication]
            """);
        var resolved = StudyResolver.resolve(study);
        String resolvedText = resolved.document().toString();
        assertTrue(resolvedText.contains("meta_analysis"));
        assertTrue(resolvedText.contains("gene_ontology"));
        assertTrue(resolvedText.contains("genetic_support"));
        assertTrue(resolvedText.contains("candidate_prioritization"));
        assertTrue(resolvedText.contains("1000g-phase3"));

        StudyPlan plan = StudyPlanner.compile(resolved);
        @SuppressWarnings("unchecked")
        List<Map<String, Object>> operations = (List<Map<String, Object>>)
            (List<?>) plan.document().get("operations");
        Map<String, Object> go = operations.stream().filter(operation ->
            operation.get("type").equals("gene_ontology")).findFirst().orElseThrow();
        Map<String, Object> pathway = operations.stream().filter(operation ->
            operation.get("type").equals("pathway_enrichment")).findFirst().orElseThrow();
        assertEquals(go.get("depends_on"), pathway.get("depends_on"));
        assertTrue(operations.stream().anyMatch(operation ->
            operation.get("type").equals("fine_mapping")));
        assertTrue(operations.stream().anyMatch(operation ->
            operation.get("type").equals("colocalization")));
        assertTrue(operations.stream().anyMatch(operation ->
            operation.get("type").equals("mendelian_randomization")));
        assertTrue(operations.stream().anyMatch(operation ->
            operation.get("type").equals("candidate_prioritization")));
    }

    @Test void oneCommandAutoRunVerifiesAndWritesLlmEvidenceContext()
            throws Exception {
        Path project = directory.resolve("one-command");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/results.tsv"),
            "feature_id\tbeta\tse\tp_value\nGENE1\t1\t.2\t.01\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project:
              id: one-command
              title: One command result-first study
              description: Exercise automatic lifecycle execution without external resources.
              objective: publication
            datasets:
              - id: results
                type: transcriptomics
                stage: association_results
                files: {results: data/results.tsv}
                columns: {feature_id: feature_id, effect: beta, standard_error: se, p_value: p_value}
            hypotheses:
              - id: H1
                description: Pool the supplied association result.
                inputs: {discovery: [results]}
                analysis: {mode: manual, goals: [meta_analysis]}
            """);
        assertEquals(0, run("study", "auto", study.toString(),
            "--project", project.toString(), "--actor-id", "scientist-1"));
        assertTrue(stdout.contains("state: completed"));
        assertTrue(Files.isRegularFile(project.resolve("index/LLM_CONTEXT.md")));
        assertTrue(Files.readString(project.resolve("index/LLM_CONTEXT.md"))
            .contains("Trace every numerical claim"));
        assertTrue(Files.isRegularFile(project.resolve(
            "index/evidence-summary.json")));
        assertEquals(0, run("study", "verify", "--project",
            project.toString()));
    }

    @Test void resultFirstWorkflowPlansRunsCachesIndexesAndVerifies()
            throws Exception {
        Path project = directory.resolve("result-first");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/discovery.tsv"), """
            feature_id	beta	se	p_value
            gene-a	0.50	0.10	0.0001
            gene-b	-0.20	0.10	0.04
            gene-c	0.10	0.20	0.80
            """);
        Files.writeString(project.resolve("data/replication.tsv"), """
            feature_id	beta	se	p_value
            gene-a	0.40	0.12	0.001
            gene-b	0.10	0.15	0.50
            gene-c	0.30	0.10	0.001
            """);
        Path study = project.resolve("study.yaml");
        Files.writeString(study, resultStudy(true));

        assertEquals(0, run("study", "resolve", study.toString(),
            "--actor-id", "scientist-1"));
        assertTrue(Files.isRegularFile(project.resolve("resolved-study.yaml")));
        assertTrue(Files.isRegularFile(project.resolve("SESSION_HISTORY.md")));
        assertTrue(Files.isRegularFile(project.resolve("AGENTS.md")));
        assertEquals(0, run("study", "plan", "--project", project.toString(),
            "--actor-id", "scientist-1"));
        assertFalse(stdout.contains("type: association"));
        String planId = activePlan(project);
        assertEquals(0, run("study", "approve", "--project", project.toString(),
            "--plan", planId, "--approved-by", "scientist-1"));

        assertEquals(0, run("study", "run", "--project", project.toString(),
            "--actor-id", "scientist-1"));
        assertTrue(stdout.contains("state: completed"));
        Path run = project.resolve("runs/run-" + planId);
        Path replication = find(run, "replication.tsv");
        assertTrue(Files.readString(replication).contains(
            "gene-a\tdiscovery-data\treplication-data"));
        assertTrue(Files.readString(replication).contains("\ttrue\ttrue\ttrue"));
        assertTrue(Files.readString(replication).contains(
            "gene-c\tdiscovery-data\treplication-data\t0.1\t0.8"));
        assertTrue(Files.readString(replication).contains("\tfalse\ttrue\tfalse"));
        Path comparisons = find(run, "effect-comparisons.tsv");
        assertTrue(Files.readString(comparisons).contains("pearson_r_display"));
        assertTrue(Files.readString(comparisons).matches("(?s).*\\t-?\\d+\\.\\d{4}\\R.*"));
        assertTrue(Files.isRegularFile(find(run, "replication-venn.tsv")));
        try (var paths = Files.walk(run)) {
            assertTrue(paths.anyMatch(path -> path.toString().contains("effect-scatterplots")
                && path.toString().endsWith(".pdf")));
        }
        try (var paths = Files.walk(run)) {
            assertTrue(paths.anyMatch(path -> path.toString().contains("replication-venn")
                && path.toString().endsWith(".png")));
        }
        Path meta = find(run, "meta-analysis.tsv");
        assertTrue(Files.readString(meta).contains("gene-a"));
        String diagnostics = Files.readString(find(run, "diagnostics.yaml"));
        assertTrue(diagnostics.contains("discovery_replication_summary:"));
        assertTrue(diagnostics.contains("replicated_in_at_least_one:"));
        assertTrue(diagnostics.contains("discovery_replication:"));
        assertTrue(diagnostics.contains("ols_slope:"));
        assertTrue(diagnostics.contains("same_sign_fraction:"));
        assertTrue(Files.isRegularFile(project.resolve("index/index.json")));
        assertTrue(Files.readString(project.resolve("index/index.html"))
            .contains("Project evidence index"));
        String llmContext = Files.readString(
            project.resolve("index/LLM_CONTEXT.md"));
        assertTrue(llmContext.contains("Key findings"));
        assertTrue(llmContext.contains("ols_slope"));
        assertEquals(0, run("study", "verify", "--project", project.toString()));

        assertEquals(0, run("study", "run", "--project", project.toString(),
            "--resume", "--actor-id", "scientist-1"));
        assertTrue(stdout.contains("reused: 6"));
        assertEquals("completed", StudyProject.status(project).get("state"));
    }

    @Test void resultFirstSupportsRRowNamesStrictIntersectionAndReducedFdr()
            throws Exception {
        Path project = directory.resolve("reduced-fdr");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/discovery.csv"), """
            ,Fx_BMI,SE_BMI,P_BMI,N_BMI,Fx_Other
            a,1,.1,.001,100,.5
            b,1,.1,.01,100,.5
            c,1,.1,.9,100,.5
            d,1,.1,.02,100,NA
            """);
        Files.writeString(project.resolve("data/replication.csv"), """
            ,Fx_BMI,SE_BMI,P_BMI,N_BMI,Fx_Other
            a,1,.1,.04,100,.5
            b,-1,.1,.04,100,.5
            c,1,.1,.9,100,.5
            d,1,.1,.001,100,.5
            """);
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project:
              id: reduced-fdr
              title: Reduced-FDR result workflow
              description: Exercise strict common-universe replication.
            cohorts:
              - {id: discovery}
              - {id: replication}
            cohort_relationships:
              - between: [discovery, replication]
                type: independent
            datasets:
              - id: discovery-data
                cohort: discovery
                type: transcriptomics
                stage: association_results
                files: {results: data/discovery.csv}
                columns:
                  feature_id: __row_name__
                  effect: Fx_BMI
                  standard_error: SE_BMI
                  p_value: P_BMI
                  validation:
                    finite: [Fx_BMI, SE_BMI, P_BMI, N_BMI, Fx_Other]
                    positive: [N_BMI]
                    nonnegative: [SE_BMI]
                    probability: [P_BMI]
              - id: replication-data
                cohort: replication
                type: transcriptomics
                stage: association_results
                files: {results: data/replication.csv}
                columns:
                  feature_id: __row_name__
                  effect: Fx_BMI
                  standard_error: SE_BMI
                  p_value: P_BMI
                  validation:
                    finite: [Fx_BMI, SE_BMI, P_BMI, N_BMI, Fx_Other]
                    positive: [N_BMI]
                    nonnegative: [SE_BMI]
                    probability: [P_BMI]
            hypotheses:
              - id: H1
                description: Discovery results replicate after reduced-FDR correction.
                inputs:
                  discovery: [discovery-data]
                  replication: [replication-data]
                analysis:
                  mode: manual
                  goals: [discovery_replication]
                  harmonization:
                    invalid_rows: exclude
                    common_universe: intersection
                  replication:
                    adjustment_scope: discovery_subset
                    require_same_sign: false
                    threshold_comparison: less_than
            """);

        assertEquals(0, run("study", "resolve", study.toString(),
            "--actor-id", "scientist-1"));
        assertEquals(0, run("study", "plan", "--project", project.toString(),
            "--actor-id", "scientist-1"));
        String planId = activePlan(project);
        assertEquals(0, run("study", "approve", "--project", project.toString(),
            "--plan", planId, "--approved-by", "scientist-1"));
        assertEquals(0, run("study", "run", "--project", project.toString(),
            "--actor-id", "scientist-1"));

        Path run = project.resolve("runs/run-" + planId);
        String harmonization = Files.readString(
            find(run, "harmonization-summary.tsv"));
        assertTrue(harmonization.contains("discovery-data\t4\t3\t1\t0\t3"));
        assertTrue(harmonization.contains("# exact_common_features\t4"));
        assertTrue(harmonization.contains("# strict_common_features\t3"));
        String replication = Files.readString(find(run, "replication.tsv"));
        assertTrue(replication.contains(
            "a\tdiscovery-data\treplication-data\t1.0\t0.001"));
        assertTrue(replication.contains("\ttrue\ttrue\ttrue"));
        assertTrue(replication.contains(
            "b\tdiscovery-data\treplication-data\t1.0\t0.01"));
        assertTrue(replication.contains("\ttrue\tfalse\ttrue"));
        assertFalse(replication.contains("d\tdiscovery-data"));
        assertEquals(0, run("study", "verify", "--project", project.toString()));
    }

    @Test void matrixInputsPlanAssociationWhileResultInputsDoNot()
            throws Exception {
        Path matrixProject = directory.resolve("matrix");
        Files.createDirectories(matrixProject);
        Path matrixStudy = matrixProject.resolve("study.yaml");
        Files.writeString(matrixStudy, """
            schema_version: 1
            project:
              id: matrix-study
              title: Matrix study
              description: Start from a processed matrix.
            datasets:
              - id: expression
                cohort: cohort-a
                type: transcriptomics
                stage: processed_matrix
                files:
                  matrix: data/matrix.tsv
                  phenotype: data/pheno.tsv
                columns:
                  sample_id: id
                  phenotypes: {trait: y}
            hypotheses:
              - id: H1
                description: Expression is associated with the trait.
                outcome: trait
                inputs: {discovery: [expression]}
            """);
        assertEquals(0, run("study", "resolve", matrixStudy.toString()));
        StudyPlan matrixPlan = StudyPlanner.compile(
            org.jlinalg.study.StudyResolver.resolve(matrixStudy));
        assertTrue(matrixPlan.document().toString().contains("association"));

        Path resultProject = directory.resolve("result-plan");
        Files.createDirectories(resultProject.resolve("data"));
        Files.writeString(resultProject.resolve("data/discovery.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Files.writeString(resultProject.resolve("data/replication.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t.8\t.3\t.02\n");
        Path resultStudy = resultProject.resolve("study.yaml");
        Files.writeString(resultStudy, resultStudy(true));
        StudyPlan resultPlan = StudyPlanner.compile(
            org.jlinalg.study.StudyResolver.resolve(resultStudy));
        List<Object> operations = list(resultPlan.document().get("operations"));
        assertTrue(operations.stream().map(Object::toString)
            .noneMatch(value -> value.contains("type=association,")));
    }

    @Test void processedMatrixExecutesAssociationAndProducesVerifiedResults()
            throws Exception {
        Path project = directory.resolve("matrix-run");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/expression.tsv"), """
            feature_id\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8
            gene-a\t0\t1\t2\t3\t4\t5\t6\t7
            gene-b\t1\t0\t1\t0\t1\t0\t1\t0
            """);
        Files.writeString(project.resolve("data/phenotypes.tsv"), """
            sample_id\ttrait\tage
            S3\t3.2\t32
            S1\t1.0\t20
            S8\t8.2\t58
            S4\t4.1\t36
            S2\t2.1\t25
            S7\t7.0\t51
            S6\t6.1\t46
            S5\t5.0\t41
            """);
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project:
              id: matrix-run
              title: Matrix execution study
              description: Execute a processed-matrix association.
            cohorts:
              - {id: cohort-a}
            datasets:
              - id: expression
                cohort: cohort-a
                type: transcriptomics
                stage: processed_matrix
                files:
                  matrix: data/expression.tsv
                  phenotype: data/phenotypes.tsv
                columns:
                  sample_id: sample_id
                  feature_id: feature_id
                  phenotypes: {trait: trait}
                  covariates: {age: age}
            hypotheses:
              - id: H1
                description: Expression is associated with the trait.
                outcome: trait
                covariates: [age]
                inputs: {discovery: [expression]}
                analysis:
                  mode: manual
                  goals: [association]
            """);
        assertEquals(0, run("study", "resolve", study.toString()));
        assertEquals(0, run("study", "plan", "--project", project.toString()));
        assertEquals(0, run("study", "approve", "--project", project.toString(),
            "--approved-by", "scientist-1"));
        assertEquals(0, run("study", "run", "--project", project.toString(),
            "--actor-id", "scientist-1"));
        assertTrue(stdout.contains("state: completed"));
        Path standard = find(project.resolve("runs"), "expression.tsv");
        String results = Files.readString(standard);
        assertTrue(results.startsWith("feature_id\tbeta\tse\tp_value"));
        assertTrue(results.contains("gene-a"));
        assertEquals(0, run("study", "verify", "--project", project.toString()));
    }

    @Test void resultFirstFollowupDependsOnHarmonizationNotVirtualDataset()
            throws Exception {
        Path project = directory.resolve("result-followup-plan");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/results.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project:
              id: followup-plan
              title: Follow-up plan
              description: Plan enrichment from supplied results.
            datasets:
              - id: supplied-results
                type: transcriptomics
                stage: association_results
                files: {results: data/results.tsv}
                columns:
                  feature_id: feature_id
                  effect: beta
                  standard_error: se
                  p_value: p_value
            hypotheses:
              - id: H1
                description: Enrich supplied association results.
                inputs: {discovery: [supplied-results]}
                analysis: {mode: manual, goals: [gene_ontology]}
            """);
        StudyPlan plan = StudyPlanner.compile(
            org.jlinalg.study.StudyResolver.resolve(study));
        assertFalse(plan.document().toString().contains("dataset:supplied-results"));
        assertTrue(plan.document().toString().contains("validate_harmonize"));
    }

    @Test void specializedRawMatrixWaitsInsteadOfRunningBulkAssociation()
            throws Exception {
        Path project = directory.resolve("single-cell-wait");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/counts.tsv"),
            "feature_id\tS1\tS2\ngene-a\t1\t2\n");
        Files.writeString(project.resolve("data/phenotypes.tsv"),
            "sample_id\ttrait\nS1\t0\nS2\t1\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project:
              id: single-cell-wait
              title: Single-cell contract
              description: Do not route a specialized raw matrix to bulk OLS.
            datasets:
              - id: cells
                type: single_cell_rna
                stage: raw_quantified
                files:
                  matrix: data/counts.tsv
                  phenotype: data/phenotypes.tsv
                columns: {sample_id: sample_id}
            hypotheses:
              - id: H1
                description: Test a cell-resolved association.
                outcome: trait
                inputs: {discovery: [cells]}
                analysis: {mode: manual, goals: [association]}
            """);
        assertEquals(0, run("study", "resolve", study.toString()));
        assertEquals(0, run("study", "plan", "--project", project.toString()));
        assertEquals(0, run("study", "approve", "--project", project.toString(),
            "--approved-by", "scientist-1"));
        assertEquals(0, run("study", "run", "--project", project.toString(),
            "--allow-partial"));
        assertTrue(stdout.contains("state: waiting_for_input"));
        assertTrue(stdout.contains("failed: 0"));
        Path request = find(project.resolve("runs"), "request.yaml");
        assertTrue(Files.readString(request).contains(
            "No automatic association adapter is registered"));
        assertEquals(0, run("study", "run", "--project", project.toString(),
            "--resume", "--allow-partial"));
        assertTrue(stdout.contains("state: waiting_for_input"));
        assertEquals(0, run("study", "verify", "--project", project.toString()));
    }

    @Test void amendmentsArePreviewedAppliedOnceAndDoNotRewriteBaseStudy()
            throws Exception {
        Path project = directory.resolve("amended");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/discovery.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, singleResultStudy());
        assertEquals(0, run("study", "resolve", study.toString(),
            "--actor-id", "scientist-1"));
        String base = String.valueOf(
            StudyProject.status(project).get("resolved_study_hash"));
        Path amendment = project.resolve("amendment-add-replication.yaml");
        Files.writeString(amendment, """
            schema_version: 1
            amendment:
              id: add-replication
              title: Add replication data
              base_study_hash: %s
              proposed_by: {type: human, id: scientist-1}
              rationale: Replication data became available.
            changes:
              add:
                cohorts:
                  - {id: cohort-b}
                datasets:
                  - id: replication-data
                    cohort: cohort-b
                    type: transcriptomics
                    stage: association_results
                    files: {results: data/replication.tsv}
                    columns:
                      feature_id: feature_id
                      effect: beta
                      standard_error: se
                      p_value: p_value
                analysis_requests:
                  - id: sensitivity-1
                    type: sensitivity_analysis
                    for_hypothesis: H1
                    description: Compare fixed and random effects.
              update:
                hypotheses:
                  - id: H1
                    inputs:
                      replication:
                        append: [replication-data]
            """.formatted(base));
        Files.writeString(project.resolve("data/replication.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t.9\t.3\t.02\n");
        String original = Files.readString(study);
        String event = String.valueOf(
            StudyProject.status(project).get("last_event_id"));
        assertEquals(0, run("study", "amend", amendment.toString(),
            "--project", project.toString(), "--check"));
        assertEquals(event, StudyProject.status(project).get("last_event_id"));
        assertEquals(0, run("study", "amend", amendment.toString(),
            "--project", project.toString(), "--approved-by", "scientist-1"));
        assertEquals(original, Files.readString(study));
        assertTrue(Files.readString(project.resolve("resolved-study.yaml"))
            .contains("replication-data"));
        assertEquals(2, run("study", "amend", amendment.toString(),
            "--project", project.toString(), "--approved-by", "scientist-1"));
        assertTrue(stderr.contains("already been applied"));
        assertEquals(0, run("study", "resolve", study.toString(),
            "--project", project.toString()));
        assertTrue(Files.readString(project.resolve("resolved-study.yaml"))
            .contains("replication-data"));
    }

    @Test void resolvedBaseIsImmutableAndAmendmentSupersedesApprovedPlan()
            throws Exception {
        Path project = directory.resolve("immutable");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/discovery.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, singleResultStudy());
        assertEquals(0, run("study", "resolve", study.toString()));
        assertEquals(0, run("study", "plan", "--project", project.toString()));
        String plan = activePlan(project);
        assertEquals(0, run("study", "approve", "--project", project.toString(),
            "--plan", plan, "--approved-by", "scientist-1"));

        String base = String.valueOf(
            StudyProject.status(project).get("resolved_study_hash"));
        Path amendment = project.resolve("amendment-add-request.yaml");
        Files.writeString(amendment, """
            schema_version: 1
            amendment:
              id: add-request
              title: Add a sensitivity request
              base_study_hash: %s
              rationale: Prespecify an additional analysis.
            changes:
              add:
                analysis_requests:
                  - id: sensitivity-2
                    type: sensitivity_analysis
                    for_hypothesis: H1
                    description: Refit after excluding influential observations.
            """.formatted(base));
        assertEquals(0, run("study", "amend", amendment.toString(),
            "--project", project.toString(), "--approved-by", "scientist-1"));
        @SuppressWarnings("unchecked")
        Map<String, Object> active = (Map<String, Object>)
            StudyProject.status(project).get("active_plan");
        assertEquals("superseded", active.get("status"));
        assertEquals(2, run("study", "run", "--project", project.toString()));
        assertTrue(stderr.contains("not approved"));

        String original = Files.readString(study);
        Files.writeString(study, original.replace(
            "Exercise immutable project amendments.", "Changed without amendment."));
        assertEquals(2, run("study", "resolve", study.toString(),
            "--project", project.toString()));
        assertTrue(stderr.contains("cannot replace its base study"));
        assertEquals(2, run("study", "verify", "--project", project.toString()));
        assertTrue(stderr.contains("study.yaml hash"));
    }

    @Test void assumptionsRequireExplicitApprovalAndTamperingIsDetected()
            throws Exception {
        Path project = directory.resolve("assumptions");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/discovery.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Files.writeString(project.resolve("data/replication.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t.8\t.3\t.02\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, resultStudy(false));
        assertEquals(0, run("study", "resolve", study.toString()));
        assertEquals(0, run("study", "plan", "--project", project.toString()));
        String plan = activePlan(project);
        assertEquals(2, run("study", "approve", "--project", project.toString(),
            "--plan", plan, "--approved-by", "scientist-1"));
        assertTrue(stderr.contains("pending assumptions"));
        assertEquals(0, run("study", "approve", "--project", project.toString(),
            "--plan", plan, "--approved-by", "scientist-1",
            "--accept-assumptions"));
        assertEquals(0L, StudyProject.status(project).get("pending_assumptions"));
        assertEquals(0, run("study", "run", "--project", project.toString()));
        Path artifact = find(project.resolve("runs"), "replication.tsv");
        Files.writeString(artifact, Files.readString(artifact) + "tamper\n");
        assertEquals(2, run("study", "verify", "--project", project.toString()));
        assertTrue(stderr.contains("Artifact hash mismatch"));
    }

    @Test void verificationAnchorsManifestsPlansAndDerivedNotebook()
            throws Exception {
        Path project = directory.resolve("audit-chain");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/discovery.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, singleResultStudy());
        assertEquals(0, run("study", "resolve", study.toString()));
        assertEquals(0, run("study", "plan", "--project", project.toString()));
        assertEquals(0, run("study", "approve", "--project", project.toString(),
            "--approved-by", "scientist-1"));
        assertEquals(0, run("study", "run", "--project", project.toString()));
        assertEquals(0, run("study", "verify", "--project", project.toString()));

        Path manifest = find(project.resolve("runs"), "manifest.yaml");
        Files.writeString(manifest, Files.readString(manifest) + "# altered\n");
        assertEquals(2, run("study", "verify", "--project", project.toString()));
        assertTrue(stderr.contains("not anchored in the ledger"));
    }

    @Test void llmObservationRecordsModelAndDraftsAnAmendment()
            throws Exception {
        Path project = directory.resolve("observation");
        Files.createDirectories(project.resolve("data"));
        Files.writeString(project.resolve("data/discovery.tsv"),
            "feature_id\tbeta\tse\tp_value\na\t1\t.2\t.01\n");
        Path study = project.resolve("study.yaml");
        Files.writeString(study, singleResultStudy());
        assertEquals(0, run("study", "resolve", study.toString()));
        Path plot = project.resolve("plot.png");
        Files.writeString(plot, "synthetic plot fixture");
        assertEquals(0, run("study", "observe",
            "--project", project.toString(), "--id", "off-diagonal-1",
            "--type", "off_diagonal_streak", "--artifact", plot.toString(),
            "--hypothesis", "H1", "--description",
            "Effect comparison contains an off-diagonal streak.",
            "--explanation", "processing batch",
            "--explanation", "biological subgroup",
            "--proposed-direction", "sensitivity_analysis",
            "--amendment-id", "check-batch",
            "--actor-type", "llm", "--actor-id", "local-reviewer",
            "--model", "sha256:model-fixture"));
        String observation = Files.readString(project.resolve(
            "observations/observation-off-diagonal-1.yaml"));
        assertTrue(observation.contains("sha256:model-fixture"));
        Path amendment = project.resolve(
            "amendments/proposed/amendment-check-batch.yaml");
        assertTrue(Files.isRegularFile(amendment));
        assertEquals(0, run("study", "amend", amendment.toString(),
            "--project", project.toString(), "--check"));

        assertEquals(2, run("study", "observe",
            "--project", project.toString(), "--id", "invalid-reference",
            "--type", "cluster", "--artifact", plot.toString(),
            "--hypothesis", "UNKNOWN", "--description", "Unknown reference.",
            "--proposed-direction", "sensitivity_analysis"));
        assertTrue(stderr.contains("unknown hypothesis"));
        assertFalse(Files.exists(project.resolve(
            "observations/observation-invalid-reference.yaml")));
    }

    private int run(String... arguments) {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        ByteArrayOutputStream err = new ByteArrayOutputStream();
        int status = JLinAlgCli.run(arguments,
            new PrintStream(out), new PrintStream(err));
        stdout = out.toString(StandardCharsets.UTF_8);
        stderr = err.toString(StandardCharsets.UTF_8);
        return status;
    }

    private static String activePlan(Path project) throws IOException {
        Map<String, Object> status = StudyProject.status(project);
        @SuppressWarnings("unchecked")
        Map<String, Object> active =
            (Map<String, Object>) status.get("active_plan");
        return String.valueOf(active.get("plan_id"));
    }

    private static Path find(Path root, String name) throws IOException {
        try (var paths = Files.walk(root)) {
            return paths.filter(path -> path.getFileName().toString().equals(name))
                .findFirst().orElseThrow();
        }
    }

    @SuppressWarnings("unchecked")
    private static List<Object> list(Object value) {
        return (List<Object>) value;
    }

    private static String resultStudy(boolean explicitRelationship) {
        return """
            schema_version: 1
            project:
              id: result-first
              title: Result-first study
              description: Validate and pool supplied association results.
            cohorts:
              - {id: cohort-a}
              - {id: cohort-b}
            %s
            datasets:
              - id: discovery-data
                cohort: cohort-a
                type: transcriptomics
                stage: association_results
                files: {results: data/discovery.tsv}
                columns:
                  feature_id: feature_id
                  effect: beta
                  standard_error: se
                  p_value: p_value
              - id: replication-data
                cohort: cohort-b
                type: transcriptomics
                stage: association_results
                files: {results: data/replication.tsv}
                columns:
                  feature_id: feature_id
                  effect: beta
                  standard_error: se
                  p_value: p_value
            hypotheses:
              - id: H1
                description: The association replicates across cohorts.
                inputs:
                  discovery: [discovery-data]
                  replication: [replication-data]
                analysis:
                  mode: manual
                  goals: [discovery_replication, meta_analysis]
            """.formatted(explicitRelationship
                ? "cohort_relationships:\n  - between: [cohort-a, cohort-b]\n    type: independent"
                : "");
    }

    private static String singleResultStudy() {
        return """
            schema_version: 1
            project:
              id: amendment-study
              title: Amendment study
              description: Exercise immutable project amendments.
            cohorts:
              - {id: cohort-a}
            datasets:
              - id: discovery-data
                cohort: cohort-a
                type: transcriptomics
                stage: association_results
                files: {results: data/discovery.tsv}
                columns:
                  feature_id: feature_id
                  effect: beta
                  standard_error: se
                  p_value: p_value
            hypotheses:
              - id: H1
                description: Test the discovery association.
                inputs: {discovery: [discovery-data]}
                analysis:
                  mode: manual
                  goals: [meta_analysis]
            """;
    }
}
