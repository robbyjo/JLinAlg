/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyProject;
import org.jlinalg.study.StudyQuestions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.yaml.snakeyaml.Yaml;

class StudyQuestionCliTest {
    @TempDir Path root;
    private String stdout, stderr;

    @Test void resultFirstQuestionsAreAuditableAndResumeDoesNotDuplicateThem() throws Exception {
        prepare(false, false);
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        List<Map<String, Object>> first = questions();
        assertFalse(first.isEmpty());
        assertTrue(first.stream().anyMatch(q -> String.valueOf(q.get("semantic_key")).contains("files.annotation")));
        assertTrue(first.stream().anyMatch(q -> String.valueOf(q.get("semantic_key")).contains("genetic-inputs")));
        assertFalse(first.stream().anyMatch(q -> String.valueOf(q.get("semantic_key")).contains("variant_scoring")));
        assertTrue(first.stream().allMatch(q -> !((List<?>) q.get("evidence")).isEmpty()));
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        assertEquals(first, questions());
        assertTrue(Files.readString(root.resolve("index/index.html")).contains("Scientist questions and decisions"));
        assertEquals(0, run("study", "questions", "--project", root.toString()), stderr);
        assertTrue(stdout.contains("question_hash"));
        assertEquals(0, run("study", "verify", "--project", root.toString()), stderr);
    }

    @Test void annotationAnswerCliRequiresExplicitModeAndScientistAndChangesOnlyByAmendment() throws Exception {
        prepare(false, false);
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        Path answers = root.resolve("answers-for-meeting.yaml");
        assertEquals(0, run("study", "questions", "--project", root.toString(), "--out", answers.toString()), stderr);
        assertEquals(2, run("study", "questions", "--project", root.toString(), "--out", answers.toString(), "--force"));
        Map<String, Object> document = new Yaml().load(Files.readString(answers));
        document.put("id", "annotation-answer");
        document.put("amendment_id", "annotation-reviewed");
        List<Object> responses = new ArrayList<>();
        for (Object raw : (List<?>) document.get("responses")) {
            Map<String, Object> response = new LinkedHashMap<>(cast(raw));
            Map<String, Object> q = questions().stream().filter(item -> item.get("question_id").equals(response.get("question_id"))).findFirst().orElseThrow();
            if (!String.valueOf(q.get("semantic_key")).contains("files.annotation")) continue;
            response.put("disposition", "answer");
            response.put("value", "annotation.tsv");
            response.put("rationale", "Reviewed platform annotation.");
            response.put("evidence", List.of(Map.of("note", "Synthetic fixture source release one; 1-based-closed GRCh38 locations.")));
            responses.add(response);
        }
        document.put("responses", responses);
        Files.writeString(root.resolve("annotation.tsv"), "feature_id\tgene_symbol\tchr\tposition\tgenome_build\np1\tGENE1\t1\t100\tGRCh38\n");
        Files.writeString(answers, new Yaml().dump(document));
        String ledger = Files.readString(root.resolve("ledger/events.jsonl"));
        assertEquals(2, run("study", "answer", answers.toString(), "--project", root.toString(), "--actor-id", "scientist"));
        assertEquals(2, run("study", "answer", answers.toString(), "--project", root.toString(), "--check"));
        assertEquals(0, run("study", "answer", answers.toString(), "--project", root.toString(), "--actor-id", "scientist", "--check"), stderr);
        assertEquals(ledger, Files.readString(root.resolve("ledger/events.jsonl")));
        assertEquals(0, run("study", "answer", answers.toString(), "--project", root.toString(), "--actor-id", "scientist", "--record"), stderr);
        assertFalse(Files.readString(root.resolve("resolved-study.yaml")).contains("annotation.tsv"));
        Path proposal = root.resolve("amendments/proposed/amendment-annotation-reviewed.yaml");
        assertEquals(0, run("study", "amend", proposal.toString(), "--project", root.toString(), "--check"), stderr);
        assertEquals(0, run("study", "amend", proposal.toString(), "--project", root.toString(), "--approved-by", "scientist"), stderr);
        assertTrue(Files.readString(root.resolve("resolved-study.yaml")).contains("annotation.tsv"));
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        assertEquals(0, run("study", "verify", "--project", root.toString()), stderr);
    }

    @Test void relationshipAnswerCannotBeReplacedByBlanketAssumptionAcceptance() throws Exception {
        prepare(true, false);
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        assertTrue(((Number) StudyProject.status(root).get("pending_assumptions")).intValue() > 0);
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline", "--accept-assumptions"), stderr);
        assertTrue(questions().stream().anyMatch(q -> "cohorts:a:b".equals(q.get("semantic_key"))));
        assertTrue(Files.readString(root.resolve("index/index.html")).contains("Scientific prerequisites require scientist input"));
        try (var files = Files.walk(root.resolve("runs"))) {
            assertFalse(files.anyMatch(p -> p.getFileName().toString().equals("meta-analysis.tsv")));
        }
        StudyProject.verify(root);
    }

    @Test void declaredDependenceRequiresSupportedAnalysisNotMereAcknowledgement() throws Exception {
        prepare(true, true);
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        assertTrue(questions().stream().anyMatch(q -> String.valueOf(q.get("semantic_key")).startsWith("dependent-inference:")));
        assertFalse(questions().stream().anyMatch(q -> "cohorts:a:b".equals(q.get("semantic_key"))));
    }

    @Test void reviewedRelationshipAmendmentUnblocksComparisonNotMissingGeneticInputs() throws Exception {
        prepare(true, false);
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        Map<String, Object> q = questions().stream().filter(item -> "cohorts:a:b".equals(item.get("semantic_key"))).findFirst().orElseThrow();
        Map<String, Object> response = Map.of("question_id", q.get("question_id"), "question_revision", q.get("revision"),
            "question_hash", q.get("question_hash"), "disposition", "answer", "value", "independent",
            "rationale", "Synthetic fixture cohorts are independently sampled despite a common study label.",
            "evidence", List.of(Map.of("note", "Reviewed synthetic sample design; no shared or related participants.")));
        Path answer = root.resolve("cohort-answer.yaml");
        Files.writeString(answer, new Yaml().dump(Map.of("schema_version", 1, "id", "cohort-answer",
            "project_id", "cli-questions", "amendment_id", "cohort-design", "responses", List.of(response))));
        assertEquals(0, run("study", "answer", answer.toString(), "--project", root.toString(), "--actor-id", "scientist", "--record"), stderr);
        assertEquals(0, run("study", "amend", root.resolve("amendments/proposed/amendment-cohort-design.yaml").toString(),
            "--project", root.toString(), "--approved-by", "scientist"), stderr);
        assertEquals(0, run("study", "auto", root.resolve("study.yaml").toString(), "--offline"), stderr);
        assertEquals(0, ((Number) StudyProject.status(root).get("pending_assumptions")).intValue());
        assertTrue(questions().stream().anyMatch(item -> "cohorts:a:b".equals(item.get("semantic_key")) && "resolved".equals(item.get("state"))));
        assertTrue(questions().stream().anyMatch(item -> String.valueOf(item.get("semantic_key")).contains("genetic-inputs")));
        try (var files = Files.walk(root.resolve("runs"))) {
            assertTrue(files.anyMatch(p -> p.getFileName().toString().equals("meta-analysis.tsv")));
        }
        StudyProject.verify(root);
    }

    @Test void questionInspectionOnNewProjectDoesNotCreateProject() throws Exception {
        Path absent = root.resolve("unresolved");
        assertEquals(0, run("study", "questions", "--project", absent.toString()), stderr);
        assertFalse(Files.exists(absent));
        assertTrue(stdout.contains("questions: []"));
    }

    @Test void diagnosticsNeverConfusesDatasetNameWithProducedReplicationArtifact() throws Exception {
        prepare(false, false);
        Path study = root.resolve("study.yaml");
        String original = Files.readString(study);
        String collision = original.replace("- id: b\n    cohort: b", "- id: replication\n    cohort: b")
            .replace("replication: [b]", "replication: [replication]")
            .replace("[discovery_replication, meta_analysis, genetic_support]", "[meta_analysis]");
        assertNotEquals(original, collision);
        assertTrue(collision.contains("- id: replication"));
        Files.writeString(study, collision);
        assertEquals(0, run("study", "auto", study.toString(), "--offline"), stderr + stdout);
        assertTrue(stdout.contains("failed: 0"), stdout);
        assertEquals("completed", StudyProject.status(root).get("state"));
    }

    private void prepare(boolean sharedParent, boolean declared) throws Exception {
        Files.writeString(root.resolve("a.tsv"), "feature_id\tbeta\tse\tp_value\np1\t1\t.1\t.001\n");
        Files.writeString(root.resolve("b.tsv"), "feature_id\tbeta\tse\tp_value\np1\t1.1\t.2\t.002\n");
        String cohorts = sharedParent ? "cohorts: [{id: a, parent: family}, {id: b, parent: family}]\n" : "cohorts: [{id: a}, {id: b}]\n";
        String relationships = declared ? "cohort_relationships: [{id: overlap, between: [a, b], type: overlapping}]\n" : "";
        Files.writeString(root.resolve("study.yaml"), """
            schema_version: 1
            project: {id: cli-questions, title: Questions, description: Result-first question tests.}
            resources: {downloads: deny, network: offline}
            %s%sdatasets:
              - id: a
                cohort: a
                type: transcriptomics
                stage: association_results
                feature_id_type: probe
                files: {results: a.tsv}
              - id: b
                cohort: b
                type: transcriptomics
                stage: association_results
                feature_id_type: probe
                files: {results: b.tsv}
            hypotheses:
              - id: H1
                description: Test association evidence.
                inputs: {discovery: [a], replication: [b]}
                analysis: {mode: manual, goals: [discovery_replication, meta_analysis, genetic_support]}
            """.formatted(cohorts, relationships));
    }

    private List<Map<String, Object>> questions() throws Exception {
        return ((List<?>) StudyQuestions.read(root).get("questions")).stream().map(StudyQuestionCliTest::cast).toList();
    }
    @SuppressWarnings("unchecked") private static Map<String, Object> cast(Object raw) { return (Map<String, Object>) raw; }
    private int run(String... args) {
        var out = new ByteArrayOutputStream();
        var err = new ByteArrayOutputStream();
        int code = JLinAlgCli.run(args, new PrintStream(out), new PrintStream(err));
        stdout = out.toString(); stderr = err.toString();
        return code;
    }
}
