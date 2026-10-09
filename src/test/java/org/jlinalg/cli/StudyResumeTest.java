/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;
import org.jlinalg.study.StudyExecutionEngine;
import org.jlinalg.study.StudyOperationExecutor;
import org.jlinalg.study.StudyPlanner;
import org.jlinalg.study.StudyAmendments;
import org.jlinalg.study.StudyProject;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyResumeTest {
    @TempDir Path root;
    private final StudyProject.Actor actor = StudyProject.Actor.human("fixture-scientist");

    @Test void retriesOnlyUnfinishedWorkAndAffectedDependantsPreservingEveryAttempt() throws Exception {
        prepare("");
        AtomicBoolean unavailable = new AtomicBoolean(true);
        StudyOperationExecutor executor = context -> {
            if ("meta_analysis".equals(context.operation().get("type")) && unavailable.getAndSet(false))
                return StudyOperationExecutor.Result.waiting("Connection fixture unavailable", Map.of("connection_status", "connection_timeout"));
            return new StudyCliExecutor().execute(context);
        };
        var first = StudyExecutionEngine.run(root, actor, executor, true);
        assertEquals(1, first.waiting());
        Path original = root.resolve("runs").resolve(first.runId());
        byte[] manifest = Files.readAllBytes(find(original, "manifest.yaml"));
        var resumed = StudyExecutionEngine.run(root, actor, executor, true);
        assertEquals("completed", resumed.state());
        assertEquals(1, resumed.reused());
        assertNotEquals(first.runId(), resumed.runId());
        assertArrayEquals(manifest, Files.readAllBytes(find(original, "manifest.yaml")));
        var repeated = StudyExecutionEngine.run(root, actor, executor, true);
        assertEquals(4, repeated.reused());
        assertEquals(0, repeated.completed());
        var redo = StudyExecutionEngine.run(root, actor, executor, false);
        assertEquals(4, redo.completed());
        assertEquals(0, redo.reused());
        StudyProject.verify(root);
        Map<?, ?> index = (Map<?, ?>) SimpleJson.parse(Files.readString(root.resolve("index/index.json")));
        assertTrue(((List<?>) index.get("blockers")).isEmpty(), "Historical waiting attempts are not current blockers");
        assertTrue(Files.readString(root.resolve("index/index.html")).contains("Show historical attempts"));
    }

    @Test void changedInputsInvalidateReuseAndTamperedEvidenceIsNeverOverwritten() throws Exception {
        prepare("");
        var first = StudyExecutionEngine.run(root, actor, new StudyCliExecutor(), true);
        Files.writeString(root.resolve("results.tsv"), "feature_id\tbeta\tse\tp_value\nGENE1\t2\t.1\t.001\n");
        var changed = StudyExecutionEngine.run(root, actor, new StudyCliExecutor(), true);
        assertEquals(0, changed.reused());
        assertEquals(4, changed.completed());
        Path original = find(root.resolve("runs").resolve(first.runId()), "meta-analysis.tsv");
        Files.writeString(original, "tampered evidence");
        assertThrows(IOException.class, () -> StudyExecutionEngine.run(root, actor, new StudyCliExecutor(), false));
        assertEquals("tampered evidence", Files.readString(original));
    }

    @Test void dictionaryRolesAreResolvedAndChangedDictionaryRequiresAmendment() throws Exception {
        Files.writeString(root.resolve("dictionary.yaml"), "columns: {effect: beta, standard_error: se, p_value: p_value, feature_id: feature_id}\ndescriptions: {beta: Association effect}\n");
        prepare(", data_dictionary: dictionary.yaml");
        String resolved = Files.readString(root.resolve("resolved-study.yaml"));
        assertTrue(resolved.contains("resolved_data_dictionary"));
        assertTrue(resolved.contains("Association effect"));
        StudyExecutionEngine.run(root, actor, new StudyCliExecutor(), true);
        Files.writeString(root.resolve("dictionary.yaml"), "columns: {effect: wrong}\n");
        IOException changed = assertThrows(IOException.class, () -> StudyExecutionEngine.run(root, actor, new StudyCliExecutor(), true));
        assertTrue(changed.getMessage().contains("Data dictionary changed"));
    }

    @Test void assistanceAndOfflinePolicyCanBeAmendedWithoutEditingBaseStudy() throws Exception {
        prepare("");
        String base = Files.readString(root.resolve("study.yaml"));
        Path amendment = root.resolve("amendment-assist.yaml");
        Files.writeString(amendment, """
            schema_version: 1
            amendment:
              id: assist
              base_study_hash: %s
            changes:
              configure:
                assistance: {enabled: true, provider: openai, model: test-model, next_steps: true}
                resources: {network: offline}
                project:
                  data_policy: {llm: {remote_transmission: deny_all}}
            """.formatted(StudyProject.status(root).get("resolved_study_hash")));
        StudyAmendments.apply(amendment, root, actor);
        assertEquals(base, Files.readString(root.resolve("study.yaml")));
        var plan = StudyPlanner.create(root, actor);
        StudyPlanner.approve(root, plan.id(), actor, true);
        var run = StudyExecutionEngine.run(root, actor, new StudyCliExecutor(), true);
        assertEquals(1, run.waiting());
        assertTrue(Files.readString(root.resolve("index/index.html")).contains("policy_denied"));
        StudyProject.verify(root);
    }

    private void prepare(String dictionary) throws IOException {
        Files.writeString(root.resolve("results.tsv"), "feature_id\tbeta\tse\tp_value\nGENE1\t1\t.1\t.001\n");
        Path study = root.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project: {id: fixture, title: Resume fixture, description: Immutable execution attempts.}
            resources: {downloads: deny}
            datasets:
              - id: data
                type: transcriptomics
                stage: association_results
                feature_id_type: gene_symbol
                files: {results: results.tsv%s}
            hypotheses:
              - id: H1
                description: Test the result-first path.
                outcome: BMI
                analysis: {mode: manual, goals: [meta_analysis]}
            """.formatted(dictionary));
        StudyProject.resolve(study, root, actor);
        var plan = StudyPlanner.create(root, actor);
        StudyPlanner.approve(root, plan.id(), actor, true);
    }
    private static Path find(Path folder, String file) throws IOException {
        try (var paths = Files.walk(folder)) { return paths.filter(path -> path.getFileName().toString().equals(file)).findFirst().orElseThrow(); }
    }
}
