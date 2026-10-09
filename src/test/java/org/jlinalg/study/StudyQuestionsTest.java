/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.study;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyQuestionsTest {
    @TempDir Path root;
    private final StudyProject.Actor actor = StudyProject.Actor.human("scientist-fixture");

    @Test void questionsAreDocumentedBeforePresentationAndPreviewIsReadOnly() throws Exception {
        prepareMatrix();
        assertEquals(2, questions().size());
        String history = Files.readString(root.resolve("SESSION_HISTORY.md"));
        for (Map<String, Object> q : questions()) {
            assertEquals("open", q.get("state"));
            assertTrue(history.contains(String.valueOf(map(q.get("need")).get("prompt"))));
            assertNotNull(q.get("documented_at"));
            assertTrue(Files.isRegularFile(root.resolve(String.valueOf(q.get("artifact")))));
        }
        Path exported = root.resolve("blank.yaml");
        StudyQuestions.export(root, exported, false);
        for (Object raw : list(StudyIO.load(exported).get("responses"))) {
            assertNull(map(raw).get("disposition"));
            assertNull(map(raw).get("value"));
        }
        assertThrows(IOException.class, () -> StudyQuestions.export(root, exported, false));
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.preview(root, exported, actor));
        Path answer = matrixAnswers("initial", "metadata");
        Map<String, String> before = snapshot();
        Map<String, Object> preview = StudyQuestions.preview(root, answer, actor);
        assertTrue(preview.containsKey("amendment"));
        assertEquals(before, snapshot(), "--check must not write, regenerate, or execute");
        StudyProject.verify(root);
    }

    @Test void answersProduceAmendmentOnlyAndIndependentBranchesContinue() throws Exception {
        prepareMatrix();
        AtomicInteger associationCalls = new AtomicInteger();
        StudyOperationExecutor executor = context -> {
            if (context.operation().get("type").equals("association")) associationCalls.incrementAndGet();
            return StudyOperationExecutor.Result.complete("fixture", List.of());
        };
        var first = StudyExecutionEngine.run(root, actor, executor, true, true);
        assertTrue(first.waiting() > 0);
        assertTrue(first.completed() > 0, "Nonblocked validation/reporting can continue");
        assertEquals(0, associationCalls.get());
        assertEquals(2, questions().size(), "Planning and runtime should not duplicate the same need");
        String base = Files.readString(root.resolve("resolved-study.yaml"));
        Path answers = matrixAnswers("initial", "metadata");
        Map<String, Object> receipt = StudyQuestions.record(root, answers, actor);
        assertEquals(base, Files.readString(root.resolve("resolved-study.yaml")));
        assertEquals(0, associationCalls.get());
        assertTrue(questions().stream().allMatch(q -> "awaiting_approval".equals(q.get("state"))));
        Path proposal = root.resolve(String.valueOf(map(receipt.get("amendment")).get("path")));
        StudyAmendments.preview(proposal, root);
        StudyAmendments.apply(proposal, root, actor);
        assertTrue(questions().stream().allMatch(q -> "resolved".equals(q.get("state"))));
        var plan = StudyPlanner.create(root, actor);
        StudyPlanner.approve(root, plan.id(), actor, true);
        var second = StudyExecutionEngine.run(root, actor, executor, true, true);
        assertEquals(0, second.waiting());
        assertEquals(1, associationCalls.get());
        StudyProject.verify(root);
    }

    @Test void duplicateBatchIsIdempotentButConflictingIdentityFails() throws Exception {
        prepareMatrix();
        Path answers = matrixAnswers("initial", "metadata");
        StudyQuestions.record(root, answers, actor);
        Map<String, String> before = snapshot();
        StudyQuestions.record(root, answers, actor);
        assertEquals(before, snapshot());
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, answers, StudyProject.Actor.human("someone-else")));
        Map<String, Object> altered = StudyIO.load(answers);
        altered.put("meeting_id", "different-meeting");
        StudyIO.writeYaml(answers, altered);
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, answers, actor));
    }

    @Test void invalidBatchCommitsNothing() throws Exception {
        prepareMatrix();
        Path answers = matrixAnswers("initial", "metadata");
        Map<String, Object> document = StudyIO.load(answers);
        List<Object> responses = list(document.get("responses"));
        Map<String, Object> invalid = map(responses.get(1));
        invalid.put("question_hash", "wrong");
        responses.set(1, invalid);
        document.put("responses", responses);
        StudyIO.writeYaml(answers, document);
        Map<String, String> before = snapshot();
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, answers, actor));
        assertEquals(before, snapshot());
        assertFalse(Files.exists(root.resolve("answers")));
    }

    @Test void changedDatasetInputsRejectAnswersAndCreateNewRevisions() throws Exception {
        prepareMatrix();
        Path answer = matrixAnswers("initial", "metadata");
        Files.writeString(root.resolve("matrix.tsv"), "sample_id\tp1\ns1\t9\n");
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, answer, actor));
        StudyExecutionEngine.run(root, actor, context -> StudyOperationExecutor.Result.complete("fixture", List.of()), true, true);
        assertEquals(4, questions().size());
        assertEquals(2, questions().stream().filter(q -> "superseded".equals(q.get("state"))).count());
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, answer, actor));
    }

    @Test void unknownAndDeferredNeedsStayBlockedWithoutRepeatedPrompts() throws Exception {
        prepareMatrix();
        List<Object> responses = new ArrayList<>();
        for (Map<String, Object> q : questions()) responses.add(response(q,
            responses.isEmpty() ? "unknown" : "defer", null));
        Path answer = saveBatch("uncertain", null, responses);
        StudyQuestions.record(root, answer, actor);
        int events = StudyLedger.read(root).size();
        StudyQuestions.reconcile(root, StudyPlanner.loadApproved(root), null);
        assertEquals(events, StudyLedger.read(root).size());
        Path next = root.resolve("next.yaml");
        StudyQuestions.export(root, next, false);
        assertTrue(list(StudyIO.load(next).get("responses")).isEmpty());
        var run = StudyExecutionEngine.run(root, actor, context -> StudyOperationExecutor.Result.complete("fixture", List.of()), true, true);
        assertTrue(run.waiting() > 0);
        assertEquals(2, questions().size());
    }

    @Test void mandatoryNeedsCannotBeDeclinedAndModelCannotAnswerAsScientist() throws Exception {
        prepareMatrix();
        Path declined = saveBatch("declined", null, List.of(response(questions().get(0), "decline", null)));
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, declined, actor));
        Path factual = matrixAnswers("initial", "metadata");
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, factual, new StudyProject.Actor("llm", "advisor", "model-fixture")));
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, factual, StudyProject.Actor.human("unspecified")));
    }

    @Test void optionalReviewCanBeDeclinedWithoutClaimingNegativeEvidence() throws Exception {
        prepareMatrix();
        StudyExecutionEngine.run(root, actor, context -> {
            if (!"diagnostic_review".equals(context.operation().get("type"))) return StudyOperationExecutor.Result.complete("fixture", List.of());
            return new StudyOperationExecutor.Result(StudyOperationExecutor.State.COMPLETE, "fixture", List.of(),
                Map.of("needs", List.of(StudyNeed.review("review-fixture", "Review <script>finding</script>?", "A review is optional, not a biological conclusion.", false).document())));
        }, true, true);
        Map<String, Object> question = questions().stream().filter(q -> "review-fixture".equals(q.get("semantic_key"))).findFirst().orElseThrow();
        Path answer = saveBatch("declined", null, List.of(response(question, "decline", null)));
        StudyQuestions.record(root, answer, actor);
        assertTrue(questions().stream().anyMatch(q -> "declined".equals(q.get("state"))));
        assertTrue(Files.readString(root.resolve("QUESTIONS.md")).contains("not a biological negative"));
        String html = Files.readString(root.resolve("index/index.html"));
        assertTrue(html.contains("&lt;script&gt;finding&lt;/script&gt;"));
        assertFalse(html.contains("<script>finding</script>"));
        StudyProject.verify(root);
    }

    @Test void correctionSupersedesAnswerAndInvalidatesOldProposal() throws Exception {
        prepareMatrix();
        Path answers = matrixAnswers("initial", "metadata");
        var receipt = StudyQuestions.record(root, answers, actor);
        List<Object> responses = new ArrayList<>();
        for (Map<String, Object> q : questions()) {
            var response = response(q, "unknown", null);
            response.put("supersedes", "initial");
            responses.add(response);
        }
        Path correction = saveBatch("correction", null, responses);
        StudyQuestions.record(root, correction, actor);
        Path proposal = root.resolve(String.valueOf(map(receipt.get("amendment")).get("path")));
        assertThrows(IllegalArgumentException.class, () -> StudyAmendments.apply(proposal, root, actor));
        assertTrue(questions().stream().allMatch(q -> list(q.get("answers")).size() == 2));
        StudyProject.verify(root);
    }

    @Test void evidenceIsSnapshottedTamperDetectedAndLocalAnswersExcludedFromLlmContext() throws Exception {
        prepareMatrix();
        Files.writeString(root.resolve("metadata.txt"), "LOCAL-PRIVATE-NOTE-84931");
        Map<String, Object> r = response(questions().get(0), "unknown", null);
        r.put("rationale", "LOCAL-PRIVATE-ANSWER-62734");
        r.put("evidence", List.of(Map.of("file", "metadata.txt", "selector", "whole note")));
        Path answers = saveBatch("private", null, List.of(r));
        StudyQuestions.record(root, answers, actor);
        Files.writeString(root.resolve("metadata.txt"), "changed source; original evidence remains captured");
        StudyProject.verify(root);
        assertFalse(Files.readString(root.resolve("index/LLM_CONTEXT.md")).contains("LOCAL-PRIVATE-ANSWER"));
        assertFalse(Files.readString(root.resolve("index/LLM_CONTEXT.md")).contains("LOCAL-PRIVATE-NOTE"));
        assertTrue(Files.readString(root.resolve("QUESTIONS.md")).contains("LOCAL-PRIVATE-ANSWER"));
        Path snapshot;
        try (var files = Files.list(root.resolve("answers/evidence"))) { snapshot = files.findFirst().orElseThrow(); }
        assertEquals("LOCAL-PRIVATE-NOTE-84931", Files.readString(snapshot));
        Files.writeString(snapshot, "tampered");
        assertThrows(IOException.class, () -> StudyProject.verify(root));
        assertThrows(IOException.class, () -> StudyProject.repairViews(root));
    }

    @Test void missingDerivedViewsCanBeRepairedButImmutableQuestionCannot() throws Exception {
        prepareMatrix();
        Files.delete(root.resolve("QUESTIONS.md"));
        assertThrows(IOException.class, () -> StudyProject.verify(root));
        String ledger = Files.readString(root.resolve("ledger/events.jsonl"));
        StudyProject.repairViews(root);
        assertEquals(ledger, Files.readString(root.resolve("ledger/events.jsonl")));
        StudyProject.verify(root);
        Files.writeString(root.resolve(String.valueOf(questions().get(0).get("artifact"))), "tampered");
        assertThrows(IOException.class, () -> StudyProject.repairViews(root));
    }

    @Test void changingAnsweredResourceBeforeApprovalIsRejected() throws Exception {
        prepareMatrix();
        var receipt = StudyQuestions.record(root, matrixAnswers("initial", "metadata"), actor);
        Files.writeString(root.resolve("phenotypes.tsv"), "changed");
        Path proposal = root.resolve(String.valueOf(map(receipt.get("amendment")).get("path")));
        assertThrows(IllegalArgumentException.class, () -> StudyAmendments.apply(proposal, root, actor));
    }

    @Test void uncommittedArtifactsHaveNoAuthorityAndIdenticalRetryCanFinish() throws Exception {
        prepareMatrix();
        Path answers = matrixAnswers("initial", "metadata");
        var preview = StudyQuestions.preview(root, answers, actor);
        Map<String, Object> stored = StudyIO.copyMap(preview);
        stored.remove("answer_file");
        StudyIO.writeYaml(root.resolve("answers/initial.yaml"), stored);
        assertTrue(questions().stream().allMatch(q -> "open".equals(q.get("state"))));
        StudyQuestions.record(root, answers, actor);
        StudyProject.verify(root);
    }

    @Test void refusesTraversalUnsupportedFieldsAndConcurrentMutation() throws Exception {
        prepareMatrix();
        Map<String, Object> q = questions().stream().filter(item -> String.valueOf(item.get("semantic_key")).endsWith("files.phenotype")).findFirst().orElseThrow();
        var traversal = response(q, "answer", "../outside.tsv");
        Path answers = saveBatch("escape", "metadata", List.of(traversal));
        assertThrows(IOException.class, () -> StudyQuestions.record(root, answers, actor));
        Path good = matrixAnswers("initial", "metadata");
        try (var channel = FileChannel.open(root.resolve(".study-run.lock"), StandardOpenOption.WRITE);
                var lock = channel.lock()) {
            assertTrue(lock.isValid());
            assertThrows(IOException.class, () -> StudyQuestions.record(root, good, actor));
            assertThrows(IOException.class, () -> StudyPlanner.create(root, actor));
        }
        var doc = StudyIO.load(good);
        doc.put("command", "arbitrary-code");
        StudyIO.writeYaml(good, doc);
        assertThrows(IllegalArgumentException.class, () -> StudyQuestions.record(root, good, actor));
    }

    private void prepareMatrix() throws Exception {
        Files.writeString(root.resolve("matrix.tsv"), "sample_id\tp1\ns1\t1\n");
        Files.writeString(root.resolve("phenotypes.tsv"), "sample_id\tBMI\ns1\t25\n");
        Files.writeString(root.resolve("study.yaml"), """
            schema_version: 1
            project: {id: question-fixture, title: Questions, description: Offline scientist question tests.}
            resources: {downloads: deny, network: offline}
            datasets:
              - id: data
                type: transcriptomics
                stage: processed_matrix
                files: {matrix: matrix.tsv}
            hypotheses:
              - id: H1
                description: Test associations.
                analysis: {mode: manual, goals: [association]}
            """);
        StudyProject.resolve(root.resolve("study.yaml"), root, actor);
        var plan = StudyPlanner.create(root, actor);
        StudyPlanner.approve(root, plan.id(), actor, true);
    }

    private Path matrixAnswers(String id, String amendment) throws Exception {
        List<Object> responses = new ArrayList<>();
        for (Map<String, Object> q : questions()) {
            String key = String.valueOf(q.get("semantic_key"));
            responses.add(response(q, "answer", key.endsWith("outcome") ? "BMI" : "phenotypes.tsv"));
        }
        return saveBatch(id, amendment, responses);
    }

    private Map<String, Object> response(Map<String, Object> q, String disposition, String value) {
        Map<String, Object> response = new LinkedHashMap<>();
        response.put("question_id", q.get("question_id"));
        response.put("question_revision", q.get("revision"));
        response.put("question_hash", q.get("question_hash"));
        response.put("disposition", disposition);
        if (value != null) response.put("value", value);
        response.put("rationale", "Fixture scientist reviewed metadata.");
        response.put("evidence", List.of(Map.of("note", "Source is the fixture data dictionary, reviewed locally.")));
        return response;
    }

    private Path saveBatch(String id, String amendment, List<Object> responses) throws IOException {
        Map<String, Object> document = new LinkedHashMap<>();
        document.put("schema_version", 1);
        document.put("id", id);
        document.put("project_id", "question-fixture");
        if (amendment != null) document.put("amendment_id", amendment);
        document.put("responses", responses);
        Path path = root.resolve(id + "-input.yaml");
        StudyIO.writeYaml(path, document);
        return path;
    }

    private List<Map<String, Object>> questions() throws IOException {
        return list(StudyQuestions.read(root).get("questions")).stream().map(StudyQuestionsTest::map).toList();
    }
    private Map<String, String> snapshot() throws IOException {
        Map<String, String> result = new LinkedHashMap<>();
        try (var paths = Files.walk(root)) {
            for (Path file : paths.filter(Files::isRegularFile).sorted().toList())
                result.put(root.relativize(file).toString(), StudyIO.hash(file));
        }
        return result;
    }
    private static Map<String, Object> map(Object value) { return StudyIO.map(value, "fixture"); }
    private static List<Object> list(Object value) { return StudyIO.list(value, "fixture"); }
}
