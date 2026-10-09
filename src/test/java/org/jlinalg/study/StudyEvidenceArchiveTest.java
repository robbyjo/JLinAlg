package org.jlinalg.study;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.ZipFile;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyEvidenceArchiveTest {
    @TempDir Path root;
    @Test void exportsVerifiedLocalSnapshotWithoutRawInputsOrSilentOverwrite() throws Exception {
        Path study = root.resolve("study.yaml");
        Files.writeString(study, """
            schema_version: 1
            project: {id: archive-fixture, title: Archive fixture, description: Synthetic archive test, objective: publication}
            datasets:
              - id: results
                type: transcriptomics
                stage: association_results
                files: {results: private-input.tsv}
            hypotheses:
              - id: H1
                description: Test archive bookkeeping
                inputs: {discovery: [results]}
                analysis: {mode: manual, goals: []}
            """);
        Files.writeString(root.resolve("private-input.tsv"), "do not include original data automatically");
        var actor = new StudyProject.Actor("human", "fixture-scientist", null);
        StudyProject.resolve(study, root, actor);
        StudyEvidenceIndex.rebuild(root);
        Path output = root.resolve("evidence.zip");
        assertThrows(IllegalArgumentException.class, () -> StudyEvidenceArchive.export(root, output, actor, false));
        String sha = StudyEvidenceArchive.export(root, output, actor, true);
        assertEquals(64, sha.length());
        try (ZipFile zip = new ZipFile(output.toFile())) {
            assertNotNull(zip.getEntry("SESSION_HISTORY.md"));
            assertNotNull(zip.getEntry("ledger/events.jsonl"));
            assertNotNull(zip.getEntry("EVIDENCE_BUNDLE.json"));
            assertNull(zip.getEntry("private-input.tsv"));
        }
        StudyProject.verify(root);
        assertTrue(Files.readString(root.resolve("SESSION_HISTORY.md")).contains("Exported local evidence snapshot"));
        assertThrows(java.io.IOException.class, () -> StudyEvidenceArchive.export(root, output, actor, true));
    }
}
