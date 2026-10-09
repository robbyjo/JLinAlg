package org.jlinalg.study;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyEvidenceIndexTest {
    @TempDir Path directory;

    @Test void olderAliasedEvidenceUsesAuthoritativeManifestWithoutModifyingArtifacts() throws Exception {
        Path summary = directory.resolve("bioinformatics-investigation.yaml");
        String original = "evidence:\n- metadata: &shared {gene: PCBP2}\n- metadata: *shared\n";
        Files.writeString(summary, original);
        Map<String, Object> details = Map.of("evidence", List.of(
            Map.of("metadata", Map.of("gene", "PCBP2")),
            Map.of("metadata", Map.of("gene", "PCBP2"))));
        StudyIO.writeYaml(directory.resolve("manifest.yaml"), Map.of(
            "operation_type", "bioinformatics_investigation", "details", details));
        assertThrows(IllegalArgumentException.class, () -> StudyIO.load(summary));
        assertEquals(details, StudyEvidenceIndex.keyFindings(directory, "bioinformatics_investigation"));
        assertEquals(original, Files.readString(summary));
    }

    @Test void standaloneSummaryAndContractSafetyRemainStrict() throws Exception {
        Path summary = directory.resolve("literature-search.yaml");
        Files.writeString(summary, "evidence: []\n");
        assertEquals(Map.of("evidence", List.of()), StudyEvidenceIndex.keyFindings(directory, "literature_search"));
        Files.writeString(summary, "evidence:\n- metadata: &shared {gene: PCBP2}\n- metadata: *shared\n");
        assertThrows(IllegalArgumentException.class, () -> StudyEvidenceIndex.keyFindings(directory, "literature_search"));
    }
}
