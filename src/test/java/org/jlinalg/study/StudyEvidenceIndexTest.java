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

    @Test void displayPreviewIsBoundedAndDoesNotMutateAuthoritativeEvidence() {
        List<Object> variants = new java.util.ArrayList<>();
        for (int i = 0; i < 2000; i++) variants.add(Map.of("variant_id", "rs" + i, "audit", "x".repeat(4000)));
        Map<String, Object> evidence = Map.of("tested_features", 2000, "variants", variants);
        String displayed = StudyIO.canonicalJson(StudyEvidenceIndex.displayPreview(evidence));
        assertTrue(displayed.length() < 60000);
        assertTrue(displayed.contains("20 of 2000"));
        assertFalse(displayed.contains("rs1999"));
        assertEquals(2000, variants.size());
        assertEquals(4000, ((Map<?, ?>) variants.get(1999)).get("audit").toString().length());
        assertEquals(2000, evidence.get("tested_features"));
    }

    @Test void htmlDoesNotExpandFullAuditsIntoSearchAttributesOrTables() {
        List<Object> variants = new java.util.ArrayList<>();
        for (int i = 0; i < 2000; i++) variants.add(Map.of("variant_id", "rs" + i, "audit", "x".repeat(4000)));
        Map<String, Object> entry = Map.of("result_id", "run/MR", "operation_type", "mendelian_randomization",
            "hypothesis", "H-BMI", "state", "waiting", "summary", "Review input", "manifest", "runs/run/manifest.yaml",
            "current", true, "key_findings", Map.of("variants", variants), "details", Map.of("variants", variants), "artifacts", List.of());
        String html = StudyEvidenceIndex.html(Map.of("project", Map.of("id", "example")), List.of(entry));
        assertTrue(html.length() < 130000);
        assertFalse(html.contains("rs1999"));
        assertTrue(html.contains("20 of 2000"));
        assertTrue(html.contains("runs/run/manifest.yaml"));
        assertTrue(html.contains("data-search=\"run/MR mendelian_randomization H-BMI waiting Review input\""));
    }

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
