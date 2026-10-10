package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyPublicationTest {
    @TempDir Path directory;
    @Test void draftStartsResultsWithTableOneAndRetainsMissingMethodsAndHypothesisText() throws Exception {
        var dataset = Map.<String, Object>of("id", "d", "cohort", "c", "analysis_description", Map.of("model", "linear regression"));
        Context c = new Context(directory, directory, Map.of("hypothesis", "H-secret", "datasets", List.of("d")),
            Map.of("datasets", List.of(dataset), "hypotheses", List.of(Map.of("id", "H-secret", "description", "Isoform ratios associate with obesity."))), Map.of(), Map.of());
        var outputs = StudyPublication.write(c);
        assertEquals(14, outputs.size());
        String draft = Files.readString(directory.resolve("publication-draft.html"));
        assertTrue(draft.contains("Isoform ratios associate with obesity.")); assertFalse(draft.contains("H-secret"));
        assertTrue(draft.contains("<h2>Results</h2><h3>Table 1."));
        assertTrue(draft.contains("linear regression")); assertTrue(draft.contains("SCIENTIST INPUT REQUIRED"));
        assertTrue(draft.contains("most impactful finding"));
        assertFalse(draft.contains("Ranked publication points"));
        assertTrue(Files.readString(directory.resolve("results-discussion-review.yaml")).contains("needs_input"));
        assertTrue(Files.readString(directory.resolve("results-discussion-template.yaml")).contains("impact_reason"));
        assertTrue(draft.indexOf("most impactful finding") < draft.indexOf("second most impactful finding"));
        assertTrue(draft.indexOf("third most impactful finding") < draft.indexOf("conditional scientific"));
        assertTrue(draft.indexOf("conditional scientific") < draft.indexOf("Propose prioritized validation"));
        assertTrue(draft.contains("publication-supplement.html"));
        assertTrue(draft.contains("reference-catalog.html"));
        assertTrue(Files.readString(directory.resolve("introduction-template.yaml")).contains("importance_and_gap"));
        assertTrue(Files.readString(directory.resolve("introduction-review.yaml")).contains("needs_input"));
        String methods = draft.substring(draft.indexOf("<h2>Methods</h2>"), draft.indexOf("<h2>Results</h2>"));
        assertTrue(methods.indexOf("describe cohort design") < methods.indexOf("linear regression"));
        assertTrue(Files.isRegularFile(directory.resolve("methods-inventory.yaml")));
        assertTrue(Files.readString(directory.resolve("methods-template.yaml")).contains("methods_role: cohorts"));
    }
    @Test void nativeResultsUseMethodsOrderNotDependencyOrImpactOrder() throws Exception {
        var dependencies = new java.util.LinkedHashMap<String, Path>();
        for (String type : List.of("colocalization", "meta_analysis", "discovery_replication")) {
            Path folder = directory.resolve(type); Files.createDirectory(folder); dependencies.put(type, folder);
            StudyPublicEvidence.write(folder.resolve("manifest.yaml"), Map.of("operation_type", type, "state", "complete", "summary", "Summary of " + type));
            if (type.equals("meta_analysis")) Files.writeString(folder.resolve("meta-analysis.tsv"), "feature_id\tbeta\tse\tp_value\nGENE1\t.5\t.1\t.001\n");
            else if (type.equals("discovery_replication")) Files.writeString(folder.resolve("replication.tsv"), "feature_id\tdiscovery_selected\nGENE1\ttrue\n");
            else StudyPublicEvidence.write(folder.resolve("genetic-evidence.yaml"), Map.of("pairs", List.of(Map.of("status", "complete", "posterior", .8))));
        }
        var c = new Context(directory, directory, Map.of("hypothesis", "H", "datasets", List.of()),
            Map.of("hypotheses", List.of(Map.of("id", "H", "description", "Synthetic hypothesis."))), Map.of(), dependencies);
        StudyPublication.write(c);
        String draft = Files.readString(directory.resolve("publication-draft.html"));
        String results = draft.substring(draft.indexOf("<h2>Results</h2>"), draft.indexOf("<h2>Discussion</h2>"));
        assertTrue(results.indexOf("Table 1") < results.indexOf("discovery replication"));
        assertTrue(results.indexOf("discovery replication") < results.indexOf("meta analysis"));
        assertTrue(results.indexOf("meta analysis") < results.indexOf("colocalization"));
        assertTrue(draft.contains("GENE1"));
        assertFalse(Files.readString(directory.resolve("publication-claims.yaml")).contains("main_text_rank"));
    }
    @Test void publicationReadsLargeGeneratedOperationManifest() throws Exception {
        Path dependency = directory.resolve("large-operation");
        Path output = directory.resolve("large-report");
        Files.createDirectories(dependency); Files.createDirectories(output);
        Files.writeString(dependency.resolve("manifest.yaml"),
            "operation_type: custom_evidence\nstate: waiting_for_input\nsummary: Synthetic limitation\ndetails:\n  padding: "
                + "x".repeat(3_300_000) + "\n");
        Context context = new Context(directory, output,
            Map.of("hypothesis", "H", "datasets", List.of()),
            Map.of("hypotheses", List.of(Map.of("id", "H",
                "description", "Synthetic hypothesis."))), Map.of(),
            Map.of("large", dependency));
        StudyPublication.write(context);
        assertTrue(Files.readString(output.resolve("publication-draft.html"))
            .contains("Synthetic limitation"));
    }
    @Test void abbreviationsExpandAtFirstTextUseWithoutChangingLinksOrAlreadyExpandedText() {
        String result = StudyPublication.expandAbbreviations("<html><body><p>BMI predicts BMI.</p><a href=\"BMI.html\">BMI</a><p>Mendelian randomization (MR)</p><p>MR</p></body></html>", Map.of());
        assertTrue(result.contains("body mass index (BMI) predicts BMI"));
        assertTrue(result.contains("href=\"BMI.html\""));
        assertFalse(result.contains("Mendelian randomization (Mendelian"));
    }
    @Test void displayRoundingPreservesTinyValuesCountsAndIdentifiers() {
        assertEquals("8.01e-21", StudyReportNumbers.format(8.012839414774448E-21));
        assertEquals("8.01e-21", StudyReportNumbers.field("FDR", "8.012839414774448E-21"));
        assertEquals("8.02e-21", StudyReportNumbers.format("8.015e-21"));
        assertEquals("1.00e-3", StudyReportNumbers.format("0.001"));
        assertEquals("0.01", StudyReportNumbers.format("0.01"));
        assertEquals("223947", StudyReportNumbers.format(223947));
        assertEquals("30.12", StudyReportNumbers.format(30.123456));
        assertEquals("0.00", StudyReportNumbers.format(0.0));
        assertEquals("ENST000001.14", StudyReportNumbers.field("feature_id", "ENST000001.14"));
        assertEquals("123.456", StudyReportNumbers.field("feature_id", "123.456"));
        assertEquals("10.1000/182", StudyReportNumbers.field("doi", "10.1000/182"));
    }
    private Map<String, Object> summary(String dataset, String selection) {
        var table = new java.util.LinkedHashMap<String, Object>();
        table.put("dataset", dataset); table.put("cohort", "fhs"); table.put("N", 3);
        table.put("phenotype_sha256", "same-input"); table.put("selection_sha256", selection);
        table.put("continuous", Map.of("BMI", Map.of("mean", 26.123456789, "sd", 2.234567, "n", 2, "missing", 1)));
        table.put("binary", Map.of("Obesity", Map.of("case_n", 1, "percent", 50.0, "denominator", 2, "missing", 1, "definition", "Obesity is BMI > 30 kg/m².")));
        return table;
    }
    @Test void tableOneMergesSamePopulationAcrossPhenotypesWithoutSummingCounts() {
        var tables = List.<Object>of(summary("fhs-bmi", "same"), summary("fhs-obesity", "same"));
        String html = StudyPublication.table1(tables);
        assertEquals(1, html.split("scope=\"col\">FHS", -1).length - 1);
        assertTrue(html.contains("<td>3</td>")); assertFalse(html.contains("<td>6</td>"));
        assertTrue(html.contains("26.12 ± 2.23 [n=2]"));
        assertTrue(html.contains("Obesity, n (%)")); assertTrue(html.contains("1 (50.00%) [n=2]"));
        assertFalse(html.contains("missing=")); assertFalse(html.contains("26.123456789"));
        assertTrue(StudyPublication.table1Csv(tables).contains("26.123456789"));
        assertTrue(StudyPublication.detailedTable1(tables).contains("missing"));
    }
    @Test void differentSelectionsOrUnprovenPopulationsAreNotSilentlyMerged() {
        var a = summary("fhs-bmi", "a"); var b = summary("fhs-obesity", "b");
        assertTrue(StudyPublication.table1(List.of(a, b)).contains("FHS (fhs-obesity)"));
        b.put("selection_sha256", "a"); b.put("N", 2);
        assertTrue(StudyPublication.table1(List.of(a, b)).contains("FHS (fhs-obesity)"));
        a.remove("selection_sha256"); b.remove("selection_sha256");
        assertTrue(StudyPublication.table1(List.of(a, b)).contains("FHS (fhs-obesity)"));
    }
}
