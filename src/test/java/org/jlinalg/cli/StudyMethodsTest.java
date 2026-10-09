package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;

class StudyMethodsTest {
    private static Map<String, Object> study() {
        return Map.of("cohorts", List.of(Map.of("id", "fhs", "description", "Declared cohort design.", "references", List.of("pmid:123")),
            Map.of("id", "cardia"), Map.of("id", "unused")), "datasets", List.of(
            Map.of("id", "a", "cohort", "fhs", "stage", "association_results", "analysis_description", Map.of("model", "linear mixed model", "relatedness", "pedigree matrix")),
            Map.of("id", "b", "cohort", "fhs", "stage", "association_results"), Map.of("id", "c", "cohort", "cardia", "stage", "association_results")));
    }
    @Test void cohortInventoryDeduplicatesSamplesAndNeverGuessesPrivateNames() {
        var cohorts = StudyMethods.cohorts(study());
        assertEquals(2, cohorts.size()); assertEquals("Framingham Heart Study", cohorts.get(0).get("name"));
        assertEquals("", StudyMethods.literatureName(Map.of("id", "PRIVATE_ID", "description", "Confidential recruitment notes")));
        assertEquals("New public cohort", StudyMethods.literatureName(Map.of("id", "c", "literature_name", "New public cohort")));
    }
    @Test void templateStartsWithEveryCohortAndDoesNotInventCitationsOrBlockedMethods() throws Exception {
        var evidence = Map.of("H.meta_analysis", Map.<String, Object>of("operation_type", "meta_analysis", "state", "complete", "manifest", "source.yaml"),
            "H.colocalization", Map.<String, Object>of("operation_type", "colocalization", "state", "waiting_for_input"),
            "H.gene_ontology", Map.<String, Object>of("operation_type", "gene_ontology", "state", "complete", "details", Map.of("method_report", Map.of("performed", false))));
        var inventory = StudyMethods.inventory(study(), evidence);
        var template = StudyMethods.template("Fixture", inventory, StudyReferenceCatalog.build(List.of()));
        var sections = StudyPublicEvidence.list(template.get("sections"));
        assertEquals("cohorts", StudyPublicEvidence.map(sections.get(0)).get("methods_role"));
        assertEquals("cohorts", StudyPublicEvidence.map(sections.get(1)).get("methods_role"));
        assertTrue(sections.toString().contains("not refitted"));
        assertFalse(sections.toString().contains("methods_role=colocalization"));
        assertFalse(sections.toString().contains("methods_role=enrichment_ontology"));
        assertFalse(sections.toString().contains("has been described previously"));
        assertTrue(StudyMethods.assess(sections, inventory).toString().contains("recover missing historical"));
    }
    @Test void onlyRetrievedAssignedCohortReferenceIsCitedNotAnArbitrarySearchHit() {
        var catalog = StudyReferenceCatalog.build(List.of(Map.of("articles", List.of(Map.of("pmid", "123", "title", "Fixture cohort design")))));
        var template = StudyMethods.template("Fixture", StudyMethods.inventory(study(), Map.of()), catalog);
        var first = StudyPublicEvidence.map(StudyPublicEvidence.list(template.get("sections")).get(0));
        assertEquals(List.of("pmid:123"), first.get("references"));
        assertTrue(first.get("text").toString().contains("has been described previously"));
        assertTrue(StudyMethods.scaffold(template, catalog).contains("reference-catalog.html#ref-"));
    }
    @Test void typedMethodsEnforceOrderAndRequireAllCohortsParametersAndJustifications() throws Exception {
        var inventory = StudyMethods.inventory(study(), Map.of());
        var paragraphs = new ArrayList<Object>();
        for (String cohort : List.of("fhs", "cardia")) paragraphs.add(Map.of("section", "methods", "methods_role", "cohorts",
            "cohort_ids", List.of(cohort), "text", "Cohort description.", "references", List.of("pmid:123")));
        paragraphs.add(Map.of("section", "methods", "methods_role", "association", "text", "Association method.", "references", List.of("pmid:456"),
            "parameters", Map.of("model", "linear"), "parameter_justifications", Map.of("model", "Not recorded in the upstream source.")));
        assertEquals("needs_input", StudyMethods.assess(paragraphs, inventory).get("status")); // Missing upstream method for dataset b/c.
        java.util.Collections.swap(paragraphs, 0, 2);
        assertThrows(IOException.class, () -> StudyMethods.assess(paragraphs, inventory));
        for (var invalid : List.of(Map.of("section", "results", "methods_role", "cohorts"), Map.of("section", "methods", "methods_role", "invented"),
                Map.of("section", "methods", "methods_role", "cohorts", "cohort_ids", List.of("wrong"))))
            assertThrows(IOException.class, () -> StudyMethods.assess(List.of(invalid), inventory));
    }
    @Test void legacyProseIsPreservedAndUnsupportedExecutionRemainsAReviewNeed() throws Exception {
        var legacy = Map.of("section", "methods", "text", "Previous Methods prose.");
        var review = StudyMethods.assess(List.of(legacy), StudyMethods.inventory(study(), Map.of()));
        assertTrue(review.toString().contains("prose is preserved")); assertTrue(review.toString().contains("Missing Methods description for cohort fhs"));
        var network = Map.of("section", "methods", "methods_role", "network", "text", "Network methods.", "references", List.of("pmid:1"));
        assertTrue(StudyMethods.assess(List.of(network), StudyMethods.inventory(study(), Map.of())).toString().contains("no completed or explicitly supplied execution evidence"));
    }
    @Test void methodRecordDistinguishesDefaultsFromScientificRationale() {
        var record = StudyMethods.record("Example", Map.of("alpha", 0.05), List.of(), null, List.of());
        assertTrue(record.get("parameter_justification").toString().contains("Not recorded"));
        assertTrue(StudyMethods.instructions().contains("pedigree/relatedness"));
        assertTrue(StudyMethods.instructions().contains("network source/version"));
        assertTrue(StudyMethods.instructions().contains("all relevant effective parameters"));
    }
    @Test void customaryParametersNeedDisclosureNotInventedJustificationAndTissuesArePreserved() throws Exception {
        var section = Map.of("section", "methods", "methods_role", "association", "text", "FDR < 0.05.", "references", List.of("pmid:1"), "parameters", Map.of("FDR", .05));
        var inventory = StudyMethods.inventory(Map.of("datasets", List.of(Map.of("id", "jhs-bmi", "cohort", "jhs", "tissue", "peripheral_blood_mononuclear_cells", "stage", "association_results"))), Map.of());
        assertTrue(inventory.toString().contains("peripheral_blood_mononuclear_cells"));
        assertFalse(StudyMethods.assess(List.of(section), inventory).get("needs").toString().contains("justif"));
        var unusual = new LinkedHashMap<String, Object>(section); unusual.put("parameters", Map.of("window", 12345)); unusual.put("uncommon_parameters", List.of("window"));
        assertTrue(StudyMethods.assess(List.of(unusual), inventory).toString().contains("justify uncommon parameter window"));
        unusual.put("parameter_justifications", Map.of("window", "Approved study-specific locus definition."));
        assertFalse(StudyMethods.assess(List.of(unusual), inventory).get("needs").toString().contains("justify uncommon"));
    }
}
