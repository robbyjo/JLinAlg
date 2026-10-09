package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;

class StudyResultsDiscussionTest {
    private Map<String, Object> evidence(String type) {
        return new LinkedHashMap<>(Map.of("operation_type", type, "state", "complete", "key_findings", Map.of("tested_features", 2)));
    }
    private Map<String, Object> result(String id, String role, String source) {
        return new LinkedHashMap<>(Map.of("id", id, "section", "results", "results_role", role, "text", "Observed association with uncertainty.", "sources", List.of(source)));
    }
    private Map<String, Object> finding(int rank, List<String> ids) {
        return new LinkedHashMap<>(Map.of("id", "discussion-" + rank, "section", "discussion", "discussion_role", "principal_finding",
            "impact_rank", rank, "impact_reason", "Replication with compatible effect direction motivates this priority, subject to uncertainty.",
            "result_ids", ids, "text", "This finding matters because it was independently reproduced, although its causal interpretation remains unresolved."));
    }
    private List<Object> tail() {
        List<Object> out = new ArrayList<>();
        for (String role : StudyResultsDiscussion.DISCUSSION.subList(1, 5)) out.add(Map.of("id", role,
            "section", role.equals("overall_conclusion") ? "conclusion" : "discussion", "discussion_role", role, "text", "Scope-appropriate interpretation."));
        return out;
    }
    @Test void fullModelAcceptsReasonedFindingsAndTailWithoutScientificCertification() throws Exception {
        var evidence = Map.of("H.association", evidence("association"), "H.meta", evidence("meta_analysis"), "H.coloc", evidence("colocalization"));
        var sections = new ArrayList<Object>(List.of(result("r1", "association", "H.association"), result("r2", "meta_analysis", "H.meta"), result("r3", "colocalization", "H.coloc"),
            finding(1, List.of("r2")), finding(2, List.of("r1")), finding(3, List.of("r3"))));
        sections.addAll(tail());
        var review = StudyResultsDiscussion.assess(Map.of("sections", sections), evidence);
        assertEquals("structure_complete_scientific_review_required", review.get("status"));
        assertTrue(review.get("scientific_review").toString().contains("do not validate scientific claims"));
    }
    @Test void resultsFollowMethodsWithNoFiveParagraphCapAndRejectMisclassifiedSources() throws Exception {
        Map<String, Map<String, Object>> evidence = new LinkedHashMap<>(); List<Object> sections = new ArrayList<>();
        for (String role : StudyResultsDiscussion.RESULTS.subList(0, 6)) {
            evidence.put(role, evidence(role)); sections.add(result("r-" + role, role, role));
        }
        assertDoesNotThrow(() -> StudyResultsDiscussion.assess(Map.of("sections", sections), evidence));
        java.util.Collections.swap(sections, 0, 1);
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", sections), evidence));
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(result("r", "colocalization", "association"))), evidence));
        var blocked = evidence("association"); blocked.put("state", "waiting_for_input");
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(result("r", "association", "a"))), Map.of("a", blocked)));
    }
    @Test void templatesNeverAutoRankAnalysesOrPromoteUnperformedOperations() throws Exception {
        var unperformed = evidence("gene_ontology"); unperformed.put("details", Map.of("method_report", Map.of("performed", false)));
        var blocked = evidence("colocalization"); blocked.put("state", "waiting_for_input");
        var evidence = Map.of("H.meta", evidence("meta_analysis"), "H.coloc", blocked, "H.go", unperformed);
        var template = StudyResultsDiscussion.template("Example", evidence);
        var sections = list(template.get("sections"));
        assertEquals(1, sections.stream().map(StudyPublicEvidence::map).filter(p -> "results".equals(p.get("section"))).count());
        assertEquals(List.of(), map(sections.get(1)).get("result_ids"));
        assertTrue(StudyResultsDiscussion.assess(template, evidence).toString().contains("link the discussed finding"));
        assertTrue(StudyResultsDiscussion.inventory(evidence).toString().contains("not an impact ranking"));
    }
    @Test void supportingEvidenceNeedsExplicitRelationshipAndFlagsPartialSources() throws Exception {
        var p = result("r", "meta_analysis", "H.meta");
        var bio = evidence("bioinformatics_investigation"); bio.put("state", "waiting_for_input");
        var evidence = Map.of("H.meta", evidence("meta_analysis"), "H.bio", bio);
        p.put("supporting_evidence", List.of(Map.of("source", "unknown", "relationship", "Same gene only.")));
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(p)), evidence));
        p.put("supporting_evidence", List.of(Map.of("source", "H.bio", "relationship", "")));
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(p)), evidence));
        p.put("supporting_evidence", List.of(Map.of("source", "H.bio", "relationship", "The retrieved protein annotation concerns the same gene, not the assayed isoform; sequence retrieval failed.")));
        var review = StudyResultsDiscussion.assess(Map.of("sections", List.of(p)), evidence);
        assertTrue(review.toString().contains("incomplete or empty"));
        assertTrue(review.toString().contains("usable source-level result"));
    }
    @Test void impactRanksLinksReasonsAndOneParagraphAreChecked() throws Exception {
        var evidence = Map.of("H.meta", evidence("meta_analysis")); var r = result("r", "meta_analysis", "H.meta");
        var p = finding(2, List.of("r"));
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(r, p)), evidence));
        p.put("impact_rank", 1); p.put("result_ids", List.of("missing"));
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(r, p)), evidence));
        p.put("result_ids", List.of("r")); p.put("text", "First paragraph.\n\nSecond paragraph.");
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(r, p)), evidence));
        p.put("text", "One paragraph."); p.remove("impact_reason");
        assertTrue(StudyResultsDiscussion.assess(Map.of("sections", List.of(r, p)), evidence).toString().contains("explain why"));
    }
    @Test void fewerThanThreeFindingsRequireReasonAndDoNotRequireInventedResults() throws Exception {
        var evidence = Map.of("H.meta", evidence("meta_analysis"));
        var sections = new ArrayList<Object>(List.of(result("r", "meta_analysis", "H.meta"), finding(1, List.of("r")))); sections.addAll(tail());
        assertEquals("needs_input", StudyResultsDiscussion.assess(Map.of("sections", sections), evidence).get("status"));
        assertEquals("structure_complete_scientific_review_required", StudyResultsDiscussion.assess(Map.of("sections", sections,
            "discussion_selection_note", "Only one finding has adequate evidence; the remaining analyses are unresolved."), evidence).get("status"));
    }
    @Test void misplacedRolesDuplicateClosingAndWrongTailOrderFail() {
        var wrong = new LinkedHashMap<>(map(tail().get(0))); wrong.put("section", "results");
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", List.of(wrong)), Map.of()));
        var sections = tail(); java.util.Collections.swap(sections, 0, 1);
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", sections), Map.of()));
        sections.clear(); var closing = new LinkedHashMap<>(map(tail().get(3))); closing.put("section", "discussion");
        sections.addAll(tail()); sections.add(closing);
        assertThrows(IOException.class, () -> StudyResultsDiscussion.assess(Map.of("sections", sections), Map.of()));
    }
    @Test void legacyTextAndRepeatedSourcePrioritiesRemainExplicitReviewNeeds() throws Exception {
        var legacy = Map.of("id", "old", "section", "results", "text", "Keep this prose unchanged.");
        var review = StudyResultsDiscussion.assess(Map.of("sections", List.of(legacy, finding(1, List.of("old")), finding(2, List.of("old")))), Map.of());
        assertTrue(review.toString().contains("older prose is preserved"));
        assertTrue(review.toString().contains("distinct findings, not repeated claims"));
        assertEquals("Keep this prose unchanged.", legacy.get("text"));
    }
}
