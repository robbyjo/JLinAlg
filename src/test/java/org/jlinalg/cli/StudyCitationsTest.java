package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;

class StudyCitationsTest {
    private Map<String, Object> catalog() {
        return StudyReferenceCatalog.build(List.of(Map.of("articles", List.of(
            Map.of("doi", "10.1234/first", "pmid", "123", "title", "First <study>", "authorString", "Smith A", "pubYear", "2020", "journalInfo", Map.of("volume", "7", "journal", Map.of("title", "Journal")), "pageInfo", "1-9"),
            Map.of("doi", "10.1234/second", "title", "Second study", "pubYear", "2021")))));
    }
    @Test void numbersAreFirstAppearanceDeduplicatedAndLinkToAnInPaperBibliography() throws Exception {
        var c = new StudyCitations(catalog());
        String body = c.render("Claim[[cite:10.1234/second]]. Other[[cite:10.1234/first;pmid:123]]. Again[[cite:10.1234/second]].", List.of());
        assertTrue(body.startsWith("Claim<sup class=\"citation\"><a href=\"#reference-1\""));
        assertEquals(2, body.split("href=\"#reference-1\"", -1).length - 1);
        assertEquals(1, body.split("href=\"#reference-2\"", -1).length - 1);
        String refs = c.bibliography(); assertTrue(refs.indexOf("Second study") < refs.indexOf("First &lt;study&gt;"));
        assertTrue(refs.contains("<b>7</b>, 1-9 (2020)")); assertTrue(refs.contains("id=\"references\""));
    }
    @Test void strictModeEnforcesLiteratureCitationsInEverySectionAndSupportingContext() throws Exception {
        for (String section : List.of("introduction", "methods", "results", "discussion", "conclusion")) {
            var p = Map.of("id", "p", "section", section, "citation_basis", "literature", "text", "Borrowed claim.");
            assertThrows(IOException.class, () -> StudyCitations.assess(Map.of("citation_policy", "strict", "sections", List.of(p)), catalog()));
        }
        var support = Map.of("source", "H.literature_search", "citation_basis", "literature", "relationship", "Claim.");
        var main = Map.of("id", "p", "citation_basis", "study_results", "text", "Our result.", "supporting_evidence", List.of(support));
        assertThrows(IOException.class, () -> StudyCitations.assess(Map.of("citation_policy", "strict", "sections", List.of(main)), catalog()));
        assertTrue(StudyCitations.assess(Map.of("citation_policy", "strict", "sections", List.of(Map.of("id", "p", "citation_basis", "study_design", "text", "Obesity was defined as BMI > 30."))), catalog()).get("needs").toString().equals("[]"));
    }
    @Test void unresolvedMalformedAndUnplacedCitationsCannotSilentlyPass() throws Exception {
        for (String text : List.of("[[cite:10.1234/missing]]", "[[cite:]]", "[[cite:10.1234/first;]]"))
            assertThrows(IOException.class, () -> StudyCitations.assess(Map.of("sections", List.of(Map.of("text", text))), catalog()));
        assertThrows(IOException.class, () -> StudyCitations.assess(Map.of("citation_policy", "strict", "sections", List.of(Map.of("text", "Claim.", "citation_basis", "literature", "references", List.of("10.1234/first")))), catalog()));
        assertTrue(StudyCitations.assess(Map.of("sections", List.of(Map.of("text", "Older prose."))), catalog()).get("needs").toString().contains("classify"));
        assertTrue(new StudyCitations(catalog()).render("<script>no</script>", List.of("pmid:123")).contains("&lt;script&gt;"));
    }
    @Test void narrativeChecksRecognizeInlineCitationsAndPayloadDoesNotRepeatPolicies() throws Exception {
        var intro = Map.of("id", "i", "section", "introduction", "introduction_role", "field_background", "text", "Field context[[cite:10.1234/first]].");
        assertFalse(StudyIntroduction.assess(List.of(intro)).get("needs").toString().contains("field_background: cite"));
        var methods = Map.of("id", "m", "section", "methods", "methods_role", "cohorts", "cohort_ids", List.of("A"), "text", "Cohort[[cite:10.1234/first]].");
        var inventory = Map.<String, Object>of("cohorts", List.of(Map.of("id", "A")));
        assertFalse(StudyMethods.assess(List.of(methods), inventory).get("needs").toString().contains("supply supporting cataloged citations"));
        var payload = StudyEditorialReview.payload(Map.of("title", "Review", "sections", List.of(intro, methods)), Map.of(), inventory, catalog());
        assertTrue(payload.containsKey("methods_policy")); assertTrue(payload.containsKey("citation_policy"));
        for (String field : List.of("methods_review", "introduction_review", "results_discussion_review", "citation_review")) {
            var review = (Map<?, ?>) payload.get(field);
            assertFalse(review.containsKey("policy")); assertFalse(review.containsKey("instructions"));
            assertTrue(review.containsKey("needs")); assertTrue(review.containsKey("status"));
        }
    }
}
