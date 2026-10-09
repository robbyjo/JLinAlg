package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;

class StudyReferenceCatalogTest {
    private static Map<String, Object> query(String strategy, Map<String, Object> article) {
        return Map.of("search_strategy", strategy, "status", "candidate_matches", "articles", List.of(article),
            "query", strategy + " public terms", "retrieval", Map.of("response_sha256", "source-hash", "started_at", "2026-10-07T00:00:00Z"));
    }
    @Test void aliasBridgeDeduplicatesWithoutDiscardingQueriesOrMetadataConflicts() {
        var catalog = StudyReferenceCatalog.build(List.of(
            query("field_background", Map.of("doi", "https://doi.org/10.1000/ABC", "title", "First title", "authorString", "A Author", "pubYear", "2025")),
            query("similar_studies", Map.of("pmid", "123", "title", "Updated title")),
            query("hypothesis_relevance", Map.of("doi", "10.1000/abc", "pmid", "123", "pmcid", "PMC456"))));
        var refs = StudyPublicEvidence.list(catalog.get("references")); assertEquals(1, refs.size());
        var ref = StudyPublicEvidence.map(refs.get(0));
        assertEquals(3, StudyPublicEvidence.list(ref.get("retrieval_context")).size());
        assertTrue(ref.get("metadata_conflicts").toString().contains("Updated title"));
        for (String id : List.of("doi:10.1000/ABC", "https://dx.doi.org/10.1000/abc", "pmid:123", "pmcid:pmc456", ref.get("reference_id").toString()))
            assertEquals(ref, StudyReferenceCatalog.lookup(catalog, id));
        assertEquals("A Author", ref.get("authorString"));
        assertEquals("unreviewed_candidate_not_hypothesis_support", ref.get("review_status"));
        var reversed = StudyReferenceCatalog.build(List.of(query("x", Map.of("pmid", "123", "doi", "10.1000/abc", "pmcid", "PMC456"))));
        assertEquals(ref.get("reference_id"), StudyPublicEvidence.map(StudyPublicEvidence.list(reversed.get("references")).get(0)).get("reference_id"));
    }
    @Test void sameTitleIsNotIdentityAndUnidentifiedItemsCannotBeCited() {
        var catalog = StudyReferenceCatalog.build(List.of(
            query("field_background", Map.of("doi", "10.1000/a", "title", "Same title")),
            query("similar_studies", Map.of("doi", "10.1000/b", "title", "Same title")),
            query("similar_studies", Map.of("title", "<script>Unidentified</script>"))));
        assertEquals(3, StudyPublicEvidence.list(catalog.get("references")).size());
        var unidentified = StudyPublicEvidence.list(catalog.get("references")).stream().map(StudyPublicEvidence::map)
            .filter(r -> Boolean.FALSE.equals(r.get("citation_eligible"))).findFirst().orElseThrow();
        assertTrue(StudyReferenceCatalog.lookup(catalog, unidentified.get("reference_id").toString()).isEmpty());
        assertFalse(StudyReferenceCatalog.html(catalog).contains("<script>"));
        assertTrue(StudyReferenceCatalog.html(catalog).contains("&lt;script&gt;"));
    }
    @Test void completedEmptySearchStillNeedsScientificReviewAndMissingBlockedTruncatedRemainDistinct() {
        List<Object> searches = new ArrayList<>();
        for (String strategy : StudyIntroduction.SEARCH_STRATEGIES) searches.add(Map.of("search_strategy", strategy, "status", "no_matches"));
        var coverage = StudyPublicEvidence.map(StudyReferenceCatalog.build(searches).get("search_coverage"));
        assertEquals("bounded_coverage_complete_scientific_review_required", coverage.get("status"));
        assertTrue(coverage.get("sufficiency").toString().contains("Not established"));
        searches.set(0, Map.of("search_strategy", "field_background", "status", "candidate_matches", "search_limited", true));
        searches.set(1, Map.of("search_strategy", "hypothesis_relevance", "status", "offline_policy"));
        searches.remove(3);
        coverage = StudyPublicEvidence.map(StudyReferenceCatalog.build(searches).get("search_coverage"));
        assertEquals("additional_search_or_input_required", coverage.get("status"));
        for (String status : List.of("bounded_results_truncated", "incomplete_connection_or_policy", "not_searched")) assertTrue(coverage.toString().contains(status));
    }
    @Test void localCatalogRetainsAuditButAdviceWhitelistOmitsPathsAndFullAbstracts() {
        var article = Map.<String, Object>of("doi", "10.1000/a", "abstractText", "PRIVATE_FULL_ABSTRACT", "title", "Public title");
        var evidence = Map.of("H1.literature_search", Map.<String, Object>of("operation_type", "literature_search", "manifest", "LOCAL_PATH",
            "key_findings", Map.of("evidence", List.of(query("field_background", article)))));
        var catalog = StudyReferenceCatalog.fromEvidence(evidence);
        var ref = StudyPublicEvidence.map(StudyPublicEvidence.list(catalog.get("references")).get(0));
        assertTrue(ref.toString().contains("LOCAL_PATH"));
        assertTrue(ref.toString().contains("source-hash"));
        String advice = StudyReferenceCatalog.publicMetadata(ref).toString();
        assertFalse(advice.contains("LOCAL_PATH")); assertFalse(advice.contains("PRIVATE_FULL_ABSTRACT"));
        assertTrue(advice.contains("Public title"));
    }
    @Test void completeSearchForOnePhenotypeCannotHideMissingCoverageForAnother() {
        List<Object> queries = new ArrayList<>();
        for (String strategy : StudyIntroduction.SEARCH_STRATEGIES)
            queries.add(Map.of("search_strategy", strategy, "status", "no_matches", "phenotype", "sleep"));
        queries.add(Map.of("search_strategy", "field_background", "status", "no_matches", "phenotype", "BMI"));
        var coverage = StudyPublicEvidence.map(StudyReferenceCatalog.build(queries).get("search_coverage"));
        assertEquals("additional_search_or_input_required", coverage.get("status"));
        assertTrue(coverage.toString().contains("not_searched"));
    }
    @Test void editorialPayloadIncludesCitedPmcAndBoundedFieldCandidatesWithoutLocalAudit() throws Exception {
        List<Object> queries = new ArrayList<>();
        for (int i = 0; i < 25; i++) queries.add(query("field_background", Map.of("pmid", Integer.toString(100 + i), "title", "Field paper " + i)));
        queries.add(query("gene_followup", Map.of("pmcid", "PMC900", "title", "Cited gene paper")));
        var evidence = Map.of("H1.literature_search", Map.<String, Object>of("operation_type", "literature_search", "state", "complete", "manifest", "LOCAL_PRIVATE_PATH",
            "key_findings", Map.of("evidence", queries)));
        var revision = Map.<String, Object>of("title", "Fixture", "sections", List.of(Map.of("id", "field", "section", "introduction",
            "introduction_role", "field_background", "text", "Reviewed prose.", "references", List.of("pmcid:PMC900"))));
        var payload = StudyEditorialReview.payload(revision, evidence);
        assertEquals(21, StudyPublicEvidence.list(payload.get("introduction_reference_candidates")).size());
        assertTrue(payload.toString().contains("Cited gene paper"));
        assertTrue(payload.toString().contains("https://pmc.ncbi.nlm.nih.gov/articles/PMC900/"));
        assertTrue(payload.toString().contains("importance_and_gap"));
        assertFalse(payload.toString().contains("LOCAL_PRIVATE_PATH"));
        assertFalse(payload.toString().contains("source-hash"));
    }
    @Test void laterAvailableAbstractKeepsItsOwnExcerptStatusOffsetsAndAudit() {
        var missing = StudyLiterature.abstractEvidence(Map.of("doi", "10.1000/a"), "", "sleep", "queries/0001");
        var available = StudyLiterature.abstractEvidence(Map.of("doi", "10.1000/a", "abstractText", "Later source text about sleep."), "", "sleep", "queries/0002");
        var catalog = StudyReferenceCatalog.build(List.of(query("field_background", missing), query("similar_studies", available)));
        var ref = StudyPublicEvidence.map(StudyPublicEvidence.list(catalog.get("references")).get(0));
        assertEquals(available.get("excerpt"), ref.get("excerpt"));
        assertEquals(available.get("excerpt_status"), ref.get("excerpt_status"));
        assertEquals(available.get("excerpt_start"), ref.get("excerpt_start"));
        assertEquals("queries/0002", ref.get("source_audit"));
        assertTrue(ref.get("excerpt_retrieval_context").toString().contains("similar_studies"));
        assertEquals(2, StudyPublicEvidence.list(ref.get("retrieval_context")).size());
    }
}
