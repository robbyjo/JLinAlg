/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** Deduplicated, query-linked references. Retrieval never implies scientific endorsement. */
final class StudyReferenceCatalog {
    private StudyReferenceCatalog() { }
    static Map<String, Object> fromEvidence(Map<String, Map<String, Object>> evidence) {
        return fromEvidence(evidence, List.of());
    }
    static Map<String, Object> fromEvidence(Map<String, Map<String, Object>> evidence, List<?> importedQueries) {
        List<Object> queries = new ArrayList<>();
        evidence.forEach((id, e) -> {
            if (!"literature_search".equals(e.get("operation_type"))) return;
            for (Object raw : list(map(e.get("key_findings")).get("evidence"))) {
                Map<String, Object> q = new LinkedHashMap<>(map(raw));
                q.put("source_operation", id); q.put("source_manifest", e.get("manifest")); queries.add(q);
            }
        });
        queries.addAll(importedQueries);
        return build(queries);
    }
    static Map<String, Object> build(List<?> queries) {
        List<Map<String, Object>> records = new ArrayList<>();
        Map<String, Map<String, Object>> aliases = new LinkedHashMap<>();
        List<Object> searches = new ArrayList<>();
        for (Object raw : queries) {
            Map<String, Object> query = map(raw);
            if (query.containsKey("search_strategy")) searches.add(select(query, List.of("source_operation", "source_manifest", "query", "search_strategy",
                "phenotype", "gene_symbol", "cohort_id", "cohort_name", "methods_role", "requested_reference", "page", "status", "hit_count", "search_limited", "audit", "retrieval", "scientific_rationale")));
            for (Object item : list(query.get("articles"))) {
                Map<String, Object> article = map(item); Set<String> ids = identifiers(article);
                if (ids.isEmpty()) ids.add("unresolved:" + digest(StudyAssistance.json(article)));
                List<Map<String, Object>> matches = new ArrayList<>();
                for (String id : ids) if (aliases.containsKey(id) && matches.stream().noneMatch(m -> m == aliases.get(id))) matches.add(aliases.get(id));
                Map<String, Object> record;
                if (matches.isEmpty()) { record = new LinkedHashMap<>(); records.add(record); }
                else {
                    record = matches.get(0);
                    for (int i = 1; i < matches.size(); i++) {
                        Map<String, Object> other = matches.get(i); merge(record, other);
                        records.removeIf(r -> r == other);
                        aliases.replaceAll((key, value) -> value == other ? record : value);
                    }
                }
                Map<String, Object> occurrence = select(query, List.of("source_operation", "source_manifest", "search_strategy", "query", "scientific_rationale", "phenotype", "gene_symbol", "feature_id", "cohort_id", "cohort_name", "methods_role", "requested_reference", "page", "audit", "retrieval"));
                article.put("identifiers", new ArrayList<>(ids)); article.put("retrieval_context", List.of(occurrence));
                merge(record, article);
                for (Object id : list(record.get("identifiers"))) aliases.put(id.toString(), record);
            }
        }
        for (Map<String, Object> r : records) {
            List<String> ids = list(r.get("identifiers")).stream().map(Object::toString).sorted().toList();
            String identity = ids.stream().filter(s -> s.startsWith("doi:")).findFirst().orElse(ids.get(0));
            r.put("reference_id", "ref-" + digest(identity).substring(0, 20));
            r.put("citation_eligible", ids.stream().anyMatch(s -> !s.startsWith("unresolved:")));
            r.put("review_status", "unreviewed_candidate_not_hypothesis_support");
            String doi = ids.stream().filter(s -> s.startsWith("doi:")).map(s -> s.substring(4)).findFirst().orElse("");
            if (!doi.isEmpty()) { r.put("doi", doi); r.put("source_url", "https://doi.org/" + doi); }
            else ids.stream().filter(s -> s.startsWith("pmid:")).findFirst().ifPresent(id -> r.put("source_url", "https://pubmed.ncbi.nlm.nih.gov/" + id.substring(5) + "/"));
            if (!r.containsKey("source_url")) ids.stream().filter(s -> s.startsWith("pmcid:")).findFirst()
                .ifPresent(id -> r.put("source_url", "https://pmc.ncbi.nlm.nih.gov/articles/" + id.substring(6) + "/"));
            Map<String, Object> journalInfo = map(r.get("journalInfo"));
            r.put("journal", map(journalInfo.get("journal")).getOrDefault("title", "not reported"));
        }
        records.sort(java.util.Comparator.comparing(r -> r.get("reference_id").toString()));
        Map<String, Object> coverage = coverage(searches.stream().filter(q -> !map(q).containsKey("methods_role")).toList());
        return Map.of("schema_version", 1, "references", records, "searches", searches, "search_coverage", coverage,
            "methods_searches", searches.stream().filter(q -> map(q).containsKey("methods_role")).toList(),
            "interpretation", "Reference catalog, not a sufficient or systematic review certification. Deduplicated by normalized DOI, PMID, PMCID or provider ID; no title-only fuzzy merges. Review full available sources, metadata conflicts, similar and contradictory studies before claims of contribution or novelty.");
    }
    private static void merge(Map<String, Object> target, Map<String, Object> source) {
        // Excerpt, offsets, status and audit must come from the same occurrence.
        // An earlier record may have no abstract, or a duplicate may lack its quote.
        if (!target.containsKey("excerpt") && source.containsKey("excerpt")) {
            target.keySet().removeIf(key -> key.startsWith("excerpt"));
            source.forEach((key, value) -> { if (key.startsWith("excerpt") || key.equals("source_audit")) target.put(key, value); });
            target.put("excerpt_retrieval_context", source.getOrDefault("excerpt_retrieval_context", source.getOrDefault("retrieval_context", List.of())));
        }
        for (var e : source.entrySet()) {
            String key = e.getKey(); Object value = e.getValue();
            if (target.containsKey("excerpt") && (key.startsWith("excerpt") || key.equals("source_audit"))) continue;
            if (List.of("identifiers", "retrieval_context", "metadata_conflicts").contains(key)) {
                List<Object> merged = new ArrayList<>(list(target.get(key))); for (Object v : list(value)) if (!merged.contains(v)) merged.add(v); target.put(key, merged);
            } else if (!target.containsKey(key) || text(target.get(key)).isBlank()) target.put(key, value);
            else if (List.of("doi", "pmid", "pmcid", "title", "authorString", "pubYear", "journalInfo").contains(key)
                    && !text(value).isBlank() && !normalize(text(target.get(key))).equals(normalize(text(value)))) {
                List<Object> conflicts = new ArrayList<>(list(target.get("metadata_conflicts")));
                conflicts.add(Map.of("field", key, "retained", target.get(key), "alternative", value)); target.put("metadata_conflicts", conflicts);
            }
        }
    }
    private static Map<String, Object> coverage(List<Object> searches) {
        Map<String, List<Object>> scopes = new LinkedHashMap<>();
        for (Object raw : searches) {
            Map<String, Object> q = map(raw);
            String scope = text(q.get("source_operation")) + " / " + text(q.get("phenotype"));
            scopes.computeIfAbsent(scope, ignored -> new ArrayList<>()).add(raw);
        }
        Map<String, Object> result = new LinkedHashMap<>(strategyCoverage(searches));
        List<Object> byScope = new ArrayList<>(); boolean incomplete = false;
        for (var entry : scopes.entrySet()) {
            Map<String, Object> item = new LinkedHashMap<>(strategyCoverage(entry.getValue())); item.put("scope", entry.getKey());
            byScope.add(item); incomplete |= "additional_search_or_input_required".equals(item.get("status"));
        }
        if (incomplete) result.put("status", "additional_search_or_input_required");
        result.put("by_operation_and_phenotype", byScope);
        return result;
    }
    private static Map<String, Object> strategyCoverage(List<Object> searches) {
        List<Object> strategies = new ArrayList<>(); boolean incomplete = false;
        for (String strategy : StudyIntroduction.SEARCH_STRATEGIES) {
            List<Map<String, Object>> rows = searches.stream().map(StudyPublicEvidence::map).filter(q -> strategy.equals(q.get("search_strategy"))).toList();
            boolean blocked = rows.stream().anyMatch(q -> !Set.of("candidate_matches", "no_matches").contains(q.get("status")));
            boolean limited = rows.stream().anyMatch(q -> Boolean.TRUE.equals(q.get("search_limited")));
            String status = rows.isEmpty() ? "not_searched" : blocked ? "incomplete_connection_or_policy" : limited ? "bounded_results_truncated" : "bounded_search_completed";
            strategies.add(Map.of("strategy", strategy, "status", status, "pages", rows.size()));
            incomplete |= rows.isEmpty() || blocked || limited;
        }
        return Map.of("status", incomplete ? "additional_search_or_input_required" : "bounded_coverage_complete_scientific_review_required",
            "strategies", strategies, "sufficiency", "Not established by article counts or query completion. Scientist/advisor must assess topic and synonym coverage, source quality, relevance, contrary findings and database/date limits.");
    }
    static Set<String> identifiers(Map<String, Object> article) {
        Set<String> result = new LinkedHashSet<>(); String doi = normalize(text(article.get("doi")));
        if (doi.startsWith("10.") && doi.contains("/") && doi.chars().noneMatch(Character::isWhitespace)) result.add("doi:" + doi);
        String pmid = text(article.get("pmid"));
        if (pmid.isBlank() && "MED".equals(article.get("source"))) pmid = text(article.get("id"));
        if (pmid.matches("[1-9][0-9]*")) result.add("pmid:" + pmid);
        String pmcid = text(article.get("pmcid")).toUpperCase(Locale.ROOT);
        if (pmcid.matches("PMC[0-9]+")) result.add("pmcid:" + pmcid);
        if (!text(article.get("source")).isBlank() && !text(article.get("id")).isBlank()) result.add("provider:" + article.get("source") + ":" + article.get("id"));
        return result;
    }
    static String normalize(String id) {
        String value = id.strip().replaceFirst("(?i)^https?://(?:dx\\.)?doi\\.org/", "").replaceFirst("(?i)^doi:", "");
        if (value.matches("(?i)pmcid:PMC[0-9]+")) return "pmcid:" + value.substring(6).toUpperCase(Locale.ROOT);
        if (value.matches("(?i)pmid:[1-9][0-9]*")) return value.toLowerCase(Locale.ROOT);
        return value.startsWith("10.") ? value.toLowerCase(Locale.ROOT) : value;
    }
    static Map<String, Object> publicMetadata(Map<String, Object> reference) {
        var result = select(reference, List.of("reference_id", "identifiers", "doi", "pmid", "pmcid", "title", "authorString", "pubYear", "journal",
            "source_url", "excerpt", "excerpt_scope", "review_status", "citation_eligible", "metadata_conflicts", "commentCorrectionList"));
        result.put("citation_purposes", list(reference.get("retrieval_context")).stream().map(StudyPublicEvidence::map)
            .filter(q -> q.containsKey("methods_role")).map(q -> select(q, List.of("methods_role", "cohort_name"))).distinct().toList());
        return result;
    }
    static Map<String, Object> lookup(Map<String, Object> catalog, String id) {
        String normalized = normalize(id);
        for (Object raw : list(catalog.get("references"))) {
            Map<String, Object> r = map(raw);
            if (!Boolean.TRUE.equals(r.get("citation_eligible"))) continue;
            if (id.equals(r.get("reference_id")) || list(r.get("identifiers")).contains(normalized) || list(r.get("identifiers")).contains("doi:" + normalized)) return r;
        }
        return Map.of();
    }
    static String html(Map<String, Object> catalog) {
        StringBuilder h = new StringBuilder("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Reference catalog</title><style>body{font:16px/1.6 system-ui;max-width:1050px;margin:2rem auto;padding:1rem}article{border-top:1px solid #ccc}dt{font-weight:bold}a{overflow-wrap:anywhere}</style><h1>Reference catalog</h1><p>");
        h.append(esc(text(catalog.get("interpretation")))).append("</p><h2>Search coverage</h2>").append(StudyPublication.render(catalog.get("search_coverage")));
        h.append("<h2>Cohort and method reference searches</h2>").append(StudyPublication.render(catalog.get("methods_searches")));
        for (Object raw : list(catalog.get("references"))) {
            Map<String, Object> r = map(raw);
            h.append("<article id=\"").append(esc(text(r.get("reference_id")))).append("\"><h2>").append(esc(text(r.getOrDefault("title", "Title unavailable")))).append("</h2><p>")
                .append(esc(text(r.get("authorString")))).append(". ").append(esc(text(r.get("pubYear")))).append(". ").append(esc(text(r.get("journal")))).append(".</p>");
            if (r.containsKey("source_url")) h.append("<p><a href=\"").append(esc(text(r.get("source_url")))).append("\">Original source</a></p>");
            h.append("<p>").append(esc(text(r.get("reference_id")))).append(" — unreviewed candidate</p><details><summary>Identifiers, source evidence and retrieval history</summary>").append(StudyPublication.render(r)).append("</details></article>");
        }
        return h.append("</html>").toString();
    }
    private static String digest(String value) {
        try { return java.util.HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(value.getBytes(StandardCharsets.UTF_8))); }
        catch (java.security.NoSuchAlgorithmException e) { throw new IllegalStateException(e); }
    }
    private static Map<String, Object> select(Map<String, Object> source, List<String> fields) { Map<String, Object> r = new LinkedHashMap<>(); for (String f : fields) if (source.containsKey(f)) r.put(f, source.get(f)); return r; }
    private static String esc(String s) { return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;"); }
}
