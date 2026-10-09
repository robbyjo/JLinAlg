/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/** Claim-adjacent citations, numbered on first appearance, with an in-paper bibliography. */
final class StudyCitations {
    private static final Pattern TOKEN = Pattern.compile("\\[\\[cite:([^\\[\\]\\r\\n]+)]]");
    private static final Set<String> BASES = Set.of("literature", "study_results", "study_design", "interpretation", "not_applicable");
    private final Map<String, Object> catalog;
    private final Map<String, Map<String, Object>> used = new LinkedHashMap<>();
    StudyCitations(Map<String, Object> catalog) { this.catalog = catalog; }

    static String instructions() {
        return "Cite every journal-derived claim wherever it occurs, including Methods, Results context and Discussion. "
            + "Place [[cite:DOI;pmid:ID]] immediately beside the supported claim; identifiers must resolve to retrieved catalog records. "
            + "Declare citation_basis: literature for prose relying on publications; study_results, study_design, interpretation or not_applicable distinguish original material. "
            + "Use citation_policy: strict to reject missing classifications, missing inline citations for literature prose and unplaced references. "
            + "Original results and declared study definitions do not require invented journal citations. Number references by first appearance, reuse numbers, "
            + "and link each superscript bracket to the final References section. This follows Nature's numerical ordering with user-requested bracketed hyperlinks. "
            + "A mechanically valid citation does not prove that it supports every sentence; source support, completeness and citation_basis require scientific review.";
    }
    static List<String> identifiers(Map<String, Object> paragraph, String field) {
        Set<String> ids = new LinkedHashSet<>();
        var matcher = TOKEN.matcher(text(paragraph.get(field)));
        while (matcher.find()) for (String id : matcher.group(1).split(";", -1)) ids.add(id.strip());
        for (Object id : list(paragraph.get("references"))) ids.add(text(id));
        return new ArrayList<>(ids);
    }
    static Map<String, Object> assess(Map<String, Object> revision, Map<String, Object> catalog) throws IOException {
        String policy = text(revision.getOrDefault("citation_policy", "review"));
        if (!Set.of("strict", "review").contains(policy)) throw new IOException("citation_policy must be strict or review");
        List<String> needs = new ArrayList<>();
        for (Object raw : list(revision.get("sections"))) {
            var p = map(raw); check(p, "text", text(p.get("id")), catalog, needs);
            int i = 0;
            for (Object s : list(p.get("supporting_evidence"))) check(map(s), "relationship", p.get("id") + ".support-" + (++i), catalog, needs);
        }
        if (policy.equals("strict") && !needs.isEmpty()) throw new IOException("Citation enforcement: " + String.join("; ", needs));
        return Map.of("policy", policy, "instructions", instructions(), "needs", needs,
            "status", needs.isEmpty() ? "citation_structure_complete_scientific_review_required" : "needs_input",
            "limitation", "Declared literature dependence and citation syntax are checked, not the truth or completeness of source support. Unclassified older prose is preserved only in review mode.");
    }
    private static void check(Map<String, Object> p, String field, String id, Map<String, Object> catalog, List<String> needs) throws IOException {
        String basis = text(p.get("citation_basis")), prose = text(p.get(field));
        if (basis.isBlank()) needs.add(id + ": classify citation_basis, including journal-derived statements in any section");
        else if (!BASES.contains(basis)) throw new IOException("Unknown citation_basis: " + basis);
        if (p.containsKey("references") && !(p.get("references") instanceof List<?>)) throw new IOException("references must be a list");
        var inline = new LinkedHashSet<String>(); var matcher = TOKEN.matcher(prose);
        while (matcher.find()) for (String key : matcher.group(1).split(";", -1)) inline.add(resolve(catalog, key.strip()).get("reference_id").toString());
        if (TOKEN.matcher(prose).replaceAll("").contains("[[cite:")) throw new IOException("Malformed inline citation in " + id);
        for (String key : identifiers(p, field)) {
            String ref = resolve(catalog, key).get("reference_id").toString();
            if (!inline.contains(ref)) needs.add(id + ": place reference " + key + " beside the supported claim with [[cite:...]]");
        }
        if (basis.equals("literature") && inline.isEmpty()) needs.add(id + ": journal-derived prose requires an inline citation");
    }
    private static Map<String, Object> resolve(Map<String, Object> catalog, String id) throws IOException {
        var r = StudyReferenceCatalog.lookup(catalog, id);
        if (id.isBlank() || r.isEmpty()) throw new IOException("Reference is not present in current retrieved evidence: " + id);
        return r;
    }
    String render(String prose, List<?> fallback) throws IOException {
        var matcher = TOKEN.matcher(prose); StringBuilder h = new StringBuilder(); int offset = 0;
        Set<String> inline = new LinkedHashSet<>();
        while (matcher.find()) {
            h.append(esc(prose.substring(offset, matcher.start())));
            List<String> keys = List.of(matcher.group(1).split(";")); h.append(cite(keys));
            for (String key : keys) inline.add(resolve(catalog, key.strip()).get("reference_id").toString());
            offset = matcher.end();
        }
        h.append(esc(prose.substring(offset)));
        // Review-mode compatibility: keep old paragraph references visible, without pretending they are claim-local.
        List<String> rest = new ArrayList<>();
        for (Object key : fallback) if (!inline.contains(resolve(catalog, key.toString()).get("reference_id"))) rest.add(key.toString());
        if (!rest.isEmpty()) h.append(cite(rest));
        return h.toString();
    }
    private String cite(List<String> keys) throws IOException {
        StringBuilder h = new StringBuilder("<sup class=\"citation\">"); Set<String> seen = new LinkedHashSet<>();
        for (String key : keys) {
            var r = resolve(catalog, key.strip()); String id = r.get("reference_id").toString();
            if (!seen.add(id)) continue;
            used.putIfAbsent(id, r); int number = new ArrayList<>(used.keySet()).indexOf(id) + 1;
            h.append("<a href=\"#reference-").append(number).append("\" aria-label=\"Reference ").append(number).append("\">[").append(number).append("]</a>");
        }
        return h.append("</sup>").toString();
    }
    String bibliography() {
        StringBuilder h = new StringBuilder("<section id=\"references\"><h2>References</h2>");
        if (used.isEmpty()) return h.append("<p>No journal references have been assigned; see the citation review before publication.</p></section>").toString();
        h.append("<ol class=\"references\">"); int n = 0;
        for (var r : used.values()) {
            var info = map(r.get("journalInfo")); var journal = map(info.get("journal"));
            String authors = authors(r), volume = text(info.getOrDefault("volume", r.get("volume"))), pages = text(r.getOrDefault("pageInfo", r.get("pages")));
            String name = text(journal.getOrDefault("medlineAbbreviation", journal.getOrDefault("title", r.get("journal"))));
            h.append("<li id=\"reference-").append(++n).append("\">");
            if (!authors.isBlank()) h.append(esc(authors)).append(" ");
            h.append(esc(text(r.get("title")))).append(" ").append("<i>").append(esc(name)).append("</i>");
            if (!volume.isBlank()) h.append(" <b>").append(esc(volume)).append("</b>");
            if (!pages.isBlank()) h.append(", ").append(esc(pages));
            h.append(" (").append(esc(text(r.get("pubYear")))).append(").");
            String url = text(r.get("source_url"));
            if (url.startsWith("https://")) h.append(" <a href=\"").append(esc(url)).append("\">").append(esc(r.containsKey("doi") ? "doi:" + r.get("doi") : "Original source")).append("</a>");
            h.append(" <a class=\"source\" href=\"reference-catalog.html#").append(esc(text(r.get("reference_id")))).append("\">Retrieval record</a></li>");
        }
        return h.append("</ol></section>").toString();
    }
    private static String authors(Map<String, Object> r) {
        List<String> names = new ArrayList<>();
        for (Object raw : list(map(r.get("authorList")).get("author"))) {
            var a = map(raw); String family = text(a.get("lastName")), initials = text(a.get("initials"));
            if (family.isBlank()) names.add(text(a.get("fullName")));
            else names.add(family + (initials.isBlank() ? "" : ", " + initials.replaceAll("([A-Z])", "$1.").replace("..", ".")));
        }
        names.removeIf(String::isBlank);
        if (names.size() > 5) return names.get(0) + " et al.";
        if (!names.isEmpty()) return String.join(", ", names) + ".";
        return text(r.get("authorString"));
    }
    private static String esc(String s) { return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;"); }
}
