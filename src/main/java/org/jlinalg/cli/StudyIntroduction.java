/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Shared, source-grounded introduction policy for native drafts and LLM advice. */
final class StudyIntroduction {
    static final List<String> ROLES = List.of("field_background", "hypothesis_relevance",
        "importance_and_gap", "hypotheses_and_motivation");
    static final List<String> SEARCH_STRATEGIES = List.of("field_background", "hypothesis_relevance",
        "similar_studies", "field_conflicting_evidence");
    private StudyIntroduction() { }

    static String instructions() {
        return "Write the Introduction as four or five connected prose paragraphs, not a list or a set of mini-sections. "
            + "Start with one or two paragraphs summarizing the relevant field and established findings, grounded in relevant literature. "
            + "Next write one paragraph explaining why the study hypothesis is relevant to that field. "
            + "Next write one paragraph explaining why the question matters, what closely related studies have already done, "
            + "and specifically how this study differs (for example molecular measurement, population, tissue, design or validation); "
            + "cite the comparisons and do not equate a difference with proven novelty. "
            + "End with one paragraph stating the overall hypotheses and their scientific motivation, using the declared study hypotheses "
            + "rather than retrospectively inventing motivations from significant findings. "
            + "Use paragraph roles field_background, hypothesis_relevance, importance_and_gap, hypotheses_and_motivation in that order. "
            + "Search the field and similar studies, not only top-hit genes. Include contradictory findings. "
            + "Use only retrieved, cataloged references for factual and prior-work claims. A bounded, failed or empty search does not establish "
            + "that no similar research exists; do not assert first-ever novelty. Identify missing coverage as a research need. "
            + "Citation presence and paragraph-role checks do not verify that a source actually supports a claim. "
            + "Keep drafting instructions, source-review notes and software/audit details out of finished introduction prose; "
            + "put unresolved writing/research needs in the review checklist instead of inventing background.";
    }

    static Map<String, Object> policy() {
        return Map.of("version", 1, "paragraphs_min", 4, "paragraphs_max", 5,
            "ordered_roles", ROLES, "background_paragraphs", "one or two; every other role exactly one",
            "literature_search_strategies", SEARCH_STRATEGIES, "instructions", instructions());
    }

    static Map<String, Object> assess(List<?> sections) throws IOException {
        List<String> needs = new ArrayList<>(); Map<String, Integer> counts = new LinkedHashMap<>();
        ROLES.forEach(role -> counts.put(role, 0));
        int previous = -1, unclassified = 0;
        for (Object raw : sections) {
            Map<String, Object> p = map(raw);
            if (!p.containsKey("introduction_role")) {
                if ("introduction".equals(p.get("section"))) unclassified++;
                continue;
            }
            if (!"introduction".equals(p.get("section"))) throw new IOException("introduction_role is only valid for introduction paragraphs");
            String role = text(p.get("introduction_role")); int position = ROLES.indexOf(role);
            if (position < 0) throw new IOException("Unknown introduction_role: " + role);
            if (position < previous) throw new IOException("Introduction paragraphs must follow field_background, hypothesis_relevance, importance_and_gap, hypotheses_and_motivation order");
            previous = position;
            String prose = text(p.get("text")).strip();
            int paragraphs = prose.isEmpty() ? 0 : prose.split("\\R\\s*\\R").length;
            counts.merge(role, paragraphs, Integer::sum);
            if (counts.get(role) > (role.equals("field_background") ? 2 : 1)) throw new IOException("Too many introduction paragraphs for " + role);
            if (prose.isBlank() || prose.contains("[SCIENTIST INPUT REQUIRED")) needs.add(role + ": supply source-grounded prose");
            if (!role.equals("hypotheses_and_motivation") && StudyCitations.identifiers(p, "text").isEmpty())
                needs.add(role + ": cite cataloged literature supporting the field, rationale or prior-work comparison");
        }
        counts.forEach((role, count) -> { if (count == 0) needs.add(role + ": missing paragraph"); });
        if (unclassified > 0) needs.add("Assign introduction_role to " + unclassified + " existing introduction block(s); older prose is preserved, not silently classified or replaced");
        return Map.of("policy", policy(), "status", needs.isEmpty() ? "structure_complete_scientific_review_required" : "needs_input",
            "paragraph_counts", counts, "needs", needs, "scientific_review", List.of(
                "Assess literature coverage, relevance and source support, including contradictory studies and retrieval limits.",
                "Verify the actual comparisons and contribution; citation presence is not novelty verification.",
                "Confirm motivations and hypotheses against the study history; record changed hypotheses in an amendment."));
    }

    static Map<String, Object> template(String title, List<String> hypotheses) {
        List<Object> sections = new ArrayList<>();
        List<String> tasks = List.of("Summarize the relevant field in one or two paragraphs with supporting references.",
            "Explain the hypothesis's relevance to the field, citing the scientific rationale.",
            "Explain importance, compare the most similar studies, and state this study's specific contribution with references.",
            "Explain the scientific motivation for these hypotheses without inventing a prespecified rationale.");
        for (int i = 0; i < ROLES.size(); i++) {
            Map<String, Object> p = new LinkedHashMap<>(); String role = ROLES.get(i);
            p.put("id", "intro-" + role.replace('_', '-')); p.put("section", "introduction"); p.put("introduction_role", role);
            String statement = i == 3 && !hypotheses.isEmpty() ? "The study hypotheses are: " + String.join(" ", hypotheses) + " " : "";
            p.put("text", statement + "[SCIENTIST INPUT REQUIRED: " + tasks.get(i) + "]");
            p.put("rationale", "Develop the " + role.replace('_', ' ') + " paragraph from reviewed study and literature evidence.");
            p.put("sources", List.of()); p.put("references", List.of()); sections.add(p);
        }
        return Map.of("schema_version", 1, "id", "introduction-revision", "title", title,
            "rationale", "Develop a field-to-hypothesis introduction using reviewed literature and the declared study hypotheses.", "sections", sections);
    }

    static String scaffold(List<?> sections) {
        StringBuilder html = new StringBuilder();
        for (Object raw : sections) html.append("<p>").append(escape(text(map(raw).get("text")))).append("</p>");
        return html.toString();
    }
    private static String escape(String s) { return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;"); }
}
