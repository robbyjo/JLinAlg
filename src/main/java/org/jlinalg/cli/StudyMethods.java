/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** Cohort-first Methods model. Reporting declarations never change analysis settings. */
final class StudyMethods {
    static final List<String> ROLES = List.of("cohorts", "association", "discovery_replication", "meta_analysis",
        "enrichment_ontology", "mendelian_randomization", "network", "fine_mapping", "colocalization", "sensitivity_analysis", "reproducibility");
    private static final Map<String, String> PUBLIC_COHORTS = Map.ofEntries(
        Map.entry("fhs", "Framingham Heart Study"), Map.entry("fhs-offspring", "Framingham Offspring Study"),
        Map.entry("fhs-gen3", "Framingham Heart Study Third Generation"),
        Map.entry("cardia", "Coronary Artery Risk Development in Young Adults"),
        Map.entry("whi", "Women's Health Initiative"), Map.entry("jhs", "Jackson Heart Study"));
    private StudyMethods() { }

    static String instructions() {
        return "Always start Methods with a description of every participating cohort, before any statistical or software description. "
            + "For established cohorts, retrieve and cite the appropriate design/recruitment paper for each cohort and subcohort; "
            + "for example, the Framingham Heart Study (FHS), Coronary Artery Risk Development in Young Adults (CARDIA), "
            + "Women's Health Initiative (WHI), and Jackson Heart Study (JHS). Expand names at first use. "
            + "A cohort paper describes the historical cohort, not this study's analysis N, visit, exclusions, consent or phenotype definition. "
            + "Then describe association, discovery/replication if performed, meta-analysis, enrichment/ontology, Mendelian randomization (MR), "
            + "network analysis, fine mapping if performed, and colocalization, followed by sensitivity analyses and reproducibility. "
            + "Describe only analyses actually performed; explicitly distinguish supplied external association results, planned work and blocked analyses. "
            + "For each performed method give its exact estimator/algorithm, an appropriate retrieved method citation, software/version, inputs, "
            + "all relevant effective parameters and defaults used, including exact phenotype definitions and comparison operators (for example obesity defined as BMI > 30 kg/m2). "
            + "State customary values such as FDR < 0.05 without requiring a special justification. Require recorded scientific justifications only for uncommon, study-specific or deviating parameter choices; never invent a rationale. "
            + "Cite both statistical methods and annotation/database versions where relevant. "
            + "Association needs model family/link and formula orientation, biological covariates, technical factors, transformations, "
            + "sample selection/missingness, pedigree/relatedness adjustment and multiple-testing rule. Report any executed effect-outlier QC, "
            + "its effective thresholds, robust-scale fallback or skipped screening, per-cohort exclusions and resulting testing universe; "
            + "identify post hoc QC changes and require sensitivity review rather than treating excluded large effects as proven artifacts. Meta-analysis needs effect scale, "
            + "weighting, fixed/random-effects estimator, heterogeneity and cohort-overlap handling. Enrichment needs selection threshold, "
            + "tested-gene universe, feature-to-gene mapping, database/release, test and multiplicity adjustment. MR needs molecular exposure/outcome "
            + "and tissue, GWAS/QTL sources, instrument P/F thresholds, LD clumping reference/r2/window, allele harmonization, overlap assumptions, estimator and sensitivity settings. "
            + "Network analysis needs network source/version, nodes/edges, weights, filtering, algorithm and null model. Colocalization needs "
            + "region/build, complete summary statistics, sample sizes, LD, signal model, priors, posterior threshold and prior sensitivity. "
            + "Do not infer settings from current software defaults for historical runs, substitute unrelated citations, invent a parameter rationale, "
            + "or present a query match as verified support. Missing source, parameter or required uncommon-parameter justification remains a review need. "
            + "Report the dataset tissue/sample compartment for each cohort and consider declared differences when interpreting heterogeneity; do not infer causation from a cohort-confounded tissue contrast. "
            + "Keep audit/editorial commentary in supporting records rather than finished Methods prose.";
    }
    static Map<String, Object> policy() { return Map.of("version", 1, "ordered_roles", ROLES, "instructions", instructions()); }
    static Map<String, Object> record(String method, Map<String, Object> parameters, List<?> references, Object justification, List<String> literatureTerms) {
        return Map.of("method", method, "parameters", parameters, "references", references,
            "parameter_justification", justification == null || text(justification).isBlank() ? "Not recorded; customary values need no special justification. Review uncommon or deviating choices separately." : justification,
            "literature_terms", literatureTerms, "software", "JLinAlg " + JLinAlgCli.version(),
            "provenance", "Recorded by the executing adapter for this run; inspect command and source artifacts for backend settings.");
    }
    static String role(String type) {
        return switch (type) {
            case "validate_harmonize" -> "association";
            case "gene_ontology", "pathway_enrichment" -> "enrichment_ontology";
            case "network_analysis", "network_followup" -> "network";
            case "leave_one_cohort_out", "fixed_random_comparison" -> "sensitivity_analysis";
            default -> ROLES.contains(type) ? type : "";
        };
    }
    static List<Map<String, Object>> cohorts(Map<String, Object> study) {
        Map<String, Map<String, Object>> declared = new LinkedHashMap<>();
        for (Object raw : list(study.get("cohorts"))) { var c = map(raw); declared.put(text(c.get("id")), c); }
        Set<String> used = new LinkedHashSet<>();
        for (Object raw : list(study.get("datasets"))) {
            var d = map(raw); if ("retired".equals(d.get("status"))) continue;
            if (!text(d.get("cohort")).isBlank()) used.add(text(d.get("cohort")));
            for (Object id : list(d.get("cohorts"))) used.add(text(id));
        }
        List<Map<String, Object>> result = new ArrayList<>();
        for (String id : used) {
            var c = new LinkedHashMap<>(declared.getOrDefault(id, Map.of("id", id)));
            if (c.containsKey("name")) c.putIfAbsent("literature_name", c.get("name"));
            c.putIfAbsent("name", PUBLIC_COHORTS.getOrDefault(id.toLowerCase(Locale.ROOT), id));
            List<Object> tissues = new ArrayList<>();
            for (Object raw : list(study.get("datasets"))) {
                var d = map(raw);
                if (!"retired".equals(d.get("status")) && (id.equals(d.get("cohort")) || list(d.get("cohorts")).contains(id)) && !text(d.get("tissue")).isBlank())
                    tissues.add(Map.of("dataset", d.get("id"), "tissue", d.get("tissue")));
            }
            c.put("dataset_tissues", tissues); result.add(c);
        }
        return result;
    }
    static String literatureName(Map<String, Object> cohort) {
        if (cohort.containsKey("literature_name")) return text(cohort.get("literature_name"));
        String explicit = text(cohort.getOrDefault("literature_name", cohort.get("name")));
        if (!explicit.isBlank() && !explicit.equals(cohort.get("id"))) return explicit;
        return PUBLIC_COHORTS.getOrDefault(text(cohort.get("id")).toLowerCase(Locale.ROOT), "");
    }
    static Map<String, Object> inventory(Map<String, Object> study, Map<String, Map<String, Object>> evidence) {
        List<Object> analyses = new ArrayList<>();
        for (Object raw : list(study.get("datasets"))) {
            var d = map(raw); if ("retired".equals(d.get("status"))) continue;
            boolean supplied = Set.of("association_results", "annotated_results").contains(text(d.get("stage")));
            if (supplied || !map(d.get("analysis_description")).isEmpty()) {
                var a = map(d.get("analysis_description"));
                analyses.add(Map.of("role", "association", "source", "dataset:" + d.get("id"), "execution_status", supplied ? "external_results_supplied_not_refitted" : "upstream_execution_unverified",
                    "method_record", a, "parameter_provenance", "Scientist-declared upstream analysis_description, not inferred from this JLinAlg run."));
            }
        }
        evidence.forEach((id, e) -> {
            String role = role(text(e.get("operation_type"))); if (role.isEmpty() || role.equals("cohorts")) return;
            Map<String, Object> a = new LinkedHashMap<>(); a.put("role", role); a.put("source", id);
            var method = map(map(e.get("details")).get("method_report"));
            a.put("execution_status", Boolean.FALSE.equals(method.get("performed")) ? "not_performed" : e.getOrDefault("state", "unknown")); a.put("manifest", e.getOrDefault("manifest", "unknown"));
            a.put("method_record", method);
            a.put("parameter_provenance", "Effective settings only when recorded by the executed adapter. Missing historical settings must be recovered from original logs, not current defaults.");
            analyses.add(a);
        });
        return Map.of("schema_version", 1, "cohorts", cohorts(study), "analyses", analyses,
            "study_contract_context", StudyContextReview.publicContext(StudyContextReview.inventory(study, null)),
            "policy", policy(), "interpretation", "Local reporting inventory. Completion alone does not prove a test was performed; inspect method_record, pair-level states and source artifacts. No patient rows or model refits are used.");
    }
    static Map<String, Object> template(String title, Map<String, Object> inventory, Map<String, Object> catalog) {
        List<Object> sections = new ArrayList<>(); int sequence = 0;
        for (Object raw : list(inventory.get("cohorts"))) {
            var c = map(raw); List<Object> refs = retrieved(list(c.get("references")), catalog);
            String prose = text(c.get("name")) + ": " + text(c.getOrDefault("description", "[SCIENTIST INPUT REQUIRED: describe cohort design, recruitment, population and participating subcohort/visit.]"));
            prose += refs.isEmpty() ? " [SCIENTIST INPUT REQUIRED: select and cite the appropriate cohort design paper.]" : " The cohort design has been described previously.";
            var p = paragraph("methods-cohort-" + (++sequence), "cohorts", prose, List.of(), refs, Map.of());
            p.put("title", c.get("name"));
            p.put("cohort_ids", List.of(c.get("id"))); sections.add(p);
        }
        if (sections.isEmpty()) sections.add(paragraph("methods-cohorts", "cohorts", "[SCIENTIST INPUT REQUIRED: identify and describe the participating cohorts, with design citations.]", List.of(), List.of(), Map.of()));
        for (String role : ROLES.subList(1, ROLES.size())) for (Object raw : list(inventory.get("analyses"))) {
            var a = map(raw); if (!role.equals(a.get("role"))) continue;
            String status = text(a.get("execution_status"));
            if (!Set.of("complete", "external_results_supplied_not_refitted").contains(status)) continue;
            var record = map(a.get("method_record"));
            String source = text(a.get("source"));
            String prose = "[SCIENTIST INPUT REQUIRED: describe " + role.replace('_', ' ') + " for " + source
                + ", citing the exact method and reporting effective parameters and recorded justifications. "
                + (status.startsWith("external") ? "These association results were supplied, not refitted by this reporting run." : "Verify that the recorded operation actually performed the test, including pair-level states.") + "]";
            var paragraph = paragraph("methods-" + role.replace('_', '-') + "-" + (++sequence), role, prose,
                source.startsWith("dataset:") ? List.of() : List.of(source), retrieved(list(record.get("references")), catalog), record);
            if (record.containsKey("parameter_justification")) paragraph.put("parameter_justifications", Map.of("recorded", record.get("parameter_justification")));
            else if (record.get("parameter_justifications") instanceof Map<?, ?>) paragraph.put("parameter_justifications", record.get("parameter_justifications"));
            sections.add(paragraph);
        }
        return Map.of("schema_version", 1, "id", "methods-revision", "title", title,
            "rationale", "Develop cohort-first Methods from source-linked execution records and reviewed citations.", "sections", sections);
    }
    static Map<String, Object> displayInventory(Map<String, Object> inventory) {
        var result = new LinkedHashMap<>(inventory); result.remove("study_contract_context"); return result;
    }
    private static Map<String, Object> paragraph(String id, String role, String prose, List<?> sources, List<?> refs, Map<String, Object> parameters) {
        var p = new LinkedHashMap<String, Object>(); p.put("id", id); p.put("section", "methods"); p.put("methods_role", role);
        p.put("title", role.replace('_', ' '));
        p.put("text", prose); p.put("sources", sources); p.put("references", refs); p.put("parameters", parameters);
        p.put("parameter_justifications", Map.of()); p.put("rationale", "Source-grounded Methods reporting; no analysis settings changed."); return p;
    }
    private static List<Object> retrieved(List<?> ids, Map<String, Object> catalog) {
        List<Object> result = new ArrayList<>(); for (Object id : ids) if (!StudyReferenceCatalog.lookup(catalog, id.toString()).isEmpty()) result.add(id); return result;
    }
    static Map<String, Object> assess(List<?> sections, Map<String, Object> inventory) throws IOException {
        List<String> needs = new ArrayList<>(); Set<String> described = new LinkedHashSet<>(), roles = new LinkedHashSet<>();
        int previous = -1, unclassified = 0;
        Set<String> knownCohorts = new LinkedHashSet<>(); for (Object raw : list(inventory.get("cohorts"))) knownCohorts.add(text(map(raw).get("id")));
        for (Object raw : sections) {
            var p = map(raw);
            if (!p.containsKey("methods_role")) { if ("methods".equals(p.get("section"))) unclassified++; continue; }
            if (!"methods".equals(p.get("section"))) throw new IOException("methods_role is only valid for Methods paragraphs");
            String role = text(p.get("methods_role")); int order = ROLES.indexOf(role);
            if (order < 0 || order < previous) throw new IOException("Methods must use the cohort-first ordered methods_role model: " + ROLES);
            previous = order; roles.add(role);
            if (text(p.get("text")).isBlank() || text(p.get("text")).contains("[SCIENTIST INPUT REQUIRED")) needs.add(role + ": supply reviewed Methods prose");
            if (!role.equals("reproducibility") && StudyCitations.identifiers(p, "text").isEmpty()) needs.add(role + ": supply supporting cataloged citations");
            for (String field : List.of("parameters", "parameter_justifications")) if (p.containsKey(field) && !(p.get(field) instanceof Map<?, ?>)) throw new IOException(field + " must be a mapping");
            if (role.equals("cohorts")) {
                if (p.containsKey("cohort_ids") && !(p.get("cohort_ids") instanceof List<?>)) throw new IOException("cohort_ids must be a list");
                if (!p.containsKey("cohort_ids")) needs.add("cohorts: explicitly identify the cohort_ids described in each paragraph");
                for (Object id : list(p.get("cohort_ids"))) {
                    if (!knownCohorts.contains(id.toString())) throw new IOException("Unknown Methods cohort: " + id);
                    described.add(id.toString());
                }
                if (list(p.get("cohort_ids")).size() > 1) needs.add("cohorts: verify each cohort-to-reference correspondence in combined prose");
            } else {
                if (map(p.get("parameters")).isEmpty()) needs.add(role + ": record effective parameters and their provenance, or explain why none apply");
                if (p.containsKey("uncommon_parameters") && !(p.get("uncommon_parameters") instanceof List<?>)) throw new IOException("uncommon_parameters must be a list of parameter keys");
                for (Object key : list(p.get("uncommon_parameters"))) {
                    if (!map(p.get("parameters")).containsKey(key.toString())) throw new IOException("Unknown uncommon parameter: " + key);
                    if (text(map(p.get("parameter_justifications")).get(key.toString())).isBlank()) needs.add(role + ": justify uncommon parameter " + key + " using a recorded reason; never invent a rationale");
                }
            }
        }
        if (!roles.contains("cohorts")) needs.add("Methods must begin with cohort descriptions and appropriate references");
        for (String id : knownCohorts) if (!described.contains(id)) needs.add("Missing Methods description for cohort " + id);
        for (Object raw : list(inventory.get("analyses"))) {
            var a = map(raw); String status = text(a.get("execution_status"));
            if (Set.of("complete", "external_results_supplied_not_refitted").contains(status)) {
                if (!roles.contains(a.get("role"))) needs.add(a.get("source") + ": missing " + a.get("role") + " Methods paragraph");
                if (map(a.get("method_record")).isEmpty()) needs.add(a.get("source") + ": recover missing historical method/parameter provenance");
            }
        }
        for (String role : roles) if (!Set.of("cohorts", "reproducibility").contains(role) && list(inventory.get("analyses")).stream().map(StudyPublicEvidence::map)
                .noneMatch(a -> role.equals(a.get("role")) && Set.of("complete", "external_results_supplied_not_refitted").contains(a.get("execution_status"))))
            needs.add(role + ": no completed or explicitly supplied execution evidence; do not present planned/blocked work as performed");
        if (unclassified > 0) needs.add("Classify " + unclassified + " existing Methods block(s); prose is preserved, not silently replaced");
        return Map.of("policy", policy(), "status", needs.isEmpty() ? "structure_complete_scientific_review_required" : "needs_input", "needs", needs,
            "scientific_review", "Verify cohort/subcohort identity, actual execution, citation support, parameters and provenance. Structural checks cannot certify these scientific facts.");
    }
    static String scaffold(Map<String, Object> template, Map<String, Object> catalog) {
        StringBuilder html = new StringBuilder();
        for (Object raw : list(template.get("sections"))) {
            var p = map(raw); html.append(StudyIntroduction.scaffold(List.of(p)));
            for (Object id : list(p.get("references"))) { var ref = StudyReferenceCatalog.lookup(catalog, id.toString());
                html.append(" <a href=\"reference-catalog.html#").append(ref.get("reference_id")).append("\">[Reference]</a>"); }
        }
        return html.toString();
    }
    static String scaffold(Map<String, Object> template, StudyCitations citations) throws IOException {
        StringBuilder html = new StringBuilder();
        for (Object raw : list(template.get("sections"))) {
            var p = map(raw); html.append("<p>").append(citations.render(text(p.get("text")), list(p.get("references")))).append("</p>");
        }
        return html.toString();
    }
}
