/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Declared study context for result interpretation, never causal attribution. */
final class StudyContextReview {
    private StudyContextReview() { }
    private static final List<String> DATA = List.of("id", "cohort", "cohorts", "type", "modality", "tissue", "stage", "anatomical_site", "compartment",
        "preparation", "assay", "species", "condition", "visit", "time", "batch", "cell_type", "genome_build", "annotation_release", "feature_id_type",
        "coordinate_system", "units", "detection_limit", "feature_universe", "spatial_resolution");
    private static final List<String> ANALYSIS = List.of("model", "family", "link", "formula", "outcome", "response", "predictor", "covariates", "technical_factors",
        "relatedness", "pedigree_adjustment", "pedigree", "transformations", "normalization", "preprocessing", "missingness", "effect_scale", "weights");
    static String instructions() {
        return "Review the study_contract_context against ALL observed results, not only high heterogeneity. "
            + "Ask whether any declared tissue/sample compartment, assay, preparation, population, model, covariate, technical factor, effect scale, "
            + "feature definition, genome/annotation version, cohort relationship, design, selection rule or assumption could explain the pattern. "
            + "For each proposed explanation name the observed result and exact metadata field/dataset IDs, separate declared facts from inference, "
            + "describe a plausible mechanism and alternatives, identify confounding and missing metadata, and propose a discriminating sensitivity check. "
            + "A PBMC-versus-whole-blood contrast may explain some heterogeneity but is confounded with cohort when only one cohort uses PBMC; "
            + "do not declare it causal or dismiss technical scale errors. Rank explanations by supporting and contrary evidence, not merely plausibility. "
            + "Return explanations as scientist_question or supported sensitivity proposals, with the explanation and metadata pointers in rationale. "
            + "Unreviewed free-text/private fields are explicitly omitted from automatic disclosure: ask the scientist for a reviewed metadata summary rather than infer their contents or request patient rows. "
            + "Do not change hypotheses, exclude cohorts, rescale effects, run code or add analyses without the existing approval process. No explanation is an established finding until tested.";
    }
    static Map<String, Object> inventory(Map<String, Object> study, String hypothesisId) {
        List<Object> omitted = new ArrayList<>(), rows = new ArrayList<>(), hypotheses = new ArrayList<>(), cohorts = new ArrayList<>(), relationships = new ArrayList<>();
        // Keep unprojected top-level settings visible to the local reviewer too:
        // a private description or an unfamiliar extension can contain relevant context.
        for (String field : study.keySet()) if (!Set.of("datasets", "hypotheses", "cohorts", "cohort_relationships", "schema_version").contains(field))
            omitted.add(Map.of("field", field, "reason", "Not automatically disclosed; inspect locally for relevant study context and supply a scientist-reviewed summary if needed."));
        Set<String> selected = new LinkedHashSet<>();
        for (Object raw : list(study.get("hypotheses"))) {
            var h = map(raw); if ("retired".equals(h.get("status")) || (hypothesisId != null && !hypothesisId.equals(h.get("id")))) continue;
            collectIds(h.get("inputs"), selected);
            var projected = project(h, List.of("id", "outcome", "exposure", "covariates", "inputs", "estimand", "contrasts", "population", "timepoint", "validation_requirement"), "hypotheses/" + h.get("id"), omitted);
            // Scientific settings only, never resource paths or assistance credentials nested in the analysis plan.
            var analysis = project(map(h.get("analysis")), List.of("mode", "goals", "multiple_testing", "replication", "meta_analysis", "transformations", "model", "covariates", "technical_factors"), "hypotheses/" + h.get("id") + "/analysis", omitted);
            if (!analysis.isEmpty()) projected.put("analysis", analysis); hypotheses.add(projected);
        }
        Set<String> usedCohorts = new LinkedHashSet<>();
        for (Object raw : list(study.get("datasets"))) {
            var d = map(raw); if ("retired".equals(d.get("status")) || (!selected.isEmpty() && !selected.contains(text(d.get("id"))))) continue;
            var row = project(d, DATA, "datasets/" + d.get("id"), omitted);
            var analysis = project(map(d.get("analysis_description")), ANALYSIS, "datasets/" + d.get("id") + "/analysis_description", omitted);
            if (!analysis.isEmpty()) row.put("analysis", analysis);
            var metadata = project(map(d.get("metadata")), List.of("sex", "sex_composition", "ancestry", "population", "age_range", "fasting", "platform", "library_preparation", "normalization", "batch_adjustment"), "datasets/" + d.get("id") + "/metadata", omitted);
            if (!metadata.isEmpty()) row.put("metadata", metadata);
            var columns = map(d.get("columns"));
            var variables = project(columns, List.of("phenotypes", "covariates", "technical_factors", "effect", "standard_error", "p_value"), "datasets/" + d.get("id") + "/columns", omitted);
            var definitions = project(map(columns.get("characteristics")), List.of("sex_constant", "male_value", "female_value", "main_phenotype"), "datasets/" + d.get("id") + "/columns/characteristics", omitted);
            if (!variables.isEmpty()) row.put("variables", variables);
            if (!definitions.isEmpty()) row.put("declared_definitions", definitions);
            usedCohorts.add(text(d.get("cohort"))); for (Object c : list(d.get("cohorts"))) usedCohorts.add(text(c));
            if (!row.isEmpty()) rows.add(row);
        }
        for (Object raw : list(study.get("cohorts"))) {
            var c = map(raw); if (!"retired".equals(c.get("status")) && usedCohorts.contains(text(c.get("id"))))
                cohorts.add(project(c, List.of("id", "parent", "species", "population"), "cohorts/" + c.get("id"), omitted));
        }
        for (Object raw : list(study.get("cohort_relationships"))) {
            var r = map(raw); if (!"retired".equals(r.get("status")) && list(r.get("between")).stream().allMatch(c -> usedCohorts.contains(text(c))))
                relationships.add(project(r, List.of("id", "between", "type", "overlap_fraction"), "cohort_relationships", omitted));
        }
        List<Object> contrasts = new ArrayList<>(); Map<String, Map<String, List<String>>> groups = new LinkedHashMap<>();
        for (Object raw : rows) {
            var row = map(raw); var flat = new LinkedHashMap<String, Object>(); flatten(row, "", flat);
            for (var field : flat.entrySet()) if (!Set.of("id", "cohort", "cohorts").contains(field.getKey()))
                groups.computeIfAbsent(field.getKey(), ignored -> new LinkedHashMap<>()).computeIfAbsent(text(field.getValue()), ignored -> new ArrayList<>()).add(text(row.get("id")));
        }
        groups.forEach((field, values) -> { if (values.size() > 1) contrasts.add(Map.of("field", field, "groups", values,
            "assessment", "Declared difference; possible explanation to investigate, not a demonstrated cause. Missing declarations are not treated as equal.")); });
        Map<String, Object> context = new LinkedHashMap<>();
        context.put("hypotheses", hypotheses); context.put("datasets", rows); context.put("cohorts", cohorts); context.put("cohort_relationships", relationships);
        context.put("declared_contrasts", contrasts); context.put("unreviewed_field_count", omitted.size());
        context.put("coverage", "Bounded structured scientific metadata only. Absence is unknown, not equality. Free text, dictionaries, donor/specimen identifiers, file paths, credentials and participant records are not automatically disclosed; review excluded fields locally and request a scientist-approved summary if relevant.");
        context.put("advisor_instructions", instructions());
        return Map.of("schema_version", 1, "context", context, "local_omitted_fields", omitted,
            "status", "candidate_explanations_require_review", "source", "Scientist-declared resolved study contract; no file contents read and no analysis performed.");
    }
    static Map<String, Object> publicContext(Map<String, Object> inventory) { return map(inventory.get("context")); }
    private static Map<String, Object> project(Map<String, Object> input, List<String> fields, String path, List<Object> omitted) {
        var output = new LinkedHashMap<String, Object>();
        for (var e : input.entrySet()) {
            Object value = fields.contains(e.getKey()) ? safe(e.getValue(), 0) : null;
            if (value != null) output.put(e.getKey(), value);
            else omitted.add(Map.of("field", path + "/" + e.getKey(), "reason", "Outside the automatic structured-metadata projection; review locally, do not infer contents."));
        }
        return output;
    }
    private static Object safe(Object value, int depth) {
        if (depth > 4) return null;
        if (value instanceof Number n) return Double.isFinite(n.doubleValue()) ? n : null;
        if (value instanceof Boolean) return value;
        if (value instanceof String s) return s.length() <= 200 && !s.contains("\\") && !s.matches("(?s).*(?:https?://|file:|[A-Za-z]:/|/home/|/Users/|@).*") ? s : null;
        if (value instanceof List<?> items) {
            if (items.size() > 40) return null; List<Object> result = new ArrayList<>();
            for (Object item : items) { Object v = safe(item, depth + 1); if (v == null) return null; result.add(v); } return result;
        }
        if (value instanceof Map<?, ?> items) {
            if (items.size() > 40) return null; var result = new LinkedHashMap<String, Object>();
            for (var item : items.entrySet()) {
                String key = text(item.getKey());
                if (key.matches("(?i).*(?:file|path|sample|patient|subject|donor|specimen|secret|token|password|dictionary|description|notes).*")) return null;
                Object v = safe(item.getValue(), depth + 1); if (v == null) return null; result.put(key, v);
            } return result;
        }
        return null;
    }
    private static void collectIds(Object value, Set<String> ids) {
        if (value instanceof String s) ids.add(s);
        else if (value instanceof List<?> l) for (Object v : l) collectIds(v, ids);
        else if (value instanceof Map<?, ?> m) for (Object v : m.values()) collectIds(v, ids);
    }
    private static void flatten(Map<String, Object> value, String prefix, Map<String, Object> result) {
        value.forEach((k, v) -> { String key = prefix + k;
            if (v instanceof Map<?, ?>) flatten(map(v), key + "/", result); else result.put(key, v); });
    }
}
