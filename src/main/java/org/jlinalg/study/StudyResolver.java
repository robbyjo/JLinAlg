/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/** Expands a concise study contract into an explicit deterministic contract. */
public final class StudyResolver {
    private StudyResolver() { }

    public static StudyContract resolve(Path source) throws IOException {
        return resolve(StudyValidator.requireValid(source, false), List.of());
    }

    public static StudyContract resolve(StudyContract source,
            List<String> appliedAmendments) {
        Map<String, Object> resolved = source.document();
        Map<String, Object> project = StudyIO.map(resolved.get("project"), "project");
        project.putIfAbsent("objective", "publication");
        resolved.put("project", project);

        Map<String, Object> resources = StudyIO.optionalMap(
            resolved.get("resources"), "resources");
        resources.putIfAbsent("downloads", "automatic");
        Map<String, Object> ld = StudyIO.optionalMap(resources.get("ld"),
            "resources.ld");
        ld.putIfAbsent("database", "1000g-phase3");
        ld.putIfAbsent("population", "auto");
        ld.putIfAbsent("auto_download", Boolean.TRUE);
        resources.put("ld", ld);
        resolved.put("resources", resources);

        if (resolved.containsKey("assistance")) {
            Map<String, Object> assistance = StudyIO.optionalMap(resolved.get("assistance"), "assistance");
            assistance.putIfAbsent("provider", "browser");
            resolved.put("assistance", assistance);
        }

        List<Object> datasets = StudyIO.list(resolved.get("datasets"), "datasets");
        List<String> activeDatasetIds = new ArrayList<>();
        Map<String, Map<String, Object>> datasetById = new LinkedHashMap<>();
        for (Object value : datasets) {
            Map<String, Object> dataset = StudyIO.map(value, "dataset");
            StudyDataDictionary.apply(source.source() == null ? Path.of("").toAbsolutePath()
                : source.source().getParent(), dataset);
            dataset.putIfAbsent("stage", "processed_matrix");
            dataset.putIfAbsent("status", "active");
            String id = StudyIO.requireString(dataset, "id", "dataset");
            if (!"retired".equals(dataset.get("status"))) activeDatasetIds.add(id);
            datasetById.put(id, dataset);
        }
        resolved.put("datasets", new ArrayList<>(datasetById.values()));

        List<Object> hypotheses = StudyIO.list(resolved.get("hypotheses"), "hypotheses");
        List<Object> resolvedHypotheses = new ArrayList<>();
        for (Object value : hypotheses) {
            Map<String, Object> hypothesis = StudyIO.map(value, "hypothesis");
            hypothesis.putIfAbsent("status", "approved");
            Map<String, Object> inputs = StudyIO.optionalMap(
                hypothesis.get("inputs"), "hypothesis.inputs");
            if (inputs.isEmpty()) inputs.put("discovery", List.copyOf(activeDatasetIds));
            normalizeRoles(inputs);
            hypothesis.put("inputs", inputs);
            Map<String, Object> analysis = StudyIO.optionalMap(
                hypothesis.get("analysis"), "hypothesis.analysis");
            analysis.putIfAbsent("mode", "auto");
            analysis.putIfAbsent("direction", "two_sided");
            analysis.putIfAbsent("multiple_testing",
                new LinkedHashMap<>(Map.of("method", "BH", "alpha", 0.05)));
            List<Object> inferred = defaultGoals(inputs, datasetById,
                project.get("objective"));
            if (!analysis.containsKey("goals")) analysis.put("goals", inferred);
            else if ("auto".equalsIgnoreCase(String.valueOf(analysis.get("mode")))) {
                LinkedHashSet<Object> goals = new LinkedHashSet<>(
                    StudyIO.list(analysis.get("goals"), "analysis.goals"));
                goals.addAll(inferred);
                analysis.put("goals", new ArrayList<>(goals));
            }
            hypothesis.put("analysis", analysis);
            resolvedHypotheses.add(hypothesis);
        }
        resolved.put("hypotheses", resolvedHypotheses);

        List<Object> assumptions = assumptions(resolved, datasetById);
        resolved.put("assumptions", assumptions);

        Map<String, Object> resolution = new LinkedHashMap<>();
        resolution.put("schema_version", 1);
        resolution.put("defaults_version", 1);
        resolution.put("source_study_hash", source.hash());
        resolution.put("applied_amendments", List.copyOf(appliedAmendments));
        resolution.put("default_seed_policy", "derived_from_approved_plan_hash");
        resolved.put("_resolution", resolution);
        String resolvedHash = hashWithoutResolvedHash(resolved);
        resolution.put("resolved_study_hash", resolvedHash);

        StudyValidation validation = StudyValidator.validateDocument(
            resolved, source.source() == null ? Path.of("").toAbsolutePath()
                : source.source().getParent(), false);
        if (!validation.valid())
            throw new IllegalArgumentException(
                "Resolved study is invalid:\n" + validation.report());
        return new StudyContract(source.source(), resolved);
    }

    static String hashWithoutResolvedHash(Map<String, Object> document) {
        Map<String, Object> copy = StudyIO.copyMap(document);
        Map<String, Object> resolution = StudyIO.optionalMap(
            copy.get("_resolution"), "_resolution");
        resolution.remove("resolved_study_hash");
        copy.put("_resolution", resolution);
        return StudyIO.hash(copy);
    }

    private static void normalizeRoles(Map<String, Object> inputs) {
        for (String role : new ArrayList<>(inputs.keySet())) {
            Object value = inputs.get(role);
            if (value instanceof String text) inputs.put(role, List.of(text));
            else if (value instanceof Collection<?>)
                inputs.put(role, StudyIO.list(value, "inputs." + role));
        }
    }

    private static List<Object> defaultGoals(Map<String, Object> inputs,
            Map<String, Map<String, Object>> datasets, Object objective) {
        Set<String> selected = new LinkedHashSet<>();
        inputs.values().forEach(value -> {
            if (value instanceof Collection<?> items)
                items.forEach(item -> selected.add(String.valueOf(item)));
        });
        Set<String> types = new HashSet<>();
        boolean matrix = false;
        for (String id : selected) {
            Map<String, Object> dataset = datasets.get(id);
            if (dataset == null) continue;
            types.add(String.valueOf(dataset.get("type")).toLowerCase(Locale.ROOT));
            String stage = String.valueOf(dataset.get("stage"));
            matrix |= Set.of("raw_quantified", "quantified_matrix",
                "processed_matrix").contains(stage);
        }
        List<Object> goals = new ArrayList<>();
        if (matrix) goals.add("association");
        boolean hasReplication = inputs.containsKey("replication")
            && !StudyIO.list(inputs.get("replication"), "replication").isEmpty();
        if (hasReplication || selected.size() > 1) goals.add("discovery_replication");
        if (selected.size() > 1) goals.add("meta_analysis");
        if (types.stream().anyMatch(type -> Set.of("transcriptomics", "expression",
                "mirna", "lincrna", "lncrna", "extracellular_rna",
                "epigenomics", "ewas", "proteomics", "single_cell_rna",
                "spatial_transcriptomics").contains(type))) {
            goals.add("gene_ontology");
            goals.add("pathway_enrichment");
        }
        if (types.stream().anyMatch(type ->
                Set.of("gwas", "genomics").contains(type)))
            goals.add("variant_scoring");
        Set<String> objectives = new HashSet<>();
        if (objective instanceof Collection<?> values)
            values.forEach(value -> objectives.add(normalize(value)));
        else objectives.add(normalize(objective));
        if (objectives.stream().anyMatch(value -> value.contains("publication")
                || value.contains("drug") || value.contains("causal")))
            goals.add("genetic_support");
        goals.add("literature_support");
        goals.add("literature_search");
        goals.add("bioinformatics_investigation");
        if (goals.contains("genetic_support")) goals.add("source_discovery");
        if (objectives.stream().anyMatch(value -> value.contains("drug")))
            goals.add("candidate_prioritization");
        return goals;
    }

    private static String normalize(Object value) {
        return String.valueOf(value).toLowerCase(Locale.ROOT)
            .replace('-', '_').replace(' ', '_');
    }

    private static List<Object> assumptions(Map<String, Object> document,
            Map<String, Map<String, Object>> datasets) {
        Map<String, Map<String, Object>> cohorts = new LinkedHashMap<>();
        for (Object value : StudyIO.optionalList(document.get("cohorts"), "cohorts")) {
            Map<String, Object> cohort = StudyIO.map(value, "cohort");
            cohorts.put(StudyIO.requireString(cohort, "id", "cohort"), cohort);
        }
        Set<String> explicit = new HashSet<>();
        for (Object value : StudyIO.optionalList(
                document.get("cohort_relationships"), "cohort_relationships")) {
            Map<String, Object> relationship = StudyIO.map(value, "relationship");
            List<Object> between = StudyIO.list(relationship.get("between"), "between");
            for (int left = 0; left < between.size(); left++)
                for (int right = left + 1; right < between.size(); right++)
                    explicit.add(pair(String.valueOf(between.get(left)),
                        String.valueOf(between.get(right))));
        }
        Set<String> used = new TreeSet<>();
        for (Map<String, Object> dataset : datasets.values()) {
            if ("retired".equals(dataset.get("status"))) continue;
            Object cohort = dataset.get("cohort");
            if (cohort instanceof String text) used.add(text);
            for (Object value : StudyIO.optionalList(
                    dataset.get("cohorts"), "dataset.cohorts"))
                used.add(String.valueOf(value));
        }
        List<String> ids = new ArrayList<>(used);
        List<Object> result = new ArrayList<>();
        for (Object value : StudyIO.optionalList(
                document.get("assumptions"), "assumptions")) {
            Map<String, Object> assumption = StudyIO.map(value, "assumption");
            assumption.putIfAbsent("review_status", "pending");
            result.add(assumption);
        }
        for (int left = 0; left < ids.size(); left++) {
            for (int right = left + 1; right < ids.size(); right++) {
                String a = ids.get(left), b = ids.get(right);
                if (explicit.contains(pair(a, b))) continue;
                String parentA = parent(cohorts, a);
                String parentB = parent(cohorts, b);
                boolean sameParent = parentA != null && parentA.equals(parentB);
                Map<String, Object> assumption = new LinkedHashMap<>();
                assumption.put("id", "relationship-" + safe(a) + "-" + safe(b));
                assumption.put("between", List.of(a, b));
                assumption.put("type", sameParent
                    ? "potentially_dependent" : "independent");
                assumption.put("statement", sameParent
                    ? a + " and " + b + " may be dependent because they share parent " + parentA
                    : a + " and " + b + " are assumed independent");
                assumption.put("source", sameParent
                    ? "shared_parent_requires_confirmation"
                    : "default_due_to_undeclared_relationship");
                assumption.put("review_status", "pending");
                result.add(assumption);
            }
        }
        return result;
    }

    private static String parent(Map<String, Map<String, Object>> cohorts,
            String id) {
        Map<String, Object> cohort = cohorts.get(id);
        return cohort == null || !(cohort.get("parent") instanceof String parent)
            ? null : parent;
    }

    private static String pair(String left, String right) {
        return left.compareTo(right) < 0 ? left + "\u0000" + right
            : right + "\u0000" + left;
    }

    private static String safe(String value) {
        return value.toLowerCase(Locale.ROOT).replaceAll("[^a-z0-9._-]+", "-");
    }
}
