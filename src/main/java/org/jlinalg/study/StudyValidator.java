/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import static org.jlinalg.study.StudyValidation.Severity.ERROR;
import static org.jlinalg.study.StudyValidation.Severity.WARNING;

/** Strict validation for versioned study contracts. */
public final class StudyValidator {
    private StudyValidator() { }

    private static final Pattern ID =
        Pattern.compile("[A-Za-z][A-Za-z0-9._-]{0,127}");
    private static final Set<String> TOP = Set.of(
        "schema_version", "project", "cohorts", "cohort_relationships",
        "datasets", "hypotheses", "analysis_requests",
        "assumptions", "resources", "assistance", "_resolution");
    private static final Set<String> PROJECT = Set.of(
        "id", "title", "description", "objective", "investigators",
        "organization", "tags", "data_policy");
    private static final Set<String> COHORT = Set.of(
        "id", "parent", "name", "literature_name", "references", "description", "species", "population", "metadata",
        "status", "retired_reason");
    private static final Set<String> RELATIONSHIP = Set.of(
        "id", "between", "type", "description", "overlap_fraction",
        "shared_participants", "metadata", "status", "retired_reason");
    private static final Set<String> DATASET = Set.of(
        "id", "cohort", "cohorts", "type", "modality", "tissue",
        "stage", "description", "files", "columns", "analysis_description",
        "provenance", "source", "anatomical_site", "compartment",
        "preparation", "assay", "species", "donor", "specimen",
        "condition", "visit", "time", "batch", "cell_type", "genome_build",
        "annotation_release", "feature_id_type", "coordinate_system", "units", "detection_limit", "feature_universe",
        "spatial_resolution", "metadata", "status", "retired_reason");
    private static final Set<String> HYPOTHESIS = Set.of(
        "id", "description", "outcome", "exposure", "covariates", "inputs",
        "analysis", "depends_on", "origin", "validation_requirement",
        "status", "estimand", "contrasts", "population", "timepoint",
        "metadata", "retired_reason");
    private static final Set<String> DEPENDENCY = Set.of(
        "hypothesis", "criterion", "on_failure");
    private static final Set<String> ANALYSIS_REQUEST = Set.of(
        "id", "type", "for_hypothesis", "for_hypotheses", "description",
        "reason", "origin_observation", "status", "metadata",
        "retired_reason");
    private static final Set<String> ASSUMPTION = Set.of(
        "id", "between", "type", "statement", "source", "review_status",
        "metadata");
    private static final Set<String> GOALS = Set.of(
        "association", "effect_comparison", "discovery_replication", "meta_analysis",
        "gene_ontology", "pathway_enrichment", "variant_scoring",
        "fine_mapping", "colocalization", "mendelian_randomization",
        "genetic_support", "literature_support", "candidate_prioritization",
        "source_discovery", "bioinformatics_investigation", "literature_search");
    private static final Set<String> STAGES = Set.of(
        "raw_quantified", "quantified_matrix", "processed_matrix",
        "association_results", "annotated_results", "meta_analysis_results",
        "enrichment_results", "variant_scores", "fine_mapping_results",
        "colocalization_results", "mr_results", "evidence_results");
    private static final Set<String> TYPES = Set.of(
        "transcriptomics", "expression", "mirna", "lincrna", "lncrna",
        "extracellular_rna", "epigenomics", "ewas", "gwas", "genomics",
        "proteomics", "metabolomics", "single_cell_rna", "single_cell_protein",
        "spatial_transcriptomics", "spatial_proteomics", "microbiome",
        "multimodal", "phenotype", "generic");

    public static StudyValidation validate(Path path, boolean checkFiles)
            throws IOException {
        List<StudyValidation.Issue> issues = new ArrayList<>();
        Map<String, Object> document;
        try {
            document = StudyIO.load(path);
        } catch (IllegalArgumentException failure) {
            issues.add(issue(ERROR, "$", failure.getMessage()));
            return new StudyValidation(issues);
        }
        validate(document, path.toAbsolutePath().normalize().getParent(),
            checkFiles, issues);
        return new StudyValidation(issues);
    }

    public static StudyContract requireValid(Path path, boolean checkFiles)
            throws IOException {
        Map<String, Object> document = StudyIO.load(path);
        List<StudyValidation.Issue> issues = new ArrayList<>();
        validate(document, path.toAbsolutePath().normalize().getParent(),
            checkFiles, issues);
        StudyValidation validation = new StudyValidation(issues);
        if (!validation.valid())
            throw new IllegalArgumentException(
                "Study validation failed:\n" + validation.report());
        return new StudyContract(path, document);
    }

    static StudyValidation validateDocument(Map<String, Object> document,
            Path root, boolean checkFiles) {
        List<StudyValidation.Issue> issues = new ArrayList<>();
        validate(document, root, checkFiles, issues);
        return new StudyValidation(issues);
    }

    private static void validate(Map<String, Object> document, Path root,
            boolean checkFiles, List<StudyValidation.Issue> issues) {
        unknown(document, TOP, "$", issues);
        if (!Integer.valueOf(1).equals(document.get("schema_version")))
            issues.add(issue(ERROR, "schema_version",
                "schema_version must be the integer 1"));

        Map<String, Object> project = mapping(document.get("project"),
            "project", true, issues);
        unknown(project, PROJECT, "project", issues);
        identifier(project.get("id"), "project.id", issues);
        requiredText(project.get("title"), "project.title", issues);
        requiredText(project.get("description"), "project.description", issues);
        scalarOrList(project.get("objective"), "project.objective", issues);
        Map<String, Object> privacy = mapping(project.get("data_policy"), "project.data_policy", false, issues);
        Map<String, Object> llmPrivacy = mapping(privacy.get("llm"), "project.data_policy.llm", false, issues);
        enumValue(privacy.get("remote_transmission"), "project.data_policy.remote_transmission", Set.of("aggregate_results", "deny_all"), issues);
        enumValue(llmPrivacy.get("remote_transmission"), "project.data_policy.llm.remote_transmission", Set.of("aggregate_results", "deny_all"), issues);
        if (llmPrivacy.containsKey("remote_patient_level_data"))
            enumValue(llmPrivacy.get("remote_patient_level_data"), "project.data_policy.llm.remote_patient_level_data", Set.of("deny"), issues);
        if (llmPrivacy.containsKey("require_human_approval") && !(llmPrivacy.get("require_human_approval") instanceof Boolean))
            issues.add(issue(ERROR, "project.data_policy.llm.require_human_approval", "must be true or false"));
        resources(document.get("resources"), issues);
        assistance(document.get("assistance"), issues);

        Map<String, Map<String, Object>> cohorts = objects(
            document.get("cohorts"), "cohorts", COHORT, false, issues);
        for (Map.Entry<String, Map<String, Object>> entry : cohorts.entrySet()) {
            for (String field : List.of("name", "literature_name")) if (entry.getValue().containsKey(field)) {
                Object value = entry.getValue().get(field);
                if (!(value instanceof String s) || s.isBlank() || s.length() > 150 || s.chars().anyMatch(Character::isISOControl))
                    issues.add(issue(ERROR, "cohorts." + entry.getKey() + "." + field, "must be nonblank public text of at most 150 characters without control characters"));
            }
            if (entry.getValue().containsKey("references")) for (Object ref : sequence(entry.getValue().get("references"), "cohorts." + entry.getKey() + ".references", true, issues))
                requiredText(ref, "cohorts." + entry.getKey() + ".references", issues);
            Object parent = entry.getValue().get("parent");
            if (parent != null && !(parent instanceof String text && ID.matcher(text).matches()))
                issues.add(issue(ERROR, "cohorts." + entry.getKey() + ".parent",
                    "parent must be a valid cohort identifier"));
        }

        List<Map<String, Object>> relationships = objectList(
            document.get("cohort_relationships"), "cohort_relationships",
            RELATIONSHIP, false, issues);
        for (int index = 0; index < relationships.size(); index++) {
            Map<String, Object> relationship = relationships.get(index);
            String path = "cohort_relationships[" + index + "]";
            List<Object> between = sequence(relationship.get("between"),
                path + ".between", true, issues);
            if (between.size() < 2)
                issues.add(issue(ERROR, path + ".between",
                    "at least two cohort identifiers are required"));
            for (Object value : between)
                if (!(value instanceof String id) || !ID.matcher(id).matches()
                        || !cohorts.isEmpty() && !cohorts.containsKey(id))
                    issues.add(issue(ERROR, path + ".between",
                        "unknown or invalid cohort identifier: " + value));
            requiredText(relationship.get("type"), path + ".type", issues);
        }

        Map<String, Map<String, Object>> datasets = objects(
            document.get("datasets"), "datasets", DATASET, true, issues);
        for (Map.Entry<String, Map<String, Object>> entry : datasets.entrySet())
            dataset(entry.getKey(), entry.getValue(), cohorts, root,
                checkFiles, issues);

        Map<String, Map<String, Object>> hypotheses = objects(
            document.get("hypotheses"), "hypotheses", HYPOTHESIS, true, issues);
        for (Map.Entry<String, Map<String, Object>> entry : hypotheses.entrySet())
            hypothesis(entry.getKey(), entry.getValue(), datasets, issues);
        for (Map.Entry<String, Map<String, Object>> entry : hypotheses.entrySet()) {
            List<Object> dependencies = sequence(entry.getValue().get("depends_on"),
                "hypotheses." + entry.getKey() + ".depends_on", false, issues);
            for (Object dependency : dependencies) {
                Map<String, Object> value = mapping(dependency,
                    "hypotheses." + entry.getKey() + ".depends_on", true, issues);
                unknown(value, DEPENDENCY,
                    "hypotheses." + entry.getKey() + ".depends_on", issues);
                Object id = value.get("hypothesis");
                if (!(id instanceof String text) || !hypotheses.containsKey(text))
                    issues.add(issue(ERROR,
                        "hypotheses." + entry.getKey() + ".depends_on",
                        "dependency references an unknown hypothesis: " + id));
            }
        }
        Map<String, Map<String, Object>> requests = objects(
            document.get("analysis_requests"), "analysis_requests",
            ANALYSIS_REQUEST, false, issues);
        for (Map.Entry<String, Map<String, Object>> entry : requests.entrySet()) {
            String path = "analysis_requests." + entry.getKey();
            Map<String, Object> request = entry.getValue();
            requiredText(request.get("type"), path + ".type", issues);
            requiredText(request.get("description"), path + ".description", issues);
            Object one = request.get("for_hypothesis");
            if (one != null && (!(one instanceof String id)
                    || !hypotheses.containsKey(id)))
                issues.add(issue(ERROR, path + ".for_hypothesis",
                    "unknown hypothesis: " + one));
            Object many = request.get("for_hypotheses");
            if (many != null)
                for (Object value : sequence(many,
                        path + ".for_hypotheses", true, issues))
                    if (!(value instanceof String id)
                            || !hypotheses.containsKey(id))
                        issues.add(issue(ERROR, path + ".for_hypotheses",
                            "unknown hypothesis: " + value));
        }
        Map<String, Map<String, Object>> assumptions = objects(
            document.get("assumptions"), "assumptions", ASSUMPTION,
            false, issues);
        for (Map.Entry<String, Map<String, Object>> entry : assumptions.entrySet())
            requiredText(entry.getValue().get("statement"),
                "assumptions." + entry.getKey() + ".statement", issues);
    }

    private static void resources(Object value,
            List<StudyValidation.Issue> issues) {
        if (value == null) return;
        Map<String, Object> resources = mapping(value, "resources", true, issues);
        unknown(resources, Set.of("downloads", "directory", "ld", "network"),
            "resources", issues);
        enumValue(resources.get("network"), "resources.network", Set.of("online", "offline"), issues);
        enumValue(resources.get("downloads"), "resources.downloads",
            Set.of("automatic", "manual", "deny"), issues);
        Object directory = resources.get("directory");
        if (directory != null) requiredText(directory,
            "resources.directory", issues);
        if (resources.containsKey("ld")) {
            Map<String, Object> ld = mapping(resources.get("ld"),
                "resources.ld", true, issues);
            unknown(ld, Set.of("database", "population", "location",
                "auto_download"), "resources.ld", issues);
            for (String key : List.of("database", "population", "location"))
                if (ld.containsKey(key)) requiredText(ld.get(key),
                    "resources.ld." + key, issues);
            Object automatic = ld.get("auto_download");
            if (automatic != null && !(automatic instanceof Boolean))
                issues.add(issue(ERROR, "resources.ld.auto_download",
                    "must be true or false"));
        }
    }

    private static void assistance(Object value, List<StudyValidation.Issue> issues) {
        if (value == null) return;
        Map<String, Object> config = mapping(value, "assistance", true, issues);
        unknown(config, Set.of("enabled", "provider", "model", "api_key_env", "endpoint", "auth", "oauth_token_env", "response_files",
            "next_steps", "related_phenotypes", "literature", "manuscript", "allow_analyses", "timeout_seconds"),
            "assistance", issues);
        enumValue(config.get("provider"), "assistance.provider", Set.of("browser", "openai", "gemini", "openai_compatible"), issues);
        for (String key : List.of("enabled", "next_steps", "related_phenotypes", "literature", "manuscript"))
            if (config.containsKey(key) && !(config.get(key) instanceof Boolean))
                issues.add(issue(ERROR, "assistance." + key, "must be true or false"));
        if (Boolean.TRUE.equals(config.get("enabled")) && !"browser".equals(config.getOrDefault("provider", "browser"))) {
            requiredText(config.get("model"), "assistance.model", issues);
        }
        for (String key : List.of("api_key_env", "oauth_token_env")) {
        Object env = config.get(key);
        if (env != null && (!(env instanceof String text) || !text.matches("[A-Za-z_][A-Za-z0-9_]*")))
            issues.add(issue(ERROR, "assistance." + key, "must be an environment variable name, not a credential"));
        }
        enumValue(config.get("auth"), "assistance.auth", Set.of("api_key", "oauth"), issues);
        if ("oauth".equals(config.get("auth"))) {
            if (!"openai_compatible".equals(config.get("provider"))) issues.add(issue(ERROR, "assistance.auth", "OAuth requires a provider-authorized OpenAI-compatible endpoint"));
            requiredText(config.get("oauth_token_env"), "assistance.oauth_token_env", issues);
            if (config.containsKey("api_key_env")) issues.add(issue(ERROR, "assistance.api_key_env", "Do not combine API key and OAuth authentication"));
        } else if (config.containsKey("oauth_token_env")) issues.add(issue(ERROR, "assistance.auth", "oauth_token_env requires auth: oauth"));
        if ("openai_compatible".equals(config.get("provider"))) {
            try {
                java.net.URI uri = java.net.URI.create(String.valueOf(config.get("endpoint")));
                if (!"https".equals(uri.getScheme()) || uri.getHost() == null || uri.getUserInfo() != null || uri.getQuery() != null || uri.getFragment() != null || !uri.getPath().endsWith("/chat/completions")) throw new IllegalArgumentException();
            } catch (IllegalArgumentException bad) { issues.add(issue(ERROR, "assistance.endpoint", "requires explicit HTTPS /chat/completions URL without credentials, query or fragment")); }
        } else if (config.containsKey("endpoint")) issues.add(issue(ERROR, "assistance.endpoint", "custom endpoints require provider: openai_compatible"));
        if ("browser".equals(config.getOrDefault("provider", "browser")) && (config.containsKey("api_key_env") || config.containsKey("auth") || config.containsKey("model")))
            issues.add(issue(ERROR, "assistance", "browser handoff does not use model or credential settings; model label is recorded on import"));
        List<Object> responseFiles = sequence(config.get("response_files"), "assistance.response_files", false, issues);
        if (!responseFiles.isEmpty() && !"browser".equals(config.getOrDefault("provider", "browser"))) issues.add(issue(ERROR, "assistance.response_files", "browser receipts require provider: browser"));
        Set<Object> seenResponses = new HashSet<>();
        for (Object file : responseFiles) {
            requiredText(file, "assistance.response_files", issues);
            if (file instanceof String path && (!Path.of(path).getFileName().toString().equals("receipt.yaml") || !seenResponses.add(path)))
                issues.add(issue(ERROR, "assistance.response_files", "use unique paths ending in receipt.yaml"));
        }
        Object timeout = config.get("timeout_seconds");
        if (timeout != null && (!(timeout instanceof Number n) || n.doubleValue() < 1
                || n.doubleValue() > 300 || n.doubleValue() != n.intValue()))
            issues.add(issue(ERROR, "assistance.timeout_seconds", "must be an integer from 1 to 300"));
        Set<Object> unique = new HashSet<>();
        for (Object analysis : sequence(config.get("allow_analyses"), "assistance.allow_analyses", false, issues)) {
            requiredText(analysis, "assistance.allow_analyses", issues);
            enumValue(analysis, "assistance.allow_analyses", Set.of("leave_one_cohort_out", "fixed_random_comparison"), issues);
            if (!unique.add(analysis)) issues.add(issue(ERROR, "assistance.allow_analyses", "duplicate analysis"));
        }
    }

    private static void dataset(String id, Map<String, Object> dataset,
            Map<String, Map<String, Object>> cohorts, Path root,
            boolean checkFiles, List<StudyValidation.Issue> issues) {
        String path = "datasets." + id;
        Object cohort = dataset.get("cohort");
        if (cohort != null && (!(cohort instanceof String text)
                || !ID.matcher(text).matches()
                || !cohorts.isEmpty() && !cohorts.containsKey(text)))
            issues.add(issue(ERROR, path + ".cohort",
                "cohort is invalid or is not declared"));
        if (dataset.containsKey("cohorts"))
            for (Object value : sequence(dataset.get("cohorts"),
                    path + ".cohorts", true, issues))
                if (!(value instanceof String text) || !cohorts.isEmpty()
                        && !cohorts.containsKey(text))
                    issues.add(issue(ERROR, path + ".cohorts",
                        "unknown cohort: " + value));
        String type = text(dataset.get("type"));
        if (type == null) issues.add(issue(ERROR, path + ".type", "type is required"));
        else if (!TYPES.contains(type.toLowerCase()))
            issues.add(issue(WARNING, path + ".type",
                "unregistered modality; planning will use the generic profile"));
        Object stageValue = dataset.get("stage");
        if (stageValue != null && (!(stageValue instanceof String stage)
                || !STAGES.contains(stage.toLowerCase())))
            issues.add(issue(ERROR, path + ".stage",
                "unsupported processing stage: " + stageValue));
        Map<String, Object> files = mapping(dataset.get("files"),
            path + ".files", true, issues);
        for (Map.Entry<String, Object> file : files.entrySet()) {
            if (!(file.getValue() instanceof String name) || name.isBlank()) {
                issues.add(issue(ERROR, path + ".files." + file.getKey(),
                    "file path must be a nonblank string"));
                continue;
            }
            if (checkFiles) try {
                Path resolved = root.resolve(name).normalize();
                if (!Files.isRegularFile(resolved))
                    issues.add(issue(ERROR, path + ".files." + file.getKey(),
                        "file is absent: " + resolved));
            } catch (InvalidPathException invalid) {
                issues.add(issue(ERROR, path + ".files." + file.getKey(),
                    "invalid path: " + name));
            }
        }
        String stage = text(dataset.getOrDefault("stage", "processed_matrix"));
        if (stage != null && Set.of("raw_quantified", "quantified_matrix",
                "processed_matrix").contains(stage.toLowerCase())
                && !files.containsKey("matrix"))
            issues.add(issue(ERROR, path + ".files",
                "matrix input stages require files.matrix"));
        if (stage != null && stage.toLowerCase().endsWith("_results")
                && !files.containsKey("results"))
            issues.add(issue(ERROR, path + ".files",
                "result input stages require files.results"));
        Map<String, Object> columns = mapping(dataset.get("columns"),
            path + ".columns", false, issues);
        enumValue(dataset.get("feature_id_type"), path + ".feature_id_type",
            Set.of("auto", "probe", "ensembl_gene", "ensembl_transcript",
                "gene_symbol", "variant", "protein", "metabolite"), issues);
        enumValue(dataset.get("coordinate_system"), path + ".coordinate_system",
            Set.of("1-based-closed", "0-based-half-open"), issues);
        if (columns.containsKey("validation")) {
            Map<String, Object> validation = mapping(columns.get("validation"),
                path + ".columns.validation", true, issues);
            unknown(validation, Set.of("finite", "positive", "nonnegative",
                "probability"), path + ".columns.validation", issues);
            for (String rule : List.of("finite", "positive", "nonnegative",
                    "probability"))
                if (validation.containsKey(rule))
                    for (Object column : sequence(validation.get(rule),
                            path + ".columns.validation." + rule, true, issues))
                        requiredText(column,
                            path + ".columns.validation." + rule, issues);
        }
    }

    private static void hypothesis(String id, Map<String, Object> hypothesis,
            Map<String, Map<String, Object>> datasets,
            List<StudyValidation.Issue> issues) {
        String path = "hypotheses." + id;
        requiredText(hypothesis.get("description"), path + ".description", issues);
        Map<String, Object> inputs = mapping(hypothesis.get("inputs"),
            path + ".inputs", false, issues);
        if (!inputs.isEmpty()) for (Map.Entry<String, Object> role : inputs.entrySet()) {
            List<Object> values = role.getValue() instanceof String
                ? List.of(role.getValue())
                : sequence(role.getValue(), path + ".inputs." + role.getKey(),
                    true, issues);
            for (Object value : values)
                if (!(value instanceof String dataset) || !datasets.containsKey(dataset))
                    issues.add(issue(ERROR, path + ".inputs." + role.getKey(),
                        "unknown dataset: " + value));
        }
        Object covariates = hypothesis.get("covariates");
        if (covariates != null)
            sequence(covariates, path + ".covariates", true, issues);
        Map<String, Object> analysis = mapping(hypothesis.get("analysis"),
            path + ".analysis", false, issues);
        enumValue(analysis.get("mode"), path + ".analysis.mode",
            Set.of("auto", "manual"), issues);
        Object goals = analysis.get("goals");
        if (analysis.containsKey("bioinformatics")) {
            String where = path + ".analysis.bioinformatics";
            Map<String, Object> bio = mapping(analysis.get("bioinformatics"), where, true, issues);
            unknown(bio, Set.of("max_features", "max_pages", "max_sequence_bp", "max_regional_snps", "window_bp", "fdr", "phenotype_terms", "relevant_tissues", "literature_terms", "remote_feature_queries", "genetic_pairs", "mr", "colocalization"), where, issues);
            if (bio.containsKey("remote_feature_queries")
                    && !(bio.get("remote_feature_queries") instanceof Boolean))
                issues.add(issue(ERROR, where + ".remote_feature_queries",
                    "must be true or false"));
            for (String field : List.of("max_features", "max_pages", "max_sequence_bp", "max_regional_snps", "window_bp")) if (bio.containsKey(field)) {
                Object number = bio.get(field);
                int max = Set.of("max_features", "max_regional_snps").contains(field) ? 100 : field.equals("max_pages") ? 20 : field.equals("window_bp") ? 5000000 : 1000000;
                if (!(number instanceof Number value) || value.doubleValue() != value.intValue() || value.intValue() < (field.equals("window_bp") ? 0 : 1) || value.intValue() > max)
                    issues.add(issue(ERROR, where + "." + field, "must be an integer from 1 to " + max));
            }
            if (bio.containsKey("fdr") && (!(bio.get("fdr") instanceof Number value) || !(value.doubleValue() > 0 && value.doubleValue() <= 1)))
                issues.add(issue(ERROR, where + ".fdr", "must be in (0,1]"));
            for (String field : List.of("phenotype_terms", "relevant_tissues", "literature_terms")) if (bio.containsKey(field)) {
                List<Object> values = sequence(bio.get(field), where + "." + field, true, issues);
                int limit = field.equals("phenotype_terms") ? 20 : 10;
                if (values.size() > limit) issues.add(issue(ERROR, where + "." + field, "at most " + limit + " labels are allowed"));
                for (Object value : values) {
                    requiredText(value, where + "." + field, issues);
                    if (value instanceof String label && (label.length() > 150 || label.chars().anyMatch(Character::isISOControl)))
                        issues.add(issue(ERROR, where + "." + field, "labels must be bounded text without control characters"));
                }
            }
            if (bio.containsKey("genetic_pairs")) {
                List<Object> pairs = sequence(bio.get("genetic_pairs"), where + ".genetic_pairs", true, issues);
                if (pairs.size() > 100) issues.add(issue(ERROR, where + ".genetic_pairs", "at most 100 reviewed regional pairs are allowed"));
                Set<String> seen = new HashSet<>();
                for (Object value : pairs) {
                    Map<String, Object> pair = mapping(value, where + ".genetic_pairs", true, issues);
                    Object key = pair.get("id");
                    if (!(key instanceof String name) || !name.matches("[A-Za-z0-9][A-Za-z0-9_-]{0,63}") || !seen.add(name))
                        issues.add(issue(ERROR, where + ".genetic_pairs.id", "IDs must be unique safe identifiers"));
                    for (String field : List.of("feature_id", "rationale")) requiredText(pair.get(field), where + ".genetic_pairs." + field, issues);
                }
            }
            if (bio.containsKey("mr")) {
                String mrWhere = where + ".mr";
                Map<String, Object> mr = mapping(bio.get("mr"), mrWhere, true, issues);
                unknown(mr, Set.of("exposure_dataset", "exposure_ancestry", "exposure_source_id", "exposure_citation",
                    "outcome_provider", "outcome_trait", "outcome_ancestry", "ld_provider", "ld_population",
                    "sample_overlap", "sample_overlap_basis", "instrument_validity_rationale", "instrument_p",
                    "minimum_f", "clump_r2", "max_features"), mrWhere, issues);
                requiredText(mr.get("exposure_ancestry"), mrWhere + ".exposure_ancestry", issues);
                for (String field : List.of("exposure_dataset", "exposure_source_id", "exposure_citation", "outcome_trait",
                        "outcome_ancestry", "ld_population", "sample_overlap_basis", "instrument_validity_rationale"))
                    if (mr.containsKey(field)) requiredText(mr.get(field), mrWhere + "." + field, issues);
                enumValue(mr.get("outcome_provider"), mrWhere + ".outcome_provider", Set.of("auto", "pan-ukb"), issues);
                enumValue(mr.get("ld_provider"), mrWhere + ".ld_provider", Set.of("1000g-30x-grch38"), issues);
                enumValue(mr.get("sample_overlap"), mrWhere + ".sample_overlap", Set.of("none"), issues);
                if (mr.containsKey("max_features") && (!(mr.get("max_features") instanceof Number value)
                        || value.doubleValue() != value.intValue() || value.intValue() < 1 || value.intValue() > 10000))
                    issues.add(issue(ERROR, mrWhere + ".max_features", "must be an integer from 1 to 10000"));
                if (mr.containsKey("instrument_p") && (!(mr.get("instrument_p") instanceof Number value)
                        || !(value.doubleValue() > 0 && value.doubleValue() <= 5e-8)))
                    issues.add(issue(ERROR, mrWhere + ".instrument_p", "must be in (0,5e-8]"));
                if (mr.containsKey("minimum_f") && (!(mr.get("minimum_f") instanceof Number value) || value.doubleValue() < 10))
                    issues.add(issue(ERROR, mrWhere + ".minimum_f", "must be at least 10"));
                if (mr.containsKey("clump_r2") && (!(mr.get("clump_r2") instanceof Number value)
                        || !(value.doubleValue() > 0 && value.doubleValue() <= 0.01)))
                    issues.add(issue(ERROR, mrWhere + ".clump_r2", "must be in (0,0.01]"));
            }
            if (bio.containsKey("colocalization")) {
                String colocWhere = where + ".colocalization";
                Map<String, Object> coloc = mapping(bio.get("colocalization"), colocWhere, true, issues);
                unknown(coloc, Set.of("exposure_dataset", "exposure_ancestry", "exposure_source_id", "exposure_citation",
                    "outcome_ancestry", "outcome_provider", "ld_provider", "ld_population", "related_phenotypes",
                    "max_related_phenotypes", "max_colocalization_traits", "max_exposure_variants", "region_window_bp",
                    "allow_partial_regions", "minimum_exposure_overlap", "sample_overlap", "sample_overlap_basis"), colocWhere, issues);
                requiredText(coloc.get("exposure_ancestry"), colocWhere + ".exposure_ancestry", issues);
                for (String field : List.of("exposure_dataset", "exposure_source_id", "exposure_citation", "outcome_ancestry",
                        "ld_population", "sample_overlap_basis"))
                    if (coloc.containsKey(field)) requiredText(coloc.get(field), colocWhere + "." + field, issues);
                enumValue(coloc.get("outcome_provider"), colocWhere + ".outcome_provider", Set.of("auto", "gwas-catalog"), issues);
                enumValue(coloc.get("ld_provider"), colocWhere + ".ld_provider", Set.of("1000g-30x-grch38"), issues);
                enumValue(coloc.get("sample_overlap"), colocWhere + ".sample_overlap", Set.of("none", "unknown"), issues);
                if (coloc.containsKey("allow_partial_regions") && !(coloc.get("allow_partial_regions") instanceof Boolean))
                    issues.add(issue(ERROR, colocWhere + ".allow_partial_regions", "must be true or false"));
                for (String field : List.of("max_related_phenotypes", "max_colocalization_traits", "max_exposure_variants", "region_window_bp"))
                    if (coloc.containsKey(field) && (!(coloc.get(field) instanceof Number value) || value.doubleValue() != value.intValue()
                            || value.intValue() < (field.equals("region_window_bp") ? 0 : 1)
                            || value.intValue() > (field.equals("region_window_bp") ? 5000000 : field.equals("max_exposure_variants") ? 5000 : 20)))
                        issues.add(issue(ERROR, colocWhere + "." + field, "is outside the supported bounded range"));
                if (coloc.containsKey("minimum_exposure_overlap") && (!(coloc.get("minimum_exposure_overlap") instanceof Number value)
                        || !(value.doubleValue() > 0 && value.doubleValue() <= 1)))
                    issues.add(issue(ERROR, colocWhere + ".minimum_exposure_overlap", "must be in (0,1]"));
                if (coloc.containsKey("related_phenotypes")) {
                    List<Object> terms = sequence(coloc.get("related_phenotypes"), colocWhere + ".related_phenotypes", true, issues);
                    if (terms.size() > 20) issues.add(issue(ERROR, colocWhere + ".related_phenotypes", "at most 20 related phenotypes are allowed"));
                    for (Object term : terms) requiredText(term, colocWhere + ".related_phenotypes", issues);
                }
            }
        }
        if (goals != null)
            for (Object goal : sequence(goals, path + ".analysis.goals", true, issues)) {
                requiredText(goal, path + ".analysis.goals", issues);
                if (goal instanceof String text && !GOALS.contains(
                        text.toLowerCase().replace('-', '_')))
                    issues.add(issue(ERROR, path + ".analysis.goals",
                        "unregistered analysis goal: " + text));
            }
        if (analysis.containsKey("harmonization")) {
            Map<String, Object> harmonization = mapping(
                analysis.get("harmonization"), path + ".analysis.harmonization",
                true, issues);
            unknown(harmonization, Set.of("invalid_rows", "common_universe"),
                path + ".analysis.harmonization", issues);
            enumValue(harmonization.get("invalid_rows"),
                path + ".analysis.harmonization.invalid_rows",
                Set.of("error", "exclude"), issues);
            enumValue(harmonization.get("common_universe"),
                path + ".analysis.harmonization.common_universe",
                Set.of("none", "intersection"), issues);
        }
        if (analysis.containsKey("replication")) {
            Map<String, Object> replication = mapping(
                analysis.get("replication"), path + ".analysis.replication",
                true, issues);
            unknown(replication, Set.of("adjustment_scope", "require_same_sign",
                "threshold_comparison"), path + ".analysis.replication", issues);
            enumValue(replication.get("adjustment_scope"),
                path + ".analysis.replication.adjustment_scope",
                Set.of("all_results", "discovery_subset"), issues);
            enumValue(replication.get("threshold_comparison"),
                path + ".analysis.replication.threshold_comparison",
                Set.of("less_than", "less_than_or_equal"), issues);
            Object sameSign = replication.get("require_same_sign");
            if (sameSign != null && !(sameSign instanceof Boolean))
                issues.add(issue(ERROR,
                    path + ".analysis.replication.require_same_sign",
                    "must be true or false"));
        }
    }

    private static Map<String, Map<String, Object>> objects(Object value,
            String path, Set<String> allowed, boolean required,
            List<StudyValidation.Issue> issues) {
        List<Map<String, Object>> list = objectList(value, path, allowed,
            required, issues);
        Map<String, Map<String, Object>> result = new LinkedHashMap<>();
        for (int index = 0; index < list.size(); index++) {
            Map<String, Object> item = list.get(index);
            String itemPath = path + "[" + index + "]";
            Object id = item.get("id");
            identifier(id, itemPath + ".id", issues);
            if (id instanceof String text && result.putIfAbsent(text, item) != null)
                issues.add(issue(ERROR, itemPath + ".id",
                    "duplicate identifier: " + text));
        }
        return result;
    }

    private static List<Map<String, Object>> objectList(Object value,
            String path, Set<String> allowed, boolean required,
            List<StudyValidation.Issue> issues) {
        List<Object> values = sequence(value, path, required, issues);
        List<Map<String, Object>> result = new ArrayList<>();
        for (int index = 0; index < values.size(); index++) {
            Map<String, Object> item = mapping(values.get(index),
                path + "[" + index + "]", true, issues);
            unknown(item, allowed, path + "[" + index + "]", issues);
            result.add(item);
        }
        return result;
    }

    private static Map<String, Object> mapping(Object value, String path,
            boolean required, List<StudyValidation.Issue> issues) {
        if (value == null) {
            if (required) issues.add(issue(ERROR, path, "mapping is required"));
            return new LinkedHashMap<>();
        }
        try {
            return StudyIO.map(value, path);
        } catch (IllegalArgumentException failure) {
            issues.add(issue(ERROR, path, failure.getMessage()));
            return new LinkedHashMap<>();
        }
    }

    private static List<Object> sequence(Object value, String path,
            boolean required, List<StudyValidation.Issue> issues) {
        if (value == null) {
            if (required) issues.add(issue(ERROR, path, "list is required"));
            return List.of();
        }
        try {
            return StudyIO.list(value, path);
        } catch (IllegalArgumentException failure) {
            issues.add(issue(ERROR, path, failure.getMessage()));
            return List.of();
        }
    }

    private static void unknown(Map<String, Object> value, Set<String> allowed,
            String path, List<StudyValidation.Issue> issues) {
        for (String key : value.keySet())
            if (!allowed.contains(key))
                issues.add(issue(ERROR, path + "." + key, "unknown field"));
    }

    private static void identifier(Object value, String path,
            List<StudyValidation.Issue> issues) {
        if (!(value instanceof String text) || !ID.matcher(text).matches())
            issues.add(issue(ERROR, path,
                "identifier must match " + ID.pattern()));
    }

    private static void requiredText(Object value, String path,
            List<StudyValidation.Issue> issues) {
        if (!(value instanceof String text) || text.isBlank())
            issues.add(issue(ERROR, path, "nonblank text is required"));
    }

    private static void scalarOrList(Object value, String path,
            List<StudyValidation.Issue> issues) {
        if (value == null) return;
        if (value instanceof String text && !text.isBlank()) return;
        if (value instanceof Collection<?> items && !items.isEmpty()
                && items.stream().allMatch(item ->
                    item instanceof String text && !text.isBlank())) return;
        issues.add(issue(ERROR, path,
            "must be a nonblank string or nonempty list of strings"));
    }

    private static void enumValue(Object value, String path,
            Set<String> allowed, List<StudyValidation.Issue> issues) {
        if (value == null) return;
        if (!(value instanceof String text)
                || !allowed.contains(text.toLowerCase()))
            issues.add(issue(ERROR, path,
                "must be one of " + allowed));
    }

    private static String text(Object value) {
        return value instanceof String text && !text.isBlank() ? text : null;
    }

    private static StudyValidation.Issue issue(
            StudyValidation.Severity severity, String path, String message) {
        return new StudyValidation.Issue(severity, path, message);
    }
}
