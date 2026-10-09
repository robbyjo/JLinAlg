/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Typed capabilities available to the deterministic study planner. */
public final class StudyOperationRegistry {
    private StudyOperationRegistry() { }

    public record Descriptor(String type, String executor,
            Set<String> acceptedStages, Set<String> modalities,
            List<String> outputs, List<String> diagnostics) {
        public Descriptor {
            acceptedStages = Set.copyOf(acceptedStages);
            modalities = Set.copyOf(modalities);
            outputs = List.copyOf(outputs);
            diagnostics = List.copyOf(diagnostics);
        }
    }

    private static final Set<String> ALL_STAGES = Set.of(
        "raw_quantified", "quantified_matrix", "processed_matrix",
        "association_results", "annotated_results", "meta_analysis_results",
        "enrichment_results", "variant_scores", "fine_mapping_results",
        "colocalization_results", "mr_results", "evidence_results");
    private static final Set<String> ALL_MODALITIES = Set.of(
        "generic", "transcriptomics", "expression", "mirna", "lincrna",
        "lncrna", "extracellular_rna", "epigenomics", "ewas", "gwas",
        "genomics", "proteomics", "metabolomics", "single_cell_rna",
        "single_cell_protein", "spatial_transcriptomics",
        "spatial_proteomics", "microbiome", "multimodal", "phenotype");
    private static final Set<String> BULK_ASSOCIATION_MODALITIES = Set.of(
        "generic", "transcriptomics", "expression", "mirna", "lincrna",
        "lncrna", "extracellular_rna", "epigenomics", "ewas", "gwas",
        "genomics", "proteomics", "metabolomics", "phenotype");
    private static final Map<String, Descriptor> OPERATIONS = operations();

    public static Descriptor descriptor(String type) {
        Descriptor result = OPERATIONS.get(type);
        if (result == null)
            throw new IllegalArgumentException("Unregistered study operation: " + type);
        return result;
    }

    public static Map<String, Descriptor> descriptors() {
        return Map.copyOf(OPERATIONS);
    }

    private static Map<String, Descriptor> operations() {
        Map<String, Descriptor> result = new LinkedHashMap<>();
        add(result, "validate_harmonize", "jlinalg", ALL_STAGES, ALL_MODALITIES,
            List.of("harmonized inputs", "input manifest"),
            List.of("identity alignment", "missing columns", "input hashes"));
        add(result, "association", "jlinalg",
            Set.of("raw_quantified", "quantified_matrix", "processed_matrix"),
            BULK_ASSOCIATION_MODALITIES, List.of("association results"),
            List.of("sample alignment", "model convergence", "tested family"));
        add(result, "effect_comparison", "jlinalg",
            Set.of("association_results", "annotated_results"),
            ALL_MODALITIES, List.of("pairwise effect scatterplots", "effect correlations"),
            List.of("Pearson correlation", "effect-scale agreement", "matched features"));
        add(result, "discovery_replication", "jlinalg",
            Set.of("association_results", "annotated_results"),
            ALL_MODALITIES, List.of("replication comparison"),
            List.of("sign concordance", "replication counts"));
        add(result, "meta_analysis", "jlinalg",
            Set.of("association_results", "annotated_results"),
            ALL_MODALITIES, List.of("pooled effects"),
            List.of("heterogeneity", "cohort directions"));
        add(result, "gene_ontology", "jlinalg",
            Set.of("association_results", "annotated_results",
                "meta_analysis_results"),
            ALL_MODALITIES, List.of("ontology enrichment"),
            List.of("tested background", "set overlap"));
        add(result, "pathway_enrichment", "jlinalg",
            Set.of("association_results", "annotated_results",
                "meta_analysis_results"),
            ALL_MODALITIES, List.of("pathway enrichment"),
            List.of("tested background", "set overlap"));
        add(result, "variant_scoring", "jlinalg",
            Set.of("association_results", "annotated_results"),
            Set.of("gwas", "genomics"), List.of("variant scores"),
            List.of("allele alignment", "annotation coverage"));
        add(result, "fine_mapping", "jlinalg",
            Set.of("association_results", "meta_analysis_results"),
            Set.of("gwas", "genomics"), List.of("credible sets", "PIPs"),
            List.of("LD alignment", "credible-set coverage"));
        add(result, "colocalization", "jlinalg",
            Set.of("association_results", "fine_mapping_results"),
            ALL_MODALITIES, List.of("colocalization posteriors"),
            List.of("allele alignment", "overlap", "prior sensitivity"));
        add(result, "mendelian_randomization", "jlinalg",
            Set.of("association_results", "meta_analysis_results"),
            ALL_MODALITIES, List.of("causal estimates"),
            List.of("instrument strength", "heterogeneity", "pleiotropy"));
        add(result, "literature_support", "evidence_request",
            ALL_STAGES, ALL_MODALITIES, List.of("literature evidence request"),
            List.of("query provenance", "manual or service disposition"));
        add(result, "candidate_prioritization", "jlinalg",
            ALL_STAGES, ALL_MODALITIES,
            List.of("transparent candidate evidence table"),
            List.of("evidence availability", "ranking limitations"));
        add(result, "diagnostic_review", "jlinalg",
            ALL_STAGES, ALL_MODALITIES, List.of("diagnostic observations"),
            List.of("unusual patterns", "technical alternatives"));
        add(result, "report", "jlinalg", ALL_STAGES, ALL_MODALITIES,
            List.of("traceable report"), List.of("result lineage completeness"));
        add(result, "participant_characteristics", "jlinalg", ALL_STAGES, ALL_MODALITIES,
            List.of("local aggregate participant characteristics"), List.of("analysis population", "missingness", "sex coding", "visit alignment"));
        add(result, "irqtl_candidate_review", "jlinalg", ALL_STAGES, ALL_MODALITIES,
            List.of("local irQTL candidate review"), List.of("exact transcript matching", "regional completeness", "instrument review flags"));
        for (String type : List.of("advisor_review", "related_phenotype_advice", "manuscript_draft", "sensitivity_analysis"))
            add(result, type, "jlinalg", ALL_STAGES, ALL_MODALITIES,
                List.of("decision record", "audited outputs"), List.of("permissions", "prerequisites", "lineage"));
        for (String type : List.of("source_discovery", "bioinformatics_investigation", "literature_search"))
            add(result, type, "jlinalg", ALL_STAGES, ALL_MODALITIES,
                List.of("cited public evidence", "request and response audit"),
                List.of("scientific relevance", "disclosure policy", "search completeness", "source availability"));
        return result;
    }

    private static void add(Map<String, Descriptor> target, String type,
            String executor, Set<String> stages, Set<String> modalities,
            List<String> outputs, List<String> diagnostics) {
        target.put(type, new Descriptor(type, executor, stages, modalities,
            outputs, diagnostics));
    }
}
