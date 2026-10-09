/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyNeed;
import org.jlinalg.study.StudyOperationExecutor;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.jlinalg.study.StudyOperationExecutor.Result;
import org.jlinalg.study.StudyOperationExecutor.State;
import org.yaml.snakeyaml.Yaml;

/** Hypothesis-scoped public evidence, with explicit search limits and no causal inference. */
final class StudyBioinformatics implements StudyOperationExecutor {
    static final String GWAS = "https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/";
    static final String EQTL = "https://www.ebi.ac.uk/eqtl/Data_access/";
    static final String GTEX = "https://doi.org/10.1126/science.aaz1776";
    static final String UNIPROT = "https://www.uniprot.org/help/api_queries";
    static final String ALPHAFOLD = "https://doi.org/10.1093/nar/gkab1061";
    static final String UCSC = "https://genome.ucsc.edu/goldenPath/help/api.html";
    static final String LITERATURE = "https://europepmc.org/RestfulWebService";
    private final StudyPublicEvidence client;
    private final Set<String> excerptedArticles = new LinkedHashSet<>();
    StudyBioinformatics() { this(new StudyPublicEvidence()); }
    StudyBioinformatics(StudyPublicEvidence client) { this.client = client; }

    @Override public Result execute(Context context) throws IOException {
        Files.createDirectories(context.outputDirectory());
        List<Object> evidence = new ArrayList<>();
        List<StudyNeed> needs = new ArrayList<>();
        Map<String, Object> hypothesis = hypothesis(context);
        Map<String, Object> config = config(context);
        String operation = text(context.operation().get("type"));
        List<Feature> features = features(context);
        Map<String, String> termOrigins = phenotypeTerms(context, hypothesis, config,
            operation.equals("source_discovery"));
        List<String> terms = new ArrayList<>(termOrigins.keySet());
        if (operation.equals("source_discovery")) {
            String taxon = species(context);
            if (!taxon.equals("9606")) {
                evidence.add(Map.of("kind", "genetic_source_discovery", "status", "species_required", "reason",
                    "Built-in GWAS Catalog, eQTL Catalogue and GTEx discovery currently supports explicitly human studies only. Supply a species-matched source adapter for other organisms; human sources are not substituted."));
            } else {
            if (terms.isEmpty()) needs.add(review(context, "phenotype", "Which phenotype label or ontology term should the GWAS search use?",
                "Specify hypotheses[].outcome or analysis.bioinformatics.phenotype_terms. Abbreviations are not expanded silently."));
            String primary = text(hypothesis.get("outcome"));
            for (String term : terms) catalog(context, term,
                term.equalsIgnoreCase(primary) ? features : List.of(), evidence, termOrigins.get(term));
            qtl(context, evidence);
            }
            needs.add(review(context, "genetic-source-selection", "Which GWAS and molecular QTL sources are scientifically appropriate for this hypothesis?",
                "Inspect source-discovery.yaml. Automatic colocalization may rank indexed harmonized GWAS candidates under analysis.bioinformatics.colocalization; custom sources use genetic_pairs. Advisor labels and catalogue matches are search candidates, not automatically verified biological relationships or causal evidence."));
        } else if (operation.equals("bioinformatics_investigation")) {
            for (Feature feature : features) {
                protein(context, feature, evidence);
                sequence(context, feature, evidence);
            }
        } else if (operation.equals("literature_search")) {
            methodsLiterature(context, evidence);
            if (terms.isEmpty()) needs.add(review(context, "literature-phenotype", "Which phenotype should the literature search address?",
                "Specify a phenotype label; do not substitute the entire project description or patient records."));
            List<String> molecularTerms = literatureTerms(context);
            if (molecularTerms.isEmpty()) needs.add(review(context, "literature-field", "Which public molecular-field terms describe this study?",
                "Specify analysis.bioinformatics.literature_terms when the declared dataset type does not identify the field. Do not send the full project description or patient information."));
            if (molecularTerms.isEmpty() || terms.isEmpty()) for (String strategy : StudyIntroduction.SEARCH_STRATEGIES)
                for (String term : terms.isEmpty() ? List.of("unspecified") : terms) evidence.add(Map.of(
                    "search_strategy", strategy, "phenotype", term, "status", "search_inputs_required", "reason",
                    "Field-level search needs explicit public molecular-field and phenotype labels; no query was sent for this scope."));
            for (String term : terms) if (!molecularTerms.isEmpty()) introductionLiterature(context, molecularTerms, term, evidence);
            for (Feature feature : features) if (!feature.symbol().isBlank())
                for (String term : terms) literature(context, feature, term, evidence);
        } else throw new IllegalArgumentException("Unsupported bioinformatics operation: " + operation);
        if (features.isEmpty() && !operation.equals("source_discovery")) {
            evidence.add(Map.of("status", "no_eligible_features", "reason",
                "No unambiguously annotated feature passed the full-family BH threshold in an available primary/meta table. No arbitrary top hit is called significant."));
        }
        if (operation.equals("literature_search") && evidence.stream().map(StudyPublicEvidence::map).anyMatch(row -> row.containsKey("methods_role")
                && Set.of("search_inputs_required", "reference_identifier_required").contains(row.get("status"))))
            needs.add(review(context, "methods-references", "Can you supply public cohort names or valid method reference identifiers for Methods?",
                "Inspect cohort/method reference search records. Specify cohorts[].literature_name and DOI or pmid:/pmcid: references; do not substitute confidential descriptions or patient data."));
        boolean blocked = evidence.stream().map(StudyPublicEvidence::map).anyMatch(row ->
            Set.of("annotation_required", "build_required", "species_required", "interval_too_large").contains(row.get("status")));
        if (blocked) needs.add(review(context, "annotation", "Can you resolve the annotation prerequisites in the bioinformatics evidence table?",
            "Supply explicit species, gene identity, unique genome build and 1-based feature interval. Protein counterparts may be absent for noncoding RNA; do not force a protein mapping."));
        boolean retryable = client.requests().stream().map(StudyPublicEvidence::map)
            .anyMatch(row -> !"retrieved".equals(row.get("status")));
        Map<String, Object> summary = new LinkedHashMap<>();
        summary.put("schema_version", 1);
        summary.put("hypothesis", context.operation().get("hypothesis"));
        summary.put("selection", "BH across every valid tested feature per primary table, or the full meta table; deterministic p-value ordering. Selection is exploratory and not independent validation.");
        summary.put("feature_limit", integer(config, "max_features", 10, 1, 100));
        summary.put("features_selected", features.stream().map(Feature::document).toList());
        summary.put("phenotype_search_terms", termOrigins.entrySet().stream()
            .map(entry -> Map.of("phenotype", entry.getKey(), "origin", entry.getValue())).toList());
        summary.put("evidence", evidence);
        summary.put("queries", client.requests());
        summary.put("retryable", retryable);
        summary.put("needs", needs.stream().map(StudyNeed::document).toList());
        summary.put("interpretation", "Search results and reference annotations are contextual evidence, not proof of causality, molecular mechanism, druggability or absence of an association. Unavailable/limited searches are not negative evidence.");
        String name = operation.replace('_', '-') + ".yaml";
        List<String> outputs = new ArrayList<>(List.of(name));
        if (operation.equals("literature_search")) {
            Map<String, Object> catalog = StudyReferenceCatalog.build(evidence);
            summary.put("reference_catalog", catalog);
            write(context.outputDirectory().resolve("reference-catalog.yaml"), catalog);
            Files.writeString(context.outputDirectory().resolve("reference-catalog.html"), StudyReferenceCatalog.html(catalog));
            outputs.addAll(List.of("reference-catalog.yaml", "reference-catalog.html"));
        }
        write(context.outputDirectory().resolve(name), summary);
        boolean waiting = retryable || blocked || !needs.isEmpty();
        return new Result(waiting ? State.WAITING_FOR_INPUT : State.COMPLETE,
            operation + ": " + evidence.size() + " evidence records; " + client.requests().size()
                + " audited queries" + (waiting ? "; inspect required inputs and service status" : ""),
            outputs, summary);
    }

    private void catalog(Context context, String term, List<Feature> features, List<Object> out,
            String origin) throws IOException {
        int pages = integer(config(context), "max_pages", 3, 1, 20);
        for (int page = 0; page < pages; page++) {
            String url = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies?efo_trait=" + encode(term)
                + "&show_child_trait=false&size=100&page=" + page;
            StudyPublicEvidence.Response response = client.get(context, "GWAS Catalog", url, GWAS,
                "Identify genome-wide association studies of " + term + " as candidate outcome sources for the declared hypothesis; review ancestry, sample overlap and regional completeness before inference.");
            if (!response.ok()) { out.add(record("gwas_studies", "", response, GWAS)); break; }
            List<?> studies = catalogItems(response.data(), "studies");
            if (studies.isEmpty() && page == 0) out.add(Map.of("kind", "gwas_studies", "phenotype", term, "status", "no_matches", "audit", response.audit(), "citation", GWAS));
            for (Object value : studies) {
                Map<String, Object> study = map(value);
                Map<String, Object> row = record("gwas_study", "", response, GWAS);
                row.put("phenotype_query", term);
                row.put("phenotype_origin", origin);
                row.put("relevance", "Exact catalogue trait filter; phenotype definition, cohort overlap and ancestry still require review");
                row.put("metadata", study);
                if (Boolean.TRUE.equals(study.get("full_summary_stats_available"))) {
                    row.put("regional_source_status", "full_summary_statistics_candidate");
                    String source = text(study.get("full_summary_stats")).replaceFirst("^http://", "https://");
                    if (!source.isBlank()) row.put("full_summary_statistics", source);
                } else row.put("regional_source_status", "curated_associations_only_incomplete_for_conventional_colocalization");
                if (text(study.get("pubmed_id")).matches("[1-9][0-9]*")) row.put("publication", "https://pubmed.ncbi.nlm.nih.gov/" + text(study.get("pubmed_id")) + "/");
                out.add(row);
            }
            if (!map(map(response.data()).get("_links")).containsKey("next")) break;
            if (page == pages - 1) out.add(Map.of("kind", "gwas_studies", "phenotype", term, "status", "search_limited", "reason", "Page budget reached; not an exhaustive search", "citation", GWAS));
        }
        // Gene-filtered hits are candidates; proximity is evaluated only against explicitly declared catalogue build.
        Set<String> genes = new LinkedHashSet<>();
        for (Feature feature : features) if (safeSymbol(feature.symbol()) && genes.add(feature.symbol())) {
            String url = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations?efo_trait=" + encode(term)
                + "&show_child_trait=false&mapped_gene=" + encode(feature.symbol()) + "&size=100";
            StudyPublicEvidence.Response response = client.get(context, "GWAS Catalog associations", url, GWAS,
                "Check whether " + term + " associations map to " + feature.symbol() + ", a selected omics finding. Gene proximity is contextual evidence, not evidence of an isoform mechanism.");
            Map<String, Object> row = record("gwas_gene_associations", feature.id(), response, GWAS);
            if (response.ok()) {
                List<?> hits = catalogItems(response.data(), "associations");
                row.put("associations", hits);
                row.put("status", list(hits).isEmpty() ? "no_matches" : "candidate_matches");
                row.put("search_limited", map(map(response.data()).get("_links")).containsKey("next"));
                row.put("relevance", "Same catalogue phenotype and mapped gene; mapping includes nearest genes, so this is not a coordinate-window or causal match. Full regional GWAS/QTL and build checks remain required.");
            }
            out.add(row);
        }
        nearby(context, term, features, out);
    }

    private static Map<String, String> phenotypeTerms(Context context,
            Map<String, Object> hypothesis, Map<String, Object> config,
            boolean includeColocalizationAdvice) throws IOException {
        Map<String, String> result = new LinkedHashMap<>();
        String primary = text(hypothesis.get("outcome"));
        if (!primary.isBlank()) addTerm(result, label(primary), "hypothesis_outcome");
        for (Object term : list(config.get("phenotype_terms")))
            addTerm(result, label(term), "scientist_declared_phenotype_term");
        if (!includeColocalizationAdvice) return result;
        Map<String, Object> coloc = map(config.get("colocalization"));
        for (Object term : list(coloc.get("related_phenotypes")))
            addTerm(result, label(term), "scientist_declared_related_phenotype");
        int maximum = integer(coloc, "max_related_phenotypes", 20, 1, 20);
        int added = 0;
        for (Path directory : context.dependencyDirectories().values()) {
            Path decision = directory.resolve("decision.yaml");
            if (!Files.isRegularFile(decision)) continue;
            Map<String, Object> document;
            try (var reader = Files.newBufferedReader(decision)) { document = map(new Yaml().load(reader)); }
            if (!"related_phenotypes".equals(document.get("purpose"))
                    || !"completed_unreviewed".equals(document.get("status"))) continue;
            for (Object raw : list(document.get("phenotypes"))) {
                if (added >= maximum) break;
                Map<String, Object> item = map(raw);
                String term = label(item.get("name"));
                if (addTerm(result, term, "llm_advisor_proposal_unreviewed")) added++;
            }
        }
        return result;
    }

    private static boolean addTerm(Map<String, String> terms, String term,
            String origin) {
        for (String current : terms.keySet()) if (current.equalsIgnoreCase(term)) return false;
        terms.put(term, origin); return true;
    }

    private void nearby(Context context, String term, List<Feature> features, List<Object> out) throws IOException {
        if (features.isEmpty()) return;
        StudyPublicEvidence.Response metadata = client.get(context, "GWAS Catalog release",
            "https://www.ebi.ac.uk/gwas/rest/api/v2/metadata", GWAS,
            "Verify the catalogue genome build before interpreting proximity to the study's annotated molecular features; prevent cross-build coordinate matches.");
        String catalogBuild = text(map(metadata.data()).get("gene_build")).split("\\.")[0];
        for (Feature feature : features) {
            if (!feature.status().equals("located") || !metadata.ok() || !feature.build().equals(catalogBuild)) {
                Map<String, Object> row = new LinkedHashMap<>();
                row.put("kind", "nearby_gwas"); row.put("feature_id", feature.id());
                row.put("status", metadata.ok() ? "build_required" : metadata.status());
                row.put("reason", "Proximity needs a unique feature interval matching the live catalogue annotation build; no implicit liftover");
                row.put("catalogue_metadata_audit", metadata.audit()); row.put("catalogue_annotation_build", catalogBuild); out.add(row); continue;
            }
            int window = integer(config(context), "window_bp", 500000, 0, 5000000);
            long start = Math.max(1, Long.parseLong(feature.start()) - window), end = Long.parseLong(feature.end()) + window;
            int limit = integer(config(context), "max_regional_snps", 20, 1, 100);
            String chromosome = feature.chromosome().replaceFirst("^chr", "");
            StudyPublicEvidence.Response variants = client.get(context, "GWAS Catalog regional SNPs",
                "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms?chromosome=" + encode(chromosome)
                    + "&bp_start=" + start + "&bp_end=" + end + "&size=" + limit, GWAS,
                "Find catalogue variants within the declared " + window + " bp window of " + feature.id() + " on " + catalogBuild + ", then test their relevance to " + term + ". This is candidate discovery, not colocalization.");
            Map<String, Object> row = record("nearby_gwas", feature.id(), variants, GWAS);
            row.put("interval", Map.of("chromosome", chromosome, "start", start, "end", end, "genome_build", catalogBuild));
            row.put("catalogue_metadata_audit", metadata.audit());
            if (variants.ok()) {
                List<?> snps = catalogItems(variants.data(), "snps");
                List<Object> associations = new ArrayList<>();
                boolean limited = map(map(variants.data()).get("_links")).containsKey("next"), unavailable = false;
                for (Object raw : list(snps)) {
                    Map<String, Object> snp = map(raw);
                    boolean inside = list(snp.get("locations")).stream().map(StudyPublicEvidence::map).anyMatch(location ->
                        chromosome.equals(text(location.get("chromosome_name"))) && location.get("chromosome_position") instanceof Number pos && pos.longValue() >= start && pos.longValue() <= end);
                    if (!inside) throw new IOException("GWAS regional result did not satisfy requested coordinates");
                    String id = text(snp.get("rs_id"));
                    if (!id.matches("rs[0-9]+")) continue;
                    StudyPublicEvidence.Response trait = client.get(context, "GWAS Catalog phenotype SNP",
                        "https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id=" + id + "&efo_trait=" + encode(term) + "&show_child_trait=false&size=100", GWAS,
                        "Determine whether nearby variant " + id + " has a reported association with " + term + "; do not count proximity alone as phenotype-specific support.");
                    if (!trait.ok()) { unavailable = true; associations.add(record("regional_snp_query", feature.id(), trait, GWAS)); continue; }
                    List<?> hits = catalogItems(trait.data(), "associations");
                    limited |= map(map(trait.data()).get("_links")).containsKey("next");
                    for (Object hit : list(hits)) associations.add(Map.of("snp", snp, "association", hit, "audit", trait.audit()));
                }
                row.put("phenotype", term); row.put("associations", associations); row.put("search_limited", limited);
                row.put("status", unavailable ? "partial_unavailable" : limited ? "search_limited" : associations.isEmpty() ? "no_matches" : "candidate_matches");
                row.put("interpretation", "Catalogue lead SNPs near the feature and matching the declared phenotype. Lead SNP proximity is not QTL evidence; use full regional statistics for fine-mapping/MR/colocalization.");
            }
            out.add(row);
        }
    }

    static List<?> catalogItems(Object data, String collection) throws IOException {
        Map<String, Object> root = map(data), embedded = map(root.get("_embedded"));
        if (embedded.get(collection) instanceof List<?> items) return items;
        Map<String, Object> page = map(root.get("page"));
        // Spring HAL omits _embedded for an empty collection. Explicit zero totals,
        // not the absence of a field alone, establish a valid empty search.
        if (embedded.isEmpty() && page.get("totalElements") instanceof Number elements && elements.doubleValue() == 0
                && page.get("totalPages") instanceof Number pages && pages.doubleValue() == 0
                && page.get("number") instanceof Number number && number.doubleValue() == 0) return List.of();
        throw new IOException("GWAS " + collection + " response schema changed; cannot interpret as no matches");
    }

    private void qtl(Context context, List<Object> out) throws IOException {
        Map<String, Object> row = new LinkedHashMap<>();
        row.put("kind", "qtl_catalogue"); row.put("status", "endpoint_retired");
        row.put("citation", EQTL); row.put("request_attempted", false);
        row.put("reason", "The provider's Data Access documentation states that its REST API is deprecated and unavailable (verified 2026-10-02); this is not evidence of no scientific matches. Use provider FTP/tabix files with a reviewed source manifest; this adapter does not yet acquire them.");
        row.put("fallback", "https://www.ebi.ac.uk/eqtl/Data_access/");
        row.put("relevance", "Dataset discovery only. Assay, exact molecular trait and tissue must match the proposed question; gene eQTL is not a substitute for a measured splicing QTL.");
        out.add(row);
        StudyPublicEvidence.Response gtex = client.get(context, "GTEx",
            "https://gtexportal.org/api/v2/dataset/tissueSiteDetail?datasetId=gtex_v8", GTEX,
            "Check availability of the study's declared or explicitly justified follow-up tissues in GTEx v8; tissue availability is not feature-specific genetic evidence.");
        Map<String, Object> tissues = record("gtex_tissue_candidates", "", gtex, GTEX);
        if (gtex.ok()) {
            Set<String> relevant = new LinkedHashSet<>();
            for (Object item : list(config(context).get("relevant_tissues"))) relevant.add(text(item).toLowerCase(java.util.Locale.ROOT));
            for (Map<String, Object> dataset : datasets(context))
                if (!text(dataset.get("tissue")).isBlank()) relevant.add(text(dataset.get("tissue")).toLowerCase(java.util.Locale.ROOT));
            tissues.put("requested_tissues", relevant.stream().toList());
            if (!(map(gtex.data()).get("data") instanceof List<?>)) throw new IOException("GTEx tissue response schema changed");
            tissues.put("tissues", list(map(gtex.data()).get("data")).stream().map(StudyPublicEvidence::map)
                .filter(item -> relevant.stream().anyMatch(tissue -> (text(item.get("tissueSiteDetail")) + " " + text(item.get("tissueSiteDetailId")))
                    .toLowerCase(java.util.Locale.ROOT).replace('_', ' ').contains(tissue.replace('_', ' '))))
                .toList());
            tissues.put("relevance", "Declared study or follow-up tissues only; additional tissues require a documented hypothesis. Availability counts are not feature-specific QTL results.");
        }
        out.add(tissues);
    }

    private void protein(Context context, Feature feature, List<Object> out) throws IOException {
        String taxon = species(context);
        if (taxon.isEmpty()) { out.add(Map.of("feature_id", feature.id(), "status", "species_required", "kind", "protein")); return; }
        if (!safeSymbol(feature.symbol())) { out.add(Map.of("feature_id", feature.id(), "status", "annotation_required", "kind", "protein")); return; }
        String url = "https://rest.uniprot.org/uniprotkb/search?query=" + encode("gene_exact:" + feature.symbol() + " AND organism_id:" + taxon + " AND reviewed:true") + "&format=json&size=10";
        StudyPublicEvidence.Response response = client.get(context, "UniProt", url, UNIPROT,
            "Identify reviewed protein counterparts of the selected gene " + feature.symbol() + " for functional context; gene-level mapping does not establish isoform equivalence.");
        Map<String, Object> row = record("protein_mapping", feature.id(), response, UNIPROT);
        if (!response.ok()) { out.add(row); return; }
        if (!(map(response.data()).get("results") instanceof List<?>)) throw new IOException("UniProt response schema changed");
        List<?> proteins = list(map(response.data()).get("results"));
        row.put("status", proteins.isEmpty() ? "no_reviewed_protein_matches" : proteins.size() == 1 ? "gene_level_match" : "ambiguous_protein_mapping");
        row.put("proteins", proteins.stream().map(StudyPublicEvidence::map).map(item -> select(item,
            List.of("primaryAccession", "entryType", "entryAudit", "organism", "genes", "proteinDescription"))).toList());
        row.put("interpretation", "Reviewed gene-level protein candidates, not a proven counterpart of the measured transcript or splicing isoform; unreviewed proteins were not searched.");
        out.add(row);
        for (Object protein : proteins) {
            String accession = text(map(protein).get("primaryAccession"));
            if (!accession.matches("[A-Z0-9]{6,10}")) continue;
            StudyPublicEvidence.Response structure = client.get(context, "AlphaFold DB", "https://alphafold.ebi.ac.uk/api/prediction/" + accession, ALPHAFOLD,
                "Assess availability and confidence of predicted structures for protein " + accession + " linked to " + feature.symbol() + "; structure availability alone is not drug-target validation.");
            Map<String, Object> found = record("predicted_structure", feature.id(), structure, ALPHAFOLD);
            found.put("uniprot", accession);
            if (structure.ok()) {
                if (!(structure.data() instanceof List<?>)) throw new IOException("AlphaFold response schema changed");
                found.put("models", list(structure.data()).stream().map(StudyPublicEvidence::map).map(item -> select(item,
                    List.of("entryId", "uniprotAccession", "latestVersion", "modelCreatedDate", "globalMetricValue", "fractionPlddtVeryLow", "fractionPlddtLow", "sequenceStart", "sequenceEnd", "cifUrl", "pdbUrl", "paeDocUrl"))).toList());
                found.put("status", list(structure.data()).isEmpty() ? "no_models" : "models_available");
            }
            found.put("interpretation", "Predicted structures and confidence metrics; no docking, functional validation, isoform equivalence or drug efficacy inferred");
            out.add(found);
        }
    }

    private void sequence(Context context, Feature feature, List<Object> out) throws IOException {
        if (!feature.status().equals("located")) { out.add(Map.of("feature_id", feature.id(), "kind", "reference_sequence", "status", "annotation_required")); return; }
        String genome = switch (feature.build()) { case "GRCh38", "hg38" -> "hg38"; case "GRCh37", "hg19" -> "hg19"; default -> ""; };
        if (genome.isEmpty()) { out.add(Map.of("feature_id", feature.id(), "kind", "reference_sequence", "status", "build_required")); return; }
        long start = Long.parseLong(feature.start()), end = Long.parseLong(feature.end());
        if (start < 1 || end < start) throw new IOException("Invalid annotated interval");
        if (end - start + 1 > integer(config(context), "max_sequence_bp", 100000, 1, 1000000)) {
            out.add(Map.of("feature_id", feature.id(), "kind", "reference_sequence", "status", "interval_too_large", "reason", "Explicit sequence budget exceeded; not silently truncated")); return;
        }
        String chromosome = feature.chromosome().startsWith("chr") ? feature.chromosome() : "chr" + feature.chromosome();
        String url = "https://api.genome.ucsc.edu/getData/sequence?genome=" + genome + ";chrom=" + encode(chromosome) + ";start=" + (start - 1) + ";end=" + end;
        StudyPublicEvidence.Response response = client.get(context, "UCSC reference sequence", url, UCSC,
            "Retrieve the build-matched genomic reference interval for " + feature.id() + " to support sequence inspection. This is reference DNA, not participant DNA or a spliced transcript.");
        Map<String, Object> row = record("reference_sequence", feature.id(), response, UCSC);
        row.put("browser", "https://genome.ucsc.edu/cgi-bin/hgTracks?db=" + genome + "&position=" + encode(chromosome + ":" + start + "-" + end));
        if (response.ok()) {
            String dna = text(map(response.data()).get("dna"));
            if (dna.length() != end - start + 1 || !dna.matches("[ACGTNacgtn]+")) throw new IOException("UCSC sequence length/content did not match the requested interval");
            String filename = "sequence-" + Integer.toHexString(url.hashCode()) + ".fasta";
            Files.writeString(context.outputDirectory().resolve(filename), ">" + genome + " " + chromosome + ":" + start + "-" + end + " forward_reference\n" + dna + "\n");
            row.put("sequence_file", filename);
            row.put("sequence_sha256", FollowupSupport.hash(context.outputDirectory().resolve(filename)));
            row.put("length", dna.length());
        }
        row.put("interpretation", "Forward-strand genomic reference sequence, not participant DNA, spliced transcript sequence or functional evidence");
        out.add(row);
    }

    private void literature(Context context, Feature feature, String term, List<Object> out) throws IOException {
        if (!safeSymbol(feature.symbol())) { out.add(Map.of("feature_id", feature.id(), "status", "annotation_required", "kind", "literature")); return; }
        String gene = "TITLE_ABS:\"" + feature.symbol() + "\"";
        String phenotype = "TITLE_ABS:\"" + term.replace("\"", "") + "\"";
        Map<String, String> strategies = new LinkedHashMap<>();
        strategies.put("direct_gene_phenotype", gene + " AND " + phenotype);
        strategies.put("molecular_mechanism", gene + " AND (TITLE_ABS:splicing OR TITLE_ABS:isoform OR TITLE_ABS:perturbation OR TITLE_ABS:knockout)");
        strategies.put("contradictions_and_replication", gene + " AND " + phenotype
            + " AND (TITLE_ABS:replication OR TITLE_ABS:null OR TITLE_ABS:negative OR TITLE_ABS:conflicting)");
        literatureQueries(context, feature.id(), feature.symbol(), term, strategies, out);
    }

    private static List<String> literatureTerms(Context context) {
        Set<String> terms = new LinkedHashSet<>();
        for (Object value : list(config(context).get("literature_terms"))) terms.add(label(value));
        if (!terms.isEmpty()) return List.copyOf(terms);
        for (Map<String, Object> dataset : datasets(context)) switch (text(dataset.get("type"))) {
            case "transcriptomics" -> terms.addAll(List.of("transcriptomics", "RNA sequencing", "RNA splicing", "transcript isoform"));
            case "epigenomics", "methylation" -> terms.addAll(List.of("DNA methylation", "epigenomics"));
            case "gwas", "genomics" -> terms.addAll(List.of("genome wide association", "genetic association"));
            case "proteomics" -> terms.addAll(List.of("proteomics", "protein abundance"));
            case "metabolomics" -> terms.add("metabolomics");
            default -> { /* Explicit terms are required for unknown assay types. */ }
        }
        return List.copyOf(terms);
    }
    private void introductionLiterature(Context context, List<String> molecularTerms, String term, List<Object> out) throws IOException {
        String field = "(" + String.join(" OR ", molecularTerms.stream().map(s -> "TITLE_ABS:\"" + s.replace("\"", "").replace("\\", "") + "\"").toList()) + ")";
        String base = field + " AND TITLE_ABS:\"" + term.replace("\"", "").replace("\\", "") + "\"";
        Map<String, String> strategies = new LinkedHashMap<>();
        strategies.put("field_background", base);
        strategies.put("hypothesis_relevance", base + " AND (TITLE_ABS:mechanism OR TITLE_ABS:regulation OR TITLE_ABS:biomarker)");
        strategies.put("similar_studies", base + " AND (TITLE_ABS:association OR TITLE_ABS:cohort OR TITLE_ABS:replication OR TITLE_ABS:longitudinal)");
        strategies.put("field_conflicting_evidence", base + " AND (TITLE_ABS:null OR TITLE_ABS:negative OR TITLE_ABS:inconsistent OR TITLE_ABS:replication)");
        literatureQueries(context, "", "", term, strategies, out);
    }
    private void literatureQueries(Context context, String featureId, String symbol, String term,
            Map<String, String> strategies, List<Object> out) throws IOException {
        literatureQueries(context, featureId, symbol, term, strategies, out, Map.of());
    }
    private void literatureQueries(Context context, String featureId, String symbol, String term,
            Map<String, String> strategies, List<Object> out, Map<String, Object> scope) throws IOException {
        for (Map.Entry<String, String> strategy : strategies.entrySet()) {
            String cursor = "*";
            int pages = integer(config(context), "max_pages", 3, 1, 20);
            long seen = 0;
            for (int page = 0; page < pages; page++) {
                String rationale = symbol.isBlank()
                    ? "Review the molecular field in relation to " + term + " for the introduction: " + strategy.getKey()
                        + ". Search terms come from explicit public literature_terms or declared dataset types, independently of significant genes. Compare existing studies and contradictory evidence; a bounded search cannot establish novelty."
                    : "Investigate " + symbol + " in relation to " + term + " using " + strategy.getKey()
                        + ". Retrieve abstracts for source-grounded review, including contradictory evidence; keyword matches are not confirmation of the hypothesis.";
                if (scope.containsKey("scientific_rationale")) rationale = text(scope.get("scientific_rationale"));
                StudyPublicEvidence.Response response = client.get(context, "Europe PMC",
                    "https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&pageSize=20&cursorMark="
                        + encode(cursor) + "&query=" + encode(strategy.getValue()), LITERATURE,
                    rationale);
                Map<String, Object> row = record("literature", featureId, response, LITERATURE);
                row.put("query", strategy.getValue()); row.put("search_strategy", strategy.getKey());
                row.put("phenotype", scope.isEmpty() ? term : ""); row.put("gene_symbol", symbol); row.put("page", page + 1);
                row.putAll(scope);
                row.put("scientific_rationale", rationale);
                row.put("retrieval", select(map(client.requests().get(client.requests().size() - 1)), List.of("started_at", "finished_at", "url", "response_sha256", "cached_response_audit")));
                if (!response.ok()) { out.add(row); break; }
                Map<String, Object> payload = map(response.data());
                if (!(payload.get("hitCount") instanceof Number)) throw new IOException("Europe PMC response schema changed");
                List<?> articles = list(map(payload.get("resultList")).get("result"));
                long hits = ((Number) payload.get("hitCount")).longValue();
                if (hits > 0 && articles.isEmpty()) throw new IOException("Europe PMC returned positive hitCount without article records");
                seen += articles.size();
                String next = text(payload.get("nextCursorMark"));
                boolean more = seen < hits;
                row.put("status", hits == 0 ? "no_matches" : "candidate_matches");
                row.put("hit_count", hits); row.put("search_limited", more && (page == pages - 1 || next.isBlank() || next.equals(cursor)));
                row.put("articles", articles.stream().map(StudyPublicEvidence::map)
                    .map(article -> uniqueExcerpt(StudyLiterature.abstractEvidence(article, symbol, term, response.audit()))).toList());
                row.put("interpretation", "Bounded multi-strategy search, not a systematic review. Excerpts are literal abstract text selected by lexical relevance, not LLM-validated support. Review study design, tissue, species, direction, retractions and conflicting results before ranking a story.");
                out.add(row);
                if (!more || next.isBlank() || next.equals(cursor)) break;
                cursor = next;
            }
        }
    }

    private void methodsLiterature(Context context, List<Object> out) throws IOException {
        var scoped = new LinkedHashMap<>(context.study()); scoped.put("datasets", datasets(context));
        for (var cohort : StudyMethods.cohorts(scoped)) {
            String name = StudyMethods.literatureName(cohort);
            Map<String, Object> scope = Map.of("cohort_id", cohort.get("id"), "cohort_name", name, "methods_role", "cohorts",
                "scientific_rationale", "Find the design/recruitment publication for this cohort or subcohort before describing it in Methods. Historical enrollment is not the current analysis sample; returned candidates need identity and source review.");
            if (!name.isBlank()) {
                String nameQuery = "TITLE_ABS:\"" + label(name).replace("\"", "").replace("\\", "") + "\"";
                if (name.equals("Framingham Heart Study Third Generation")) nameQuery = "TITLE_ABS:Framingham AND TITLE_ABS:\"Third Generation\"";
                literatureQueries(context, "", "", name, Map.of("cohort_description",
                    nameQuery + " AND (TITLE_ABS:design OR TITLE_ABS:recruitment OR TITLE_ABS:cohort)"), out, scope);
            }
            else { var row = new LinkedHashMap<>(scope); row.put("search_strategy", "cohort_description"); row.put("status", "search_inputs_required");
                row.put("reason", "Supply cohorts[].literature_name or an unambiguous public name; private IDs and descriptions are not sent as search queries."); out.add(row); }
            referenceQueries(context, list(cohort.get("references")), scope, out);
        }
        for (Map<String, Object> dataset : datasets(context)) referenceQueries(context, list(map(dataset.get("analysis_description")).get("references")),
            Map.of("methods_role", "association", "scientific_rationale", "Retrieve the scientist-declared method citation for supplied association results; publication identity is not verification of the actual upstream model."), out);
        for (String role : StudyMethods.ROLES) referenceQueries(context, list(map(map(hypothesis(context).get("analysis")).get(role)).get("references")),
            Map.of("methods_role", role, "scientific_rationale", "Retrieve the declared " + role + " method citation. This lookup does not establish execution or parameter settings."), out);
        for (Path directory : context.dependencyDirectories().values()) {
            Path file = directory.resolve("manifest.yaml"); if (!Files.isRegularFile(file)) continue;
            Map<String, Object> manifest; try (var reader = Files.newBufferedReader(file)) { manifest = map(new Yaml().load(reader)); }
            var report = map(map(manifest.get("details")).get("method_report"));
            if (report.isEmpty()) continue;
            var scope = Map.<String, Object>of("methods_role", StudyMethods.role(text(manifest.get("operation_type"))),
                "scientific_rationale", "Find a methods source corresponding to an executed adapter's recorded algorithm. Review estimator/version applicability before citation; query matches do not verify scientific support.");
            referenceQueries(context, list(report.get("references")), scope, out);
            var terms = list(report.get("literature_terms")).stream().map(t -> "TITLE_ABS:\"" + label(t).replace("\"", "").replace("\\", "") + "\"").toList();
            if (!terms.isEmpty()) literatureQueries(context, "", "", "", Map.of("analysis_method", String.join(" AND ", terms)), out, scope);
        }
    }
    private void referenceQueries(Context context, List<?> refs, Map<String, Object> scope, List<Object> out) throws IOException {
        for (Object id : refs) {
            String normalized = StudyReferenceCatalog.normalize(text(id)); String query;
            if (normalized.matches("pmid:[1-9][0-9]*")) query = "EXT_ID:" + normalized.substring(5) + " AND SRC:MED";
            else if (normalized.matches("pmcid:PMC[0-9]+")) query = "PMCID:" + normalized.substring(6);
            else if (normalized.matches("10\\.[0-9]{4,9}/[^\\s\"\\\\]+")) query = "DOI:\"" + normalized + "\"";
            else { var row = new LinkedHashMap<>(scope); row.put("search_strategy", "method_reference"); row.put("status", "reference_identifier_required"); row.put("reason", "Supply a DOI, pmid:NUMBER or pmcid:PMCNUMBER for an exact reference lookup."); out.add(row); continue; }
            var specific = new LinkedHashMap<>(scope); specific.put("requested_reference", normalized);
            literatureQueries(context, "", "", "", Map.of("method_reference", query), out, specific);
        }
    }

    private Map<String, Object> uniqueExcerpt(Map<String, Object> article) {
        String identity = !text(article.get("doi")).isBlank() ? text(article.get("doi")).toLowerCase(java.util.Locale.ROOT)
            : text(article.get("source")) + ":" + text(article.get("id"));
        if (!identity.equals(":") && article.containsKey("excerpt") && !excerptedArticles.add(identity)) {
            for (String key : List.of("excerpt", "excerpt_word_count", "excerpt_start", "excerpt_end")) article.remove(key);
            article.put("excerpt_status", "already_excerpted_in_this_operation");
            article.put("excerpt_reference", identity);
        }
        return article;
    }

    static Map<String, Object> hypothesis(Context context) {
        return list(context.study().get("hypotheses")).stream().map(StudyPublicEvidence::map)
            .filter(item -> item.get("id").equals(context.operation().get("hypothesis"))).findFirst().orElseThrow();
    }
    private static Map<String, Object> select(Map<String, Object> source, List<String> keys) {
        Map<String, Object> result = new LinkedHashMap<>();
        for (String key : keys) if (source.containsKey(key)) result.put(key, source.get(key));
        return result;
    }
    static Map<String, Object> config(Context context) { return map(map(hypothesis(context).get("analysis")).get("bioinformatics")); }
    static List<Map<String, Object>> datasets(Context context) {
        return list(context.study().get("datasets")).stream().map(StudyPublicEvidence::map)
            .filter(item -> list(context.operation().get("datasets")).contains(item.get("id"))).toList();
    }
    static String species(Context context) {
        Set<String> species = new LinkedHashSet<>();
        for (Map<String, Object> dataset : datasets(context)) {
            String value = text(dataset.get("species"));
            species.add(Set.of("human", "Homo sapiens", "9606").contains(value) ? "9606"
                : Set.of("mouse", "Mus musculus", "10090").contains(value) ? "10090" : "");
        }
        return species.size() == 1 ? species.iterator().next() : "";
    }
    static int integer(Map<String, Object> config, String key, int fallback, int min, int max) {
        Object value = config.getOrDefault(key, fallback);
        if (!(value instanceof Number number) || number.doubleValue() != number.intValue() || number.intValue() < min || number.intValue() > max)
            throw new IllegalArgumentException(key + " must be an integer in [" + min + "," + max + "]");
        return number.intValue();
    }
    private static String label(Object value) {
        String result = text(value);
        if (result.isBlank() || result.length() > 150 || result.chars().anyMatch(Character::isISOControl)) throw new IllegalArgumentException("Phenotype labels must be nonblank bounded text");
        return result;
    }
    private static boolean safeSymbol(String value) { return value.matches("[A-Za-z0-9][A-Za-z0-9._-]{0,63}"); }
    private static Map<String, Object> record(String kind, String feature, StudyPublicEvidence.Response response, String citation) {
        Map<String, Object> row = new LinkedHashMap<>();
        row.put("kind", kind); row.put("feature_id", feature); row.put("status", response.status());
        row.put("audit", response.audit()); row.put("citation", citation); row.put("reason", response.reason());
        return row;
    }
    private static StudyNeed review(Context context, String key, String prompt, String reason) {
        return StudyNeed.review("bioinformatics:" + context.operation().get("hypothesis") + ":" + key, prompt, reason, true);
    }

    record Feature(String id, String symbol, String gene, String chromosome, String start, String end, String build, String status, double p, double fdr) {
        Map<String, Object> document() { return Map.of("feature_id", id, "gene_symbol", symbol, "gene_id", gene, "chromosome", chromosome,
            "start", start, "end", end, "genome_build", build, "mapping_status", status, "p_value", p, "FDR", fdr); }
    }

    static List<Feature> features(Context context) throws IOException {
        Map<String, List<String>> annotations = new LinkedHashMap<>();
        Set<String> ambiguous = new LinkedHashSet<>();
        List<Path> primary = new ArrayList<>();
        Path meta = null;
        for (Path directory : context.dependencyDirectories().values()) {
            Path manifest = directory.resolve("manifest.yaml");
            if (!Files.isRegularFile(manifest)) continue;
            Map<String, Object> declared;
            try (var reader = Files.newBufferedReader(manifest)) { declared = map(new Yaml().load(reader)); }
            if (!"complete".equals(declared.get("state"))) continue;
            String type = text(declared.get("operation_type"));
            if (type.equals("meta_analysis") && Files.isRegularFile(directory.resolve("meta-analysis.tsv"))) meta = directory.resolve("meta-analysis.tsv");
            if (!Set.of("validate_harmonize", "association").contains(type)) continue;
            Path standard = directory.resolve("standard-results");
            if (Files.isDirectory(standard)) try (var files = Files.list(standard)) { primary.addAll(files.filter(Files::isRegularFile).sorted().toList()); }
            Path annotation = directory.resolve("feature-annotation.tsv");
            if (!Files.isRegularFile(annotation)) continue;
            DelimitedData table = DelimitedData.read(annotation);
            for (String[] row : table.rows()) {
                String id = row[table.column("feature_id")];
                List<String> fields = List.of("gene_symbol", "gene_id", "chromosome", "start", "end", "genome_build", "mapping_status").stream()
                    .map(key -> row[table.column(key)]).toList();
                List<String> old = annotations.putIfAbsent(id, fields);
                if (old != null && !old.equals(fields)) ambiguous.add(id);
            }
        }
        if (meta != null) primary = List.of(meta);
        Map<String, Feature> chosen = new LinkedHashMap<>();
        double threshold = ((Number) config(context).getOrDefault("fdr", 0.05)).doubleValue();
        if (!(threshold > 0 && threshold <= 1)) throw new IllegalArgumentException("bioinformatics.fdr must be in (0,1]");
        for (Path source : primary) {
            DelimitedData table = DelimitedData.read(source);
            Map<String, Double> probabilities = new LinkedHashMap<>();
            for (String[] row : table.rows()) {
                String id = row[table.column("feature_id")];
                if (table.header().contains("status") && !Set.of("ok", "single_cohort").contains(row[table.column("status")])) continue;
                try { double p = Double.parseDouble(row[table.column("p_value")]); if (Double.isFinite(p) && p >= 0 && p <= 1) probabilities.put(id, p); }
                catch (NumberFormatException ignored) { /* Invalid/not-estimable tests are not findings. */ }
            }
            List<Map.Entry<String, Double>> ranked = probabilities.entrySet().stream().sorted(Map.Entry.<String, Double>comparingByValue().thenComparing(Map.Entry::getKey)).toList();
            double adjusted = 1;
            for (int i = ranked.size() - 1; i >= 0; i--) {
                var entry = ranked.get(i);
                adjusted = Math.min(adjusted, entry.getValue() * ranked.size() / (i + 1));
                List<String> a = annotations.get(entry.getKey());
                if (adjusted > threshold || a == null || ambiguous.contains(entry.getKey()) || !Set.of("located", "missing_location").contains(a.get(6))) continue;
                Feature f = new Feature(entry.getKey(), a.get(0), a.get(1), a.get(2), a.get(3), a.get(4), a.get(5), a.get(6), entry.getValue(), adjusted);
                chosen.merge(f.id(), f, (left, right) -> left.p() <= right.p() ? left : right);
            }
        }
        return chosen.values().stream().sorted(Comparator.comparingDouble(Feature::p).thenComparing(Feature::id))
            .limit(integer(config(context), "max_features", 10, 1, 100)).toList();
    }
}
