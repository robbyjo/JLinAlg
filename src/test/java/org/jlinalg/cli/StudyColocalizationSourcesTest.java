package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyColocalizationSourcesTest {
    @TempDir Path project;

    @Test void ranksAncestryMatchedDownloadableGwasAndBuildsPartialPairs() throws Exception {
        Path upstream = project.resolve("upstream"); Files.createDirectories(upstream.resolve("standard-results"));
        Files.writeString(upstream.resolve("manifest.yaml"), "operation_type: validate_harmonize\nstate: complete\n");
        Files.writeString(upstream.resolve("standard-results/d1.tsv"), "feature_id\tp_value\nprobe1\t1e-9\n");
        Files.writeString(upstream.resolve("feature-annotation.tsv"), "feature_id\tgene_symbol\tgene_id\tchromosome\tstart\tend\tgenome_build\tmapping_status\nprobe1\tTEST\tENSG1\t1\t100\t300\tGRCh38\tlocated\n");
        Files.writeString(upstream.resolve("irqtl-1-candidates.tsv"), "feature_id\tvariant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tlog10_p\teaf\tn\treview_flags\n"
            + "probe1\tchr1:100:G:A:rs1\t1\t100\tA\tG\t.1\t.01\t-9\t.2\t4700\t.\n"
            + "probe1\tchr1:200:G:A:rs2\t1\t200\tA\tG\t.2\t.02\t-8\t.3\t4700\t.\n");
        Path discovery = project.resolve("discovery"); Files.createDirectory(discovery);
        Map<String, Object> smaller = study("GCST90000001", 100000), larger = study("GCST90000002", 400000);
        StudyPublicEvidence.write(discovery.resolve("source-discovery.yaml"), Map.of(
            "phenotype_search_terms", List.of(Map.of("phenotype", "body mass index", "origin", "hypothesis_outcome")),
            "evidence", List.of(row(smaller), row(larger))));
        Map<String, Object> config = Map.of("exposure_dataset", "d1", "exposure_ancestry", "European",
            "ld_population", "EUR", "related_phenotypes", List.of("glycated hemoglobin"));
        Map<String, Object> study = Map.of("datasets", List.of(Map.of("id", "d1", "cohort", "FHS", "type", "transcriptomics",
                "tissue", "Whole Blood", "files", Map.of("irqtl", "fixture"), "columns", Map.of("irqtl", Map.of("export_p_threshold", 1e-7)))),
            "hypotheses", List.of(Map.of("id", "h1", "outcome", "body mass index", "analysis", Map.of("bioinformatics", Map.of("colocalization", config)))));
        Context context = new Context(project, project.resolve("output"), Map.of("hypothesis", "h1", "type", "colocalization", "datasets", List.of("d1")),
            study, Map.of(), Map.of("upstream", upstream, "discovery", discovery));
        var prepared = StudyColocalizationSources.prepare(context);
        assertEquals(1, prepared.pairs().size());
        Map<String, Object> pair = StudyPublicEvidence.map(prepared.pairs().get(0));
        assertEquals(true, pair.get("partial_colocalization"));
        assertEquals("GCST90000002", StudyPublicEvidence.map(pair.get("outcome")).get("gwas_catalog_accession"));
        assertEquals("unknown", pair.get("sample_overlap"));
        assertTrue(Files.isRegularFile(project.resolve("output/automatic-colocalization-sources.yaml")));
    }

    private static Map<String, Object> study(String accession, int n) {
        return Map.of("accession_id", accession, "full_summary_stats_available", true, "snp_count", 1000000,
            "initial_sample_size", n + " European ancestry individuals", "discovery_ancestry", List.of(n + " European"),
            "disease_trait", "Body mass index", "pubmed_id", 12345);
    }

    private static Map<String, Object> row(Map<String, Object> metadata) {
        return Map.of("kind", "gwas_study", "status", "retrieved", "phenotype_query", "body mass index",
            "phenotype_origin", "hypothesis_outcome", "metadata", metadata);
    }
}
