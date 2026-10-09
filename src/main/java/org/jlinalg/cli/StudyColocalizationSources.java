/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
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
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.yaml.snakeyaml.Yaml;

/** Expands related-phenotype discovery into audited exploratory irQTL-GWAS pairs. */
final class StudyColocalizationSources {
    record Prepared(List<Object> pairs, Map<String, Object> audit) { }
    private record Candidate(String feature, String id, String chromosome,
            long position, String effect, String other, double beta, double se,
            double log10p, int n) { }
    private record Study(String trait, String origin, String accession,
            int sampleSize, Map<String, Object> metadata) { }

    private StudyColocalizationSources() { }

    static boolean configured(Context context) {
        return "colocalization".equals(context.operation().get("type"))
            && !map(StudyBioinformatics.config(context).get("colocalization")).isEmpty();
    }

    static Prepared prepare(Context context) throws IOException {
        Map<String, Object> config = map(StudyBioinformatics.config(context).get("colocalization"));
        if (config.isEmpty()) return new Prepared(List.of(), Map.of());
        String ancestry = ancestry(required(config, "exposure_ancestry"));
        String population = text(config.getOrDefault("ld_population", population(ancestry))).toUpperCase(Locale.ROOT);
        if (!ThousandGenomesGrch38Ld.ancestry(population).equals(ancestry))
            throw new IllegalArgumentException("colocalization.exposure_ancestry and ld_population must describe the same ancestry");
        String outcomeAncestry = ancestry(text(config.getOrDefault("outcome_ancestry", ancestry)));
        if (!outcomeAncestry.equals(ancestry))
            throw new IllegalArgumentException("Automatic colocalization currently requires ancestry-matched exposure, outcome, and LD");
        String provider = text(config.getOrDefault("outcome_provider", "auto"));
        if (!Set.of("auto", "gwas-catalog").contains(provider))
            throw new IllegalArgumentException("Automatic colocalization supports the GWAS Catalog harmonized provider");
        if (!text(config.getOrDefault("ld_provider", ThousandGenomesGrch38Ld.DATABASE)).equals(ThousandGenomesGrch38Ld.DATABASE))
            throw new IllegalArgumentException("Automatic GRCh38 colocalization supports ld_provider: " + ThousandGenomesGrch38Ld.DATABASE);
        boolean partial = !Boolean.FALSE.equals(config.get("allow_partial_regions"));
        Map<String, Object> dataset = exposureDataset(context, config);
        String cohort = text(dataset.get("cohort"));
        if (cohort.isBlank()) throw new IllegalArgumentException("The irQTL exposure dataset must declare its cohort for overlap assessment");
        Path candidateFile = candidateFile(context);
        int maximumVariants = StudyBioinformatics.integer(config, "max_exposure_variants", 500, 2, 5000);
        Map<String, List<Candidate>> byFeature = candidates(candidateFile);
        Path discoveryFile = discoveryFile(context);
        Map<String, Object> discovery;
        try (var reader = Files.newBufferedReader(discoveryFile)) { discovery = map(new Yaml().load(reader)); }
        List<Study> studies = studies(discovery, outcomeAncestry);
        int maximumTraits = StudyBioinformatics.integer(config, "max_colocalization_traits", 5, 1, 20);
        List<Study> selectedStudies = select(discovery, studies, maximumTraits);
        if (selectedStudies.isEmpty())
            throw new IllegalArgumentException("GWAS discovery found no ancestry-matched study with indexed harmonized full summary statistics");
        String overlap = text(config.getOrDefault("sample_overlap", "unknown"));
        String overlapBasis = text(config.getOrDefault("sample_overlap_basis",
            "Exposure cohort " + cohort + "; GWAS Catalog cohort metadata retained per selected study. Absence of a named shared cohort is not proof of zero overlap."));
        double minimumExposureOverlap = number(config.getOrDefault("minimum_exposure_overlap", 0.5));
        int window = StudyBioinformatics.integer(config, "region_window_bp", 500000, 0, 5000000);
        Path root = context.outputDirectory().resolve("automatic-colocalization-sources");
        Files.createDirectories(root);
        List<Object> pairs = new ArrayList<>(), featureAudit = new ArrayList<>();
        for (StudyBioinformatics.Feature feature : StudyBioinformatics.features(context)) {
            if (!"located".equals(feature.status()) || !"GRCh38".equals(feature.build())) {
                featureAudit.add(Map.of("feature_id", feature.id(), "status", "verified_grch38_location_required"));
                continue;
            }
            List<Candidate> available = new ArrayList<>(byFeature.getOrDefault(feature.id(), List.of()));
            available.sort(Comparator.comparingDouble(Candidate::log10p).thenComparing(Candidate::id));
            int total = available.size();
            if (available.size() > maximumVariants) available = new ArrayList<>(available.subList(0, maximumVariants));
            available.sort(Comparator.comparingLong(Candidate::position).thenComparing(Candidate::id));
            if (available.size() < 2) {
                featureAudit.add(Map.of("feature_id", feature.id(), "status", "fewer_than_two_exported_irQTL_variants", "available", total));
                continue;
            }
            long start = Math.max(1, Math.min(Long.parseLong(feature.start()) - window,
                available.stream().mapToLong(Candidate::position).min().orElseThrow()));
            long end = Math.max(Long.parseLong(feature.end()) + window,
                available.stream().mapToLong(Candidate::position).max().orElseThrow());
            if (end - start > 10_000_000L) {
                featureAudit.add(Map.of("feature_id", feature.id(), "status", "regional_window_exceeds_10mb"));
                continue;
            }
            Path featureDirectory = root.resolve(safe(feature.id())); Files.createDirectories(featureDirectory);
            Path exposureFile = featureDirectory.resolve("exposure.tsv"); writeExposure(exposureFile, available);
            int sampleSize = available.stream().mapToInt(Candidate::n).max().orElseThrow();
            for (Study study : selectedStudies) {
                Map<String, Object> pair = new LinkedHashMap<>();
                pair.put("id", safe(feature.id() + "-" + study.accession())); pair.put("feature_id", feature.id());
                pair.put("approved", true); pair.put("automatic_colocalization_source", true);
                pair.put("approval_basis", "The colocalization source policy is in the approved study contract. Trait suggestions remain advisory; the exact GWAS ranking, metadata, and source are recorded.");
                pair.put("rationale", "Exploratory colocalization of the exact measured irQTL feature with an ancestry-matched harmonized GWAS for " + study.trait());
                pair.put("sample_overlap", overlap); pair.put("sample_overlap_basis", overlapBasis);
                pair.put("partial_colocalization", partial); pair.put("minimum_exposure_overlap", minimumExposureOverlap);
                pair.put("region", Map.of("chromosome", feature.chromosome().replaceFirst("^chr", ""), "start", start, "end", end));
                Map<String, Object> declaration = map(map(dataset.get("columns")).get("irqtl"));
                Map<String, Object> exposure = new LinkedHashMap<>();
                exposure.put("prepared_path", exposureFile.toAbsolutePath().normalize().toString());
                exposure.put("source_id", text(config.getOrDefault("exposure_source_id", cohort + "-irqtl")));
                exposure.put("citation", text(config.getOrDefault("exposure_citation", "Study-investigator supplied irQTL results; source path and SHA-256 are recorded by irqtl_candidate_review")));
                exposure.put("genome_build", "GRCh38"); exposure.put("ancestry", ancestry); exposure.put("sample_size", sampleSize);
                exposure.put("complete_region", false); exposure.put("selection_basis", "significance_filtered_irQTL_export");
                if (declaration.containsKey("export_p_threshold")) exposure.put("export_p_threshold", declaration.get("export_p_threshold"));
                exposure.put("trait", "isoform ratio"); exposure.put("feature_id", feature.id()); exposure.put("assay", "irQTL");
                exposure.put("tissue", text(dataset.getOrDefault("tissue", "not_declared")));
                Map<String, Object> outcome = new LinkedHashMap<>();
                outcome.put("gwas_catalog_accession", study.accession()); outcome.put("source_id", study.accession());
                outcome.put("citation", GwasCatalogRegionalOutcome.DOCUMENTATION); outcome.put("publication", publication(study.metadata()));
                outcome.put("genome_build", "GRCh38"); outcome.put("ancestry", ancestry); outcome.put("sample_size", study.sampleSize());
                outcome.put("complete_region", true); outcome.put("selection_basis", "indexed_harmonized_full_summary_statistics");
                outcome.put("trait", study.trait()); outcome.put("phenotype_origin", study.origin()); outcome.put("catalog_metadata", study.metadata());
                Map<String, Object> ld = new LinkedHashMap<>(); ld.put("database", ThousandGenomesGrch38Ld.DATABASE);
                ld.put("population", population); ld.put("genome_build", "GRCh38"); ld.put("ancestry", ancestry);
                ld.put("citation", ThousandGenomesGrch38Ld.CITATION);
                pair.put("exposure", exposure); pair.put("outcome", outcome); pair.put("ld", ld); pairs.add(pair);
            }
            featureAudit.add(Map.of("feature_id", feature.id(), "status", "prepared", "exported_irQTL_variants", total,
                "retained_under_limit", available.size(), "exposure_sha256", FollowupSupport.hash(exposureFile)));
        }
        Map<String, Object> audit = new LinkedHashMap<>();
        audit.put("interface", "analysis.bioinformatics.colocalization"); audit.put("source_discovery", discoveryFile.toString());
        audit.put("exposure_dataset", dataset.get("id")); audit.put("exposure_cohort", cohort); audit.put("exposure_ancestry", ancestry);
        audit.put("outcome_provider", "GWAS Catalog harmonized summary statistics"); audit.put("selected_studies", selectedStudies.stream().map(StudyColocalizationSources::document).toList());
        audit.put("ld_provider", ThousandGenomesGrch38Ld.DATABASE); audit.put("ld_population", population);
        audit.put("partial_region_policy", partial ? "run_and_emit_mandatory_major_caveat" : "require_complete_regions");
        audit.put("candidate_file", candidateFile.toString()); audit.put("features", featureAudit); audit.put("generated_pairs", pairs.size());
        StudyPublicEvidence.write(context.outputDirectory().resolve("automatic-colocalization-sources.yaml"), audit);
        if (pairs.isEmpty()) throw new IllegalArgumentException("No feature had at least two exported irQTL variants for automatic colocalization");
        return new Prepared(List.copyOf(pairs), audit);
    }

    static List<Study> studies(Map<String, Object> discovery, String ancestry) {
        List<Study> result = new ArrayList<>();
        for (Object raw : list(discovery.get("evidence"))) {
            Map<String, Object> row = map(raw), metadata = map(row.get("metadata"));
            if (!"gwas_study".equals(row.get("kind")) || !"retrieved".equals(row.get("status"))
                    || !Boolean.TRUE.equals(metadata.get("full_summary_stats_available"))) continue;
            String accession = text(metadata.get("accession_id"));
            if (!accession.matches("GCST[0-9]{6,}") || !ancestryMatch(metadata, ancestry)) continue;
            if (metadata.get("snp_count") instanceof Number count && count.longValue() < 100000) continue;
            if (text(metadata.get("disease_trait")).toLowerCase(Locale.ROOT).contains("cnv")) continue;
            int sample = sampleSize(metadata.get("initial_sample_size"));
            if (sample < 1000) continue;
            result.add(new Study(text(row.get("phenotype_query")), text(row.get("phenotype_origin")), accession, sample, metadata));
        }
        return result;
    }

    private static List<Study> select(Map<String, Object> discovery, List<Study> studies, int maximum) {
        List<String> order = list(discovery.get("phenotype_search_terms")).stream().map(StudyPublicEvidence::map)
            .map(term -> text(term.get("phenotype"))).toList();
        List<Study> result = new ArrayList<>();
        for (String trait : order) {
            studies.stream().filter(study -> study.trait().equalsIgnoreCase(trait))
                .sorted(Comparator.comparingInt(Study::sampleSize).reversed().thenComparing(Study::accession))
                .findFirst().ifPresent(result::add);
            if (result.size() >= maximum) break;
        }
        return result;
    }

    private static boolean ancestryMatch(Map<String, Object> metadata, String ancestry) {
        String wanted = ancestry.toLowerCase(Locale.ROOT);
        List<String> declarations = list(metadata.get("discovery_ancestry")).stream().map(StudyPublicEvidence::text)
            .map(value -> value.toLowerCase(Locale.ROOT)).toList();
        if (declarations.isEmpty() || declarations.stream().noneMatch(value -> value.contains(wanted))) return false;
        Set<String> other = new LinkedHashSet<>(List.of("european", "african", "east asian", "south asian", "hispanic", "latino", "admixed american"));
        other.remove(wanted);
        return declarations.stream().noneMatch(value -> other.stream().anyMatch(value::contains));
    }

    private static int sampleSize(Object raw) {
        Matcher matcher = Pattern.compile("([0-9][0-9,]*)").matcher(text(raw));
        if (!matcher.find()) return -1;
        try { long value = Long.parseLong(matcher.group(1).replace(",", "")); return value > Integer.MAX_VALUE ? -1 : (int) value; }
        catch (NumberFormatException invalid) { return -1; }
    }

    private static String publication(Map<String, Object> metadata) {
        String pmid = text(metadata.get("pubmed_id"));
        return pmid.matches("[1-9][0-9]*") ? "https://pubmed.ncbi.nlm.nih.gov/" + pmid + "/" : GwasCatalogRegionalOutcome.DOCUMENTATION;
    }

    private static Map<String, Object> document(Study study) {
        return Map.of("phenotype", study.trait(), "phenotype_origin", study.origin(), "accession", study.accession(),
            "sample_size", study.sampleSize(), "publication", publication(study.metadata()), "selection",
            "Largest ancestry-matched >=100,000-variant study with downloadable summary statistics among the bounded GWAS Catalog pages retrieved for this exact trait label");
    }

    private static Map<String, Object> exposureDataset(Context context, Map<String, Object> config) {
        String requested = text(config.get("exposure_dataset"));
        List<Map<String, Object>> eligible = StudyBioinformatics.datasets(context).stream()
            .filter(dataset -> map(dataset.get("files")).containsKey("irqtl"))
            .filter(dataset -> requested.isBlank() || requested.equals(dataset.get("id"))).toList();
        if (eligible.size() != 1) throw new IllegalArgumentException("colocalization.exposure_dataset must identify the one selected dataset with an irQTL source");
        return eligible.get(0);
    }

    private static Path candidateFile(Context context) throws IOException {
        List<Path> found = new ArrayList<>();
        for (Path directory : context.dependencyDirectories().values()) if (Files.isDirectory(directory)) try (var files = Files.list(directory)) {
            found.addAll(files.filter(path -> path.getFileName().toString().matches("irqtl-[0-9]+-candidates\\.tsv")).sorted().toList());
        }
        if (found.size() != 1) throw new IOException("Automatic colocalization needs exactly one irqtl_candidate_review candidate file dependency; found " + found.size());
        return found.get(0);
    }

    private static Path discoveryFile(Context context) throws IOException {
        List<Path> found = context.dependencyDirectories().values().stream().map(path -> path.resolve("source-discovery.yaml"))
            .filter(Files::isRegularFile).sorted().toList();
        if (found.size() != 1) throw new IOException("Automatic colocalization needs exactly one source_discovery artifact; found " + found.size());
        return found.get(0);
    }

    private static Map<String, List<Candidate>> candidates(Path path) throws IOException {
        DelimitedData table = DelimitedData.read(path); Map<String, List<Candidate>> result = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            String[] tokens = row[table.column("variant_id")].split(":", -1); String id = tokens[tokens.length - 1];
            if (!id.matches("rs[1-9][0-9]*")) continue;
            double n = number(row[table.column("n")]); if (n < 3 || n != Math.rint(n)) continue;
            Candidate value = new Candidate(row[table.column("feature_id")], id,
                row[table.column("chromosome")].replaceFirst("^chr", ""), Long.parseLong(row[table.column("position")]),
                row[table.column("effect_allele")].toUpperCase(Locale.ROOT), row[table.column("other_allele")].toUpperCase(Locale.ROOT),
                number(row[table.column("beta")]), number(row[table.column("se")]), number(row[table.column("log10_p")]), (int) n);
            result.computeIfAbsent(value.feature(), ignored -> new ArrayList<>()).add(value);
        }
        return result;
    }

    private static void writeExposure(Path output, List<Candidate> variants) throws IOException {
        StringBuilder table = new StringBuilder("feature_id\tvariant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n");
        for (Candidate value : variants) table.append(value.feature()).append('\t').append(value.id()).append('\t').append(value.chromosome()).append('\t')
            .append(value.position()).append('\t').append(value.effect()).append('\t').append(value.other()).append('\t').append(value.beta()).append('\t')
            .append(value.se()).append('\t').append(Math.pow(10, value.log10p())).append('\n');
        Files.writeString(output, table);
    }

    private static String ancestry(String value) {
        String normalized = value.toLowerCase(Locale.ROOT).replaceAll("[^a-z]", "");
        return switch (normalized) {
            case "eur", "european" -> "European"; case "afr", "african" -> "African";
            case "amr", "admixedamerican" -> "Admixed American"; case "eas", "eastasian" -> "East Asian";
            case "sas", "southasian" -> "South Asian";
            default -> throw new IllegalArgumentException("colocalization ancestry must map explicitly to EUR, AFR, AMR, EAS, or SAS");
        };
    }

    private static String population(String ancestry) {
        return switch (ancestry) {
            case "European" -> "EUR"; case "African" -> "AFR"; case "Admixed American" -> "AMR";
            case "East Asian" -> "EAS"; case "South Asian" -> "SAS"; default -> throw new IllegalArgumentException("Unsupported ancestry");
        };
    }

    private static String required(Map<String, Object> source, String key) {
        String result = text(source.get(key)); if (result.isBlank()) throw new IllegalArgumentException("colocalization." + key + " is required"); return result;
    }

    private static double number(Object value) {
        double result = Double.parseDouble(text(value)); if (!Double.isFinite(result)) throw new IllegalArgumentException("Nonfinite colocalization source value"); return result;
    }

    private static String safe(String value) {
        String result = value.replaceAll("[^A-Za-z0-9_-]+", "-");
        if (result.length() > 64) result = result.substring(0, 55) + "-" + Integer.toUnsignedString(value.hashCode(), 16);
        return result;
    }
}
