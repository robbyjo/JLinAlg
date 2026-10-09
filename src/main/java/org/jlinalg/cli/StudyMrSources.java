/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;

/** Expands a concise MR source policy into per-feature reviewed pair inputs. */
final class StudyMrSources {
    record Prepared(List<Object> pairs, Map<String, Object> audit) { }
    private record Candidate(String feature, String id, String chromosome,
            long position, String effect, String other, double beta, double se,
            double log10p, int n) { }

    private StudyMrSources() { }

    static boolean configured(Context context) {
        return "mendelian_randomization".equals(context.operation().get("type"))
            && !map(StudyBioinformatics.config(context).get("mr")).isEmpty();
    }

    static Prepared prepare(Context context) throws IOException {
        Map<String, Object> config = map(StudyBioinformatics.config(context).get("mr"));
        if (config.isEmpty()) return new Prepared(List.of(), Map.of());
        String ancestry = ancestry(required(config, "exposure_ancestry"));
        String population = text(config.getOrDefault("ld_population", population(ancestry))).toUpperCase(Locale.ROOT);
        if (!ThousandGenomesGrch38Ld.ancestry(population).equals(ancestry))
            throw new IllegalArgumentException("mr.exposure_ancestry and mr.ld_population must describe the same ancestry");
        String outcomeAncestry = text(config.getOrDefault("outcome_ancestry", population)).toUpperCase(Locale.ROOT);
        String trait = text(config.getOrDefault("outcome_trait", StudyBioinformatics.hypothesis(context).get("outcome")));
        String provider = text(config.getOrDefault("outcome_provider", "auto"));
        if (provider.equals("auto") && PanUkbbGwasOutcome.supports(trait, outcomeAncestry)) provider = PanUkbbGwasOutcome.PROVIDER;
        if (!provider.equals(PanUkbbGwasOutcome.PROVIDER) || !PanUkbbGwasOutcome.supports(trait, outcomeAncestry))
            throw new IllegalArgumentException("No automatic full-summary-statistics adapter is registered for outcome '" + trait
                + "' and ancestry " + outcomeAncestry + ". Use GWAS Catalog discovery, select a full harmonized study, and declare a reviewed genetic_pair; catalog lead hits are insufficient");
        String ldProvider = text(config.getOrDefault("ld_provider", ThousandGenomesGrch38Ld.DATABASE));
        if (!ldProvider.equals(ThousandGenomesGrch38Ld.DATABASE))
            throw new IllegalArgumentException("The automatic GRCh38 MR adapter currently supports ld_provider: " + ThousandGenomesGrch38Ld.DATABASE);
        Map<String, Object> dataset = exposureDataset(context, config);
        String cohort = text(dataset.get("cohort"));
        if (cohort.isBlank()) throw new IllegalArgumentException("The irQTL exposure dataset must declare its cohort so outcome overlap can be assessed");
        if (cohort.toLowerCase(Locale.ROOT).contains("ukb") || cohort.toLowerCase(Locale.ROOT).contains("biobank"))
            throw new IllegalArgumentException("Pan-UKB cannot be marked nonoverlapping when the exposure cohort may be UK Biobank");
        String overlap = text(config.getOrDefault("sample_overlap", "none"));
        if (!overlap.equals("none")) throw new IllegalArgumentException("Automatic MR requires an explicitly independent outcome; overlapping-sample MR needs a reviewed adapter");
        String overlapBasis = text(config.getOrDefault("sample_overlap_basis",
            "Exposure cohort " + cohort + "; outcome cohort UK Biobank only; distinct named source cohorts"));
        double instrumentP = number(config.getOrDefault("instrument_p", 5e-8));
        if (!(instrumentP > 0 && instrumentP <= 5e-8)) throw new IllegalArgumentException("mr.instrument_p must be in (0,5e-8]");
        int maximum = config.containsKey("max_instruments")
            ? StudyBioinformatics.integer(config, "max_instruments", 100, 2, 500) : 100;
        Path candidates = candidateFile(context);
        Map<String, List<Candidate>> byFeature = candidates(candidates, instrumentP);
        List<Object> pairs = new ArrayList<>(), featureAudit = new ArrayList<>();
        Path root = context.outputDirectory().resolve("automatic-mr-sources");
        Files.createDirectories(root);
        for (StudyBioinformatics.Feature feature : StudyBioinformatics.features(context)) {
            if (!"located".equals(feature.status())) {
                featureAudit.add(Map.of("feature_id", feature.id(), "status", "missing_verified_grch38_location"));
                continue;
            }
            List<Candidate> available = new ArrayList<>(byFeature.getOrDefault(feature.id(), List.of()));
            available.sort(Comparator.comparingDouble(Candidate::log10p).thenComparing(Candidate::id));
            int allStrong = available.size();
            if (available.size() > maximum) available = new ArrayList<>(available.subList(0, maximum));
            Map<String, Object> oneAudit = new LinkedHashMap<>();
            oneAudit.put("feature_id", feature.id()); oneAudit.put("strong_candidate_instruments", allStrong);
            oneAudit.put("instrument_lookup_limit", maximum);
            if (available.size() < 2) {
                oneAudit.put("status", "insufficient_strong_candidates_before_outcome_or_ld_lookup");
                featureAudit.add(oneAudit); continue;
            }
            Path directory = root.resolve(safe(feature.id())); Files.createDirectories(directory);
            Path exposureFile = directory.resolve("exposure.tsv");
            writeExposure(exposureFile, available);
            List<PanUkbbGwasOutcome.Instrument> instruments = available.stream()
                .map(v -> new PanUkbbGwasOutcome.Instrument(v.id(), v.chromosome(), v.position(), v.effect(), v.other())).toList();
            Path outcomeFile = directory.resolve("outcome.tsv");
            Map<String, Object> outcomeAudit = PanUkbbGwasOutcome.extract(context, instruments, outcomeFile);
            int sampleSize = available.stream().mapToInt(Candidate::n).max().orElseThrow();
            long start = Math.min(Long.parseLong(feature.start()), available.stream().mapToLong(Candidate::position).min().orElseThrow());
            long end = Math.max(Long.parseLong(feature.end()), available.stream().mapToLong(Candidate::position).max().orElseThrow());
            if (end - start > 10_000_000L) throw new IllegalArgumentException("Automatic MR instrument span exceeds 10 Mb for " + feature.id());
            Map<String, Object> pair = new LinkedHashMap<>();
            pair.put("id", safe(feature.id())); pair.put("feature_id", feature.id()); pair.put("approved", true);
            pair.put("approval_basis", "The concise mr source policy is part of the approved study contract; provider selection and every generated field are recorded in automatic-mr-sources.yaml");
            pair.put("rationale", "Exploratory cis-irQTL MR of the exact measured transcript using an independent phenotype-matched UK Biobank outcome");
            pair.put("sample_overlap", overlap); pair.put("sample_overlap_basis", overlapBasis);
            pair.put("instrument_validity_rationale", text(config.getOrDefault("instrument_validity_rationale",
                "Cis irQTLs must pass p<=5e-8, F>=10, and ancestry-matched LD clumping at r2<=0.001; this operational screen does not prove exclusion restriction or absence of horizontal pleiotropy")));
            pair.put("instrument_p", instrumentP); pair.put("minimum_f", number(config.getOrDefault("minimum_f", 10)));
            pair.put("clump_r2", number(config.getOrDefault("clump_r2", 0.001)));
            pair.put("region", Map.of("chromosome", feature.chromosome().replaceFirst("^chr", ""), "start", start, "end", end));
            Map<String, Object> declaration = map(map(dataset.get("columns")).get("irqtl"));
            Map<String, Object> exposure = new LinkedHashMap<>();
            exposure.put("prepared_path", exposureFile.toAbsolutePath().normalize().toString());
            exposure.put("source_id", text(config.getOrDefault("exposure_source_id", cohort + "-irqtl")));
            exposure.put("citation", text(config.getOrDefault("exposure_citation", "Study-investigator supplied irQTL results; source path and SHA-256 are recorded by irqtl_candidate_review")));
            exposure.put("genome_build", "GRCh38"); exposure.put("ancestry", ancestry); exposure.put("sample_size", sampleSize);
            exposure.put("complete_region", false); exposure.put("export_p_threshold", number(declaration.get("export_p_threshold")));
            exposure.put("trait", "isoform ratio"); exposure.put("feature_id", feature.id()); exposure.put("assay", "irQTL");
            exposure.put("tissue", text(dataset.getOrDefault("tissue", "not_declared")));
            Map<String, Object> outcome = new LinkedHashMap<>();
            outcome.put("prepared_path", outcomeFile.toAbsolutePath().normalize().toString());
            outcome.put("source_id", "pan-ukb-21001-eur"); outcome.put("citation", PanUkbbGwasOutcome.DOCUMENTATION);
            outcome.put("genome_build", "GRCh38"); outcome.put("source_genome_build", "GRCh37");
            outcome.put("build_projection", "NCBI SPDI all_equivalent_contextual, audited per instrument");
            outcome.put("ancestry", ancestry); outcome.put("sample_size", PanUkbbGwasOutcome.BMI_EUR_N);
            outcome.put("complete_region", false); outcome.put("selection_basis", "exposure_instrument_lookup"); outcome.put("trait", text(StudyBioinformatics.hypothesis(context).get("outcome")));
            Map<String, Object> ld = new LinkedHashMap<>(); ld.put("database", ldProvider); ld.put("population", population);
            ld.put("genome_build", "GRCh38"); ld.put("ancestry", ancestry); ld.put("citation", ThousandGenomesGrch38Ld.CITATION);
            pair.put("exposure", exposure); pair.put("outcome", outcome); pair.put("ld", ld);
            pairs.add(pair);
            oneAudit.put("status", "prepared"); oneAudit.put("pair_id", pair.get("id"));
            oneAudit.put("exposure_sha256", FollowupSupport.hash(exposureFile)); oneAudit.put("outcome", outcomeAudit);
            featureAudit.add(oneAudit);
        }
        Map<String, Object> audit = new LinkedHashMap<>();
        audit.put("interface", "analysis.bioinformatics.mr"); audit.put("legacy_interface", "analysis.bioinformatics.genetic_pairs remains supported for custom sources");
        audit.put("exposure_dataset", dataset.get("id")); audit.put("exposure_cohort", cohort); audit.put("exposure_ancestry", ancestry);
        audit.put("outcome", PanUkbbGwasOutcome.descriptor()); audit.put("ld_provider", ldProvider); audit.put("ld_population", population);
        audit.put("sample_overlap", overlap); audit.put("sample_overlap_basis", overlapBasis); audit.put("candidate_file", candidates.toString());
        audit.put("features", featureAudit); audit.put("generated_pairs", pairs.size());
        StudyPublicEvidence.write(context.outputDirectory().resolve("automatic-mr-sources.yaml"), audit);
        if (pairs.isEmpty()) throw new IllegalArgumentException("No feature has at least two genome-wide-significant irQTL candidates for automatic MR");
        return new Prepared(List.copyOf(pairs), audit);
    }

    private static Map<String, Object> exposureDataset(Context context,
            Map<String, Object> config) {
        String requested = text(config.get("exposure_dataset"));
        List<Map<String, Object>> eligible = StudyBioinformatics.datasets(context).stream()
            .filter(dataset -> map(dataset.get("files")).containsKey("irqtl"))
            .filter(dataset -> requested.isBlank() || requested.equals(dataset.get("id"))).toList();
        if (eligible.size() != 1) throw new IllegalArgumentException("mr.exposure_dataset must identify the one selected dataset with an irQTL source");
        return eligible.get(0);
    }

    private static Path candidateFile(Context context) throws IOException {
        List<Path> found = new ArrayList<>();
        for (Path directory : context.dependencyDirectories().values()) {
            if (!Files.isDirectory(directory)) continue;
            try (var files = Files.list(directory)) {
                found.addAll(files.filter(path -> path.getFileName().toString().matches("irqtl-[0-9]+-candidates\\.tsv"))
                    .sorted().toList());
            }
        }
        if (found.size() != 1) throw new IOException("Automatic MR needs exactly one irqtl_candidate_review candidate file dependency; found " + found.size());
        return found.get(0);
    }

    private static Map<String, List<Candidate>> candidates(Path path,
            double threshold) throws IOException {
        DelimitedData table = DelimitedData.read(path);
        Map<String, List<Candidate>> result = new LinkedHashMap<>();
        double cutoff = Math.log10(threshold);
        for (String[] row : table.rows()) {
            double logp = number(row[table.column("log10_p")]);
            if (logp > cutoff) continue;
            String raw = row[table.column("variant_id")];
            String[] tokens = raw.split(":", -1); String id = tokens[tokens.length - 1];
            if (!id.matches("rs[1-9][0-9]*")) continue;
            double n = number(row[table.column("n")]);
            if (n < 3 || n != Math.rint(n)) continue;
            Candidate value = new Candidate(row[table.column("feature_id")], id,
                row[table.column("chromosome")].replaceFirst("^chr", ""),
                Long.parseLong(row[table.column("position")]),
                row[table.column("effect_allele")].toUpperCase(Locale.ROOT),
                row[table.column("other_allele")].toUpperCase(Locale.ROOT),
                number(row[table.column("beta")]), number(row[table.column("se")]),
                logp, (int) n);
            result.computeIfAbsent(value.feature(), key -> new ArrayList<>()).add(value);
        }
        return result;
    }

    private static void writeExposure(Path output, List<Candidate> variants)
            throws IOException {
        StringBuilder table = new StringBuilder("feature_id\tvariant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n");
        for (Candidate v : variants) table.append(v.feature()).append('\t').append(v.id()).append('\t')
            .append(v.chromosome()).append('\t').append(v.position()).append('\t').append(v.effect()).append('\t')
            .append(v.other()).append('\t').append(v.beta()).append('\t').append(v.se()).append('\t')
            .append(Math.pow(10, v.log10p())).append('\n');
        Files.writeString(output, table);
    }

    private static String ancestry(String value) {
        String normalized = value.toLowerCase(Locale.ROOT).replaceAll("[^a-z]", "");
        return switch (normalized) {
            case "eur", "european" -> "European";
            case "afr", "african" -> "African";
            case "amr", "admixedamerican" -> "Admixed American";
            case "eas", "eastasian" -> "East Asian";
            case "sas", "southasian" -> "South Asian";
            default -> throw new IllegalArgumentException("mr.exposure_ancestry must map explicitly to EUR, AFR, AMR, EAS, or SAS");
        };
    }

    private static String population(String ancestry) {
        return switch (ancestry) {
            case "European" -> "EUR"; case "African" -> "AFR";
            case "Admixed American" -> "AMR"; case "East Asian" -> "EAS";
            case "South Asian" -> "SAS"; default -> throw new IllegalArgumentException("Unsupported ancestry");
        };
    }

    private static String required(Map<String, Object> source, String key) {
        String result = text(source.get(key));
        if (result.isBlank()) throw new IllegalArgumentException("mr." + key + " is required");
        return result;
    }

    private static double number(Object value) {
        double result = Double.parseDouble(text(value));
        if (!Double.isFinite(result)) throw new IllegalArgumentException("Nonfinite numeric MR source declaration");
        return result;
    }

    private static String safe(String value) {
        String result = value.replaceAll("[^A-Za-z0-9_-]+", "-");
        if (result.length() > 64) result = result.substring(0, 55) + "-" + Integer.toUnsignedString(value.hashCode(), 16);
        return result;
    }
}
