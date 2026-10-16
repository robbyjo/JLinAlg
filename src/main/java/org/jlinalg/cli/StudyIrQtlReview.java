/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.BufferedReader;
import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.jlinalg.study.StudyOperationExecutor;

/** Streaming local isoform-ratio QTL review; never substitutes gene-expression QTLs. */
final class StudyIrQtlReview implements StudyOperationExecutor {
    @Override public Result execute(Context context) throws IOException {
        Files.createDirectories(context.outputDirectory());
        List<StudyBioinformatics.Feature> features = StudyBioinformatics.geneticFeatures(context);
        Set<String> selected = new HashSet<>();
        features.forEach(f -> selected.add(f.id()));
        List<Object> sources = new ArrayList<>();
        Set<Path> seen = new HashSet<>();
        for (Map<String, Object> dataset : StudyBioinformatics.datasets(context)) {
            Map<String, Object> files = map(dataset.get("files"));
            if (!files.containsKey("irqtl")) continue;
            Path source = context.project().resolve(text(files.get("irqtl"))).toAbsolutePath().normalize();
            if (!seen.add(source)) continue;
            Map<String, Object> declaration = map(map(dataset.get("columns")).get("irqtl"));
            String stem = "irqtl-" + seen.size();
            Map<String, Object> reviewed = new LinkedHashMap<>(scan(source, declaration, selected,
                context.outputDirectory().resolve(stem + "-candidates.tsv")));
            reviewed.put("dataset_id", dataset.get("id"));
            sources.add(reviewed);
        }
        Map<String, Object> details = new LinkedHashMap<>();
        details.put("molecular_trait", "isoform ratio quantitative trait locus (irQTL)");
        if (StudyBioinformatics.mrRequested(context)) details.put("mr_feature_selection", StudyBioinformatics.mrFeatureSelection(context));
        details.put("sources", sources); details.put("selected_features", features.stream().map(StudyBioinformatics.Feature::document).toList());
        details.put("interpretation", "Candidate review only, not a Mendelian randomization (MR) estimate or colocalization result. Exact transcript versions are matched; no gene-expression QTL substitution. Absence from a significant-only export is not evidence of no association.");
        details.put("downstream_requirements", Map.of(
            "mendelian_randomization", "Phenotype-matched outcome effects, compatible alleles/build/ancestry, LD, sample-overlap assessment and instrument review. A configured analysis.bioinformatics.mr provider policy can prepare the supported public inputs.",
            "colocalization", "Complete regional irQTL and outcome-GWAS statistics including nonsignificant variants; this significance-filtered irQTL export is not sufficient."));
        details.put("needs", List.of());
        write(context.outputDirectory().resolve("irqtl-review.yaml"), details);
        return new Result(State.COMPLETE, "Reviewed " + sources.size() + " local irQTL sources and wrote audited candidate inputs; downstream genetic inference is tracked separately",
            List.of("irqtl-review.yaml"), details);
    }

    static Map<String, Object> scan(Path source, Map<String, Object> declaration, Set<String> selected, Path output) throws IOException {
        if (!"GRCh38".equals(declaration.get("genome_build")) || !"ALT".equals(declaration.get("effect_allele")))
            throw new IllegalArgumentException("This irQTL adapter requires explicitly declared GRCh38 and ALT-allele effects; it does not guess from filenames");
        if (!(declaration.get("complete_region") instanceof Boolean) || text(declaration.get("rationale")).isBlank())
            throw new IllegalArgumentException("Declare irQTL complete_region and a scientific rationale");
        if (!Boolean.TRUE.equals(declaration.get("complete_region")) && !(declaration.get("export_p_threshold") instanceof Number))
            throw new IllegalArgumentException("A filtered irQTL source requires its export_p_threshold");
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("source_path", source.toString()); result.put("declaration", declaration);
        result.put("scientific_rationale", declaration.get("rationale")); result.put("request_attempted", false);
        result.put("started_at", java.time.Instant.now().toString());
        result.put("colocalization_status", Boolean.TRUE.equals(declaration.get("complete_region")) ? "other_prerequisites_required" : "blocked_significance_filtered_irQTL");
        Path audit = output.resolveSibling(output.getFileName() + ".source.yaml");
        write(audit, result);
        String hash = FollowupSupport.hash(source);
        Map<String, Long> counts = new LinkedHashMap<>(), flaggedCounts = new LinkedHashMap<>();
        selected.stream().sorted().forEach(id -> { counts.put(id, 0L); flaggedCounts.put(id, 0L); });
        long rows = 0, selectedRows = 0, matched = 0, retained = 0, invalid = 0,
            nonCis = 0;
        double largestLogP = Double.NEGATIVE_INFINITY;
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(open(source), StandardCharsets.UTF_8), 1024 * 1024);
                BufferedWriter writer = Files.newBufferedWriter(output)) {
            String first = reader.readLine();
            if (first == null) throw new IOException("Empty irQTL source");
            List<String> header = DelimitedData.parse(first, ',', 1, source);
            for (String required : List.of("SNP", "SNP_Chr", "SNP_Pos", "EAF", "N", "Tx", "Is_Cis", "RSq", "Fx", "T", "log10P"))
                if (!header.contains(required)) throw new IOException("irQTL source is missing column " + required);
            int tx = header.indexOf("Tx");
            writer.write("feature_id\tvariant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tlog10_p\teaf\tn\treview_flags\n");
            for (String line; (line = reader.readLine()) != null;) {
                if (line.isBlank()) continue;
                rows++;
                List<String> row = DelimitedData.parse(line, ',', rows + 1, source);
                if (row.size() != header.size()) throw new IOException("irQTL row width differs at row " + rows);
                double logP = Double.parseDouble(row.get(header.indexOf("log10P")));
                if (Double.isFinite(logP)) largestLogP = Math.max(largestLogP, logP);
                String feature = row.get(tx);
                if (!selected.contains(feature)) continue;
                selectedRows++;
                if (!row.get(header.indexOf("Is_Cis")).equals("1")) {
                    nonCis++;
                    continue;
                }
                matched++; counts.merge(feature, 1L, Long::sum);
                String variant = row.get(header.indexOf("SNP")); String[] allele = variant.split(":", -1);
                double beta = Double.parseDouble(row.get(header.indexOf("Fx"))), t = Double.parseDouble(row.get(header.indexOf("T")));
                double eaf = Double.parseDouble(row.get(header.indexOf("EAF"))), n = Double.parseDouble(row.get(header.indexOf("N")));
                double r2 = Double.parseDouble(row.get(header.indexOf("RSq"))), se = Math.abs(beta / t);
                List<String> flags = new ArrayList<>();
                if (allele.length < 4 || !allele[2].matches("[ACGT]+") || !allele[3].matches("[ACGT]+")
                        || !Double.isFinite(beta) || !Double.isFinite(t) || !Double.isFinite(se) || se <= 0
                        || !Double.isFinite(eaf) || eaf < 0 || eaf > 1 || !Double.isFinite(n) || n < 3 || n != Math.rint(n)
                        || !Double.isFinite(logP) || logP > 0 || !Double.isFinite(r2) || r2 < 0 || r2 > 1) { invalid++; continue; }
                if (!allele[0].replaceFirst("^chr", "").equals(row.get(header.indexOf("SNP_Chr")).replaceFirst("^chr", ""))
                        || !allele[1].equals(row.get(header.indexOf("SNP_Pos")))) throw new IOException("irQTL variant identifier and coordinate columns disagree");
                if (Math.min(eaf, 1 - eaf) < 0.01) flags.add("minor_allele_frequency_below_0.01");
                if (r2 >= 0.99 || Math.abs(t) > 10000) flags.add("near_perfect_fit_or_extreme_statistic");
                if (logP > Math.log10(5e-8)) flags.add("above_default_MR_instrument_threshold");
                if (t * t < 10) flags.add("weak_instrument_statistic");
                if (!flags.isEmpty()) flaggedCounts.merge(feature, 1L, Long::sum);
                retained++;
                writer.write(String.join("\t", feature, variant, row.get(header.indexOf("SNP_Chr")), row.get(header.indexOf("SNP_Pos")),
                    allele[3], allele[2], Double.toString(beta), Double.toString(se), Double.toString(logP), Double.toString(eaf), Double.toString(n),
                    flags.isEmpty() ? "unreviewed_no_heuristic_flags" : String.join(";", flags)) + "\n");
            }
        }
        if (!hash.equals(FollowupSupport.hash(source))) throw new IOException("irQTL source changed during scan");
        result.put("source_sha256", hash); result.put("rows_scanned", rows);
        result.put("selected_feature_rows", selectedRows);
        result.put("non_cis_rows_excluded", nonCis);
        result.put("exact_feature_matches", matched);
        result.put("matches_by_feature", counts); result.put("flagged_by_feature", flaggedCounts);
        result.put("invalid_matched_rows_excluded", invalid); result.put("candidate_rows_written", retained);
        result.put("output_limited", false); result.put("output_order", "Source order; every valid cis-irQTL row for the selected features is written");
        result.put("maximum_observed_log10_p", Double.isFinite(largestLogP) ? largestLogP : "unavailable");
        result.put("candidate_file", output.getFileName().toString()); result.put("candidate_sha256", FollowupSupport.hash(output));
        result.put("standard_error_derivation", "abs(Fx / T); requires confirmation that T is the signed coefficient statistic. log10P is preserved to avoid probability underflow.");
        result.put("review_flags", "Heuristics identify records to investigate; they do not establish artifacts or instrument validity. No flags is not approval.");
        result.put("finished_at", java.time.Instant.now().toString()); write(audit, result);
        return result;
    }

    private static InputStream open(Path source) throws IOException {
        // The bzip2 bit decoder consumes compressed bytes individually. Buffer
        // before decoding, not only after it, to avoid a disk syscall per byte.
        InputStream input = new BufferedInputStream(Files.newInputStream(source), 1024 * 1024);
        try {
            if (source.toString().endsWith(".bz2")) return new BZip2CompressorInputStream(input, true);
            if (source.toString().endsWith(".gz")) return new GZIPInputStream(input);
            return input;
        } catch (IOException | RuntimeException e) { input.close(); throw e; }
    }
}
