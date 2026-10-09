/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;

import htsjdk.tribble.readers.TabixReader;
import java.io.IOException;
import java.net.URI;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyOperationExecutor.Context;

/** Bounded Pan-UK Biobank instrument lookup with explicit NCBI build projection. */
final class PanUkbbGwasOutcome {
    static final String PROVIDER = "pan-ukb";
    static final String BMI_TRAIT = "21001";
    static final int BMI_EUR_N = 419163;
    static final String BMI_FILE = "continuous-21001-both_sexes-irnt.tsv.bgz";
    static final String BMI_MD5 = "ca9ceaa99c79d510fb0e26a780787416";
    static final String BMI_INDEX_MD5 = "15552b6a7670abdfd74834b5d3b5a3b6";
    static final String ROOT = "https://pan-ukb-us-east-1.s3.amazonaws.com/";
    static final String DATA = ROOT + "sumstats_flat_files/" + BMI_FILE;
    static final String INDEX = ROOT + "sumstats_flat_files_tabix/" + BMI_FILE + ".tbi";
    static final String DOCUMENTATION = "https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html";
    static final String NCBI = "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/";
    private static long lastNcbiRequest;

    record Instrument(String id, String chromosome, long position,
            String effect, String other) { }
    record Projected(String chromosome, long position, String effect,
            String other, String ncbiUrl) { }
    @FunctionalInterface interface Projector {
        Projected project(Instrument variant) throws IOException;
    }

    private PanUkbbGwasOutcome() { }

    static boolean supports(String trait, String ancestry) {
        String normalized = trait == null ? "" : trait.toLowerCase(Locale.ROOT)
            .replaceAll("[^a-z0-9]+", "");
        return Set.of("21001", "bmi", "bodymassindex").contains(normalized)
            && "EUR".equalsIgnoreCase(ancestry);
    }

    static Map<String, Object> descriptor() {
        return Map.ofEntries(
            Map.entry("provider", PROVIDER),
            Map.entry("trait", BMI_TRAIT),
            Map.entry("description", "Body mass index (BMI), inverse-rank normal transformed"),
            Map.entry("ancestry", "EUR"),
            Map.entry("sample_size", BMI_EUR_N),
            Map.entry("source_genome_build", "GRCh37"),
            Map.entry("effect_allele", "ALT"),
            Map.entry("summary_statistics", DATA),
            Map.entry("tabix_index", INDEX),
            Map.entry("summary_statistics_md5", BMI_MD5),
            Map.entry("tabix_index_md5", BMI_INDEX_MD5),
            Map.entry("documentation", DOCUMENTATION),
            Map.entry("selection_note", "UK Biobank-only EUR outcome avoids FHS exposure-sample overlap; do not substitute a UKB-plus-GIANT meta-analysis without reassessing overlap"));
    }

    static Map<String, Object> extract(Context context, List<Instrument> variants,
            Path output) throws IOException {
        return extract(context, variants, output, PanUkbbGwasOutcome::project);
    }

    static Map<String, Object> extract(Context context, List<Instrument> variants,
            Path output, Projector projector) throws IOException {
        if (variants.isEmpty()) throw new IllegalArgumentException("No exposure instruments were supplied for outcome lookup");
        String blocked = StudyPublicEvidence.blocked(context.study());
        if (!blocked.isEmpty()) throw new IOException("Pan-UKB outcome lookup not sent: " + blocked);
        StudyPublicEvidence.validate(URI.create(DATA));
        StudyPublicEvidence.validate(URI.create(INDEX));
        Files.createDirectories(output.getParent());
        List<Object> audit = new ArrayList<>();
        int retained = 0;
        try (TabixReader reader = new TabixReader(DATA, INDEX)) {
            Map<String, Integer> columns = header(reader);
            String population = "EUR";
            for (String required : List.of("chr", "pos", "ref", "alt",
                    "beta_" + population, "se_" + population,
                    "neglog10_pval_" + population,
                    "low_confidence_" + population))
                if (!columns.containsKey(required)) throw new IOException("Pan-UKB file is missing column " + required);
            StringBuilder table = new StringBuilder("variant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n");
            for (Instrument variant : variants) {
                Map<String, Object> rowAudit = new LinkedHashMap<>();
                rowAudit.put("variant_id", variant.id());
                rowAudit.put("source_build", "GRCh38");
                rowAudit.put("source_position", variant.chromosome() + ":" + variant.position());
                try {
                    Projected mapped = projector.project(variant);
                    rowAudit.put("ncbi_projection", mapped.ncbiUrl());
                    rowAudit.put("pan_ukb_query", mapped.chromosome() + ":" + mapped.position() + "-" + mapped.position());
                    String contig = chooseContig(reader, mapped.chromosome());
                    TabixReader.Iterator matches = reader.query(contig + ":" + mapped.position() + "-" + mapped.position());
                    String selected = null;
                    for (String line; matches != null && (line = matches.next()) != null;) {
                        String[] fields = line.split("\t", -1);
                        if (fields.length < columns.size()) continue;
                        if (!fields[columns.get("pos")].equals(Long.toString(mapped.position()))) continue;
                        if (!fields[columns.get("ref")].equalsIgnoreCase(mapped.other())
                                || !fields[columns.get("alt")].equalsIgnoreCase(mapped.effect())) continue;
                        if (selected != null) throw new IOException("Multiple Pan-UKB rows match " + variant.id());
                        selected = line;
                    }
                    if (selected == null) {
                        rowAudit.put("status", "not_found_with_projected_alleles");
                        audit.add(rowAudit); continue;
                    }
                    String[] fields = selected.split("\t", -1);
                    if (Boolean.parseBoolean(fields[columns.get("low_confidence_" + population)])) {
                        rowAudit.put("status", "excluded_low_confidence");
                        audit.add(rowAudit); continue;
                    }
                    double beta = finite(fields[columns.get("beta_" + population)], "beta", variant.id());
                    double se = finite(fields[columns.get("se_" + population)], "standard error", variant.id());
                    double logp = finite(fields[columns.get("neglog10_pval_" + population)], "-log10 p", variant.id());
                    if (!(se > 0) || logp < 0) throw new IOException("Invalid Pan-UKB statistics for " + variant.id());
                    double p = Math.pow(10, -logp);
                    table.append(variant.id()).append('\t').append(variant.chromosome()).append('\t')
                        .append(variant.position()).append('\t').append(variant.effect()).append('\t')
                        .append(variant.other()).append('\t').append(beta).append('\t').append(se).append('\t')
                        .append(p).append('\n');
                    rowAudit.put("status", "retained");
                    rowAudit.put("source_alleles", mapped.other() + ">" + mapped.effect());
                    retained++;
                } catch (IOException | RuntimeException exception) {
                    rowAudit.put("status", "projection_or_lookup_failed");
                    rowAudit.put("reason", String.valueOf(exception.getMessage()));
                }
                audit.add(rowAudit);
            }
            Files.writeString(output, table);
        }
        Map<String, Object> provenance = new LinkedHashMap<>(descriptor());
        provenance.put("retrieved_at", Instant.now().toString());
        provenance.put("requested_instruments", variants.size());
        provenance.put("retained_instruments", retained);
        provenance.put("coordinate_policy", "Pan-UKB GRCh37 alleles were matched after allele-preserving NCBI SPDI projection; output coordinates are the original verified GRCh38 exposure coordinates");
        provenance.put("instrument_audit", audit);
        provenance.put("output_sha256", FollowupSupport.hash(output));
        StudyPublicEvidence.write(output.resolveSibling(output.getFileName() + ".source.yaml"), provenance);
        if (retained == 0) throw new IOException("No Pan-UKB EUR BMI effects matched the projected exposure instruments");
        return provenance;
    }

    static Projected project(Instrument variant) throws IOException {
        String grch38 = refseq(variant.chromosome(), true);
        String grch37 = refseq(variant.chromosome(), false);
        String spdi = grch38 + ":" + (variant.position() - 1) + ":"
            + variant.other().toUpperCase(Locale.ROOT) + ":"
            + variant.effect().toUpperCase(Locale.ROOT);
        String url = NCBI + URLEncoder.encode(spdi, StandardCharsets.UTF_8)
            .replace("+", "%20") + "/all_equivalent_contextual";
        try {
            synchronized (PanUkbbGwasOutcome.class) {
                long remaining = 350_000_000L - (System.nanoTime() - lastNcbiRequest);
                if (remaining > 0) java.util.concurrent.TimeUnit.NANOSECONDS.sleep(remaining);
                lastNcbiRequest = System.nanoTime();
            }
            PublicResourceDownloader.Response response = PublicResourceDownloader.responseStrict(URI.create(url), 2 * 1024 * 1024);
            if (response.status() < 200 || response.status() >= 300)
                throw new IOException("NCBI SPDI projection returned HTTP " + response.status());
            return parseProjection(variant, grch37, url, response.body());
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IOException("NCBI SPDI projection interrupted", exception);
        }
    }

    static Projected parseProjection(Instrument variant, String grch37,
            String url, String json) throws IOException {
        Map<String, Object> root = map(SimpleJson.parse(json));
        Projected result = null;
        for (Object raw : list(map(root.get("data")).get("spdis"))) {
            Map<String, Object> value = map(raw);
            if (!grch37.equals(text(value.get("seq_id")))) continue;
            long position = ((Number) value.get("position")).longValue() + 1;
            String other = text(value.get("deleted_sequence")).toUpperCase(Locale.ROOT);
            String effect = text(value.get("inserted_sequence")).toUpperCase(Locale.ROOT);
            if (!other.matches("[ACGT]") || !effect.matches("[ACGT]") || other.equals(effect))
                throw new IOException("NCBI projection is not a biallelic SNP for " + variant.id());
            if (result != null) throw new IOException("NCBI returned multiple GRCh37 contextual alleles for " + variant.id());
            result = new Projected(variant.chromosome(), position, effect, other, url);
        }
        if (result != null) return result;
        throw new IOException("NCBI did not return a unique GRCh37 contextual allele for " + variant.id());
    }

    private static Map<String, Integer> header(TabixReader reader) throws IOException {
        for (String line; (line = reader.readLine()) != null;) {
            if (line.startsWith("##")) continue;
            String[] fields = line.replaceFirst("^#", "").split("\t", -1);
            Map<String, Integer> result = new LinkedHashMap<>();
            for (int i = 0; i < fields.length; i++) result.put(fields[i], i);
            if (result.containsKey("chr") && result.containsKey("pos")) return result;
            throw new IOException("Unexpected Pan-UKB header");
        }
        throw new IOException("Pan-UKB file has no header");
    }

    private static String chooseContig(TabixReader reader, String chromosome) throws IOException {
        if (reader.getChromosomes().contains(chromosome)) return chromosome;
        if (reader.getChromosomes().contains("chr" + chromosome)) return "chr" + chromosome;
        throw new IOException("Pan-UKB index has no chromosome " + chromosome);
    }

    private static double finite(String value, String label, String id) throws IOException {
        try {
            double result = Double.parseDouble(value);
            if (!Double.isFinite(result)) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IOException("Nonfinite Pan-UKB " + label + " for " + id, exception);
        }
    }

    private static String refseq(String chromosome, boolean grch38) {
        String chr = chromosome.replaceFirst("^chr", "").toUpperCase(Locale.ROOT);
        int number;
        try { number = Integer.parseInt(chr); }
        catch (NumberFormatException exception) {
            if (chr.equals("X")) return grch38 ? "NC_000023.11" : "NC_000023.10";
            throw new IllegalArgumentException("Only autosomes and chromosome X are supported for NCBI projection");
        }
        if (number < 1 || number > 22) throw new IllegalArgumentException("Unsupported chromosome " + chromosome);
        int[] v38 = {0,11,12,12,12,10,12,14,11,12,11,10,12,11,9,10,10,11,10,10,11,9,11};
        int[] v37 = {0,10,11,11,11,9,11,13,10,11,10,9,11,10,8,9,9,10,9,9,10,8,10};
        return String.format(Locale.ROOT, "NC_%06d.%d", number, grch38 ? v38[number] : v37[number]);
    }
}
