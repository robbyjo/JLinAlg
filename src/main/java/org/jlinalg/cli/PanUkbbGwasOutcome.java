/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.yaml.snakeyaml.Yaml;

/** Bounded Pan-UK Biobank instrument lookup with explicit NCBI build projection. */
final class PanUkbbGwasOutcome {
    static final String PROVIDER = "pan-ukb";
    static final String BMI_TRAIT = "21001";
    static final int BMI_EUR_N = 419163;
    static final String BMI_FILE = "continuous-21001-both_sexes-irnt.tsv.bgz";
    static final String BMI_MD5 = "ca9ceaa99c79d510fb0e26a780787416";
    static final String BMI_INDEX_MD5 = "15552b6a7670abdfd74834b5d3b5a3b6";
    static final String OBESITY_TRAIT = "icd10-E66";
    static final int OBESITY_EUR_N = 420531;
    static final String OBESITY_FILE = "icd10-E66-both_sexes.tsv.bgz";
    static final String OBESITY_MD5 = "e226a9283e8f7a9ee1bab459eb2791fc";
    static final String OBESITY_INDEX_MD5 = "2a83e61111eb250497ef905cabe567d3";
    static final String ROOT = "https://pan-ukb-us-east-1.s3.amazonaws.com/";
    static final String DATA = ROOT + "sumstats_flat_files/" + BMI_FILE;
    static final String INDEX = ROOT + "sumstats_flat_files_tabix/" + BMI_FILE + ".tbi";
    static final String DOCUMENTATION = "https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html";
    static final String CHAIN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz";

    record Instrument(String id, String chromosome, long position,
            String effect, String other) { }
    record Profile(String key, String trait, Set<String> aliases,
            String description, String file, String md5, String indexMd5,
            int sampleSize, int cases, int controls, String effectScale,
            String selectionNote) { }
    record Projected(String chromosome, long position, String effect,
            String other, String audit) { }
    @FunctionalInterface interface Projector {
        Projected project(Instrument variant) throws IOException;
    }

    private PanUkbbGwasOutcome() { }

    private static final Profile BMI = new Profile("21001-eur", BMI_TRAIT,
        Set.of("21001", "bmi", "bodymassindex"),
        "Body mass index (BMI), inverse-rank normal transformed", BMI_FILE,
        BMI_MD5, BMI_INDEX_MD5, BMI_EUR_N, -1, -1, "standard-deviation units",
        "UK Biobank-only EUR outcome avoids FHS exposure-sample overlap; do not substitute a UKB-plus-GIANT meta-analysis without reassessing overlap");
    private static final Profile OBESITY = new Profile("icd10-E66-eur",
        OBESITY_TRAIT, Set.of("obesity", "e66", "icd10e66"),
        "ICD-10 E66 obesity", OBESITY_FILE, OBESITY_MD5,
        OBESITY_INDEX_MD5, OBESITY_EUR_N, 15917, 404614, "log odds",
        "Pan-UKB ICD-10 E66 provides a UK Biobank-only obesity outcome independent of the FHS irQTL exposure cohort; GIANT remains a possible alternative only after explicit cohort-overlap review");

    static boolean supports(String trait, String ancestry) {
        try { profile(trait, ancestry); return true; }
        catch (IllegalArgumentException exception) { return false; }
    }

    static Map<String, Object> descriptor() {
        return descriptor(BMI);
    }

    static Profile profile(String trait, String ancestry) {
        if (!"EUR".equalsIgnoreCase(ancestry))
            throw new IllegalArgumentException("Pan-UKB automatic outcomes currently require EUR ancestry");
        String normalized = trait == null ? "" : trait.toLowerCase(Locale.ROOT)
            .replaceAll("[^a-z0-9]+", "");
        for (Profile profile : List.of(BMI, OBESITY))
            if (profile.aliases().contains(normalized)) return profile;
        throw new IllegalArgumentException("Unsupported Pan-UKB outcome trait: " + trait);
    }

    static Map<String, Object> descriptor(Profile profile) {
        Map<String, Object> result = new LinkedHashMap<>(Map.ofEntries(
            Map.entry("provider", PROVIDER),
            Map.entry("trait", profile.trait()),
            Map.entry("description", profile.description()),
            Map.entry("ancestry", "EUR"),
            Map.entry("sample_size", profile.sampleSize()),
            Map.entry("source_genome_build", "GRCh37"),
            Map.entry("effect_allele", "ALT"),
            Map.entry("effect_scale", profile.effectScale()),
            Map.entry("summary_statistics", ROOT + "sumstats_flat_files/" + profile.file()),
            Map.entry("tabix_index", ROOT + "sumstats_flat_files_tabix/" + profile.file() + ".tbi"),
            Map.entry("summary_statistics_md5", profile.md5()),
            Map.entry("tabix_index_md5", profile.indexMd5()),
            Map.entry("documentation", DOCUMENTATION),
            Map.entry("resource_policy", "complete files installed in the shared cache; analysis reads local files only"),
            Map.entry("selection_note", profile.selectionNote())));
        if (profile.cases() >= 0) {
            result.put("cases", profile.cases());
            result.put("controls", profile.controls());
        }
        return result;
    }

    static Map<String, Object> extract(Context context, List<Instrument> variants,
            Path output) throws IOException {
        return extract(context, BMI, variants, output);
    }

    static Map<String, Object> extract(Context context, Profile profile,
            List<Instrument> variants, Path output) throws IOException {
        String dataUrl = ROOT + "sumstats_flat_files/" + profile.file();
        String indexUrl = ROOT + "sumstats_flat_files_tabix/" + profile.file() + ".tbi";
        PublicGeneticResourceCache.Installed data =
            PublicGeneticResourceCache.acquire(context, "pan-ukb",
                profile.key(), profile.file(), dataUrl, 10L * 1024 * 1024 * 1024,
                new PublicGeneticResourceCache.ExpectedChecksum("MD5", profile.md5()));
        PublicGeneticResourceCache.Installed index =
            PublicGeneticResourceCache.acquire(context, "pan-ukb",
                profile.key(), profile.file() + ".tbi", indexUrl,
                1024L * 1024 * 1024,
                new PublicGeneticResourceCache.ExpectedChecksum("MD5",
                    profile.indexMd5()));
        Path chain = chain(context);
        LiftOver liftOver = new LiftOver(chain.toFile());
        return extract(context, profile, variants, output,
            variant -> project(variant, liftOver), data, index, chain);
    }

    private static Map<String, Object> extract(Context context,
            Profile profile, List<Instrument> variants, Path output, Projector projector,
            PublicGeneticResourceCache.Installed data,
            PublicGeneticResourceCache.Installed index, Path chain)
            throws IOException {
        if (variants.isEmpty()) throw new IllegalArgumentException("No exposure instruments were supplied for outcome lookup");
        Files.createDirectories(output.getParent());
        List<Object> audit = new ArrayList<>();
        int retained = 0;
        try (TabixReader reader = new TabixReader(data.path().toString(),
                index.path().toString())) {
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
                    rowAudit.put("local_build_projection", mapped.audit());
                    rowAudit.put("local_pan_ukb_lookup", mapped.chromosome()
                        + ":" + mapped.position() + "-" + mapped.position());
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
        Map<String, Object> provenance = new LinkedHashMap<>(descriptor(profile));
        provenance.put("retrieved_at", Instant.now().toString());
        provenance.put("local_summary_statistics", data.path().toString());
        provenance.put("local_summary_statistics_sha256", data.sha256());
        provenance.put("local_tabix_index", index.path().toString());
        provenance.put("local_tabix_index_sha256", index.sha256());
        provenance.put("local_liftover_chain", chain.toString());
        provenance.put("local_liftover_chain_sha256",
            FollowupSupport.hash(chain));
        provenance.put("requested_instruments", variants.size());
        provenance.put("retained_instruments", retained);
        provenance.put("coordinate_policy", "Pan-UKB GRCh37 alleles were matched after local hg38-to-hg19 chain projection with reverse-strand allele complementation; output coordinates are the original verified GRCh38 exposure coordinates");
        provenance.put("instrument_audit", audit);
        provenance.put("output_sha256", FollowupSupport.hash(output));
        StudyPublicEvidence.write(output.resolveSibling(output.getFileName() + ".source.yaml"), provenance);
        if (retained == 0) throw new IOException("No Pan-UKB EUR " + profile.description()
            + " effects matched the projected exposure instruments");
        return provenance;
    }

    static Projected project(Instrument variant, LiftOver liftOver)
            throws IOException {
        if (variant.position() < 1 || variant.position() > Integer.MAX_VALUE)
            throw new IOException("Variant position is outside local liftover range: "
                + variant.id());
        String chromosome = "chr" + variant.chromosome().replaceFirst("^chr", "");
        Interval mapped = liftOver.liftOver(new Interval(chromosome,
            (int) variant.position(), (int) variant.position(), false,
            variant.id()), 1.0);
        if (mapped == null || mapped.length() != 1)
            throw new IOException("Local hg38-to-hg19 projection failed for "
                + variant.id());
        String effect = variant.effect().toUpperCase(Locale.ROOT);
        String other = variant.other().toUpperCase(Locale.ROOT);
        if (mapped.isNegativeStrand()) {
            effect = complement(effect); other = complement(other);
        }
        return new Projected(mapped.getContig().replaceFirst("^chr", ""),
            mapped.getStart(), effect, other,
            "hg38ToHg19 chain; strand="
                + (mapped.isNegativeStrand() ? "reverse" : "forward"));
    }

    private static String complement(String allele) throws IOException {
        return switch (allele) {
            case "A" -> "T"; case "T" -> "A";
            case "C" -> "G"; case "G" -> "C";
            default -> throw new IOException(
                "Local liftover supports biallelic SNP alleles only");
        };
    }

    private static synchronized Path chain(Context context)
            throws IOException {
        PublicGeneticResourceCache.Installed compressed =
            PublicGeneticResourceCache.acquire(context, "liftover",
                "hg38-to-hg19", "hg38ToHg19.over.chain.gz", CHAIN,
                64L * 1024 * 1024, null);
        Path chain = compressed.path().resolveSibling(
            "hg38ToHg19.over.chain");
        Path manifest = chain.resolveSibling(
            "hg38ToHg19.over.chain.derived.yaml");
        if (Files.isRegularFile(chain) && Files.isRegularFile(manifest)) {
            Map<String, Object> installed;
            try (var reader = Files.newBufferedReader(manifest)) {
                installed = map(new Yaml().load(reader));
            }
            if (compressed.sha256().equals(installed.get("source_sha256"))
                    && FollowupSupport.hash(chain).equals(
                        installed.get("sha256"))) return chain;
            throw new IOException("Local liftover chain derivation hash mismatch: "
                + chain);
        }
        if (Files.exists(chain) || Files.exists(manifest))
            throw new IOException("Incomplete local liftover chain derivation: "
                + chain);
        Path stage = Files.createTempFile(chain.getParent(), ".chain-", ".part");
        try (var input = new GZIPInputStream(
                    Files.newInputStream(compressed.path()));
                var output = Files.newOutputStream(stage)) {
            input.transferTo(output);
        }
        Files.move(stage, chain);
        StudyPublicEvidence.write(manifest, Map.of(
            "schema_version", 1, "source", compressed.path().toString(),
            "source_sha256", compressed.sha256(), "sha256",
            FollowupSupport.hash(chain), "method", "gzip decompression"));
        return chain;
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

}
