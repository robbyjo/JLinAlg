/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import htsjdk.tribble.readers.TabixReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.cli.StudyGeneticEvidence.Variant;
import org.jlinalg.study.StudyOperationExecutor.Context;

/** Regional LD from the official 30x 1000 Genomes GRCh38 phased callset. */
final class ThousandGenomesGrch38Ld {
    static final String DATABASE = "1000g-30x-grch38";
    static final String RELEASE = "20220422_3202_phased_SNV_INDEL_SV";
    static final String ROOT = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
        + RELEASE + "/";
    static final String PANEL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
        + "integrated_call_samples_v3.20130502.ALL.panel";
    static final String CITATION = "https://doi.org/10.1016/j.cell.2022.08.004";
    static final Set<String> POPULATIONS = Set.of("AFR", "AMR", "EAS", "EUR", "SAS");

    record Result(Path matrix, Path alleles, Map<String, Object> provenance) { }
    private record Dosage(Variant variant, String ref, String alt,
            double[] values, double callRate, double alternateFrequency) { }

    private ThousandGenomesGrch38Ld() { }

    static String ancestry(String population) {
        return switch (population) {
            case "AFR" -> "African";
            case "AMR" -> "Admixed American";
            case "EAS" -> "East Asian";
            case "EUR" -> "European";
            case "SAS" -> "South Asian";
            default -> throw new IllegalArgumentException("Choose an explicit 1000 Genomes superpopulation: AFR, AMR, EAS, EUR, or SAS");
        };
    }

    static Result extract(Context context, String population,
            List<Variant> requested, Path directory) throws IOException {
        population = population.toUpperCase(Locale.ROOT);
        ancestry(population);
        if (requested.isEmpty()) throw new IOException("No overlapping variants for GRCh38 LD extraction");
        String chromosome = requested.get(0).chromosome().replaceFirst("^chr", "");
        if (requested.stream().anyMatch(v -> !v.chromosome().replaceFirst("^chr", "").equals(chromosome)))
            throw new IOException("One regional LD query cannot span chromosomes");
        if (!chromosome.matches("(?:[1-9]|1[0-9]|2[0-2]|X)"))
            throw new IOException("1000 Genomes GRCh38 regional LD supports chromosomes 1-22 and X");
        String filename = chromosome.equals("X")
            ? "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
            : "1kGP_high_coverage_Illumina.chr" + chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz";
        String data = ROOT + filename, index = data + ".tbi";
        PublicGeneticResourceCache.Installed installedPanel =
            PublicGeneticResourceCache.acquire(context, DATABASE, RELEASE,
                "integrated_call_samples_v3.20130502.ALL.panel", PANEL,
                64L * 1024 * 1024, null);
        PublicGeneticResourceCache.Installed installedData =
            PublicGeneticResourceCache.acquire(context, DATABASE, RELEASE,
                filename, data, 100L * 1024 * 1024 * 1024, null);
        PublicGeneticResourceCache.Installed installedIndex =
            PublicGeneticResourceCache.acquire(context, DATABASE, RELEASE,
                filename + ".tbi", index, 4L * 1024 * 1024 * 1024,
                null);
        Files.createDirectories(directory);
        String panelText = Files.readString(installedPanel.path());
        Path panelFile = directory.resolve("1000g-unrelated-samples.panel");
        Files.writeString(panelFile, panelText);
        Set<String> samples = samples(panelText, population);
        if (samples.size() < 50) throw new IOException("Fewer than 50 unrelated 1000 Genomes samples for " + population);
        List<Dosage> found = new ArrayList<>();
        List<Object> exclusions = new ArrayList<>();
        try (TabixReader reader = new TabixReader(
                installedData.path().toString(),
                installedIndex.path().toString())) {
            String[] header = vcfHeader(reader);
            Set<String> available = new LinkedHashSet<>(Arrays.asList(header).subList(9, header.length));
            if (!available.containsAll(samples)) {
                Set<String> missing = new LinkedHashSet<>(samples); missing.removeAll(available);
                throw new IOException("GRCh38 VCF is missing " + missing.size() + " declared panel samples");
            }
            Map<String, Integer> sampleColumns = new LinkedHashMap<>();
            for (int i = 9; i < header.length; i++) sampleColumns.put(header[i], i);
            String contig = reader.getChromosomes().contains("chr" + chromosome) ? "chr" + chromosome
                : reader.getChromosomes().contains(chromosome) ? chromosome : null;
            if (contig == null) throw new IOException("GRCh38 VCF index has no chromosome " + chromosome);
            for (Variant variant : requested) {
                String[] selected = null;
                TabixReader.Iterator query = reader.query(contig + ":" + variant.position() + "-" + variant.position());
                for (String line; query != null && (line = query.next()) != null;) {
                        String[] candidate = line.split("\t", -1);
                        if (candidate.length != header.length || candidate[4].contains(",")) continue;
                        String ref = candidate[3].toUpperCase(Locale.ROOT), alt = candidate[4].toUpperCase(Locale.ROOT);
                        if (!ref.matches("[ACGT]") || !alt.matches("[ACGT]")) continue;
                        boolean id = Arrays.asList(candidate[2].split(";", -1)).contains(variant.id());
                        boolean alleles = (ref.equals(variant.other()) && alt.equals(variant.effect()))
                            || (ref.equals(variant.effect()) && alt.equals(variant.other()));
                        if (!id && !alleles) continue;
                        if (selected != null) throw new IOException("Multiple GRCh38 LD records match " + variant.id());
                        selected = candidate;
                }
                if (selected == null) {
                    exclusions.add(Map.of("variant_id", variant.id(), "reason", "absent_from_1000g_30x_grch38"));
                    continue;
                }
                String[] format = selected[8].split(":", -1); int gt = -1;
                for (int i = 0; i < format.length; i++) if (format[i].equals("GT")) gt = i;
                if (gt < 0) throw new IOException("GRCh38 VCF record lacks GT for " + variant.id());
                double[] dosage = new double[samples.size()];
                int indexInArray = 0, called = 0; double alternateCount = 0;
                for (String sample : samples) {
                    String[] values = selected[sampleColumns.get(sample)].split(":", -1);
                    String genotype = gt < values.length ? values[gt] : ".";
                    String[] alleles = genotype.split("[|/]", -1);
                    if (alleles.length != 2 || !alleles[0].matches("[01]") || !alleles[1].matches("[01]"))
                        dosage[indexInArray++] = Double.NaN;
                    else { int count = Integer.parseInt(alleles[0]) + Integer.parseInt(alleles[1]);
                        dosage[indexInArray++] = count; alternateCount += count; called++; }
                }
                double callRate = called / (double) samples.size();
                double frequency = called == 0 ? Double.NaN : alternateCount / (2.0 * called);
                if (callRate < 0.95) {
                    exclusions.add(Map.of("variant_id", variant.id(), "reason", "call_rate_below_0.95", "call_rate", callRate));
                    continue;
                }
                if (!Double.isFinite(frequency) || Math.min(frequency, 1 - frequency) < 0.01) {
                    exclusions.add(Map.of("variant_id", variant.id(), "reason", "minor_allele_frequency_below_0.01", "alternate_frequency", frequency));
                    continue;
                }
                found.add(new Dosage(variant, selected[3].toUpperCase(Locale.ROOT),
                    selected[4].toUpperCase(Locale.ROOT), dosage, callRate, frequency));
            }
        }
        if (found.size() < 2) throw new IOException("Fewer than two requested variants are usable in 1000 Genomes 30x GRCh38 " + population);
        double[][] correlation = correlation(found);
        Path matrix = directory.resolve("ld-source.tsv"), alleles = directory.resolve("ld-alleles.tsv");
        StringBuilder table = new StringBuilder("variant_id");
        for (Dosage value : found) table.append('\t').append(value.variant().id());
        table.append('\n');
        StringBuilder orientation = new StringBuilder("variant_id\teffect_allele\tother_allele\n");
        for (int row = 0; row < found.size(); row++) {
            Dosage value = found.get(row);
            table.append(value.variant().id());
            for (double r : correlation[row]) table.append('\t').append(r);
            table.append('\n');
            orientation.append(value.variant().id()).append('\t').append(value.alt()).append('\t').append(value.ref()).append('\n');
        }
        Files.writeString(matrix, table); Files.writeString(alleles, orientation);
        Map<String, Object> provenance = new LinkedHashMap<>();
        provenance.put("database", DATABASE); provenance.put("release", RELEASE);
        provenance.put("genome_build", "GRCh38"); provenance.put("population", population);
        provenance.put("ancestry", ancestry(population)); provenance.put("unrelated_samples", samples.size());
        provenance.put("sample_definition", "Samples in the official 2504-sample Phase 3 panel with matching superpopulation; relatives added to the 3202-sample VCF are excluded");
        provenance.put("sample_panel", PANEL); provenance.put("sample_panel_sha256", FollowupSupport.hash(panelFile));
        provenance.put("regional_vcf", data); provenance.put("regional_vcf_index", index);
        provenance.put("local_sample_panel", installedPanel.path().toString());
        provenance.put("local_sample_panel_sha256", installedPanel.sha256());
        provenance.put("local_chromosome_vcf", installedData.path().toString());
        provenance.put("local_chromosome_vcf_sha256", installedData.sha256());
        provenance.put("local_tabix_index", installedIndex.path().toString());
        provenance.put("local_tabix_index_sha256", installedIndex.sha256());
        provenance.put("retrieval_policy", "Complete public chromosome VCF and index installed once; LD extraction performed locally");
        provenance.put("requested_variants", requested.size()); provenance.put("retained_variants", found.size());
        provenance.put("exclusions", exclusions); provenance.put("call_rate_minimum", 0.95);
        provenance.put("minor_allele_frequency_minimum", 0.01);
        provenance.put("method", "Pearson ALT-dosage correlation among unrelated ancestry-panel samples; missing genotypes mean-imputed after call-rate filtering");
        provenance.put("citation", CITATION); provenance.put("retrieved_at", Instant.now().toString());
        provenance.put("matrix_sha256", FollowupSupport.hash(matrix)); provenance.put("alleles_sha256", FollowupSupport.hash(alleles));
        StudyPublicEvidence.write(directory.resolve("ld-reference.yaml"), provenance);
        return new Result(matrix, alleles, provenance);
    }

    private static Set<String> samples(String text, String population) throws IOException {
        Set<String> result = new LinkedHashSet<>();
        String[] lines = text.split("\\R");
        if (lines.length == 0) throw new IOException("Empty 1000 Genomes sample panel");
        String[] header = lines[0].replaceFirst("^#", "").split("\\s+", -1);
        int sample = index(header, "sample"), superPop = index(header, "super_pop");
        for (int i = 1; i < lines.length; i++) {
            if (lines[i].isBlank()) continue;
            String[] row = lines[i].trim().split("\\s+", -1);
            if (row.length <= Math.max(sample, superPop)) throw new IOException("Malformed 1000 Genomes sample panel row " + (i + 1));
            if (row[superPop].equals(population)) result.add(row[sample]);
        }
        return result;
    }

    private static String[] vcfHeader(TabixReader reader) throws IOException {
        for (String line; (line = reader.readLine()) != null;)
            if (line.startsWith("#CHROM\t")) {
                String[] fields = line.substring(1).split("\t", -1);
                if (fields.length < 10 || !fields[8].equals("FORMAT"))
                    throw new IOException("Malformed 1000 Genomes VCF column header");
                return fields;
            }
        throw new IOException("1000 Genomes VCF has no #CHROM header");
    }

    private static int index(String[] header, String name) throws IOException {
        for (int i = 0; i < header.length; i++) if (header[i].equals(name)) return i;
        throw new IOException("1000 Genomes sample panel is missing column " + name);
    }

    private static double[][] correlation(List<Dosage> variants) throws IOException {
        int n = variants.size();
        double[][] centered = new double[n][];
        double[] sumSquares = new double[n];
        for (int i = 0; i < n; i++) {
            double[] source = variants.get(i).values();
            double sum = 0; int count = 0;
            for (double value : source) if (Double.isFinite(value)) { sum += value; count++; }
            double mean = sum / count;
            centered[i] = new double[source.length];
            for (int j = 0; j < source.length; j++) {
                centered[i][j] = Double.isFinite(source[j]) ? source[j] - mean : 0;
                sumSquares[i] += centered[i][j] * centered[i][j];
            }
            if (!(sumSquares[i] > 0)) throw new IOException("Zero-variance LD variant " + variants.get(i).variant().id());
        }
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
            double covariance = 0;
            for (int k = 0; k < centered[i].length; k++) covariance += centered[i][k] * centered[j][k];
            double r = covariance / Math.sqrt(sumSquares[i] * sumSquares[j]);
            if (i == j) r = 1;
            result[i][j] = result[j][i] = r;
        }
        return result;
    }
}
