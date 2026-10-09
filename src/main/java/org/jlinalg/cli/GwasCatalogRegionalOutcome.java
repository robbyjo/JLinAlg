/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;

import htsjdk.tribble.readers.TabixReader;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
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

/** Indexed regional access to GWAS Catalog harmonized GRCh38 summary statistics. */
final class GwasCatalogRegionalOutcome {
    static final String DOCUMENTATION = "https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/";
    static final String ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/";
    private static final Pattern HREF = Pattern.compile("href=\"([A-Za-z0-9._-]+)\"");

    record Result(int retainedVariants, String dataUrl, String indexUrl,
            Map<String, Object> metadata, Map<String, Object> provenance) { }
    private record Row(String id, String chromosome, long position,
            String effect, String other, double beta, double se, double p) { }

    private GwasCatalogRegionalOutcome() { }

    static Result extract(Context context, String accession, String ancestry,
            String chromosome, long start, long end, Path output) throws IOException {
        if (!accession.matches("GCST[0-9]{6,}")) throw new IllegalArgumentException("Invalid GWAS Catalog accession");
        if (start < 1 || end < start || end - start > 10_000_000L)
            throw new IllegalArgumentException("GWAS Catalog regional query must be a valid interval no wider than 10 Mb");
        String blocked = StudyPublicEvidence.blocked(context.study());
        if (!blocked.isEmpty()) throw new IOException("GWAS Catalog regional query not sent: " + blocked);
        String directory = ROOT + bin(accession) + "/" + accession + "/harmonised/";
        String listing = fetchText(directory, 2 * 1024 * 1024);
        Set<String> files = links(listing);
        List<String> dataFiles = files.stream().filter(name -> name.endsWith(".h.tsv.gz") && !name.endsWith(".f.tsv.gz"))
            .sorted().toList();
        if (dataFiles.isEmpty()) throw new IOException(accession + " has no harmonized regional summary-statistics file");
        String preferred = accession + ".h.tsv.gz";
        String dataFile = dataFiles.contains(preferred) ? preferred
            : dataFiles.size() == 1 ? dataFiles.get(0)
            : throwIo(accession + " has multiple harmonized files; an explicit source adapter is required");
        String indexFile = dataFile + ".tbi", metadataFile = dataFile + "-meta.yaml";
        if (!files.contains(indexFile)) throw new IOException(accession + " harmonized results have no tabix index");
        if (!files.contains(metadataFile)) throw new IOException(accession + " harmonized results have no metadata file");
        String dataUrl = directory + dataFile, indexUrl = directory + indexFile, metadataUrl = directory + metadataFile;
        validate(URI.create(dataUrl)); validate(URI.create(indexUrl)); validate(URI.create(metadataUrl));
        Map<String, Object> metadata = map(new Yaml().load(fetchText(metadataUrl, 2 * 1024 * 1024)));
        if (!accession.equals(text(metadata.get("gwas_id"))) || !Boolean.TRUE.equals(metadata.get("is_harmonised"))
                || !Boolean.TRUE.equals(metadata.get("is_sorted")) || !"GRCh38".equals(metadata.get("genome_assembly")))
            throw new IOException("GWAS Catalog metadata is not sorted harmonized GRCh38 for " + accession);
        verifyAncestry(metadata, ancestry, accession);
        List<Row> rows = new ArrayList<>(); Set<String> duplicates = new LinkedHashSet<>();
        try (TabixReader reader = new TabixReader(dataUrl, indexUrl)) {
            Map<String, Integer> columns = header(reader);
            String contig = reader.getChromosomes().contains(chromosome) ? chromosome
                : reader.getChromosomes().contains("chr" + chromosome) ? "chr" + chromosome : null;
            if (contig == null) throw new IOException(accession + " index has no chromosome " + chromosome);
            TabixReader.Iterator query = reader.query(contig + ":" + start + "-" + end);
            Map<String, Row> unique = new LinkedHashMap<>();
            for (String line; query != null && (line = query.next()) != null;) {
                String[] fields = line.split("\t", -1);
                Row row = parse(fields, columns);
                if (row == null || !row.chromosome().replaceFirst("^chr", "").equals(chromosome)
                        || row.position() < start || row.position() > end) continue;
                if (unique.putIfAbsent(row.id(), row) != null) duplicates.add(row.id());
            }
            duplicates.forEach(unique::remove);
            rows.addAll(unique.values());
        }
        rows.sort(Comparator.comparingLong(Row::position).thenComparing(Row::id));
        if (rows.size() < 2) throw new IOException("Fewer than two unambiguous effect estimates were available in " + accession + " for the requested region");
        Files.createDirectories(output.getParent());
        StringBuilder table = new StringBuilder("variant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n");
        for (Row row : rows) table.append(row.id()).append('\t').append(row.chromosome()).append('\t').append(row.position()).append('\t')
            .append(row.effect()).append('\t').append(row.other()).append('\t').append(row.beta()).append('\t').append(row.se()).append('\t').append(row.p()).append('\n');
        Files.writeString(output, table);
        Map<String, Object> provenance = new LinkedHashMap<>();
        provenance.put("provider", "GWAS Catalog harmonized summary statistics"); provenance.put("accession", accession);
        provenance.put("documentation", DOCUMENTATION); provenance.put("directory", directory);
        provenance.put("summary_statistics", dataUrl); provenance.put("tabix_index", indexUrl); provenance.put("metadata", metadataUrl);
        provenance.put("genome_build", "GRCh38"); provenance.put("ancestry", ancestry);
        provenance.put("region", Map.of("chromosome", chromosome, "start", start, "end", end));
        provenance.put("range_retrieval", "Indexed HTTPS byte ranges; the whole remote file was not downloaded or rehashed");
        provenance.put("catalog_data_file_md5", metadata.getOrDefault("data_file_md5sum", "not_declared"));
        provenance.put("retained_variants", rows.size()); provenance.put("duplicate_rsid_rows_excluded", duplicates.size());
        provenance.put("retrieved_at", Instant.now().toString()); provenance.put("output_sha256", FollowupSupport.hash(output));
        provenance.put("metadata_document", metadata);
        StudyPublicEvidence.write(output.resolveSibling(output.getFileName() + ".source.yaml"), provenance);
        return new Result(rows.size(), dataUrl, indexUrl, metadata, provenance);
    }

    private static Row parse(String[] fields, Map<String, Integer> columns) throws IOException {
        if (fields.length <= columns.values().stream().mapToInt(Integer::intValue).max().orElse(-1)) return null;
        String id = value(fields, columns, "rsid", "variant_id");
        if (!id.matches("rs[1-9][0-9]*")) return null;
        String chromosome = value(fields, columns, "chromosome", "hm_chrom").replaceFirst("^chr", "");
        long position;
        try { position = Long.parseLong(value(fields, columns, "base_pair_location", "hm_pos")); }
        catch (NumberFormatException invalid) { return null; }
        String effect = value(fields, columns, "effect_allele", "hm_effect_allele").toUpperCase(Locale.ROOT);
        String other = value(fields, columns, "other_allele", "hm_other_allele").toUpperCase(Locale.ROOT);
        if (!effect.matches("[ACGT]") || !other.matches("[ACGT]") || effect.equals(other)) return null;
        double beta;
        String rawBeta = value(fields, columns, "beta");
        try {
            if (!rawBeta.isBlank() && !rawBeta.equalsIgnoreCase("NA")) beta = Double.parseDouble(rawBeta);
            else {
                double odds = Double.parseDouble(value(fields, columns, "odds_ratio"));
                if (!(odds > 0)) return null;
                beta = Math.log(odds);
            }
            double se = Double.parseDouble(value(fields, columns, "standard_error"));
            double p = Double.parseDouble(value(fields, columns, "p_value"));
            if (!Double.isFinite(beta) || !Double.isFinite(se) || !(se > 0) || !Double.isFinite(p) || p < 0 || p > 1) return null;
            return new Row(id, chromosome, position, effect, other, beta, se, p);
        } catch (NumberFormatException invalid) { return null; }
    }

    private static Map<String, Integer> header(TabixReader reader) throws IOException {
        for (String line; (line = reader.readLine()) != null;) {
            if (line.startsWith("##")) continue;
            String[] fields = line.replaceFirst("^#", "").split("\t", -1);
            Map<String, Integer> result = new LinkedHashMap<>();
            for (int i = 0; i < fields.length; i++) result.put(fields[i], i);
            for (List<String> alternatives : List.of(List.of("chromosome", "hm_chrom"), List.of("base_pair_location", "hm_pos"),
                    List.of("effect_allele", "hm_effect_allele"), List.of("other_allele", "hm_other_allele"),
                    List.of("rsid", "variant_id"), List.of("standard_error"), List.of("p_value")))
                if (alternatives.stream().noneMatch(result::containsKey)) throw new IOException("GWAS harmonized file is missing " + alternatives);
            if (!result.containsKey("beta") && !result.containsKey("odds_ratio"))
                throw new IOException("GWAS harmonized file has neither beta nor odds_ratio");
            return result;
        }
        throw new IOException("GWAS harmonized file has no header");
    }

    private static String value(String[] fields, Map<String, Integer> columns,
            String... alternatives) {
        for (String name : alternatives) if (columns.containsKey(name)) return fields[columns.get(name)];
        return "";
    }

    private static void verifyAncestry(Map<String, Object> metadata,
            String requested, String accession) throws IOException {
        String normalized = requested.toLowerCase(Locale.ROOT).replaceAll("[^a-z]", "");
        boolean found = false;
        for (Object raw : list(metadata.get("samples"))) for (Object category : list(map(raw).get("sample_ancestry_category"))) {
            String value = text(category).toLowerCase(Locale.ROOT).replaceAll("[^a-z]", "");
            if (value.contains(normalized) || normalized.equals("european") && value.contains("europe")) found = true;
        }
        if (!found) throw new IOException(accession + " metadata does not declare requested ancestry " + requested);
    }

    private static Set<String> links(String html) {
        Set<String> result = new LinkedHashSet<>(); Matcher matcher = HREF.matcher(html);
        while (matcher.find()) result.add(matcher.group(1));
        return result;
    }

    private static String bin(String accession) {
        String digits = accession.substring(4); long value = Long.parseLong(digits);
        long start = ((value - 1) / 1000) * 1000 + 1, end = start + 999;
        String format = "%0" + digits.length() + "d";
        return "GCST" + String.format(Locale.ROOT, format, start) + "-GCST" + String.format(Locale.ROOT, format, end);
    }

    private static String fetchText(String url, int maximum) throws IOException {
        validate(URI.create(url));
        try { return PublicResourceDownloader.text(URI.create(url), maximum, "JLinAlg-GWAS-Catalog-region/1"); }
        catch (InterruptedException interrupted) {
            Thread.currentThread().interrupt(); throw new IOException("GWAS Catalog retrieval interrupted", interrupted);
        }
    }

    private static String throwIo(String message) throws IOException { throw new IOException(message); }
}
