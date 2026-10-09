/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.yaml.snakeyaml.Yaml;

/** Local, conservative feature annotation. Coordinates are always 1-based closed. */
final class StudyFeatureAnnotations {
    private StudyFeatureAnnotations() { }

    record Mapping(String gene, String symbol, String chromosome, String start,
            String end, String build, String strand, String kind) {
        boolean located() {
            return !chromosome.isBlank() && !start.isBlank()
                && !end.isBlank() && !build.isBlank();
        }
        boolean named() { return !gene.isBlank() || !symbol.isBlank(); }
    }

    static List<String> write(Context context, Map<String, Path> results)
            throws IOException {
        Map<String, Set<Mapping>> mappings = new LinkedHashMap<>();
        Set<String> unresolved = new LinkedHashSet<>();
        List<Object> sources = new ArrayList<>();
        List<Object> requests = new ArrayList<>();
        for (Object value : (List<?>) context.study().get("datasets")) {
            Map<String, Object> dataset = map(value);
            String id = text(dataset.get("id"));
            Path standardized = results.get(id);
            if (standardized == null) continue;
            DelimitedData rows = DelimitedData.read(standardized);
            int featureColumn = rows.column("feature_id");
            Set<String> features = new LinkedHashSet<>();
            for (String[] row : rows.rows()) features.add(row[featureColumn]);
            Map<String, Object> files = map(dataset.get("files"));
            Map<String, Object> columns = map(dataset.get("columns"));
            Map<String, Set<Mapping>> annotation = new LinkedHashMap<>();
            Path source = null;
            String release = text(dataset.get("annotation_release"));
            if (files.get("annotation") != null) {
                source = context.project().resolve(text(files.get("annotation"))).normalize();
                String name = source.getFileName().toString().toLowerCase(Locale.ROOT);
                if (name.endsWith(".gtf") || name.endsWith(".gtf.gz"))
                    readGtf(source, dataset, features, annotation);
                else readTable(source, dataset, columns, true, features, annotation);
            } else if (files.get("results") != null) {
                source = context.project().resolve(text(files.get("results"))).normalize();
                readTable(source, dataset, columns, false, features, annotation);
            }
            // Resource retrieval is keyed by declared release, never by private feature IDs.
            if (files.get("annotation") == null && features.stream().anyMatch(StudyFeatureAnnotations::ensembl)
                    && annotation.values().stream().flatMap(Set::stream).noneMatch(Mapping::located)) {
                try {
                    StudyResourceManager.Annotation reference = StudyResourceManager.gencodeCoordinates(context.project(), context.study(), dataset);
                    if (reference != null) {
                        annotation.clear();
                        source = reference.table();
                        readGtf(source, dataset, features, annotation);
                        Path copy = context.outputDirectory().resolve(id + "-coordinate-reference.properties");
                        Files.copy(reference.manifest(), copy);
                        sources.add(Map.of("dataset", id, "manifest", copy.getFileName().toString(), "sha256", FollowupSupport.hash(copy)));
                    }
                } catch (InterruptedException interrupted) {
                    Thread.currentThread().interrupt();
                    throw new IOException("Annotation download interrupted", interrupted);
                } catch (IOException unavailable) {
                    requests.add(Map.of("dataset", id, "resource_status", "unavailable", "retryable", true,
                        "reason", String.valueOf(unavailable.getMessage()), "required_input", "Supply files.annotation or resolve reference connectivity and retry"));
                }
            }
            if (files.get("annotation") == null && annotation.isEmpty()
                    && features.stream().anyMatch(StudyFeatureAnnotations::ensembl)) {
                try {
                    StudyResourceManager.Annotation reference = StudyResourceManager.gencode(
                        context.project(), context.study(), List.of(dataset));
                    if (reference != null) {
                        source = reference.table();
                        readTable(source, dataset, Map.of(), true, features, annotation);
                        Path copy = context.outputDirectory().resolve(id + "-annotation-reference.properties");
                        Files.copy(reference.manifest(), copy);
                        sources.add(Map.of("dataset", id, "manifest", copy.getFileName().toString(),
                            "sha256", FollowupSupport.hash(copy)));
                    } else if ("offline".equals(map(context.study().get("resources")).get("network"))) {
                        requests.add(Map.of("dataset", id, "resource_status", "offline_policy",
                            "retryable", true, "required_input", "Supply local annotation or repeat without --offline when network access is available"));
                    }
                } catch (InterruptedException interrupted) {
                    Thread.currentThread().interrupt();
                    throw new IOException("Annotation download interrupted", interrupted);
                } catch (IOException unavailable) {
                    requests.add(Map.of("dataset", id, "resource_status", "failed",
                        "retryable", true,
                        "reason", String.valueOf(unavailable.getMessage()), "required_input",
                        "Supply a local files.annotation table/GTF or resolve the connection and repeat the same command"));
                }
            }
            if (source != null) sources.add(Map.of("dataset", id, "path", source.toString(),
                "sha256", FollowupSupport.hash(source), "release", release));
            int missing = 0, missingLocation = 0, ambiguous = 0;
            for (String feature : features) {
                Set<Mapping> found = annotation.get(feature);
                if (found == null || found.isEmpty()) {
                    Mapping identity = identity(feature, dataset);
                    found = identity == null ? Set.of() : Set.of(identity);
                }
                mappings.computeIfAbsent(feature, ignored -> new LinkedHashSet<>()).addAll(found);
                if (found.isEmpty()) { unresolved.add(feature); missing++; }
                if (found.size() > 1) ambiguous++;
                if (found.size() != 1 || !found.iterator().next().located()) missingLocation++;
            }
            if (missing > 0 || missingLocation > 0 || ambiguous > 0)
                requests.add(Map.of("dataset", id, "unmapped_features", missing,
                    "features_without_unique_location", missingLocation,
                    "ambiguous_features", ambiguous,
                    "required_input", "For genomic follow-up, supply files.annotation with feature_id, chromosome, start/end (or position), and genome_build; gene_id/gene_symbol support enrichment. Declare coordinate_system for tables and resolve ambiguous mappings.",
                    "field_paths", List.of("datasets." + id + ".files.annotation",
                        "datasets." + id + ".genome_build",
                        "datasets." + id + ".coordinate_system")));
        }
        Map<String, Integer> counts = new LinkedHashMap<>();
        try (var writer = Files.newBufferedWriter(context.outputDirectory().resolve("feature-annotation.tsv"))) {
            writer.write("feature_id\tgene_id\tgene_symbol\tchromosome\tstart\tend\tgenome_build\tstrand\tfeature_type\tmapping_status\n");
            for (var entry : mappings.entrySet()) {
                Set<Mapping> values = entry.getValue();
                String status = values.isEmpty() ? "unmapped"
                    : values.size() > 1 ? "ambiguous"
                    : unresolved.contains(entry.getKey()) ? "incomplete_across_datasets"
                    : values.iterator().next().located() ? "located" : "missing_location";
                counts.merge(status, 1, Integer::sum);
                Set<Mapping> printed = values.isEmpty()
                    ? Set.of(new Mapping("", "", "", "", "", "", "", "unknown")) : values;
                for (Mapping mapping : printed)
                    writer.write(String.join("\t", entry.getKey(), mapping.gene(), mapping.symbol(),
                        mapping.chromosome(), mapping.start(), mapping.end(), mapping.build(),
                        mapping.strand(), mapping.kind(), status) + "\n");
            }
        }
        Map<String, Object> summary = new LinkedHashMap<>();
        summary.put("schema_version", 1);
        summary.put("coordinate_system", "1-based-closed");
        summary.put("feature_counts", counts);
        summary.put("sources", sources);
        if (counts.getOrDefault("ambiguous", 0) > 0 || counts.getOrDefault("incomplete_across_datasets", 0) > 0)
            requests.add(Map.of("scope", "pooled_features", "required_input",
                "Reconcile feature identity, gene assignment, genome build, and intervals across datasets; inspect every ambiguous/incomplete row in feature-annotation.tsv"));
        summary.put("input_requests", requests);
        summary.put("interpretation", "Gene or transcript intervals describe the annotated feature; they do not establish a causal variant or an sQTL. Ambiguous and incomplete mappings must be reviewed before genomic follow-up.");
        Files.writeString(context.outputDirectory().resolve("annotation-summary.yaml"), new Yaml().dump(summary));
        return List.of("feature-annotation.tsv", "annotation-summary.yaml");
    }

    private static Mapping identity(String feature, Map<String, Object> dataset) {
        String type = text(dataset.get("feature_id_type"));
        if (type.equals("gene_symbol"))
            return new Mapping("", feature, "", "", "", "", "", "gene");
        if (feature.matches("ENSG\\d+(?:\\.\\d+)?"))
            return new Mapping(feature.replaceFirst("\\.\\d+$", ""), "", "", "", "", "", "", "gene");
        return null;
    }

    private static void readTable(Path source, Map<String, Object> dataset,
            Map<String, Object> columns, boolean separate, Set<String> features,
            Map<String, Set<Mapping>> output) throws IOException {
        DelimitedData table = DelimitedData.readResult(source);
        String prefix = separate ? "annotation_" : "";
        int id = column(table, columns, prefix + "feature_id", "feature_id", "transcript_id", "probe_id", "gene_id", "id", "variant", "rsid", "__row_name__");
        if (id < 0) throw new IOException("Annotation source has no feature identifier column: " + source);
        int gene = column(table, columns, prefix + "gene_id", "gene_id", "ensembl_gene_id");
        int symbol = column(table, columns, prefix + "gene_symbol", "gene_symbol", "gene_name", "symbol");
        int chromosome = column(table, columns, prefix + "chromosome", "chromosome", "chr");
        int start = column(table, columns, prefix + "start", "start", "position", "pos", "bp");
        if (start < 0) start = column(table, columns, prefix + "position", "position", "pos", "bp");
        int end = column(table, columns, prefix + "end", "end");
        int strand = column(table, columns, prefix + "strand", "strand");
        int build = column(table, columns, prefix + "genome_build", "genome_build");
        int kind = column(table, columns, prefix + "feature_type", "feature_type");
        String coordinateSystem = text(dataset.get("coordinate_system"));
        for (String[] row : table.rows()) {
            String feature = row[id];
            if (!features.contains(feature)) continue;
            String g = field(row, gene), s = field(row, symbol);
            // Older GENCODE caches used a stable identifier when HGNC was absent.
            // Retain the gene identifier without presenting it as a gene symbol.
            if (s.matches("ENS[GT]\\d+(?:\\.\\d+)?")) s = "";
            String chr = field(row, chromosome), from = field(row, start), to = field(row, end);
            String assembly = build(text(dataset.get("genome_build")));
            String rowBuild = build(field(row, build));
            if (!rowBuild.isEmpty() && !assembly.isEmpty() && !rowBuild.equals(assembly))
                throw new IOException("Conflicting genome builds for " + feature + " in " + source);
            if (!rowBuild.isEmpty()) assembly = rowBuild;
            if (!from.isEmpty()) {
                long first = integer(from, feature), last = to.isEmpty() ? first : integer(to, feature);
                if (coordinateSystem.equals("0-based-half-open")) {
                    if (to.isEmpty()) throw new IOException("0-based-half-open annotation requires end: " + feature);
                    first++;
                }
                if (first < (coordinateSystem.isEmpty() ? 0 : 1) || last < first)
                    throw new IOException("Invalid annotation interval for " + feature);
                from = Long.toString(first); to = Long.toString(last);
                // Unknown coordinate conventions must not be silently assumed.
                if (coordinateSystem.isEmpty()) { chr = ""; from = ""; to = ""; }
            }
            String featureKind = field(row, kind);
            if (featureKind.isEmpty()) featureKind = text(dataset.getOrDefault("feature_id_type", "feature"));
            Mapping mapping = new Mapping(g, s, chr, from, to, assembly, field(row, strand), featureKind);
            if (mapping.named() || !chr.isBlank() || !from.isBlank())
                output.computeIfAbsent(feature, ignored -> new LinkedHashSet<>()).add(mapping);
        }
    }

    private static final Pattern ATTRIBUTE = Pattern.compile("(?:^|;)\\s*([\\w]+)\\s+\"([^\"]*)\"");

    private static void readGtf(Path source, Map<String, Object> dataset,
            Set<String> features, Map<String, Set<Mapping>> output) throws IOException {
        try (var stream = Files.newInputStream(source);
                BufferedReader reader = new BufferedReader(new InputStreamReader(
                    source.toString().endsWith(".gz") ? new GZIPInputStream(stream) : stream,
                    StandardCharsets.UTF_8))) {
            for (String line; (line = reader.readLine()) != null;) {
                if (line.startsWith("#") || line.isBlank()) continue;
                String[] fields = line.split("\t", -1);
                if (fields.length != 9) throw new IOException("Invalid GTF row in " + source);
                if (!Set.of("gene", "transcript").contains(fields[2])) continue;
                Map<String, String> attributes = new LinkedHashMap<>();
                var matcher = ATTRIBUTE.matcher(fields[8]);
                while (matcher.find()) attributes.put(matcher.group(1), matcher.group(2));
                String id = attributes.getOrDefault(fields[2].equals("gene") ? "gene_id" : "transcript_id", "");
                String bare = id.replaceFirst("\\.\\d+$", "");
                if (!features.contains(id) && !features.contains(bare)) continue;
                long first = integer(fields[3], id), last = integer(fields[4], id);
                if (first < 1 || last < first) throw new IOException("Invalid GTF interval: " + id);
                Mapping mapping = new Mapping(attributes.getOrDefault("gene_id", ""),
                    attributes.getOrDefault("gene_name", ""), fields[0], Long.toString(first),
                    Long.toString(last), build(text(dataset.get("genome_build"))), fields[6], fields[2]);
                if (features.contains(id)) output.computeIfAbsent(id, ignored -> new LinkedHashSet<>()).add(mapping);
                // An explicitly versioned study ID never falls back to a different GTF version.
                if (!bare.equals(id) && features.contains(bare))
                    output.computeIfAbsent(bare, ignored -> new LinkedHashSet<>()).add(mapping);
            }
        }
    }

    private static long integer(String value, String feature) throws IOException {
        try { return Long.parseLong(value); }
        catch (NumberFormatException invalid) { throw new IOException("Invalid coordinate for " + feature + ": " + value, invalid); }
    }

    static String build(String value) {
        return switch (value.toLowerCase(Locale.ROOT)) {
            case "hg38", "grch38" -> "GRCh38";
            case "hg19", "grch37" -> "GRCh37";
            default -> value;
        };
    }

    private static int column(DelimitedData table, Map<String, Object> columns,
            String logical, String... defaults) {
        if (columns.containsKey(logical)) return table.column(text(columns.get(logical)));
        for (String name : defaults) if (table.header().contains(name)) return table.column(name);
        return -1;
    }

    private static boolean ensembl(String feature) { return feature.matches("ENS[GT]\\d+(?:\\.\\d+)?"); }
    private static String field(String[] row, int column) {
        if (column < 0) return "";
        String value = row[column].trim();
        return Set.of("", "NA", ".", "null").contains(value) ? "" : value;
    }
    private static String text(Object value) { return value == null ? "" : String.valueOf(value); }
    private static Map<String, Object> map(Object value) {
        Map<String, Object> result = new LinkedHashMap<>();
        if (value instanceof Map<?, ?> raw) raw.forEach((key, item) -> result.put(String.valueOf(key), item));
        return result;
    }
}
