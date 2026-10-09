/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Instant;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/** Project-scoped, provenance-recorded public reference resolution. */
final class StudyResourceManager {
    private static final Pattern GENCODE_RELEASE = Pattern.compile(
        "(?i)(?:GENCODE[_-]?)?v?(\\d+)(?:lift37)?");
    private static final long MAX_GENCODE_BYTES = 2L * 1024 * 1024 * 1024;

    record Annotation(Path table, Path manifest, String release) { }
    record LdReference(Path directory, String database, String population) { }

    private StudyResourceManager() { }

    /** Fetch only a declared human GRCh38 GENCODE release; never guess a build. */
    static Annotation gencodeCoordinates(Path project, Map<String, Object> study,
            Map<String, Object> dataset) throws IOException, InterruptedException {
        String releaseText = String.valueOf(dataset.getOrDefault("annotation_release", ""));
        Matcher matcher = GENCODE_RELEASE.matcher(releaseText);
        if (!matcher.matches() || releaseText.toLowerCase(java.util.Locale.ROOT).endsWith("lift37")
                || Integer.parseInt(matcher.group(1)) < 20
                || !List.of("human", "Homo sapiens", "9606").contains(String.valueOf(dataset.get("species")))
                || !"GRCh38".equals(StudyFeatureAnnotations.build(String.valueOf(dataset.getOrDefault("genome_build", ""))))) return null;
        String release = matcher.group(1);
        Path directory = resourceRoot(project, study).resolve("gencode-coordinates").resolve("v" + release);
        Path gtf = directory.resolve("gencode.v" + release + ".annotation.gtf.gz");
        Path manifest = directory.resolve("manifest.properties");
        if (Files.isRegularFile(gtf) && Files.isRegularFile(manifest)) {
            Properties installed = new Properties();
            try (var reader = Files.newBufferedReader(manifest)) { installed.load(reader); }
            if (!sha256(gtf).equals(installed.getProperty("source_sha256"))) throw new IOException("GENCODE GTF checksum does not match its manifest");
            return new Annotation(gtf, manifest, "GENCODE_v" + release);
        }
        if (!automatic(study)) return null;
        Files.createDirectories(directory.getParent());
        if (Files.exists(directory)) throw new IOException("Partial GENCODE coordinate directory requires review: " + directory);
        Path stage = Files.createTempDirectory(directory.getParent(), ".gencode-coordinates-");
        String url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_" + release + "/gencode.v" + release + ".annotation.gtf.gz";
        String rationale = "Resolve exact transcript coordinates using the dataset-declared human GENCODE release " + release + " on GRCh38 before regional genetic or sequence follow-up; no build or release is inferred.";
        StudyPublicEvidence.write(stage.resolve("acquisition-intent.yaml"), Map.of("scientific_rationale", rationale, "source_url", url, "started_at", Instant.now().toString()));
        boolean complete = false;
        try {
            Path source = stage.resolve(gtf.getFileName());
            PublicResourceDownloader.downloadStrict(URI.create(url), source, MAX_GENCODE_BYTES);
            Properties properties = new Properties();
            properties.setProperty("release", releaseText); properties.setProperty("genome_build", "GRCh38");
            properties.setProperty("source_url", url); properties.setProperty("source_sha256", sha256(source));
            properties.setProperty("downloaded_at", Instant.now().toString());
            properties.setProperty("citation", "https://www.gencodegenes.org/human/");
            properties.setProperty("scientific_rationale", rationale);
            try (var writer = Files.newBufferedWriter(stage.resolve("manifest.properties"))) { properties.store(writer, "Public reference download; no study feature identifiers transmitted"); }
            move(stage, directory);
            complete = true;
        } finally { if (!complete && Files.exists(stage)) deleteTree(stage); }
        return new Annotation(gtf, manifest, "GENCODE_v" + release);
    }

    static Annotation gencode(Path project, Map<String, Object> study,
            List<Map<String, Object>> datasets) throws IOException,
            InterruptedException {
        String release = null;
        for (Map<String, Object> dataset : datasets) {
            Object value = dataset.get("annotation_release");
            if (value == null) continue;
            Matcher matcher = GENCODE_RELEASE.matcher(String.valueOf(value));
            if (!matcher.matches()) continue;
            String species = String.valueOf(dataset.getOrDefault("species", "human"));
            if (!List.of("human", "Homo sapiens", "9606").contains(species)
                    || String.valueOf(value).toLowerCase(java.util.Locale.ROOT).endsWith("lift37"))
                return null;
            String assembly = StudyFeatureAnnotations.build(
                String.valueOf(dataset.getOrDefault("genome_build", "GRCh38")));
            // This provider currently serves modern human GRCh38 metadata only.
            if (!assembly.equals("GRCh38") || Integer.parseInt(matcher.group(1)) < 20)
                return null;
            if (release != null && !release.equals(matcher.group(1)))
                throw new IOException("Selected datasets declare conflicting "
                    + "GENCODE releases");
            release = matcher.group(1);
        }
        if (release == null) return null;
        Path destination = resourceRoot(project, study).resolve("gencode")
            .resolve("v" + release).toAbsolutePath().normalize();
        Path table = destination.resolve("feature-to-gene.tsv");
        Path manifest = destination.resolve("manifest.properties");
        if (Files.isRegularFile(table) && Files.isRegularFile(manifest)) {
            Properties installed = new Properties();
            try (var reader = Files.newBufferedReader(manifest)) {
                installed.load(reader);
            }
            if (!sha256(table).equals(installed.getProperty("mapping_sha256")))
                throw new IOException("Installed GENCODE mapping hash does not "
                    + "match its manifest: " + table);
            for (String name : List.of("transcript_rankings", "hgnc")) {
                Path source = destination.resolve(installed.getProperty(
                    "source_file." + name, ""));
                if (!Files.isRegularFile(source)
                        || !sha256(source).equals(installed.getProperty(
                            "source_sha256." + name)))
                    throw new IOException("Installed GENCODE " + name
                        + " source hash does not match its manifest: " + source);
            }
            return new Annotation(table, manifest, "GENCODE_v" + release);
        }
        if (!automatic(study)) return null;
        if (Files.exists(destination))
            throw new IOException("GENCODE resource directory is partial or "
                + "conflicting: " + destination);
        Files.createDirectories(destination.getParent());
        Path stage = Files.createTempDirectory(destination.getParent(),
            ".gencode-v" + release + "-");
        boolean complete = false;
        String base = "https://ftp.ebi.ac.uk/pub/databases/gencode/"
            + "Gencode_human/release_" + release + "/gencode.v" + release;
        String rankingsSource = base + ".transcript_rankings.txt.gz";
        String hgncSource = base + ".metadata.HGNC.gz";
        String rationale = "Map the study's transcript identifiers to gene identifiers and symbols using its declared GENCODE release for transparent gene-set enrichment and functional evidence searches.";
        StudyPublicEvidence.write(stage.resolve("acquisition-intent.yaml"), Map.of("scientific_rationale", rationale, "sources", List.of(rankingsSource, hgncSource), "started_at", Instant.now().toString()));
        try {
            Path rankings = stage.resolve("gencode.v" + release
                + ".transcript_rankings.txt.gz");
            Path hgnc = stage.resolve("gencode.v" + release
                + ".metadata.HGNC.gz");
            download(rankingsSource, rankings);
            download(hgncSource, hgnc);
            parseMetadata(rankings, hgnc,
                stage.resolve("feature-to-gene.tsv"));
            Properties properties = new Properties();
            properties.setProperty("resource", "GENCODE");
            properties.setProperty("scientific_rationale", rationale);
            properties.setProperty("release", "GENCODE_v" + release);
            properties.setProperty("species", "human");
            properties.setProperty("genome_build", "GRCh38");
            properties.setProperty("source_url.transcript_rankings",
                rankingsSource);
            properties.setProperty("source_url.hgnc", hgncSource);
            properties.setProperty("source_file.transcript_rankings",
                rankings.getFileName().toString());
            properties.setProperty("source_file.hgnc",
                hgnc.getFileName().toString());
            properties.setProperty("downloaded_at", Instant.now().toString());
            properties.setProperty("source_sha256.transcript_rankings",
                sha256(rankings));
            properties.setProperty("source_sha256.hgnc", sha256(hgnc));
            properties.setProperty("mapping_sha256", sha256(
                stage.resolve("feature-to-gene.tsv")));
            try (BufferedWriter writer = Files.newBufferedWriter(
                    stage.resolve("manifest.properties"))) {
                properties.store(writer,
                    "JLinAlg public annotation resource; no study data transmitted");
            }
            move(stage, destination);
            complete = true;
        } finally {
            if (!complete && Files.exists(stage)) deleteTree(stage);
        }
        return new Annotation(table, manifest, "GENCODE_v" + release);
    }

    static LdReference ld(Path project, Map<String, Object> study,
            java.io.PrintStream output) throws IOException,
            InterruptedException {
        Map<String, Object> resources = map(study.get("resources"));
        Map<String, Object> ld = map(resources.get("ld"));
        String databaseId = String.valueOf(ld.getOrDefault(
            "database", "1000g-phase3"));
        LdDatabaseSpec database = LdDatabaseCatalog.find(databaseId);
        if (database == null)
            throw new IOException("Unknown LD reference database: " + databaseId);
        String population = String.valueOf(ld.getOrDefault(
            "population", "auto"));
        Path location = ld.get("location") == null
            ? resourceRoot(project, study).resolve("ld").resolve(databaseId)
            : resolve(project, String.valueOf(ld.get("location")));
        Path manifest = location.resolve("jlinalg-ld-reference.json");
        if (!Files.isRegularFile(manifest)) {
            boolean enabled = !Boolean.FALSE.equals(ld.get("auto_download"));
            if (!automatic(study) || !enabled) return null;
            Path intent = resourceRoot(project, study).resolve("ld-acquisition-intent.yaml");
            Files.createDirectories(intent.getParent());
            StudyPublicEvidence.write(intent, Map.of("scientific_rationale", "Acquire the declared linkage-disequilibrium reference for instrument clumping or regional fine-mapping; downstream inference must still verify genome-build and ancestry compatibility.", "database", databaseId, "population", population, "started_at", Instant.now().toString()));
            new LdDatabaseInstaller().install(database, location, output);
        }
        return new LdReference(location, databaseId, population);
    }

    static boolean automatic(Map<String, Object> study) {
        return !"offline".equals(map(study.get("resources")).get("network"))
            && "automatic".equalsIgnoreCase(String.valueOf(
            map(study.get("resources")).getOrDefault("downloads", "automatic")));
    }

    static Path resourceRoot(Path project, Map<String, Object> study) {
        return resourceRoot(project, study,
            Path.of(System.getProperty("user.home")));
    }

    static Path resourceRoot(Path project, Map<String, Object> study,
            Path home) {
        Object configured = map(study.get("resources")).get("directory");
        if (configured != null)
            return resolve(project, String.valueOf(configured));
        try {
            Object global = ProjectConfiguration.globalResources(home)
                .get("directory");
            if (global != null)
                return Path.of(String.valueOf(global)).toAbsolutePath()
                    .normalize();
        } catch (IOException failure) {
            throw new IllegalArgumentException(
                "Cannot read global JLinAlg resource configuration", failure);
        }
        return project.resolve("resources").toAbsolutePath().normalize();
    }

    private static Path resolve(Path project, String value) {
        Path path = Path.of(value);
        return (path.isAbsolute() ? path : project.resolve(path))
            .toAbsolutePath().normalize();
    }

    private static void download(String source, Path output)
            throws IOException, InterruptedException {
        PublicResourceDownloader.download(URI.create(source), output,
            MAX_GENCODE_BYTES, "JLinAlg-study-resources/1", null, 0);
    }

    private static void parseMetadata(Path rankings, Path hgnc, Path output)
            throws IOException {
        Map<String, String> symbols = new HashMap<>();
        try (var input = gzipReader(hgnc)) {
            for (String line; (line = input.readLine()) != null;) {
                String[] fields = line.split("\\t", -1);
                if (fields.length >= 2 && !blank(fields[0])
                        && !blank(fields[1]))
                    symbols.put(fields[0], fields[1]);
            }
        }
        Map<String, String> geneSymbols = new HashMap<>();
        try (var input = gzipReader(rankings)) {
            for (String line; (line = input.readLine()) != null;) {
                String[] fields = line.split("\\t", -1);
                if (fields.length < 5) continue;
                String symbol = symbols.get(fields[4]);
                if (!blank(fields[1]) && !blank(symbol))
                    geneSymbols.putIfAbsent(fields[1], symbol);
            }
        }
        long mappings = 0;
        var writtenGenes = new HashSet<String>();
        try (var input = gzipReader(rankings);
                BufferedWriter writer = Files.newBufferedWriter(output)) {
            writer.write("feature_id\tgene_id\tgene_symbol\tfeature_type\n");
            for (String line; (line = input.readLine()) != null;) {
                String[] fields = line.split("\\t", -1);
                if (fields.length < 5) continue;
                String gene = fields[1];
                String feature = fields[4];
                String symbol = symbols.getOrDefault(feature,
                    geneSymbols.getOrDefault(gene, ""));
                if (blank(feature) || blank(gene)) continue;
                writeMapping(writer, feature, gene, symbol, "transcript");
                mappings++;
                String unversioned = feature.replaceFirst("\\.\\d+$", "");
                if (!feature.equals(unversioned)) {
                    writeMapping(writer, unversioned, gene, symbol,
                        "transcript");
                    mappings++;
                }
                if (writtenGenes.add(gene)) {
                    writeMapping(writer, gene, gene, symbol, "gene");
                    mappings++;
                    String unversionedGene = gene.replaceFirst("\\.\\d+$", "");
                    if (!gene.equals(unversionedGene)) {
                        writeMapping(writer, unversionedGene, gene, symbol,
                            "gene");
                        mappings++;
                    }
                }
            }
        }
        if (mappings == 0)
            throw new IOException("GENCODE metadata contained no gene mappings");
    }

    private static java.io.BufferedReader gzipReader(Path path)
            throws IOException {
        return new java.io.BufferedReader(new java.io.InputStreamReader(
            new GZIPInputStream(Files.newInputStream(path), 65536),
            StandardCharsets.UTF_8));
    }

    private static void writeMapping(BufferedWriter writer, String feature,
            String gene, String symbol, String type) throws IOException {
        writer.write(tsv(feature)); writer.write('\t');
        writer.write(tsv(gene.replaceFirst("\\.\\d+$", "")));
        writer.write('\t'); writer.write(tsv(symbol)); writer.write('\t');
        writer.write(type); writer.write('\n');
    }

    private static String tsv(String value) {
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

    private static boolean blank(String value) {
        return value == null || value.isBlank();
    }

    private static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (InputStream input = Files.newInputStream(path)) {
                byte[] buffer = new byte[65536];
                for (int count; (count = input.read(buffer)) >= 0;)
                    if (count > 0) digest.update(buffer, 0, count);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException impossible) {
            throw new IllegalStateException(impossible);
        }
    }

    private static void move(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE);
        } catch (AtomicMoveNotSupportedException unsupported) {
            Files.move(source, target);
        }
    }

    private static void deleteTree(Path directory) throws IOException {
        try (var paths = Files.walk(directory)) {
            for (Path path : paths.sorted(Comparator.reverseOrder()).toList())
                Files.deleteIfExists(path);
        }
    }

    private static Map<String, Object> map(Object value) {
        Map<String, Object> result = new LinkedHashMap<>();
        if (value instanceof Map<?, ?> raw)
            raw.forEach((key, item) -> result.put(String.valueOf(key), item));
        return result;
    }
}
