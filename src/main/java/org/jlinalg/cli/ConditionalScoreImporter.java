/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.jlinalg.raremetal.RareMetalStudy;
import org.jlinalg.settest.ConditionalScoreStudy;
import org.jlinalg.settest.ScoreVariantKey;
import org.jlinalg.settest.SetTestScoreState;

/** Strict streaming reader for one completed conditional-GWAS export trio.
 * Only the requested score block is retained in memory. */
public final class ConditionalScoreImporter {
    private static final List<String> COVARIANCE_COLUMNS = List.of(
        "variant_i", "variant_j", "score_covariance", "null_model_id",
        "score_covariance_block");

    private ConditionalScoreImporter() { }

    /** Import the requested variants from a summary, covariance and completion
     * manifest written by {@code --conditional-gwas-summary}.
     *
     * <p>Every requested variant must occur in one complete covariance block.
     * A target and condition in separate blocks is rejected because their
     * covariance is unknown, never treated as zero.</p>
     */
    public static ConditionalScoreStudy read(String cohort, Path summary,
            Path covariance, Path manifest, List<String> requestedVariants)
            throws IOException {
        if (summary == null || covariance == null || manifest == null)
            throw new IllegalArgumentException(
                "summary, covariance and score manifest are required");
        summary = summary.toAbsolutePath().normalize();
        covariance = covariance.toAbsolutePath().normalize();
        manifest = manifest.toAbsolutePath().normalize();
        for (Path path : List.of(summary, covariance, manifest))
            if (!Files.isRegularFile(path))
                throw new IOException("conditional score input is absent: " + path);
        ManifestData metadata = manifest(manifest);
        List<ScoreVariantKey> requested = keys(requestedVariants);
        Map<String,Integer> requestedByIdentity = new HashMap<>();
        for (int i = 0; i < requested.size(); i++)
            requestedByIdentity.put(requested.get(i).unorientedIdentity(), i);

        Selected[] selected = new Selected[requested.size()];
        Map<Long,Integer> blockSizes = new TreeMap<>();
        long exported = scanSummary(summary, metadata, requestedByIdentity,
            selected, blockSizes);
        if (exported != metadata.exportedVariants)
            throw new IOException("summary exported-variant count " + exported
                + " differs from manifest " + metadata.exportedVariants);
        if (blockSizes.size() != metadata.covarianceBlocks)
            throw new IOException("summary covariance-block count "
                + blockSizes.size() + " differs from manifest "
                + metadata.covarianceBlocks);
        long expectedBlock = 1;
        for (Map.Entry<Long,Integer> entry : blockSizes.entrySet()) {
            if (entry.getKey() != expectedBlock++)
                throw new IOException("score covariance blocks must be contiguous from 1");
            if (entry.getValue() > metadata.maximumBlockVariants)
                throw new IOException("score covariance block exceeds manifest maximum");
        }
        for (int i = 0; i < selected.length; i++) if (selected[i] == null)
            throw new IOException("requested score variant is absent: "
                + requested.get(i));
        long selectedBlock = selected[0].block;
        for (Selected value : selected) if (value.block != selectedBlock)
            throw new IOException("requested target/conditioning covariance is unavailable: "
                + value.key + " is in block " + value.block + " while "
                + selected[0].key + " is in block " + selectedBlock
                + "; cross-block covariance is unknown, not zero");

        List<String> blockKeys = blockKeys(summary, selectedBlock);
        if (blockKeys.size() != blockSizes.get(selectedBlock))
            throw new IOException("selected covariance block size changed while reading");
        double[] matrix = covariance(covariance, metadata, blockSizes,
            selected, selectedBlock, blockKeys);
        double[] scores = new double[selected.length];
        List<String> sourceKeys = new ArrayList<>();
        for (int i = 0; i < selected.length; i++) {
            scores[i] = selected[i].score;
            sourceKeys.add(selected[i].key.toString());
            double expected = selected[i].variance, observed = matrix[i*selected.length+i];
            if (Math.abs(expected-observed) > 1e-8*Math.max(expected, Math.abs(observed)))
                throw new IOException("score variance/covariance diagonal disagree for "
                    + selected[i].key);
        }
        return new ConditionalScoreStudy(cohort, metadata.nullModelId,
            metadata.conditioningSetId, metadata.contract, sourceKeys,
            new SetTestScoreState(scores, matrix, scores.length));
    }

    private static List<ScoreVariantKey> keys(List<String> values) {
        if (values == null || values.isEmpty())
            throw new IllegalArgumentException("requested score variants are required");
        List<ScoreVariantKey> result = new ArrayList<>(); Set<String> seen = new HashSet<>();
        for (String value : values) {
            ScoreVariantKey key = ScoreVariantKey.parse(value);
            if (!seen.add(key.unorientedIdentity()))
                throw new IllegalArgumentException(
                    "duplicate or oppositely coded requested variant: " + value);
            result.add(key);
        }
        return List.copyOf(result);
    }

    private static long scanSummary(Path path, ManifestData metadata,
            Map<String,Integer> requested, Selected[] selected,
            Map<Long,Integer> blockSizes) throws IOException {
        long exported = 0; ScoreVariantKey previous = null;
        Set<String> keysAtPosition = new HashSet<>();
        try (BufferedReader input = Files.newBufferedReader(path,
                StandardCharsets.UTF_8)) {
            TableHeader header = TableHeader.read(input, path);
            header.require("status"); header.require(ConditionalGwasExport.COLUMNS);
            for (long line = 2; ; line++) {
                String text = input.readLine(); if (text == null) break;
                if (text.isBlank()) continue;
                String[] row = header.row(text, line);
                String value = header.value(row, "score_variant_key");
                if (value.isBlank()) {
                    for (String column : ConditionalGwasExport.COLUMNS)
                        if (!header.value(row, column).isBlank())
                            throw new IOException("partial score fields at line "
                                + line + " in " + path);
                    continue;
                }
                if (!header.value(row, "status").equals("ok"))
                    throw new IOException("scored row is not status=ok at line "
                        + line + " in " + path);
                ScoreVariantKey key;
                try { key = ScoreVariantKey.parse(value); }
                catch (IllegalArgumentException failure) {
                    throw new IOException("invalid score key at line " + line
                        + " in " + path + ": " + failure.getMessage(), failure);
                }
                if (previous != null) {
                    int chromosome = RareMetalStudy.compareChromosome(
                        previous.chromosome(), key.chromosome());
                    if (chromosome > 0 || chromosome == 0
                            && previous.position() > key.position())
                        throw new IOException("score rows are not chromosome/position sorted");
                    if (chromosome != 0 || previous.position() != key.position())
                        keysAtPosition.clear();
                }
                if (!keysAtPosition.add(key.toString()))
                    throw new IOException("duplicate score variant key: " + key);
                previous = key;
                requireId(header.value(row, "null_model_id"),
                    metadata.nullModelId, "summary null_model_id", line, path);
                requireId(header.value(row, "conditioning_set_id"),
                    metadata.conditioningSetId, "summary conditioning_set_id",
                    line, path);
                requireId(header.value(row, "n_analyzed"),
                    Long.toString(metadata.nAnalyzed), "summary n_analyzed",
                    line, path);
                double score = number(header.value(row, "score_u"),
                    "score_u", line, path);
                double variance = number(header.value(row, "score_variance"),
                    "score_variance", line, path);
                if (!(variance > 0))
                    throw new IOException("score variance must be positive at line "
                        + line + " in " + path);
                double beta = number(header.value(row, "beta_score"),
                    "beta_score", line, path);
                double standardError = number(header.value(row, "se_score"),
                    "se_score", line, path);
                double p = number(header.value(row, "p_score_normal"),
                    "p_score_normal", line, path);
                double expectedBeta = score/variance,
                    expectedSe = 1/Math.sqrt(variance);
                double expectedP = 2*jdistlib.Normal.cumulative(
                    -Math.abs(score/Math.sqrt(variance)), 0, 1, true, false);
                if (Math.abs(beta-expectedBeta) > 1e-12*Math.max(1,Math.abs(expectedBeta))
                        || Math.abs(standardError-expectedSe) > 1e-12*Math.max(1,expectedSe)
                        || Math.abs(p-expectedP) > Math.max(1e-15,1e-10*expectedP))
                    throw new IOException("redundant score estimates disagree at line "
                        + line + " in " + path);
                if (!header.value(row, "p_score_calibrated").isBlank()
                        || !header.value(row, "calibration_method").equals("normal-score")
                        || !header.value(row, "calibration_status").equals("normal_only"))
                    throw new IOException("unsupported score calibration at line "
                        + line + " in " + path);
                long block = positiveLong(header.value(row,
                    "score_covariance_block"), "score_covariance_block", line,
                    path);
                blockSizes.merge(block, 1, Math::addExact); exported++;
                Integer index = requested.get(key.unorientedIdentity());
                if (index != null) {
                    if (selected[index] != null)
                        throw new IOException("both allele orientations occur for requested variant: "
                            + key);
                    selected[index] = new Selected(key, score, variance, block);
                }
            }
        }
        return exported;
    }

    private static List<String> blockKeys(Path path, long block)
            throws IOException {
        List<String> result = new ArrayList<>(); Set<String> seen = new HashSet<>();
        try (BufferedReader input = Files.newBufferedReader(path,
                StandardCharsets.UTF_8)) {
            TableHeader header = TableHeader.read(input, path);
            for (long line = 2; ; line++) {
                String text = input.readLine(); if (text == null) break;
                if (text.isBlank()) continue;
                String[] row = header.row(text, line);
                String key = header.value(row, "score_variant_key");
                if (!key.isBlank() && positiveLong(header.value(row,
                        "score_covariance_block"), "score_covariance_block",
                        line, path) == block) {
                    if (!seen.add(key))
                        throw new IOException("duplicate score key in selected block: " + key);
                    result.add(key);
                }
            }
        }
        return List.copyOf(result);
    }

    private static double[] covariance(Path path, ManifestData metadata,
            Map<Long,Integer> blockSizes, Selected[] selected,
            long selectedBlock, List<String> blockKeys) throws IOException {
        int n = selected.length, b = blockKeys.size();
        double[] matrix = new double[Math.multiplyExact(n, n)];
        Arrays.fill(matrix, Double.NaN);
        boolean[] selectedPairs = new boolean[Math.multiplyExact(n, n)];
        boolean[] fullBlockPairs = new boolean[Math.multiplyExact(b, b)];
        Map<String,Integer> selectedIndex = new HashMap<>(), blockIndex = new HashMap<>();
        for (int i = 0; i < n; i++) selectedIndex.put(selected[i].key.toString(), i);
        for (int i = 0; i < b; i++) blockIndex.put(blockKeys.get(i), i);
        Map<Long,Long> observedCounts = new TreeMap<>();
        try (BufferedReader input = Files.newBufferedReader(path,
                StandardCharsets.UTF_8)) {
            TableHeader header = TableHeader.read(input, path);
            header.require(COVARIANCE_COLUMNS);
            for (long line = 2; ; line++) {
                String text = input.readLine(); if (text == null) break;
                if (text.isBlank()) continue;
                String[] row = header.row(text, line);
                requireId(header.value(row, "null_model_id"),
                    metadata.nullModelId, "covariance null_model_id", line, path);
                long block = positiveLong(header.value(row,
                    "score_covariance_block"), "score_covariance_block", line,
                    path);
                if (!blockSizes.containsKey(block))
                    throw new IOException("covariance references unknown block "
                        + block + " at line " + line + " in " + path);
                observedCounts.merge(block, 1L, Math::addExact);
                String first = header.value(row, "variant_i"),
                    second = header.value(row, "variant_j");
                try { ScoreVariantKey.parse(first); ScoreVariantKey.parse(second); }
                catch (IllegalArgumentException failure) {
                    throw new IOException("invalid covariance variant key at line "
                        + line + " in " + path + ": " + failure.getMessage(), failure);
                }
                double value = number(header.value(row, "score_covariance"),
                    "score_covariance", line, path);
                if (block == selectedBlock) {
                    Integer fi = blockIndex.get(first), fj = blockIndex.get(second);
                    if (fi == null || fj == null)
                        throw new IOException("selected covariance block references a variant absent from its score block");
                    int low = Math.min(fi, fj), high = Math.max(fi, fj);
                    int pair = low*b+high;
                    if (fullBlockPairs[pair])
                        throw new IOException("duplicate covariance pair in selected block: "
                            + first + " / " + second);
                    fullBlockPairs[pair] = true;
                } else if (selectedIndex.containsKey(first)
                        || selectedIndex.containsKey(second)) {
                    throw new IOException("selected variant occurs in the wrong covariance block");
                }
                Integer i = selectedIndex.get(first), j = selectedIndex.get(second);
                if (i != null && j != null) {
                    int low = Math.min(i, j), high = Math.max(i, j);
                    int pair = low*n+high;
                    if (selectedPairs[pair])
                        throw new IOException("duplicate requested covariance pair: "
                            + first + " / " + second);
                    selectedPairs[pair] = true;
                    matrix[i*n+j] = matrix[j*n+i] = value;
                }
            }
        }
        if (!observedCounts.keySet().equals(blockSizes.keySet()))
            throw new IOException("covariance block coverage differs from the summary");
        for (Map.Entry<Long,Integer> entry : blockSizes.entrySet()) {
            long size = entry.getValue();
            long expected = Math.multiplyExact(size, size+1L)/2L;
            if (observedCounts.get(entry.getKey()) != expected)
                throw new IOException("covariance block " + entry.getKey()
                    + " has " + observedCounts.get(entry.getKey())
                    + " rows; expected complete upper triangle " + expected);
        }
        for (int i = 0; i < b; i++) for (int j = i; j < b; j++)
            if (!fullBlockPairs[i*b+j])
                throw new IOException("missing covariance in selected block for "
                    + blockKeys.get(i) + " / " + blockKeys.get(j));
        for (int i = 0; i < n; i++) for (int j = i; j < n; j++)
            if (!selectedPairs[i*n+j] || !Double.isFinite(matrix[i*n+j]))
                throw new IOException("missing requested covariance for "
                    + selected[i].key + " / " + selected[j].key
                    + "; unknown covariance is not zero");
        return matrix;
    }

    private static ManifestData manifest(Path path) throws IOException {
        Object parsed = SimpleJson.parse(Files.readString(path,
            StandardCharsets.UTF_8));
        if (!(parsed instanceof Map<?,?> raw))
            throw new IOException("conditional score manifest is not a JSON object: " + path);
        Map<String,Object> values = new LinkedHashMap<>();
        for (Map.Entry<?,?> entry : raw.entrySet()) {
            if (!(entry.getKey() instanceof String key))
                throw new IOException("conditional score manifest has a non-string key");
            values.put(key, entry.getValue());
        }
        Object schema = values.get("schema_version");
        if (!(schema instanceof Number number)
                || Double.compare(number.doubleValue(), 1.0) != 0)
            throw new IOException("conditional score manifest requires schema_version 1");
        require(values, "status", "complete");
        require(values, "export_schema", "jlinalg-conditional-score-v1");
        require(values, "null_converged", "true");
        require(values, "effect_allele", "ALT");
        String coverage = text(values, "covariance_coverage");
        if (!coverage.contains("cross-block covariance unavailable, not zero"))
            throw new IOException("manifest does not preserve unknown cross-block covariance");
        String nullId = text(values, "null_model_id"),
            conditionId = text(values, "conditioning_set_id");
        List<String> conditions = strings(values, "conditioning_variants");
        if (conditions.isEmpty() != conditionId.equals("none"))
            throw new IOException("conditioning_set_id disagrees with conditioning_variants");
        long analyzed = positive(text(values, "n_analyzed"), "n_analyzed");
        long exported = nonnegative(text(values, "exported_variants"),
            "exported_variants");
        long blocks = nonnegative(text(values, "covariance_blocks"),
            "covariance_blocks");
        long maximum = positive(text(values, "maximum_block_variants"),
            "maximum_block_variants");
        if (exported < 1 || blocks < 1 || blocks > Integer.MAX_VALUE
                || maximum > 512)
            throw new IOException("invalid score manifest variant/block counts");
        // Presence is part of the linked-file provenance contract. Relocated
        // copies are allowed; UUIDs and row-level joins establish identity.
        text(values, "summary_file"); text(values, "covariance_file");
        ConditionalScoreStudy.ModelContract contract = new ConditionalScoreStudy.ModelContract(
            text(values, "export_schema"), text(values, "genome_build"),
            text(values, "model"), text(values, "family"),
            text(values, "formula"), strings(values, "null_covariate_columns"),
            text(values, "phenotype_transform"),
            nullableText(values, "case_value"),
            nullableText(values, "control_value"),
            text(values, "analysis_sample"),
            text(values, "genotype_coding"),
            text(values, "missing_genotypes"),
            text(values, "score_covariance_scale"),
            text(values, "variance_method"),
            text(values, "relatedness_adjustment"),
            nullableText(values, "ties"), text(values, "tail_calibration"),
            conditions);
        return new ManifestData(nullId, conditionId, contract, analyzed,
            exported, (int)blocks, (int)maximum);
    }

    private static void require(Map<String,Object> values, String key,
            String expected) throws IOException {
        String actual = text(values, key);
        if (!actual.equals(expected))
            throw new IOException("conditional score manifest " + key
                + " must be " + expected + "; found " + actual);
    }
    private static String text(Map<String,Object> values, String key)
            throws IOException {
        Object value = values.get(key);
        if (!(value instanceof String result) || result.isBlank())
            throw new IOException("conditional score manifest requires string " + key);
        return result;
    }
    private static String nullableText(Map<String,Object> values, String key)
            throws IOException {
        Object value = values.get(key);
        if (value == null) return "";
        if (!(value instanceof String result))
            throw new IOException("conditional score manifest " + key
                + " must be a string or null");
        return result;
    }
    private static List<String> strings(Map<String,Object> values, String key)
            throws IOException {
        Object value = values.get(key);
        if (!(value instanceof List<?> items))
            throw new IOException("conditional score manifest requires array " + key);
        List<String> result = new ArrayList<>();
        for (Object item : items) {
            if (!(item instanceof String text) || text.isBlank())
                throw new IOException("conditional score manifest " + key
                    + " requires nonblank strings");
            result.add(text);
        }
        return List.copyOf(result);
    }
    private static long positive(String value, String field) throws IOException {
        long result = nonnegative(value, field);
        if (result < 1) throw new IOException(field + " must be positive");
        return result;
    }
    private static long nonnegative(String value, String field)
            throws IOException {
        try {
            long result = Long.parseLong(value);
            if (result < 0) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException failure) {
            throw new IOException("manifest " + field
                + " must be a nonnegative integer", failure);
        }
    }
    private static void requireId(String actual, String expected, String field,
            long line, Path path) throws IOException {
        if (!actual.equals(expected))
            throw new IOException(field + " mismatch at line " + line
                + " in " + path);
    }
    private static double number(String value, String field, long line,
            Path path) throws IOException {
        try {
            double result = Double.parseDouble(value);
            if (!Double.isFinite(result)) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException failure) {
            throw new IOException("finite " + field + " required at line "
                + line + " in " + path, failure);
        }
    }
    private static long positiveLong(String value, String field, long line,
            Path path) throws IOException {
        try {
            long result = Long.parseLong(value);
            if (result < 1) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException failure) {
            throw new IOException("positive integer " + field
                + " required at line " + line + " in " + path, failure);
        }
    }

    private record Selected(ScoreVariantKey key, double score,
            double variance, long block) { }
    private record ManifestData(String nullModelId, String conditioningSetId,
            ConditionalScoreStudy.ModelContract contract, long nAnalyzed,
            long exportedVariants, int covarianceBlocks,
            int maximumBlockVariants) { }

    private record TableHeader(List<String> columns, Map<String,Integer> index,
            char delimiter, Path path) {
        static TableHeader read(BufferedReader input, Path path)
                throws IOException {
            String first = input.readLine();
            if (first == null) throw new IOException("table is empty: " + path);
            char delimiter = path.getFileName().toString().toLowerCase(Locale.ROOT)
                .endsWith(".csv") ? ',' : '\t';
            List<String> columns = DelimitedData.parse(
                first.replaceFirst("^\\uFEFF", ""), delimiter, 1, path);
            Map<String,Integer> index = new LinkedHashMap<>();
            for (int i = 0; i < columns.size(); i++) {
                String name = columns.get(i).trim();
                if (name.isBlank() || index.put(name, i) != null)
                    throw new IOException("table column names must be unique and nonblank: " + path);
                columns.set(i, name);
            }
            return new TableHeader(List.copyOf(columns), Map.copyOf(index),
                delimiter, path);
        }
        void require(String column) {
            if (!index.containsKey(column))
                throw new IllegalArgumentException("column is absent: " + column
                    + " in " + path);
        }
        void require(List<String> required) { for (String value : required) require(value); }
        String[] row(String text, long line) throws IOException {
            List<String> values = DelimitedData.parse(text, delimiter, line, path);
            if (values.size() != columns.size())
                throw new IOException("expected " + columns.size()
                    + " fields but found " + values.size() + " at line "
                    + line + " in " + path);
            return values.toArray(String[]::new);
        }
        String value(String[] row, String column) { return row[index.get(column)]; }
    }
}
