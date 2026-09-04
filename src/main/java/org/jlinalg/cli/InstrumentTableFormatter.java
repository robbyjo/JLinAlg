/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PushbackInputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/** Streams heterogeneous GWAS/QTL tables into MRInstruments-style columns. */
final class InstrumentTableFormatter {
    static final List<String> COLUMNS = List.of(
        "Phenotype", "SNP", "beta", "se", "eaf", "effect_allele",
        "other_allele", "pval", "units", "ncase", "ncontrol",
        "samplesize", "gene");
    private static final Set<String> REQUIRED = Set.of(
        "SNP", "beta", "se", "effect_allele", "other_allele", "pval");
    private static final Map<String, List<String>> ALIASES = aliases();

    enum EffectScale { AUTO, BETA, ODDS_RATIO, HAZARD_RATIO }

    record Options(char delimiter, Map<String, String> mappings,
            String phenotype, double pThreshold, EffectScale effectScale) {
        Options {
            mappings = Map.copyOf(mappings);
            if (delimiter != '\t' && delimiter != ',')
                throw new IllegalArgumentException("delimiter must be tab or comma");
            if (!Double.isFinite(pThreshold)
                    || pThreshold < 0.0 || pThreshold > 1.0)
                throw new IllegalArgumentException(
                    "p-value threshold must lie in [0, 1]");
        }
    }

    record Result(long inputRows, long outputRows, long aboveThreshold,
            long invalidRows, EffectScale effectScale) { }

    private InstrumentTableFormatter() { }

    static Result format(InputStream source, OutputStream destination,
            Options options) throws IOException {
        try (BufferedReader input = reader(source);
                BufferedWriter output = new BufferedWriter(
                    new OutputStreamWriter(destination, StandardCharsets.UTF_8))) {
            return format(input, output, options);
        }
    }

    static Result format(BufferedReader input, BufferedWriter output,
            Options options) throws IOException {
        String first = input.readLine();
        if (first == null) throw new IOException("instrument table is empty");
        List<String> header = parse(first, options.delimiter(), 1);
        Map<String, Integer> sourceColumns = columns(header);
        Resolved resolved = resolve(header, sourceColumns, options);
        output.write(String.join("\t", COLUMNS));
        output.newLine();

        long inputRows = 0;
        long outputRows = 0;
        long thresholdRows = 0;
        long invalidRows = 0;
        long lineNumber = 1;
        for (String line; (line = input.readLine()) != null;) {
            lineNumber++;
            if (line.isBlank()) continue;
            inputRows++;
            List<String> fields = parse(line, options.delimiter(), lineNumber);
            if (fields.size() != header.size()) throw new IOException(
                "expected " + header.size() + " fields but found "
                    + fields.size() + " at line " + lineNumber);
            Row row;
            try {
                row = row(fields, resolved, options);
            } catch (IllegalArgumentException exception) {
                invalidRows++;
                continue;
            }
            if (row.pValue() > options.pThreshold()) {
                thresholdRows++;
                continue;
            }
            for (int index = 0; index < row.values().size(); index++) {
                if (index > 0) output.write('\t');
                output.write(tsv(row.values().get(index)));
            }
            output.newLine();
            outputRows++;
        }
        output.flush();
        return new Result(inputRows, outputRows, thresholdRows, invalidRows,
            resolved.effectScale());
    }

    private static Resolved resolve(List<String> header,
            Map<String, Integer> sourceColumns, Options options) {
        Map<String, Integer> selected = new LinkedHashMap<>();
        for (String target : COLUMNS) {
            String explicit = options.mappings().get(target);
            Integer index = explicit == null ? find(sourceColumns,
                ALIASES.get(target)) : findExact(header, explicit);
            if (index == null && target.equals("Phenotype")
                    && options.phenotype() != null) index = -1;
            if (index == null && REQUIRED.contains(target)
                    && !target.equals("beta"))
                throw new IllegalArgumentException("cannot map required MR column "
                    + target + "; use --map " + target + "=SOURCE");
            if (index != null) selected.put(target, index);
        }

        EffectScale scale = options.effectScale();
        if (!selected.containsKey("beta")) {
            if (scale == EffectScale.BETA)
                throw new IllegalArgumentException(
                    "cannot map required MR column beta; use --map beta=SOURCE");
            Integer odds = find(sourceColumns,
                List.of("odds_ratio", "oddsratio", "or"));
            Integer hazard = find(sourceColumns,
                List.of("hazard_ratio", "hazardratio", "hr"));
            if (odds != null && (scale == EffectScale.AUTO
                    || scale == EffectScale.ODDS_RATIO)) {
                selected.put("beta", odds);
                scale = EffectScale.ODDS_RATIO;
            } else if (hazard != null && (scale == EffectScale.AUTO
                    || scale == EffectScale.HAZARD_RATIO)) {
                selected.put("beta", hazard);
                scale = EffectScale.HAZARD_RATIO;
            } else throw new IllegalArgumentException(
                "cannot map required MR column beta; map an effect column and "
                    + "set --effect-scale when it is a ratio");
        } else if (scale == EffectScale.AUTO) {
            String source = header.get(selected.get("beta"));
            String normalized = normalize(source);
            scale = normalized.equals("oddsratio") || normalized.equals("or")
                ? EffectScale.ODDS_RATIO
                : normalized.equals("hazardratio") || normalized.equals("hr")
                    ? EffectScale.HAZARD_RATIO : EffectScale.BETA;
        }
        return new Resolved(selected, scale);
    }

    private static Row row(List<String> fields, Resolved resolved,
            Options options) {
        List<String> output = new ArrayList<>(COLUMNS.size());
        for (String target : COLUMNS) {
            Integer index = resolved.columns().get(target);
            String value = index == null ? "" : index < 0
                ? options.phenotype() : fields.get(index).trim();
            output.add(value == null ? "" : value);
        }
        int snp = COLUMNS.indexOf("SNP");
        int beta = COLUMNS.indexOf("beta");
        int se = COLUMNS.indexOf("se");
        int eaf = COLUMNS.indexOf("eaf");
        int effectAllele = COLUMNS.indexOf("effect_allele");
        int otherAllele = COLUMNS.indexOf("other_allele");
        int pval = COLUMNS.indexOf("pval");
        require(output.get(snp), "SNP");
        output.set(effectAllele, allele(output.get(effectAllele)));
        output.set(otherAllele, allele(output.get(otherAllele)));
        if (output.get(effectAllele).equals(output.get(otherAllele)))
            throw new IllegalArgumentException("alleles are equal");
        double effect = finite(output.get(beta), "beta");
        if (resolved.effectScale() == EffectScale.ODDS_RATIO
                || resolved.effectScale() == EffectScale.HAZARD_RATIO) {
            if (!(effect > 0.0))
                throw new IllegalArgumentException("effect ratio is not positive");
            effect = Math.log(effect);
            output.set(beta, Double.toString(effect));
        }
        double standardError = finite(output.get(se), "se");
        if (!(standardError > 0.0))
            throw new IllegalArgumentException("se is not positive");
        double p = finite(output.get(pval), "pval");
        if (p < 0.0 || p > 1.0)
            throw new IllegalArgumentException("pval is outside [0, 1]");
        if (!output.get(eaf).isBlank()) {
            double frequency = finite(output.get(eaf), "eaf");
            if (frequency < 0.0 || frequency > 1.0)
                throw new IllegalArgumentException("eaf is outside [0, 1]");
        }
        return new Row(output, p);
    }

    private static String allele(String value) {
        String result = require(value, "allele").toUpperCase(Locale.ROOT);
        if (result.length() != 1 || "ACGT".indexOf(result.charAt(0)) < 0)
            throw new IllegalArgumentException("allele is not a biallelic SNP base");
        return result;
    }

    private static String require(String value, String name) {
        if (value == null || value.isBlank())
            throw new IllegalArgumentException(name + " is blank");
        return value;
    }

    private static double finite(String value, String name) {
        try {
            double result = Double.parseDouble(value);
            if (!Double.isFinite(result)) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(name + " is not finite", exception);
        }
    }

    private static BufferedReader reader(InputStream source) throws IOException {
        PushbackInputStream input = new PushbackInputStream(source, 2);
        int first = input.read();
        int second = input.read();
        if (second >= 0) input.unread(second);
        if (first >= 0) input.unread(first);
        InputStream decoded = first == 0x1f && second == 0x8b
            ? new GZIPInputStream(input) : input;
        return new BufferedReader(new InputStreamReader(
            decoded, StandardCharsets.UTF_8));
    }

    private static Map<String, Integer> columns(List<String> header) {
        Map<String, Integer> result = new LinkedHashMap<>();
        for (int index = 0; index < header.size(); index++) {
            String name = header.get(index).trim();
            if (name.isEmpty()) throw new IllegalArgumentException(
                "input column names must be nonblank");
            String normalized = normalize(name);
            if (result.put(normalized, index) != null)
                throw new IllegalArgumentException(
                    "input column names are ambiguous after normalization: " + name);
        }
        return result;
    }

    private static Integer find(Map<String, Integer> columns,
            List<String> aliases) {
        if (aliases == null) return null;
        for (String alias : aliases) {
            Integer index = columns.get(normalize(alias));
            if (index != null) return index;
        }
        return null;
    }

    private static Integer findExact(List<String> header, String name) {
        for (int index = 0; index < header.size(); index++)
            if (header.get(index).equals(name)) return index;
        for (int index = 0; index < header.size(); index++)
            if (header.get(index).equalsIgnoreCase(name)) return index;
        throw new IllegalArgumentException("mapped source column is absent: " + name);
    }

    private static String normalize(String value) {
        return value.toLowerCase(Locale.ROOT).replaceAll("[^a-z0-9]", "");
    }

    private static List<String> parse(String line, char delimiter,
            long lineNumber) throws IOException {
        List<String> result = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean quoted = false;
        for (int index = 0; index < line.length(); index++) {
            char value = line.charAt(index);
            if (value == '"') {
                if (quoted && index + 1 < line.length()
                        && line.charAt(index + 1) == '"') {
                    current.append('"');
                    index++;
                } else quoted = !quoted;
            } else if (value == delimiter && !quoted) {
                result.add(current.toString());
                current.setLength(0);
            } else current.append(value);
        }
        if (quoted) throw new IOException(
            "unterminated quote at line " + lineNumber);
        result.add(current.toString());
        return result;
    }

    private static String tsv(String value) {
        if (value.indexOf('\t') < 0 && value.indexOf('\n') < 0
                && value.indexOf('\r') < 0 && value.indexOf('"') < 0)
            return value;
        return '"' + value.replace("\"", "\"\"") + '"';
    }

    private static Map<String, List<String>> aliases() {
        Map<String, List<String>> result = new LinkedHashMap<>();
        result.put("Phenotype", List.of("Phenotype", "phenotype", "trait"));
        result.put("SNP", List.of("SNP", "rsid", "rs_id", "variant_id",
            "markername", "marker", "id"));
        result.put("beta", List.of("beta", "effect", "effect_size", "b"));
        result.put("se", List.of("se", "standard_error", "stderr", "sebeta"));
        result.put("eaf", List.of("eaf", "effect_allele_frequency",
            "effect_allele_freq", "af", "frequency"));
        result.put("effect_allele", List.of("effect_allele", "ea", "a1",
            "allele1", "tested_allele", "alt"));
        result.put("other_allele", List.of("other_allele", "non_effect_allele",
            "nea", "a2", "allele2", "reference_allele", "ref"));
        result.put("pval", List.of("pval", "p_value", "pvalue", "p"));
        result.put("units", List.of("units", "unit"));
        result.put("ncase", List.of("ncase", "n_case", "cases"));
        result.put("ncontrol", List.of("ncontrol", "n_control", "controls"));
        result.put("samplesize", List.of("samplesize", "sample_size", "n",
            "n_total", "total_n"));
        result.put("gene", List.of("gene", "gene_name", "symbol"));
        return Map.copyOf(result);
    }

    private record Resolved(Map<String, Integer> columns,
            EffectScale effectScale) { }
    private record Row(List<String> values, double pValue) { }
}
