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
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PushbackInputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import org.jlinalg.genetics.LdClumpCandidate;
import org.jlinalg.genetics.LdClumpExclusion;
import org.jlinalg.genetics.LdClumpExclusionReason;
import org.jlinalg.genetics.LdClumpOptions;
import org.jlinalg.genetics.LdClumpResult;
import org.jlinalg.genetics.LdReferenceLayout;
import org.jlinalg.genetics.PlinkBedLdClumper;

/** TwoSampleMR-compatible local clumping command over installed LD panels. */
final class LdClumpCli {
    private static final String DEFAULT_GROUP = "jlinalg-clump";

    private LdClumpCli() { }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput) {
        return run(arguments, output, errorOutput,
            Path.of("").toAbsolutePath().normalize());
    }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput, Path currentDirectory) {
        try {
            Options options = Options.parse(arguments, currentDirectory);
            if (options.help) {
                output.println(help());
                return 0;
            }
            return execute(options, output);
        } catch (IllegalArgumentException exception) {
            errorOutput.println("jlinalg: " + exception.getMessage());
            errorOutput.println("Use clump --help for usage.");
            return 2;
        } catch (IOException exception) {
            errorOutput.println("jlinalg: " + exception.getMessage());
            return 1;
        }
    }

    private static int execute(Options options, PrintStream output)
            throws IOException {
        validateDatabase(options.database, options.population);
        Table table = Table.read(options.instrument);
        int snpColumn = table.column(options.snpColumn,
            List.of("SNP", "rsid", "rs_id"), true);
        int pValueColumn = table.column(options.pValueColumn,
            List.of("pval.exposure", "pval.outcome", "pval", "p"), false);
        int groupColumn = table.column(options.groupColumn,
            List.of("id.exposure", "id"), false);
        if (options.pValueColumn != null && pValueColumn < 0)
            throw new IllegalArgumentException("p-value column is absent: "
                + options.pValueColumn);
        if (options.groupColumn != null && groupColumn < 0)
            throw new IllegalArgumentException("group column is absent: "
                + options.groupColumn);

        Map<String, LdClumpCandidate> unique = new LinkedHashMap<>();
        List<String> rowKeys = new ArrayList<>(table.rows.size());
        for (int rowIndex = 0; rowIndex < table.rows.size(); rowIndex++) {
            List<String> row = table.rows.get(rowIndex);
            String snp = row.get(snpColumn).trim();
            if (snp.isEmpty()) throw new IllegalArgumentException(
                "blank SNP at input row " + (rowIndex + 2));
            String group = groupColumn < 0 ? DEFAULT_GROUP
                : row.get(groupColumn).trim();
            if (group.isEmpty()) throw new IllegalArgumentException(
                "blank clumping group at input row " + (rowIndex + 2));
            double pValue = pValueColumn < 0 ? 0.99
                : pValue(row.get(pValueColumn), rowIndex + 2);
            String key = key(group, snp);
            rowKeys.add(key);
            LdClumpCandidate candidate = new LdClumpCandidate(snp, pValue, group);
            LdClumpCandidate previous = unique.get(key);
            if (previous == null || pValue < previous.pValue())
                unique.put(key, candidate);
        }

        Path prefix = LdReferenceLayout.panelPrefix(
            options.database, options.population);
        LdClumpResult result = PlinkBedLdClumper.clump(prefix,
            List.copyOf(unique.values()), new LdClumpOptions(
                options.windowKilobases, options.rSquaredThreshold,
                options.indexPValueThreshold));
        Set<String> retainedKeys = new HashSet<>();
        for (LdClumpCandidate candidate : result.retained())
            retainedKeys.add(key(candidate.group(), candidate.variantId()));
        List<Integer> retainedRows = new ArrayList<>();
        for (int index = 0; index < rowKeys.size(); index++)
            if (retainedKeys.contains(rowKeys.get(index))) retainedRows.add(index);
        table.write(options.output, retainedRows, options.overwrite);

        long absent = count(result.exclusions(),
            LdClumpExclusionReason.ABSENT_FROM_REFERENCE);
        long linked = count(result.exclusions(),
            LdClumpExclusionReason.IN_LINKAGE_DISEQUILIBRIUM);
        long above = count(result.exclusions(),
            LdClumpExclusionReason.ABOVE_INDEX_P_VALUE_THRESHOLD);
        output.println("Clumped " + unique.size() + " unique variants in "
            + groups(unique.values()) + " group(s) using " + options.population
            + " (" + options.windowKilobases + " kb, r2 > "
            + options.rSquaredThreshold + ")");
        output.println("Retained " + result.retained().size()
            + " unique variants (" + retainedRows.size() + " rows); removed "
            + linked + " for LD, " + absent + " absent from the reference, and "
            + above + " above the index p-value threshold");
        output.println("Wrote " + options.output);
        return 0;
    }

    private static void validateDatabase(Path database, String population)
            throws IOException {
        Path manifest = LdReferenceLayout.manifest(database);
        if (!Files.isRegularFile(manifest)) throw new IOException(
            "not a complete JLinAlg LD database (manifest is absent): "
                + database);
        Object parsed = SimpleJson.parse(Files.readString(manifest,
            StandardCharsets.UTF_8));
        if (!(parsed instanceof Map<?, ?> root)
                || !LdReferenceLayout.FORMAT_NAME.equals(root.get("format"))
                || !(root.get("format_version") instanceof Number version)
                || version.intValue() != LdReferenceLayout.FORMAT_VERSION)
            throw new IOException("unsupported JLinAlg LD database manifest: "
                + manifest);
        boolean found = false;
        Object panelValue = root.get("panels");
        if (panelValue instanceof List<?> panels) for (Object value : panels)
            if (value instanceof Map<?, ?> panel
                    && population.equalsIgnoreCase(text(panel.get("id")))) {
                found = true;
                break;
            }
        if (!found) throw new IllegalArgumentException("population " + population
            + " is not installed; available panels: " + panelIds(panelValue));
        Path prefix = LdReferenceLayout.panelPrefix(database, population);
        for (String suffix : List.of(".bed", ".bim", ".fam"))
            if (!Files.isRegularFile(Path.of(prefix + suffix)))
                throw new IOException("LD panel file is absent: " + prefix + suffix);
    }

    private static String panelIds(Object value) {
        List<String> result = new ArrayList<>();
        if (value instanceof List<?> panels) for (Object item : panels)
            if (item instanceof Map<?, ?> panel) {
                String id = text(panel.get("id"));
                if (!id.isBlank()) result.add(id);
            }
        return result.isEmpty() ? "none" : String.join(", ", result);
    }

    private static String text(Object value) {
        return value instanceof String string ? string : "";
    }

    private static long count(List<LdClumpExclusion> exclusions,
            LdClumpExclusionReason reason) {
        return exclusions.stream().filter(value -> value.reason() == reason).count();
    }

    private static long groups(java.util.Collection<LdClumpCandidate> values) {
        return values.stream().map(LdClumpCandidate::group).distinct().count();
    }

    private static String key(String group, String snp) {
        return group + '\u0000' + snp;
    }

    private static double pValue(String value, int row) {
        try {
            double result = Double.parseDouble(value.trim());
            if (!Double.isFinite(result) || result < 0.0 || result > 1.0)
                throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                "invalid p-value at input row " + row + ": " + value, exception);
        }
    }

    private static String help() {
        return """
            Usage:
              java -jar jlinalg.jar clump --database DIRECTORY
                --instrument FILE --ld-threshold 0.001 --output FILE

            TwoSampleMR-compatible defaults:
              --population EUR            AFR, AMR, EAS, EUR, or SAS
              --clump-kb 10000            Radius around each index SNP
              --ld-threshold 0.001        Remove pairs with r2 above this value
              --p-threshold 1             Maximum index-SNP p-value

            Input columns are detected in TwoSampleMR order: SNP/rsid, then
            pval.exposure, pval.outcome, pval, or p. If no p-value column exists,
            every value is set to 0.99. id.exposure or id defines independent
            clumping groups. Override with --snp-column, --pval-column, or
            --id-column. CSV, TSV, and gzip-compressed input are supported.
            """;
    }

    private static final class Options {
        Path database;
        Path instrument;
        Path output;
        String population = "EUR";
        String snpColumn;
        String pValueColumn;
        String groupColumn;
        int windowKilobases = 10_000;
        double rSquaredThreshold = 0.001;
        double indexPValueThreshold = 1.0;
        boolean overwrite;
        boolean help;

        static Options parse(String[] arguments, Path currentDirectory) {
            Options result = new Options();
            for (int index = 0; index < arguments.length; index++) {
                String option = arguments[index];
                switch (option) {
                    case "--database" -> result.database = path(currentDirectory,
                        value(arguments, ++index, option));
                    case "--instrument" -> result.instrument = path(currentDirectory,
                        value(arguments, ++index, option));
                    case "--output", "--out" -> result.output = path(currentDirectory,
                        value(arguments, ++index, option));
                    case "--population", "--pop" -> result.population =
                        value(arguments, ++index, option).toUpperCase(Locale.ROOT);
                    case "--snp-column" -> result.snpColumn =
                        value(arguments, ++index, option);
                    case "--pval-column" -> result.pValueColumn =
                        value(arguments, ++index, option);
                    case "--id-column" -> result.groupColumn =
                        value(arguments, ++index, option);
                    case "--clump-kb" -> result.windowKilobases = positiveInteger(
                        value(arguments, ++index, option), option);
                    case "--ld-threshold", "--clump-r2" ->
                        result.rSquaredThreshold = probability(
                            value(arguments, ++index, option), option);
                    case "--p-threshold", "--clump-p1" ->
                        result.indexPValueThreshold = probability(
                            value(arguments, ++index, option), option);
                    case "--overwrite" -> result.overwrite = true;
                    case "--help", "-h" -> result.help = true;
                    default -> throw new IllegalArgumentException(
                        "unknown clump option: " + option);
                }
            }
            if (!result.help && (result.database == null
                    || result.instrument == null || result.output == null))
                throw new IllegalArgumentException(
                    "--database, --instrument, and --output are required");
            return result;
        }

        private static String value(String[] arguments, int index,
                String option) {
            if (index >= arguments.length)
                throw new IllegalArgumentException(option + " requires a value");
            if (arguments[index].isBlank()) throw new IllegalArgumentException(
                option + " requires a nonblank value");
            return arguments[index];
        }

        private static Path path(Path currentDirectory, String value) {
            Path path = Path.of(value);
            return (path.isAbsolute() ? path : currentDirectory.resolve(path))
                .toAbsolutePath().normalize();
        }

        private static int positiveInteger(String value, String option) {
            try {
                int result = Integer.parseInt(value);
                if (result < 1) throw new NumberFormatException();
                return result;
            } catch (NumberFormatException exception) {
                throw new IllegalArgumentException(
                    option + " requires a positive integer: " + value, exception);
            }
        }

        private static double probability(String value, String option) {
            try {
                double result = Double.parseDouble(value);
                if (!Double.isFinite(result) || result < 0.0 || result > 1.0)
                    throw new NumberFormatException();
                return result;
            } catch (NumberFormatException exception) {
                throw new IllegalArgumentException(
                    option + " requires a value in [0, 1]: " + value, exception);
            }
        }
    }

    private static final class Table {
        final List<String> header;
        final List<List<String>> rows;
        final char inputDelimiter;

        Table(List<String> header, List<List<String>> rows, char delimiter) {
            this.header = List.copyOf(header);
            this.rows = List.copyOf(rows);
            inputDelimiter = delimiter;
        }

        static Table read(Path path) throws IOException {
            char delimiter = delimiter(path);
            try (InputStream file = Files.newInputStream(path);
                    BufferedReader input = reader(file)) {
                String first = input.readLine();
                if (first == null) throw new IOException(
                    "instrument table is empty: " + path);
                List<String> header = parse(first, delimiter, 1);
                Set<String> names = new HashSet<>();
                for (String value : header)
                    if (value.isBlank() || !names.add(normalize(value)))
                        throw new IOException(
                            "instrument column names must be unique and nonblank");
                List<List<String>> rows = new ArrayList<>();
                int lineNumber = 1;
                for (String line; (line = input.readLine()) != null;) {
                    lineNumber++;
                    if (line.isBlank()) continue;
                    List<String> row = parse(line, delimiter, lineNumber);
                    if (row.size() != header.size()) throw new IOException(
                        "expected " + header.size() + " fields but found "
                            + row.size() + " at line " + lineNumber);
                    rows.add(row);
                }
                if (rows.isEmpty()) throw new IOException(
                    "instrument table contains no rows: " + path);
                return new Table(header, rows, delimiter);
            }
        }

        int column(String explicit, List<String> aliases, boolean required) {
            if (explicit != null) {
                int index = indexOf(explicit);
                if (index < 0 && required) throw new IllegalArgumentException(
                    "column is absent: " + explicit);
                return index;
            }
            for (String alias : aliases) {
                int index = indexOf(alias);
                if (index >= 0) return index;
            }
            if (required) throw new IllegalArgumentException(
                "required SNP column is absent; use --snp-column");
            return -1;
        }

        private int indexOf(String name) {
            String target = normalize(name);
            for (int index = 0; index < header.size(); index++)
                if (normalize(header.get(index)).equals(target)) return index;
            return -1;
        }

        void write(Path destination, List<Integer> retained,
                boolean overwrite) throws IOException {
            if (Files.exists(destination) && !overwrite) throw new IOException(
                "output exists; use --overwrite: " + destination);
            Path parent = destination.getParent();
            if (parent != null) Files.createDirectories(parent);
            Path part = destination.resolveSibling("." + destination.getFileName()
                + ".part");
            Files.deleteIfExists(part);
            char delimiter = delimiter(destination);
            boolean complete = false;
            try {
                try (BufferedWriter output = new BufferedWriter(
                        new OutputStreamWriter(Files.newOutputStream(part,
                            StandardOpenOption.CREATE_NEW,
                            StandardOpenOption.WRITE), StandardCharsets.UTF_8))) {
                    writeRow(output, header, delimiter);
                    for (int index : retained)
                        writeRow(output, rows.get(index), delimiter);
                }
                move(part, destination, overwrite);
                complete = true;
            } finally {
                if (!complete) Files.deleteIfExists(part);
            }
        }

        private static void move(Path source, Path destination,
                boolean overwrite) throws IOException {
            try {
                if (overwrite) Files.move(source, destination,
                    StandardCopyOption.ATOMIC_MOVE,
                    StandardCopyOption.REPLACE_EXISTING);
                else Files.move(source, destination,
                    StandardCopyOption.ATOMIC_MOVE);
            } catch (AtomicMoveNotSupportedException exception) {
                if (overwrite) Files.move(source, destination,
                    StandardCopyOption.REPLACE_EXISTING);
                else Files.move(source, destination);
            }
        }

        private static void writeRow(BufferedWriter output,
                List<String> values, char delimiter) throws IOException {
            for (int index = 0; index < values.size(); index++) {
                if (index > 0) output.write(delimiter);
                output.write(quoted(values.get(index), delimiter));
            }
            output.newLine();
        }

        private static String quoted(String value, char delimiter) {
            if (value.indexOf(delimiter) < 0 && value.indexOf('"') < 0
                    && value.indexOf('\r') < 0 && value.indexOf('\n') < 0)
                return value;
            return '"' + value.replace("\"", "\"\"") + '"';
        }

        private static BufferedReader reader(InputStream source)
                throws IOException {
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

        private static char delimiter(Path path) {
            String name = path.getFileName().toString()
                .toLowerCase(Locale.ROOT);
            if (name.endsWith(".gz"))
                name = name.substring(0, name.length() - 3);
            return name.endsWith(".csv") ? ',' : '\t';
        }

        private static List<String> parse(String line, char delimiter,
                int lineNumber) throws IOException {
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

        private static String normalize(String value) {
            return value.trim().toLowerCase(Locale.ROOT)
                .replaceAll("[^a-z0-9]", "");
        }
    }
}
