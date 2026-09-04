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
import java.io.PrintStream;
import java.io.PushbackInputStream;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import org.jlinalg.mr.HarmonizationExclusion;
import org.jlinalg.mr.MrAnalysisResult;
import org.jlinalg.mr.MrEggerResult;
import org.jlinalg.mr.MrEstimate;
import org.jlinalg.mr.MrOptions;
import org.jlinalg.mr.SummaryAssociation;
import org.jlinalg.mr.XwasMrBatchResult;
import org.jlinalg.mr.XwasMrExposure;
import org.jlinalg.mr.XwasMrFailure;
import org.jlinalg.mr.XwasMrHit;
import org.jlinalg.mr.XwasMrOptions;
import org.jlinalg.mr.XwasMrOutcome;
import org.jlinalg.mr.XwasMrPipeline;
import org.jlinalg.mr.XwasMrScreeningMethod;
import org.jlinalg.mr.XwasMrScreeningResult;
import org.jlinalg.mr.XwasMrSignificanceFilter;

/** File-to-file parallel xWAS Mendelian-randomization command. */
final class MrXwasCli {
    private static final List<String> FDR_HEADER = List.of(
        "exposure_id", "exposure_label", "outcome_category", "outcome_id",
        "outcome_label", "screen_method", "nsnp", "beta", "se", "statistic",
        "p_value", "log10_p_value", "negative_log10_p_value", "ci_lower",
        "ci_upper", "cochran_q", "heterogeneity_df",
        "heterogeneity_p_value", "i_squared", "threshold_passed");

    private MrXwasCli() { }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput) {
        try {
            Options options = Options.parse(arguments,
                Path.of("").toAbsolutePath().normalize());
            if (options.help) {
                output.println(help());
                return 0;
            }
            execute(options, output);
            return 0;
        } catch (IllegalArgumentException exception) {
            errorOutput.println("jlinalg: " + exception.getMessage());
            errorOutput.println("Use mr-xwas --help for usage.");
            return 2;
        } catch (IOException exception) {
            errorOutput.println("jlinalg: " + exception.getMessage());
            return 1;
        }
    }

    private static void execute(Options options, PrintStream output)
            throws IOException {
        Table exposureTable = Table.read(options.exposure);
        Table outcomeTable = Table.read(options.outcome);
        List<XwasMrExposure> exposures = exposures(exposureTable, options);
        List<XwasMrOutcome> outcomes = outcomes(outcomeTable, options);
        XwasMrPipeline pipeline = XwasMrPipeline.prepare(exposures, outcomes);
        XwasMrOptions scanOptions = new XwasMrOptions(
            options.threads, options.pairBlockSize, options.screeningMethod,
            options.filter, new MrOptions(options.confidenceLevel,
                options.bootstrapReplicates, options.seed));
        Path failures = options.failures == null
            ? failurePath(options.output) : options.failures;
        validateOutputs(options, failures);

        XwasMrBatchResult result;
        long adjustedTests = 0L;
        if (options.fdrOutput == null) {
            result = pipeline.scan(scanOptions);
        } else {
            try (ExternalBh fdr = new ExternalBh(options.fdrOutput,
                    options.overwrite)) {
                fdr.writeHeader(FDR_HEADER);
                try {
                    result = pipeline.scan(scanOptions,
                        screened -> writeFdr(fdr, screened));
                } catch (UncheckedIOException exception) {
                    throw exception.getCause();
                }
                fdr.finish();
                adjustedTests = fdr.tests();
            }
        }

        writeAtomic(options.output, options.overwrite,
            writer -> writeHits(writer, delimiter(options.output), result));
        writeAtomic(failures, options.overwrite,
            writer -> writeFailures(writer, delimiter(failures), result));

        output.println("Scanned " + result.totalPairs() + " exposure-outcome "
            + "pairs (" + exposures.size() + " x " + outcomes.size() + ") with "
            + result.workersUsed() + " worker(s)");
        output.println("Retained " + result.hits().size() + "; below threshold "
            + result.belowThresholdPairs() + "; fewer than 3 harmonized "
            + "instruments " + result.insufficientInstrumentPairs()
            + "; diagnostic failures " + result.failures().size());
        output.printf(Locale.ROOT, "Analysis time %.6f seconds%n",
            result.elapsedNanoseconds() / 1e9);
        output.println("Wrote " + options.output);
        output.println("Wrote " + failures);
        if (options.fdrOutput != null) {
            output.println("BH-adjusted " + adjustedTests
                + " successfully screened pairs");
            output.println("Wrote " + options.fdrOutput);
        }
    }

    private static void validateOutputs(Options options, Path failures)
            throws IOException {
        if (options.output.equals(failures))
            throw new IllegalArgumentException(
                "--output and --failures must be different files");
        if (options.fdrOutput != null) {
            if (options.fdrOutput.equals(options.output)
                    || options.fdrOutput.equals(failures))
                throw new IllegalArgumentException(
                    "--fdr-output must differ from result and failure files");
            if (gzip(options.fdrOutput)
                    || delimiter(options.fdrOutput) != '\t')
                throw new IllegalArgumentException(
                    "--fdr-output must be an uncompressed TSV file");
        }
        if (!options.overwrite) {
            for (Path path : List.of(options.output, failures))
                if (Files.exists(path)) throw new IOException(
                    "output exists; use --overwrite: " + path);
        }
    }

    private static void writeFdr(ExternalBh output,
            XwasMrScreeningResult result) {
        MrEstimate estimate = result.estimate();
        try {
            output.write(List.of(result.exposureId(), result.exposureLabel(),
                result.outcomeCategory(), result.outcomeId(),
                result.outcomeLabel(), estimate.method().name(),
                integer(estimate.instrumentCount()), number(estimate.estimate()),
                number(estimate.standardError()), number(estimate.statistic()),
                number(estimate.pValue()), number(log10(estimate.pValue())),
                number(result.negativeLog10PValue()),
                number(estimate.confidenceLower()),
                number(estimate.confidenceUpper()), number(estimate.cochranQ()),
                integer(estimate.heterogeneityDegreesOfFreedom()),
                number(estimate.heterogeneityPValue()),
                number(estimate.iSquared()),
                Boolean.toString(result.thresholdPassed())),
                estimate.pValue());
        } catch (IOException exception) {
            throw new UncheckedIOException(exception);
        }
    }

    private static List<XwasMrExposure> exposures(Table table,
            Options options) {
        Columns columns = Columns.resolve(table, Role.EXPOSURE,
            options.exposureIdColumn, null);
        Map<String, Group> groups = groups(table, columns, Role.EXPOSURE);
        List<XwasMrExposure> result = new ArrayList<>(groups.size());
        for (Group group : groups.values()) result.add(new XwasMrExposure(
            group.id, group.label, group.associations));
        return List.copyOf(result);
    }

    private static List<XwasMrOutcome> outcomes(Table table,
            Options options) {
        Columns columns = Columns.resolve(table, Role.OUTCOME,
            options.outcomeIdColumn, options.categoryColumn);
        Map<String, Group> groups = groups(table, columns, Role.OUTCOME);
        List<XwasMrOutcome> result = new ArrayList<>(groups.size());
        for (Group group : groups.values()) result.add(new XwasMrOutcome(
            group.id, group.label, group.category, group.associations));
        return List.copyOf(result);
    }

    private static Map<String, Group> groups(Table table, Columns columns,
            Role role) {
        Map<String, Group> result = new LinkedHashMap<>();
        for (int rowIndex = 0; rowIndex < table.rows.size(); rowIndex++) {
            List<String> row = table.rows.get(rowIndex);
            int line = rowIndex + 2;
            String id = required(row.get(columns.id), role.text + " id", line);
            String label = columns.label < 0 ? id
                : required(row.get(columns.label), role.text + " label", line);
            String category = columns.category < 0 ? ""
                : row.get(columns.category).trim();
            Group group = result.computeIfAbsent(id,
                ignored -> new Group(id, label, category));
            if (!group.label.equals(label) || !group.category.equals(category))
                throw new IllegalArgumentException("inconsistent " + role.text
                    + " label/category for " + id + " at input row " + line);
            String eaf = columns.eaf < 0 ? "" : row.get(columns.eaf).trim();
            double frequency = eaf.isEmpty() ? Double.NaN
                : number(eaf, "effect-allele frequency", line);
            group.associations.add(new SummaryAssociation(
                required(row.get(columns.snp), "SNP", line),
                required(row.get(columns.effectAllele), "effect allele", line),
                required(row.get(columns.otherAllele), "other allele", line),
                number(row.get(columns.beta), "beta", line),
                positive(row.get(columns.se), "standard error", line),
                frequency));
        }
        return result;
    }

    private static String required(String value, String name, int line) {
        String result = value.trim();
        if (result.isEmpty()) throw new IllegalArgumentException(
            "blank " + name + " at input row " + line);
        return result;
    }

    private static double number(String value, String name, int line) {
        try {
            double result = Double.parseDouble(value.trim());
            if (!Double.isFinite(result)) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                "invalid " + name + " at input row " + line + ": " + value,
                exception);
        }
    }

    private static double positive(String value, String name, int line) {
        double result = number(value, name, line);
        if (!(result > 0.0)) throw new IllegalArgumentException(
            name + " must be positive at input row " + line);
        return result;
    }

    private static void writeHits(BufferedWriter output, char delimiter,
            XwasMrBatchResult batch) throws IOException {
        row(output, delimiter, List.of("exposure_id", "exposure_label",
            "outcome_category", "outcome_id", "outcome_label",
            "screen_method", "nsnp", "beta", "se", "statistic", "p_value",
            "log10_p_value", "negative_log10_p_value", "ci_lower", "ci_upper",
            "cochran_q", "heterogeneity_df", "heterogeneity_p_value",
            "i_squared", "mean_f", "egger_beta", "egger_se", "egger_p_value",
            "egger_intercept", "egger_intercept_se", "egger_intercept_p_value",
            "i_squared_gx", "weighted_median_beta", "weighted_median_se",
            "weighted_median_p_value", "harmonization_exclusions", "warnings"));
        for (XwasMrHit hit : batch.hits()) {
            MrEstimate screen = hit.screeningEstimate();
            MrAnalysisResult analysis = hit.analysis();
            MrEggerResult egger = analysis.egger();
            MrEstimate median = analysis.weightedMedian();
            row(output, delimiter, List.of(hit.exposureId(), hit.exposureLabel(),
                hit.outcomeCategory(), hit.outcomeId(), hit.outcomeLabel(),
                screen.method().name(), integer(screen.instrumentCount()),
                number(screen.estimate()), number(screen.standardError()),
                number(screen.statistic()), number(screen.pValue()),
                number(log10(screen.pValue())),
                number(hit.negativeLog10PValue()),
                number(screen.confidenceLower()), number(screen.confidenceUpper()),
                number(screen.cochranQ()),
                integer(screen.heterogeneityDegreesOfFreedom()),
                number(screen.heterogeneityPValue()), number(screen.iSquared()),
                number(analysis.meanFStatistic()),
                number(egger.slope().estimate()),
                number(egger.slope().standardError()),
                number(egger.slope().pValue()), number(egger.intercept()),
                number(egger.interceptStandardError()),
                number(egger.interceptPValue()), number(egger.iSquaredGx()),
                number(median.estimate()), number(median.standardError()),
                number(median.pValue()), exclusions(hit.harmonizationExclusions()),
                String.join(" | ", analysis.warnings())));
        }
    }

    private static void writeFailures(BufferedWriter output, char delimiter,
            XwasMrBatchResult batch) throws IOException {
        row(output, delimiter, List.of("exposure_index", "outcome_index",
            "exposure_id", "outcome_id", "exception_type", "message"));
        for (XwasMrFailure failure : batch.failures()) row(output, delimiter,
            List.of(integer(failure.exposureIndex()),
                integer(failure.outcomeIndex()), failure.exposureId(),
                failure.outcomeId(), failure.exceptionType(), failure.message()));
    }

    private static String exclusions(List<HarmonizationExclusion> values) {
        Map<String, Integer> counts = new LinkedHashMap<>();
        for (HarmonizationExclusion value : values)
            counts.merge(value.reason().name(), 1, Integer::sum);
        List<String> result = new ArrayList<>();
        for (Map.Entry<String, Integer> entry : counts.entrySet())
            result.add(entry.getKey() + "=" + entry.getValue());
        return String.join(" | ", result);
    }

    private static double log10(double pValue) {
        return pValue == 0.0 ? Double.NEGATIVE_INFINITY : Math.log10(pValue);
    }

    private static String number(double value) {
        return Double.toString(value);
    }

    private static String integer(int value) {
        return Integer.toString(value);
    }

    private static void row(BufferedWriter output, char delimiter,
            List<String> values) throws IOException {
        for (int index = 0; index < values.size(); index++) {
            if (index > 0) output.write(delimiter);
            output.write(quote(values.get(index), delimiter));
        }
        output.newLine();
    }

    private static String quote(String value, char delimiter) {
        if (value.indexOf(delimiter) < 0 && value.indexOf('"') < 0
                && value.indexOf('\r') < 0 && value.indexOf('\n') < 0)
            return value;
        return '"' + value.replace("\"", "\"\"") + '"';
    }

    private static void writeAtomic(Path destination, boolean overwrite,
            WriterOperation operation) throws IOException {
        if (Files.exists(destination) && !overwrite) throw new IOException(
            "output exists; use --overwrite: " + destination);
        Path parent = destination.getParent();
        if (parent != null) Files.createDirectories(parent);
        Path part = destination.resolveSibling("." + destination.getFileName()
            + ".part");
        Files.deleteIfExists(part);
        boolean complete = false;
        try {
            try (OutputStream file = Files.newOutputStream(part,
                        StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
                    OutputStream encoded = gzip(destination)
                        ? new GZIPOutputStream(file) : file;
                    BufferedWriter writer = new BufferedWriter(
                        new OutputStreamWriter(encoded, StandardCharsets.UTF_8))) {
                operation.write(writer);
            }
            move(part, destination, overwrite);
            complete = true;
        } finally {
            if (!complete) Files.deleteIfExists(part);
        }
    }

    private static void move(Path source, Path destination, boolean overwrite)
            throws IOException {
        try {
            if (overwrite) Files.move(source, destination,
                StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING);
            else Files.move(source, destination, StandardCopyOption.ATOMIC_MOVE);
        } catch (AtomicMoveNotSupportedException exception) {
            if (overwrite) Files.move(source, destination,
                StandardCopyOption.REPLACE_EXISTING);
            else Files.move(source, destination);
        }
    }

    private static Path failurePath(Path output) {
        String name = output.getFileName().toString();
        boolean compressed = name.toLowerCase(Locale.ROOT).endsWith(".gz");
        if (compressed) name = name.substring(0, name.length() - 3);
        int dot = name.lastIndexOf('.');
        String result = dot < 0 ? name + ".failures.tsv"
            : name.substring(0, dot) + ".failures" + name.substring(dot);
        if (compressed) result += ".gz";
        return output.resolveSibling(result);
    }

    private static char delimiter(Path path) {
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        if (name.endsWith(".gz")) name = name.substring(0, name.length() - 3);
        return name.endsWith(".csv") ? ',' : '\t';
    }

    private static boolean gzip(Path path) {
        return path.getFileName().toString().toLowerCase(Locale.ROOT)
            .endsWith(".gz");
    }

    private static String help() {
        return """
            Usage:
              java -jar jlinalg.jar mr-xwas --exposure CLUMPED_FILE
                --outcome MULTI_PHENOTYPE_FILE --output RESULTS_FILE
                (--p-threshold X | --log10-p-threshold X |
                 --negative-log10-p-threshold X)

            Required long-format columns:
              exposure: SNP, beta, se, effect_allele, other_allele, and
                        id.exposure, gene, or Phenotype
              outcome:  SNP, beta, se, effect_allele, other_allele, and
                        id.outcome or Phenotype
              eaf is optional. CSV, TSV, and gzip are supported.

            Execution:
              --threads N                 Default: available processors
              --pair-block-size N         Default: max(32, 8 * threads)
              --screen-method fixed|random  Default: random
              --bootstrap-replicates N    Default: 1000; hits only
              --confidence-level X        Default: 0.95
              --seed N                    Default: 20260831

            Column overrides:
              --exposure-id-column NAME
              --outcome-id-column NAME
              --category-column NAME
              --failures FILE             Default: RESULTS.failures.ext
              --fdr-output FILE            All screened pairs plus BH q-values;
                                           uncompressed TSV
              --overwrite

            Threshold directions are inclusive: p <= X, log10(p) <= X,
            and -log10(p) >= X. Full diagnostics run only for retained pairs.
            BH uses every successfully screened pair, without prefiltering.
            """;
    }

    private enum Role {
        EXPOSURE("exposure"), OUTCOME("outcome");
        final String text;
        Role(String text) { this.text = text; }
    }

    private static final class Group {
        final String id;
        final String label;
        final String category;
        final List<SummaryAssociation> associations = new ArrayList<>();
        Group(String id, String label, String category) {
            this.id = id;
            this.label = label;
            this.category = category;
        }
    }

    private record Columns(int id, int label, int category, int snp,
            int beta, int se, int effectAllele, int otherAllele, int eaf) {
        static Columns resolve(Table table, Role role, String explicitId,
                String explicitCategory) {
            int id = explicitId == null ? table.required(role == Role.EXPOSURE
                ? List.of("id.exposure", "Phenotype", "gene")
                : List.of("id.outcome", "Phenotype"), role.text + " id")
                : table.exact(explicitId);
            int label = table.optional(role == Role.EXPOSURE
                ? List.of("exposure", "Phenotype")
                : List.of("outcome", "Phenotype"));
            int category = explicitCategory == null
                ? table.optional(List.of("category", "phenotype_category"))
                : table.exact(explicitCategory);
            return new Columns(id, label, category,
                table.required(List.of("SNP", "rsid", "rs_id"), "SNP"),
                table.required(role == Role.EXPOSURE
                    ? List.of("beta.exposure", "beta", "estimate", "effect")
                    : List.of("beta.outcome", "beta", "estimate", "effect"),
                    "beta"),
                table.required(role == Role.EXPOSURE
                    ? List.of("se.exposure", "se", "stderr", "standard_error")
                    : List.of("se.outcome", "se", "stderr", "standard_error"),
                    "standard error"),
                table.required(role == Role.EXPOSURE
                    ? List.of("effect_allele.exposure", "effect_allele", "ea")
                    : List.of("effect_allele.outcome", "effect_allele", "ea"),
                    "effect allele"),
                table.required(role == Role.EXPOSURE
                    ? List.of("other_allele.exposure", "other_allele", "oa")
                    : List.of("other_allele.outcome", "other_allele", "oa"),
                    "other allele"),
                table.optional(role == Role.EXPOSURE
                    ? List.of("eaf.exposure", "eaf", "frequency")
                    : List.of("eaf.outcome", "eaf", "frequency")));
        }
    }

    private static final class Table {
        final List<String> header;
        final List<List<String>> rows;

        Table(List<String> header, List<List<String>> rows) {
            this.header = List.copyOf(header);
            this.rows = List.copyOf(rows);
        }

        static Table read(Path path) throws IOException {
            char delimiter = delimiter(path);
            try (InputStream file = Files.newInputStream(path);
                    BufferedReader input = reader(file)) {
                String first = input.readLine();
                if (first == null) throw new IOException("table is empty: " + path);
                List<String> header = parse(first, delimiter, 1, path);
                Set<String> names = new HashSet<>();
                for (String value : header)
                    if (value.isBlank() || !names.add(normalize(value)))
                        throw new IOException(
                            "column names must be unique and nonblank: " + path);
                List<List<String>> rows = new ArrayList<>();
                int line = 1;
                for (String text; (text = input.readLine()) != null;) {
                    line++;
                    if (text.isBlank()) continue;
                    List<String> row = parse(text, delimiter, line, path);
                    if (row.size() != header.size()) throw new IOException(
                        "expected " + header.size() + " fields but found "
                            + row.size() + " at line " + line + " in " + path);
                    rows.add(row);
                }
                if (rows.isEmpty()) throw new IOException(
                    "table contains no rows: " + path);
                return new Table(header, rows);
            }
        }

        int exact(String name) {
            int result = find(name);
            if (result < 0) throw new IllegalArgumentException(
                "column is absent: " + name);
            return result;
        }

        int required(List<String> aliases, String description) {
            int result = optional(aliases);
            if (result < 0) throw new IllegalArgumentException(
                "required " + description + " column is absent");
            return result;
        }

        int optional(List<String> aliases) {
            for (String alias : aliases) {
                int result = find(alias);
                if (result >= 0) return result;
            }
            return -1;
        }

        private int find(String name) {
            String wanted = normalize(name);
            for (int index = 0; index < header.size(); index++)
                if (normalize(header.get(index)).equals(wanted)) return index;
            return -1;
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

        private static List<String> parse(String line, char delimiter,
                int lineNumber, Path path) throws IOException {
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
            if (quoted) throw new IOException("unterminated quote at line "
                + lineNumber + " in " + path);
            result.add(current.toString());
            return result;
        }

        private static String normalize(String value) {
            return value.trim().toLowerCase(Locale.ROOT)
                .replaceAll("[^a-z0-9]", "");
        }
    }

    @FunctionalInterface
    private interface WriterOperation {
        void write(BufferedWriter output) throws IOException;
    }

    private static final class Options {
        Path exposure;
        Path outcome;
        Path output;
        Path failures;
        Path fdrOutput;
        String exposureIdColumn;
        String outcomeIdColumn;
        String categoryColumn;
        int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
        int pairBlockSize;
        int bootstrapReplicates = 1_000;
        long seed = 20260831L;
        double confidenceLevel = 0.95;
        XwasMrScreeningMethod screeningMethod =
            XwasMrScreeningMethod.IVW_MULTIPLICATIVE_RANDOM;
        XwasMrSignificanceFilter filter;
        boolean overwrite;
        boolean help;

        static Options parse(String[] arguments, Path currentDirectory) {
            Options result = new Options();
            for (int index = 0; index < arguments.length; index++) {
                String option = arguments[index];
                switch (option) {
                    case "--exposure", "--instrument" -> result.exposure = path(
                        currentDirectory, value(arguments, ++index, option));
                    case "--outcome" -> result.outcome = path(currentDirectory,
                        value(arguments, ++index, option));
                    case "--output", "--out" -> result.output = path(
                        currentDirectory, value(arguments, ++index, option));
                    case "--failures" -> result.failures = path(currentDirectory,
                        value(arguments, ++index, option));
                    case "--fdr-output" -> result.fdrOutput = path(
                        currentDirectory, value(arguments, ++index, option));
                    case "--threads" -> result.threads = positiveInteger(
                        value(arguments, ++index, option), option);
                    case "--pair-block-size" -> result.pairBlockSize =
                        positiveInteger(value(arguments, ++index, option), option);
                    case "--bootstrap-replicates" -> result.bootstrapReplicates =
                        positiveInteger(value(arguments, ++index, option), option);
                    case "--seed" -> result.seed = longValue(
                        value(arguments, ++index, option), option);
                    case "--confidence-level" -> result.confidenceLevel =
                        probability(value(arguments, ++index, option), option, true);
                    case "--screen-method" -> result.screeningMethod = method(
                        value(arguments, ++index, option));
                    case "--p-threshold" -> result.setFilter(
                        XwasMrSignificanceFilter.pValueAtMost(probability(
                            value(arguments, ++index, option), option, false)));
                    case "--log10-p-threshold" -> result.setFilter(
                        XwasMrSignificanceFilter.log10PAtMost(number(
                            value(arguments, ++index, option), option)));
                    case "--negative-log10-p-threshold",
                            "--minus-log10-p-threshold" -> result.setFilter(
                        XwasMrSignificanceFilter.negativeLog10PAtLeast(number(
                            value(arguments, ++index, option), option)));
                    case "--exposure-id-column" -> result.exposureIdColumn =
                        value(arguments, ++index, option);
                    case "--outcome-id-column" -> result.outcomeIdColumn =
                        value(arguments, ++index, option);
                    case "--category-column" -> result.categoryColumn =
                        value(arguments, ++index, option);
                    case "--overwrite" -> result.overwrite = true;
                    case "--help", "-h" -> result.help = true;
                    default -> throw new IllegalArgumentException(
                        "unknown mr-xwas option: " + option);
                }
            }
            if (result.pairBlockSize == 0)
                result.pairBlockSize = Math.max(32, result.threads * 8);
            if (!result.help && (result.exposure == null
                    || result.outcome == null || result.output == null))
                throw new IllegalArgumentException(
                    "--exposure, --outcome, and --output are required");
            if (!result.help && result.filter == null)
                throw new IllegalArgumentException(
                    "one p-value threshold option is required");
            if (result.bootstrapReplicates < 2)
                throw new IllegalArgumentException(
                    "--bootstrap-replicates must be at least 2");
            return result;
        }

        void setFilter(XwasMrSignificanceFilter value) {
            if (filter != null) throw new IllegalArgumentException(
                "specify exactly one p-value threshold scale");
            filter = value;
        }

        private static String value(String[] arguments, int index,
                String option) {
            if (index >= arguments.length || arguments[index].isBlank())
                throw new IllegalArgumentException(option + " requires a value");
            return arguments[index];
        }

        private static Path path(Path directory, String value) {
            Path path = Path.of(value);
            return (path.isAbsolute() ? path : directory.resolve(path))
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

        private static long longValue(String value, String option) {
            try {
                return Long.parseLong(value);
            } catch (NumberFormatException exception) {
                throw new IllegalArgumentException(
                    option + " requires an integer: " + value, exception);
            }
        }

        private static double probability(String value, String option,
                boolean open) {
            double result = number(value, option);
            boolean valid = open ? result > 0.0 && result < 1.0
                : result >= 0.0 && result <= 1.0;
            if (!valid) throw new IllegalArgumentException(option
                + (open ? " must lie in (0, 1)" : " must lie in [0, 1]"));
            return result;
        }

        private static double number(String value, String option) {
            try {
                double result = Double.parseDouble(value);
                if (!Double.isFinite(result)) throw new NumberFormatException();
                return result;
            } catch (NumberFormatException exception) {
                throw new IllegalArgumentException(
                    option + " requires a finite number: " + value, exception);
            }
        }

        private static XwasMrScreeningMethod method(String value) {
            return switch (value.toLowerCase(Locale.ROOT)) {
                case "fixed", "ivw-fixed" -> XwasMrScreeningMethod.IVW_FIXED;
                case "random", "ivw-random", "mre" ->
                    XwasMrScreeningMethod.IVW_MULTIPLICATIVE_RANDOM;
                default -> throw new IllegalArgumentException(
                    "--screen-method must be fixed or random: " + value);
            };
        }
    }
}
