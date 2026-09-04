/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Search, download, and normalize candidate instruments for MR. */
final class MrInstrumentCli {
    private static final Map<String, String> TRAIT_EXPANSIONS = Map.ofEntries(
        Map.entry("bmi", "body mass index"),
        Map.entry("whr", "waist-hip ratio"),
        Map.entry("ldl", "low density lipoprotein cholesterol"),
        Map.entry("hdl", "high density lipoprotein cholesterol"),
        Map.entry("t2d", "type 2 diabetes mellitus"),
        Map.entry("cad", "coronary artery disease"));
    private static final double DEFAULT_INSTRUMENT_P = 5e-8;

    private MrInstrumentCli() { }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput) {
        return run(arguments, output, errorOutput,
            Path.of("").toAbsolutePath().normalize(),
            new GwasCatalogClient());
    }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput, Path currentDirectory,
            InstrumentCatalog catalog) {
        if (arguments.length == 0) {
            error(errorOutput, "an mr-instruments command is required");
            errorOutput.println(help());
            return 2;
        }
        try {
            return switch (arguments[0].toLowerCase(Locale.ROOT)) {
                case "search" -> search(
                    Arrays.copyOfRange(arguments, 1, arguments.length),
                    output, catalog);
                case "download" -> download(
                    Arrays.copyOfRange(arguments, 1, arguments.length),
                    output, currentDirectory, catalog);
                case "format" -> format(
                    Arrays.copyOfRange(arguments, 1, arguments.length),
                    output, currentDirectory);
                case "--help", "-h", "help" -> {
                    output.println(help());
                    yield 0;
                }
                default -> throw new IllegalArgumentException(
                    "unknown mr-instruments command: " + arguments[0]);
            };
        } catch (HelpRequested request) {
            output.println(request.getMessage());
            return 0;
        } catch (IllegalArgumentException exception) {
            error(errorOutput, exception.getMessage());
            errorOutput.println("Use mr-instruments --help for usage.");
            return 2;
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            error(errorOutput, "MR instrument request was interrupted");
            return 1;
        } catch (IOException exception) {
            error(errorOutput, exception.getMessage());
            return 1;
        }
    }

    private static int search(String[] arguments, PrintStream output,
            InstrumentCatalog catalog)
            throws IOException, InterruptedException {
        String trait = null;
        int limit = 20;
        for (int index = 0; index < arguments.length; index++) {
            String option = arguments[index];
            switch (option) {
                case "--trait" -> trait = value(arguments, ++index, option);
                case "--limit" -> limit = positiveInteger(
                    value(arguments, ++index, option), option, 200);
                case "--help", "-h" -> {
                    output.println(searchHelp());
                    return 0;
                }
                default -> throw new IllegalArgumentException(
                    "unknown mr-instruments search option: " + option);
            }
        }
        if (trait == null || trait.isBlank())
            throw new IllegalArgumentException("--trait is required");
        String resolvedTrait = expandTrait(trait);
        List<InstrumentStudy> studies = catalog.search(resolvedTrait, limit);
        output.println("accession\ttrait\tontology_traits\tancestry\t"
            + "sample_description\tvariant_count\tlicense\tsummary_statistics");
        for (InstrumentStudy study : studies) output.println(String.join("\t",
            tsv(study.accession()), tsv(study.trait()),
            tsv(String.join("|", study.ontologyTraits())),
            tsv(String.join("|", study.ancestry())),
            tsv(study.sampleDescription()),
            study.variantCount() < 0 ? "" : Long.toString(study.variantCount()),
            study.license() == null ? "" : tsv(study.license().toString()),
            tsv(study.summaryStatistics().toString())));
        if (!resolvedTrait.equalsIgnoreCase(trait))
            output.println("# expanded trait: " + trait + " -> " + resolvedTrait);
        output.println("# " + studies.size()
            + " downloadable GWAS Catalog studies shown");
        return 0;
    }

    private static int download(String[] arguments, PrintStream output,
            Path currentDirectory, InstrumentCatalog catalog)
            throws IOException, InterruptedException {
        CommonOptions options = parseCommon(arguments, true);
        if (options.study == null)
            throw new IllegalArgumentException("--study is required");
        if (options.output == null)
            throw new IllegalArgumentException("--out is required");
        Path destination = resolve(currentDirectory, options.output);
        InstrumentStudy study = catalog.study(options.study);
        output.println("Study: " + study.accession() + " - " + study.trait());
        if (study.license() != null)
            output.println("License: " + study.license());
        try (InstrumentCatalog.RemoteSummary remote =
                catalog.openSummaryStatistics(study)) {
            output.println("Source: " + remote.source());
            InstrumentTableFormatter.Result result = write(destination,
                remote.input(), new InstrumentTableFormatter.Options('\t',
                    options.mappings, options.phenotype == null
                        ? study.trait() : options.phenotype,
                    options.pThreshold, options.effectScale), options.overwrite);
            report(output, destination, result);
        }
        output.println("Reminder: p-value filtering selects candidates; perform "
            + "ancestry-matched LD clumping before independent-instrument MR.");
        return 0;
    }

    private static int format(String[] arguments, PrintStream output,
            Path currentDirectory) throws IOException {
        CommonOptions options = parseCommon(arguments, false);
        if (options.input == null)
            throw new IllegalArgumentException("--input is required");
        if (options.output == null)
            throw new IllegalArgumentException("--out is required");
        Path source = resolve(currentDirectory, options.input);
        Path destination = resolve(currentDirectory, options.output);
        char delimiter = options.delimiter == null
            ? inferredDelimiter(source) : options.delimiter;
        try (InputStream input = Files.newInputStream(source)) {
            InstrumentTableFormatter.Result result = write(destination, input,
                new InstrumentTableFormatter.Options(delimiter,
                    options.mappings, options.phenotype,
                    options.pThreshold, options.effectScale), options.overwrite);
            report(output, destination, result);
        }
        return 0;
    }

    private static CommonOptions parseCommon(String[] arguments,
            boolean download) {
        CommonOptions result = new CommonOptions();
        result.pThreshold = download ? DEFAULT_INSTRUMENT_P : 1.0;
        for (int index = 0; index < arguments.length; index++) {
            String option = arguments[index];
            switch (option) {
                case "--study" -> result.study =
                    value(arguments, ++index, option).toUpperCase(Locale.ROOT);
                case "--input" -> result.input = value(arguments, ++index, option);
                case "--out" -> result.output = value(arguments, ++index, option);
                case "--trait" -> result.phenotype = value(arguments, ++index, option);
                case "--p-threshold", "--pval-threshold" -> result.pThreshold =
                    probability(value(arguments, ++index, option), option);
                case "--delimiter" -> result.delimiter = delimiter(
                    value(arguments, ++index, option));
                case "--effect-scale" -> result.effectScale = effectScale(
                    value(arguments, ++index, option));
                case "--map" -> mappings(result.mappings,
                    value(arguments, ++index, option));
                case "--overwrite" -> result.overwrite = true;
                case "--help", "-h" -> throw new HelpRequested(
                    download ? downloadHelp() : formatHelp());
                default -> throw new IllegalArgumentException("unknown mr-instruments "
                    + (download ? "download" : "format") + " option: " + option);
            }
        }
        return result;
    }

    private static InstrumentTableFormatter.Result write(Path destination,
            InputStream input, InstrumentTableFormatter.Options options,
            boolean overwrite) throws IOException {
        if (Files.exists(destination) && !overwrite)
            throw new IOException("output exists; use --overwrite: " + destination);
        Path parent = destination.toAbsolutePath().normalize().getParent();
        if (parent != null) Files.createDirectories(parent);
        Path part = destination.resolveSibling("." + destination.getFileName()
            + ".part");
        Files.deleteIfExists(part);
        boolean complete = false;
        try {
            InstrumentTableFormatter.Result result;
            try (OutputStream output = Files.newOutputStream(part,
                    StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE)) {
                result = InstrumentTableFormatter.format(input, output, options);
            }
            move(part, destination, overwrite);
            complete = true;
            return result;
        } finally {
            if (!complete) Files.deleteIfExists(part);
        }
    }

    private static void move(Path source, Path destination, boolean overwrite)
            throws IOException {
        StandardCopyOption[] options = overwrite
            ? new StandardCopyOption[] {StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING}
            : new StandardCopyOption[] {StandardCopyOption.ATOMIC_MOVE};
        try {
            Files.move(source, destination, options);
        } catch (AtomicMoveNotSupportedException exception) {
            if (overwrite) Files.move(source, destination,
                StandardCopyOption.REPLACE_EXISTING);
            else Files.move(source, destination);
        }
    }

    private static void report(PrintStream output, Path destination,
            InstrumentTableFormatter.Result result) {
        output.println("Wrote " + result.outputRows() + " MR-ready rows to "
            + destination.toAbsolutePath().normalize());
        output.println("Input rows: " + result.inputRows()
            + "; above p-value threshold: " + result.aboveThreshold()
            + "; invalid/non-biallelic rows: " + result.invalidRows()
            + "; effect scale: " + scaleName(result.effectScale()));
    }

    private static void mappings(Map<String, String> target, String value) {
        for (String assignment : value.split(",", -1)) {
            int equals = assignment.indexOf('=');
            if (equals < 1 || equals == assignment.length() - 1)
                throw new IllegalArgumentException(
                    "--map requires TARGET=SOURCE assignments");
            String canonical = canonical(assignment.substring(0, equals).trim());
            String source = assignment.substring(equals + 1).trim();
            if (target.put(canonical, source) != null)
                throw new IllegalArgumentException(
                    "duplicate mapping for " + canonical);
        }
    }

    private static String canonical(String value) {
        String normalized = value.replaceAll("[^A-Za-z0-9]", "")
            .toLowerCase(Locale.ROOT);
        for (String column : InstrumentTableFormatter.COLUMNS)
            if (column.replaceAll("[^A-Za-z0-9]", "")
                    .toLowerCase(Locale.ROOT).equals(normalized)) return column;
        throw new IllegalArgumentException("unknown MR target column: " + value);
    }

    private static char inferredDelimiter(Path path) {
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        if (name.endsWith(".gz")) name = name.substring(0, name.length() - 3);
        return name.endsWith(".csv") ? ',' : '\t';
    }

    private static char delimiter(String value) {
        return switch (value.toLowerCase(Locale.ROOT)) {
            case "tab", "tsv", "\\t" -> '\t';
            case "comma", "csv", "," -> ',';
            default -> throw new IllegalArgumentException(
                "--delimiter must be tab or comma");
        };
    }

    private static InstrumentTableFormatter.EffectScale effectScale(String value) {
        return switch (value.toLowerCase(Locale.ROOT)) {
            case "auto" -> InstrumentTableFormatter.EffectScale.AUTO;
            case "beta" -> InstrumentTableFormatter.EffectScale.BETA;
            case "odds-ratio", "or" ->
                InstrumentTableFormatter.EffectScale.ODDS_RATIO;
            case "hazard-ratio", "hr" ->
                InstrumentTableFormatter.EffectScale.HAZARD_RATIO;
            default -> throw new IllegalArgumentException(
                "--effect-scale must be auto, beta, odds-ratio, or hazard-ratio");
        };
    }

    private static double probability(String value, String option) {
        try {
            double result = Double.parseDouble(value);
            if (!Double.isFinite(result) || result < 0.0 || result > 1.0)
                throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                option + " requires a number in [0, 1]: " + value, exception);
        }
    }

    private static int positiveInteger(String value, String option, int maximum) {
        try {
            int result = Integer.parseInt(value);
            if (result < 1 || result > maximum) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(option + " requires an integer from 1 to "
                + maximum + ": " + value, exception);
        }
    }

    private static Path resolve(Path currentDirectory, String value) {
        Path path = Path.of(value);
        return (path.isAbsolute() ? path : currentDirectory.resolve(path))
            .toAbsolutePath().normalize();
    }

    private static String expandTrait(String value) {
        return TRAIT_EXPANSIONS.getOrDefault(
            value.strip().toLowerCase(Locale.ROOT), value.strip());
    }

    private static String value(String[] values, int index, String option) {
        if (index >= values.length)
            throw new IllegalArgumentException(option + " requires a value");
        if (values[index].isBlank())
            throw new IllegalArgumentException(option + " requires a nonblank value");
        return values[index];
    }

    private static String tsv(String value) {
        if (value == null) return "";
        if (value.indexOf('\t') < 0 && value.indexOf('\n') < 0
                && value.indexOf('\r') < 0 && value.indexOf('"') < 0)
            return value;
        return '"' + value.replace("\"", "\"\"") + '"';
    }

    private static String scaleName(InstrumentTableFormatter.EffectScale scale) {
        return scale.name().toLowerCase(Locale.ROOT).replace('_', '-');
    }

    private static void error(PrintStream output, String message) {
        output.println("jlinalg: " + message);
    }

    private static String help() {
        return """
            Usage:
              java -jar jlinalg-<version>.jar mr-instruments search
                --trait TEXT [--limit N]
              java -jar jlinalg-<version>.jar mr-instruments download
                --study GCST... --out FILE [--p-threshold 5e-8]
              java -jar jlinalg-<version>.jar mr-instruments format
                --input FILE --out FILE [--map TARGET=SOURCE,...]

            search queries downloadable studies in the public NHGRI-EBI GWAS
            Catalog. Common acronyms such as BMI are expanded before searching.

            download streams a study's standardized summary statistics and emits
            genome-wide-significant, MR-ready candidate rows. It does not perform
            LD clumping.

            format accepts CSV/TSV (optionally gzip-compressed), auto-detects common
            GWAS/QTL column names, and writes MRInstruments-style canonical columns.
            Use --effect-scale odds-ratio or hazard-ratio to log-transform ratios.
            """;
    }

    private static String searchHelp() {
        return "Usage: mr-instruments search --trait TEXT [--limit 20]";
    }

    private static String downloadHelp() {
        return "Usage: mr-instruments download --study GCST... --out FILE "
            + "[--p-threshold 5e-8] [--trait TEXT] [--overwrite]";
    }

    private static String formatHelp() {
        return "Usage: mr-instruments format --input FILE --out FILE "
            + "[--map SNP=rsid,beta=effect,...] [--delimiter tab|comma] "
            + "[--effect-scale auto|beta|odds-ratio|hazard-ratio] "
            + "[--trait TEXT] [--p-threshold X] [--overwrite]";
    }

    private static final class CommonOptions {
        String study;
        String input;
        String output;
        String phenotype;
        Character delimiter;
        double pThreshold;
        InstrumentTableFormatter.EffectScale effectScale =
            InstrumentTableFormatter.EffectScale.AUTO;
        final Map<String, String> mappings = new LinkedHashMap<>();
        boolean overwrite;
    }

    private static final class HelpRequested extends IllegalArgumentException {
        private static final long serialVersionUID = 1L;
        HelpRequested(String message) { super(message); }
    }
}
