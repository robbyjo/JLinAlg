/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import org.jlinalg.coloc.ColocOptions;
import org.jlinalg.coloc.ColocSignalPair;
import org.jlinalg.coloc.ColocSusie;
import org.jlinalg.coloc.ColocSusieInput;
import org.jlinalg.coloc.ColocSusieResult;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.susie.CredibleSet;
import org.jlinalg.susie.Susie;
import org.jlinalg.susie.SusieOptions;
import org.jlinalg.susie.SusieResult;

/** File-oriented summary-statistic SuSiE and colocalization commands. */
final class FineMappingCli {
    private FineMappingCli() { }

    static int run(String command, String[] arguments,
            PrintStream output, PrintStream error) {
        try {
            if (command.equals("susie")) {
                SusieCliOptions options = SusieCliOptions.parse(arguments);
                if (options.help) {
                    output.println(susieHelp());
                    return 0;
                }
                runSusie(options, output);
            } else {
                ColocCliOptions options = ColocCliOptions.parse(arguments);
                if (options.help) {
                    output.println(colocHelp());
                    return 0;
                }
                runColoc(options, output);
            }
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static void runSusie(SusieCliOptions options, PrintStream console)
            throws IOException {
        Path effectsOutput = options.effectsOutput != null
            ? options.effectsOutput : Path.of(options.output + ".effects.tsv");
        Path log = options.log != null ? options.log
            : Path.of(options.output + ".log");
        prepareOutputs(options.overwrite,
            List.of(options.output, effectsOutput, log),
            List.of(options.summary, options.ld));
        SummaryInput summary = readSummary(options);
        double[][] ld = readLd(options.ld, summary.variantNames());
        SusieOptions fitOptions = new SusieOptions(options.effects,
            options.maximumIterations, options.tolerance,
            options.priorVariance, options.estimateResidualVariance,
            options.coverage, options.minimumPurity);
        SusieResult result = Susie.fitSummary(summary.zScores(), ld,
            options.sampleSize, summary.variantNames(), fitOptions,
            options.backend);
        write(options.output, variantTable(result, delimiter(options.output)),
            options.overwrite);
        write(effectsOutput, effectTable(result, delimiter(effectsOutput)),
            options.overwrite);
        String logText = "command=susie\n"
            + "summary=" + options.summary.toAbsolutePath() + "\n"
            + "ld=" + options.ld.toAbsolutePath() + "\n"
            + "variants=" + result.variableNames().size() + "\n"
            + "sample_size=" + (long) options.sampleSize + "\n"
            + "effects=" + result.effects() + "\n"
            + "credible_sets=" + result.credibleSets().size() + "\n"
            + "iterations=" + result.iterations() + "\n"
            + "converged=" + result.converged() + "\n"
            + "objective=" + result.objective() + "\n"
            + "residual_variance=" + result.residualVariance() + "\n"
            + "requested_backend=" + result.backend().requested() + "\n"
            + "selected_backend=" + result.backend().selectedBackend() + "\n"
            + "output=" + options.output.toAbsolutePath() + "\n"
            + "effects_output=" + effectsOutput.toAbsolutePath() + "\n";
        write(log, logText, options.overwrite);
        console.println("SuSiE variants: " + result.variableNames().size()
            + "; credible sets: " + result.credibleSets().size()
            + "; converged: " + result.converged());
        console.println("SuSiE output: " + options.output);
        console.println("Colocalization input: " + effectsOutput);
        console.println("Log: " + log);
    }

    private static void runColoc(ColocCliOptions options, PrintStream console)
            throws IOException {
        Path variantOutput = options.variantOutput != null
            ? options.variantOutput : Path.of(options.output + ".variants.tsv");
        Path log = options.log != null ? options.log
            : Path.of(options.output + ".log");
        prepareOutputs(options.overwrite,
            List.of(options.output, variantOutput, log),
            List.of(options.trait1, options.trait2));
        ColocSusieInput trait1 = readColocInput(options.trait1);
        ColocSusieInput trait2 = readColocInput(options.trait2);
        ColocOptions colocOptions = new ColocOptions(
            options.trait1Prior, options.trait2Prior, options.sharedPrior,
            options.minimumOverlap, options.trim, null, null);
        ColocSusieResult result = ColocSusie.analyze(
            trait1, trait2, colocOptions);
        if (result.commonVariants().isEmpty()) {
            throw new IllegalArgumentException(
                "trait inputs share no variant identifiers");
        }
        write(options.output, colocTable(result, delimiter(options.output)),
            options.overwrite);
        write(variantOutput,
            sharedVariantTable(result, delimiter(variantOutput)),
            options.overwrite);
        String logText = "command=coloc\n"
            + "trait1=" + options.trait1.toAbsolutePath() + "\n"
            + "trait2=" + options.trait2.toAbsolutePath() + "\n"
            + "trait1_signals=" + trait1.signals() + "\n"
            + "trait2_signals=" + trait2.signals() + "\n"
            + "common_variants=" + result.commonVariants().size() + "\n"
            + "retained_signal_pairs=" + result.signalPairs().size() + "\n"
            + "skipped_signal_pairs=" + result.skippedSignalPairs() + "\n"
            + "trait1_prior=" + options.trait1Prior + "\n"
            + "trait2_prior=" + options.trait2Prior + "\n"
            + "shared_prior=" + options.sharedPrior + "\n"
            + "minimum_posterior_overlap=" + options.minimumOverlap + "\n"
            + "trim_by_posterior=" + options.trim + "\n"
            + "output=" + options.output.toAbsolutePath() + "\n"
            + "variant_output=" + variantOutput.toAbsolutePath() + "\n";
        write(log, logText, options.overwrite);
        console.println("Colocalization common variants: "
            + result.commonVariants().size() + "; retained signal pairs: "
            + result.signalPairs().size() + "; skipped: "
            + result.skippedSignalPairs());
        console.println("Colocalization output: " + options.output);
        console.println("Shared-variant output: " + variantOutput);
        console.println("Log: " + log);
    }

    private static SummaryInput readSummary(SusieCliOptions options)
            throws IOException {
        DelimitedData table = DelimitedData.read(options.summary);
        int variant = column(table.header(), options.variantColumn,
            "variant_id", "SNP", "rsid", "id");
        int z = options.betaColumn == null
            ? optionalColumn(table.header(), options.zColumn,
                "z", "zscore", "z_score")
            : -1;
        int beta = optionalColumn(table.header(), options.betaColumn,
            "beta", "effect", "estimate");
        int se = optionalColumn(table.header(), options.seColumn,
            "se", "standard_error", "stderr");
        if (z < 0 && (beta < 0 || se < 0)) {
            throw new IllegalArgumentException(
                "summary input needs a z column or both beta and se columns; "
                    + "use --z-column or --beta-column/--se-column");
        }
        List<String> names = new ArrayList<>();
        double[] zScores = new double[table.rows().size()];
        Set<String> unique = new HashSet<>();
        for (int row = 0; row < table.rows().size(); row++) {
            String[] fields = table.rows().get(row);
            String name = fields[variant].trim();
            if (name.isEmpty() || !unique.add(name)) {
                throw new IllegalArgumentException(
                    "variant identifiers must be nonblank and unique: " + name);
            }
            names.add(name);
            if (z >= 0) {
                zScores[row] = finite(fields[z], "z", row + 2);
            } else {
                double effect = finite(fields[beta], "beta", row + 2);
                double error = finite(fields[se], "se", row + 2);
                if (!(error > 0.0)) {
                    throw new IllegalArgumentException(
                        "se must be positive at row " + (row + 2));
                }
                zScores[row] = effect / error;
            }
        }
        return new SummaryInput(List.copyOf(names), zScores);
    }

    private static double[][] readLd(Path path, List<String> variants)
            throws IOException {
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8)
            .stream().filter(line -> !line.isBlank()).toList();
        if (lines.isEmpty()) throw new IOException("LD matrix is empty: " + path);
        char delimiter = delimiter(path);
        String[] first = split(lines.get(0), delimiter);
        if (first.length == 0) throw new IOException("LD matrix is empty: " + path);
        if (isNumber(first[0])) {
            int size = variants.size();
            if (lines.size() != size) {
                throw new IOException("unlabeled LD matrix row count "
                    + lines.size() + " does not equal summary variant count "
                    + size);
            }
            double[][] result = new double[size][size];
            for (int row = 0; row < size; row++) {
                String[] fields = split(lines.get(row), delimiter);
                if (fields.length != size) {
                    throw new IOException("unlabeled LD matrix must be square; "
                        + "row " + (row + 1) + " has " + fields.length
                        + " fields");
                }
                for (int column = 0; column < size; column++) {
                    result[row][column] = finite(fields[column], "LD",
                        row + 1);
                }
            }
            return result;
        }
        if (first.length < 2) {
            throw new IOException(
                "labeled LD matrix needs a row-ID column and variant columns");
        }
        List<String> columns = Arrays.stream(first).skip(1)
            .map(String::trim).toList();
        int size = columns.size();
        if (lines.size() != size + 1) {
            throw new IOException("labeled LD matrix must have one data row "
                + "per variant column");
        }
        Map<String, double[]> rows = new LinkedHashMap<>();
        for (int row = 0; row < size; row++) {
            String[] fields = split(lines.get(row + 1), delimiter);
            if (fields.length != size + 1) {
                throw new IOException("labeled LD row " + (row + 2)
                    + " has the wrong field count");
            }
            String name = fields[0].trim();
            double[] values = new double[size];
            for (int column = 0; column < size; column++) {
                values[column] = finite(fields[column + 1], "LD", row + 2);
            }
            if (name.isEmpty() || rows.put(name, values) != null) {
                throw new IOException(
                    "LD row identifiers must be nonblank and unique");
            }
        }
        if (new HashSet<>(columns).size() != columns.size()
                || !rows.keySet().equals(new LinkedHashSet<>(columns))) {
            throw new IOException(
                "labeled LD row and column variant identifiers must match");
        }
        if (!rows.keySet().containsAll(variants)
                || rows.size() != variants.size()) {
            throw new IOException(
                "labeled LD and summary variant identifier sets must match");
        }
        Map<String, Integer> columnIndex = new HashMap<>();
        for (int index = 0; index < columns.size(); index++) {
            columnIndex.put(columns.get(index), index);
        }
        double[][] result = new double[variants.size()][variants.size()];
        for (int row = 0; row < variants.size(); row++) {
            double[] source = rows.get(variants.get(row));
            for (int column = 0; column < variants.size(); column++) {
                result[row][column] = source[columnIndex.get(
                    variants.get(column))];
            }
        }
        return result;
    }

    private static ColocSusieInput readColocInput(Path path)
            throws IOException {
        DelimitedData table = DelimitedData.read(path);
        if (table.rows().isEmpty()) {
            throw new IllegalArgumentException(
                "colocalization input has no signal rows: " + path);
        }
        int effect = column(table.header(), null, "effect_index", "effect");
        int variant = column(table.header(), null,
            "variant_id", "SNP", "rsid", "id");
        int logBf = column(table.header(), null,
            "log_bayes_factor", "lbf", "logbf");
        LinkedHashMap<Integer, LinkedHashMap<String, Double>> groups =
            new LinkedHashMap<>();
        for (int row = 0; row < table.rows().size(); row++) {
            String[] fields = table.rows().get(row);
            int effectIndex;
            try {
                effectIndex = Integer.parseInt(fields[effect].trim());
            } catch (NumberFormatException exception) {
                throw new IllegalArgumentException(
                    "invalid effect_index at row " + (row + 2), exception);
            }
            if (effectIndex < 0) {
                throw new IllegalArgumentException(
                    "effect_index must be nonnegative at row " + (row + 2));
            }
            String variantId = fields[variant].trim();
            if (variantId.isEmpty()) {
                throw new IllegalArgumentException(
                    "variant_id is blank at row " + (row + 2));
            }
            double value = logBayesFactor(fields[logBf], row + 2);
            LinkedHashMap<String, Double> values = groups.computeIfAbsent(
                effectIndex, ignored -> new LinkedHashMap<>());
            if (values.put(variantId, value) != null) {
                throw new IllegalArgumentException(
                    "duplicate effect/variant row at row " + (row + 2));
            }
        }
        List<String> variants = List.copyOf(
            groups.values().iterator().next().keySet());
        double[] values = new double[groups.size() * variants.size()];
        int[] effects = new int[groups.size()];
        int group = 0;
        Set<String> expected = new LinkedHashSet<>(variants);
        for (Map.Entry<Integer, LinkedHashMap<String, Double>> entry
                : groups.entrySet()) {
            if (!entry.getValue().keySet().equals(expected)) {
                throw new IllegalArgumentException(
                    "every effect in " + path
                        + " must contain the same variant identifiers");
            }
            effects[group] = entry.getKey();
            for (int variantIndex = 0;
                    variantIndex < variants.size(); variantIndex++) {
                values[group * variants.size() + variantIndex] =
                    entry.getValue().get(variants.get(variantIndex));
            }
            group++;
        }
        return new ColocSusieInput(variants, values, effects);
    }

    private static String variantTable(SusieResult result, char delimiter) {
        StringBuilder text = new StringBuilder();
        row(text, delimiter, "variant_id", "pip", "posterior_mean");
        double[] pip = result.pip();
        double[] mean = result.posteriorMean();
        for (int index = 0; index < pip.length; index++) {
            row(text, delimiter, result.variableNames().get(index),
                Double.toString(pip[index]), Double.toString(mean[index]));
        }
        return text.toString();
    }

    private static String effectTable(SusieResult result, char delimiter) {
        StringBuilder text = new StringBuilder();
        row(text, delimiter, "effect_index", "variant_id", "alpha",
            "effect_posterior_mean", "log_bayes_factor",
            "in_credible_set", "credible_set_coverage",
            "credible_set_purity");
        double[] alpha = result.alpha();
        double[] mean = result.effectPosteriorMean();
        double[] logBf = result.logBayesFactors();
        int variants = result.variableNames().size();
        for (CredibleSet credibleSet : result.credibleSets()) {
            int effect = credibleSet.effectIndex();
            Set<String> members = new HashSet<>(credibleSet.variables());
            for (int variant = 0; variant < variants; variant++) {
                int index = effect * variants + variant;
                String name = result.variableNames().get(variant);
                row(text, delimiter, Integer.toString(effect + 1), name,
                    Double.toString(alpha[index]),
                    Double.toString(mean[index]),
                    Double.toString(logBf[index]),
                    Boolean.toString(members.contains(name)),
                    Double.toString(credibleSet.posteriorCoverage()),
                    Double.toString(credibleSet.purity()));
            }
        }
        return text.toString();
    }

    private static String colocTable(ColocSusieResult result, char delimiter) {
        StringBuilder text = new StringBuilder();
        row(text, delimiter, "signal_pair", "trait1_effect_index",
            "trait2_effect_index", "variants", "trait1_lead_variant",
            "trait2_lead_variant", "posterior_h0", "posterior_h1",
            "posterior_h2", "posterior_h3", "posterior_h4");
        for (int index = 0; index < result.signalPairs().size(); index++) {
            ColocSignalPair pair = result.signalPairs().get(index);
            row(text, delimiter, Integer.toString(index + 1),
                Integer.toString(pair.trait1EffectIndex()),
                Integer.toString(pair.trait2EffectIndex()),
                Integer.toString(pair.variants()),
                pair.trait1LeadVariant(), pair.trait2LeadVariant(),
                Double.toString(pair.posteriorH0()),
                Double.toString(pair.posteriorH1()),
                Double.toString(pair.posteriorH2()),
                Double.toString(pair.posteriorH3()),
                Double.toString(pair.posteriorH4()));
        }
        return text.toString();
    }

    private static String sharedVariantTable(
            ColocSusieResult result, char delimiter) {
        StringBuilder text = new StringBuilder();
        row(text, delimiter, "signal_pair", "trait1_effect_index",
            "trait2_effect_index", "variant_id",
            "posterior_shared_given_h4");
        for (int pairIndex = 0;
                pairIndex < result.signalPairs().size(); pairIndex++) {
            ColocSignalPair pair = result.signalPairs().get(pairIndex);
            double[] posterior = result.sharedVariantPosterior(pairIndex);
            for (int variant = 0;
                    variant < result.commonVariants().size(); variant++) {
                row(text, delimiter, Integer.toString(pairIndex + 1),
                    Integer.toString(pair.trait1EffectIndex()),
                    Integer.toString(pair.trait2EffectIndex()),
                    result.commonVariants().get(variant),
                    Double.toString(posterior[variant]));
            }
        }
        return text.toString();
    }

    private static void prepareOutputs(boolean overwrite,
            List<Path> outputs, List<Path> inputs) throws IOException {
        for (int index = 0; index < outputs.size(); index++) {
            Path target = outputs.get(index).toAbsolutePath().normalize();
            for (Path input : inputs) {
                if (target.equals(input.toAbsolutePath().normalize())) {
                    throw new IllegalArgumentException(
                        "output path aliases an input: " + outputs.get(index));
                }
            }
            for (int previous = 0; previous < index; previous++) {
                if (target.equals(
                        outputs.get(previous).toAbsolutePath().normalize())) {
                    throw new IllegalArgumentException(
                        "output paths must be distinct: " + outputs.get(index));
                }
            }
        }
        if (!overwrite) {
            PipelinePaths.requireFreshOutputs(outputs.toArray(Path[]::new));
        }
    }

    private static void write(Path path, String text, boolean overwrite)
            throws IOException {
        if (overwrite) {
            Files.writeString(path, text, StandardCharsets.UTF_8,
                StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
        } else {
            Files.writeString(path, text, StandardCharsets.UTF_8,
                StandardOpenOption.CREATE_NEW);
        }
    }

    private static int column(List<String> header, String requested,
            String... aliases) {
        int result = optionalColumn(header, requested, aliases);
        if (result < 0) {
            throw new IllegalArgumentException("input is missing column "
                + (requested == null ? Arrays.toString(aliases) : requested));
        }
        return result;
    }

    private static int optionalColumn(List<String> header, String requested,
            String... aliases) {
        if (requested != null) {
            for (int index = 0; index < header.size(); index++) {
                if (header.get(index).equalsIgnoreCase(requested)) return index;
            }
            throw new IllegalArgumentException(
                "column is absent: " + requested);
        }
        for (String alias : aliases) {
            for (int index = 0; index < header.size(); index++) {
                if (header.get(index).equalsIgnoreCase(alias)) return index;
            }
        }
        return -1;
    }

    private static double finite(String value, String label, int row) {
        try {
            double result = Double.parseDouble(value.trim());
            if (!Double.isFinite(result)) {
                throw new NumberFormatException("not finite");
            }
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                label + " is not finite at row " + row + ": " + value,
                exception);
        }
    }

    private static double logBayesFactor(String value, int row) {
        try {
            double result = Double.parseDouble(value.trim());
            if (Double.isNaN(result) || result == Double.POSITIVE_INFINITY) {
                throw new NumberFormatException("invalid log BF");
            }
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                "log_bayes_factor is invalid at row " + row + ": " + value,
                exception);
        }
    }

    private static boolean isNumber(String value) {
        try {
            Double.parseDouble(value.trim());
            return true;
        } catch (NumberFormatException ignored) {
            return false;
        }
    }

    private static char delimiter(Path path) {
        return path.toString().toLowerCase(Locale.ROOT).endsWith(".csv")
            ? ',' : '\t';
    }

    private static String[] split(String line, char delimiter) {
        return line.split(Pattern.quote(String.valueOf(delimiter)), -1);
    }

    private static void row(StringBuilder output, char delimiter,
            String... values) {
        for (int index = 0; index < values.length; index++) {
            if (index > 0) output.append(delimiter);
            output.append(escape(values[index], delimiter));
        }
        output.append('\n');
    }

    private static String escape(String value, char delimiter) {
        if (value.indexOf(delimiter) < 0 && value.indexOf('"') < 0
                && value.indexOf('\n') < 0 && value.indexOf('\r') < 0) {
            return value;
        }
        return "\"" + value.replace("\"", "\"\"") + "\"";
    }

    private static BackendPolicy backend(String value) {
        try {
            return BackendPolicy.valueOf(value.toUpperCase(Locale.ROOT));
        } catch (IllegalArgumentException failure) {
            throw new IllegalArgumentException(
                "--backend must be preferred, cholmod, gpu, cuda, opencl, "
                    + "vulkan, onemkl, openblas, auto, or cpu: " + value,
                failure);
        }
    }

    private static double number(String value, String option) {
        try {
            double result = Double.parseDouble(value);
            if (!Double.isFinite(result)) {
                throw new NumberFormatException("not finite");
            }
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                option + " requires a finite number: " + value, exception);
        }
    }

    private static double positive(String value, String option) {
        double result = number(value, option);
        if (!(result > 0.0)) {
            throw new IllegalArgumentException(option + " must be positive");
        }
        return result;
    }

    private static double probability(String value, String option,
            boolean allowZero) {
        double result = number(value, option);
        if ((allowZero ? result < 0.0 : result <= 0.0)
                || (allowZero ? result > 1.0 : result >= 1.0)) {
            throw new IllegalArgumentException(option
                + (allowZero ? " must be in [0,1]" : " must be in (0,1)"));
        }
        return result;
    }

    private static int positiveInteger(String value, String option) {
        try {
            int result = Integer.parseInt(value);
            if (result < 1) throw new NumberFormatException("not positive");
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                option + " requires a positive integer: " + value,
                exception);
        }
    }

    private static String value(String[] values, int index, String option) {
        if (index >= values.length || values[index].isBlank()) {
            throw new IllegalArgumentException(option + " requires a value");
        }
        return values[index];
    }

    private static String susieHelp() {
        return """
            Usage: susie --summary FILE --ld FILE --sample-size N --out FILE
              [--variant-column COLUMN] [--z-column COLUMN]
              [--beta-column COLUMN --se-column COLUMN]
              [--effects 10] [--max-iterations 200] [--tolerance 1e-6]
              [--prior-variance 0.2] [--fixed-residual-variance]
              [--coverage 0.95] [--min-purity 0.5]
              [--backend preferred|cpu|...] [--effects-output FILE]
              [--log FILE] [--overwrite]

            Summary input is CSV/TSV with a unique variant ID and either z, or
            beta and se. Unlabeled LD must be a square CSV/TSV matrix in exact
            summary-row order. A labeled LD matrix may instead use variant IDs
            in its first row and first column; it is aligned automatically.
            The main output contains variant PIPs and posterior means. The
            effects output contains credible-effect log Bayes factors and is
            accepted directly by the coloc command.
            """;
    }

    private static String colocHelp() {
        return """
            Usage: coloc --trait1 FILE --trait2 FILE --out FILE
              [--variant-output FILE]
              [--trait1-prior 1e-4] [--trait2-prior 1e-4]
              [--shared-prior 5e-6] [--minimum-overlap 0.5]
              [--no-trim] [--log FILE] [--overwrite]

            Trait inputs are the *.effects.tsv files produced by the susie
            command, or CSV/TSV tables with effect_index, variant_id, and
            log_bayes_factor. Variant IDs are intersected in trait-1 order.
            The main output reports H0-H4 for each retained signal pair. The
            variant output reports the variant posterior conditional on H4.
            """;
    }

    private record SummaryInput(List<String> variantNames, double[] zScores) {
        private SummaryInput {
            zScores = zScores.clone();
        }
        @Override public double[] zScores() { return zScores.clone(); }
    }

    private static final class SusieCliOptions {
        Path summary, ld, output, effectsOutput, log;
        String variantColumn, zColumn, betaColumn, seColumn;
        double sampleSize;
        int effects = 10;
        int maximumIterations = 200;
        double tolerance = 1e-6;
        double priorVariance = 0.2;
        boolean estimateResidualVariance = true;
        double coverage = 0.95;
        double minimumPurity = 0.5;
        BackendPolicy backend = BackendPolicy.PREFERRED;
        boolean overwrite, help;

        static SusieCliOptions parse(String[] arguments) {
            SusieCliOptions result = new SusieCliOptions();
            for (int index = 0; index < arguments.length; index++) {
                String option = arguments[index];
                switch (option) {
                    case "--summary", "--input", "--in" ->
                        result.summary = Path.of(value(
                            arguments, ++index, option));
                    case "--ld", "--ld-correlation" ->
                        result.ld = Path.of(value(arguments, ++index, option));
                    case "--sample-size" ->
                        result.sampleSize = positive(value(
                            arguments, ++index, option), option);
                    case "--out", "--output" ->
                        result.output = Path.of(value(
                            arguments, ++index, option));
                    case "--effects-output" ->
                        result.effectsOutput = Path.of(value(
                            arguments, ++index, option));
                    case "--log" ->
                        result.log = Path.of(value(arguments, ++index, option));
                    case "--variant-column" ->
                        result.variantColumn = value(
                            arguments, ++index, option);
                    case "--z-column" ->
                        result.zColumn = value(arguments, ++index, option);
                    case "--beta-column" ->
                        result.betaColumn = value(arguments, ++index, option);
                    case "--se-column" ->
                        result.seColumn = value(arguments, ++index, option);
                    case "--effects" ->
                        result.effects = positiveInteger(value(
                            arguments, ++index, option), option);
                    case "--max-iterations" ->
                        result.maximumIterations = positiveInteger(value(
                            arguments, ++index, option), option);
                    case "--tolerance" ->
                        result.tolerance = positive(value(
                            arguments, ++index, option), option);
                    case "--prior-variance" ->
                        result.priorVariance = positive(value(
                            arguments, ++index, option), option);
                    case "--fixed-residual-variance" ->
                        result.estimateResidualVariance = false;
                    case "--coverage" ->
                        result.coverage = probability(value(
                            arguments, ++index, option), option, false);
                    case "--min-purity" ->
                        result.minimumPurity = probability(value(
                            arguments, ++index, option), option, true);
                    case "--backend" ->
                        result.backend = backend(value(
                            arguments, ++index, option));
                    case "--overwrite" -> result.overwrite = true;
                    case "--help", "-h" -> result.help = true;
                    default -> throw new IllegalArgumentException(
                        "unknown susie option: " + option);
                }
            }
            if (!result.help && (result.summary == null || result.ld == null
                    || result.sampleSize == 0.0 || result.output == null)) {
                throw new IllegalArgumentException(
                    "--summary, --ld, --sample-size, and --out are required");
            }
            if (result.sampleSize != Math.rint(result.sampleSize)) {
                throw new IllegalArgumentException(
                    "--sample-size must be an integer");
            }
            if ((result.betaColumn == null) != (result.seColumn == null)) {
                throw new IllegalArgumentException(
                    "--beta-column and --se-column must be specified together");
            }
            if (result.zColumn != null && result.betaColumn != null) {
                throw new IllegalArgumentException(
                    "choose --z-column or --beta-column/--se-column, not both");
            }
            return result;
        }
    }

    private static final class ColocCliOptions {
        Path trait1, trait2, output, variantOutput, log;
        double trait1Prior = 1e-4;
        double trait2Prior = 1e-4;
        double sharedPrior = 5e-6;
        double minimumOverlap = 0.5;
        boolean trim = true;
        boolean overwrite, help;

        static ColocCliOptions parse(String[] arguments) {
            ColocCliOptions result = new ColocCliOptions();
            for (int index = 0; index < arguments.length; index++) {
                String option = arguments[index];
                switch (option) {
                    case "--trait1" -> result.trait1 = Path.of(value(
                        arguments, ++index, option));
                    case "--trait2" -> result.trait2 = Path.of(value(
                        arguments, ++index, option));
                    case "--out", "--output" -> result.output = Path.of(value(
                        arguments, ++index, option));
                    case "--variant-output" ->
                        result.variantOutput = Path.of(value(
                            arguments, ++index, option));
                    case "--log" -> result.log = Path.of(value(
                        arguments, ++index, option));
                    case "--trait1-prior" ->
                        result.trait1Prior = probability(value(
                            arguments, ++index, option), option, false);
                    case "--trait2-prior" ->
                        result.trait2Prior = probability(value(
                            arguments, ++index, option), option, false);
                    case "--shared-prior" ->
                        result.sharedPrior = probability(value(
                            arguments, ++index, option), option, false);
                    case "--minimum-overlap" ->
                        result.minimumOverlap = probability(value(
                            arguments, ++index, option), option, true);
                    case "--no-trim" -> result.trim = false;
                    case "--overwrite" -> result.overwrite = true;
                    case "--help", "-h" -> result.help = true;
                    default -> throw new IllegalArgumentException(
                        "unknown coloc option: " + option);
                }
            }
            if (!result.help && (result.trait1 == null
                    || result.trait2 == null || result.output == null)) {
                throw new IllegalArgumentException(
                    "--trait1, --trait2, and --out are required");
            }
            return result;
        }
    }
}
