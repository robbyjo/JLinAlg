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
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mediation.MediationAnalysis;
import org.jlinalg.mediation.MediationEffect;
import org.jlinalg.mediation.MediationMixedResult;
import org.jlinalg.mediation.MediationResult;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.RankDeficiencyStrategy;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.reml.RemlOptions;

/** File-oriented Gaussian OLS, grouped-REML, and pedigree mediation. */
final class MediationCli {
    private MediationCli() { }

    static int run(String[] arguments, PrintStream output, PrintStream error) {
        try {
            Options options = Options.parse(arguments);
            if (options.help) {
                output.println(help());
                return 0;
            }
            Path log = options.log == null
                ? Path.of(options.output + ".log") : options.log;
            prepareOutputs(options, log);
            Input input = Input.read(options);
            Fit fit = fit(input, options);
            write(options.output,
                effects(fit, delimiter(options.output)), options.overwrite);
            write(log, runLog(input, fit, options), options.overwrite);
            output.println("Analysis samples: " + input.rows()
                + " (input=" + input.originalRows()
                + ", missing omitted="
                + (input.originalRows() - input.rows()) + ")");
            if (fit.singletonFamilies() > 0) {
                output.println("Pedigree singletons: "
                    + fit.singletonFamilies() + " families across "
                    + fit.singletonObservations() + " observations");
            }
            output.println("Mediation model: " + fit.model()
                + "; converged: " + fit.converged());
            output.println("Wrote " + options.output);
            output.println("Wrote " + log);
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static Fit fit(Input input, Options options) throws IOException {
        if (options.pedigree != null) {
            PedigreeReader.Loaded loaded = PedigreeReader.read(
                options.pedigree, options.pedigreeId, options.sireId,
                options.damId, options.pedigreeFamilyId);
            List<PedigreeIndividual> members =
                new ArrayList<>(loaded.individuals());
            List<String> individualIds = new ArrayList<>(
                input.individualIds().size());
            int aliasesResolved = 0;
            for (String id : input.individualIds()) {
                String resolved = loaded.resolveObservationId(id);
                individualIds.add(resolved);
                if (!resolved.equals(id.trim())) aliasesResolved++;
            }
            Set<String> known = new HashSet<>();
            for (PedigreeIndividual member : members) known.add(member.id());
            Set<String> singletons = new LinkedHashSet<>();
            int singletonObservations = 0;
            int matchedObservations = 0;
            for (String id : individualIds) {
                if (known.contains(id)) {
                    matchedObservations++;
                } else {
                    singletonObservations++;
                    if (singletons.add(id)) {
                        members.add(PedigreeIndividual.founder(id));
                    }
                }
            }
            if (!individualIds.isEmpty() && matchedObservations == 0)
                throw new IllegalArgumentException(
                    "no input observations match pedigree members; check "
                        + "--individual-id and pedigree ID qualification");
            double[] inbreeding = Arrays.copyOf(
                loaded.inbreedingCoefficients(), members.size());
            PedigreeRandomEffectTerm pedigree =
                PedigreeRandomEffectTerm.ofSparse("pedigree",
                    individualIds, members, inbreeding);
            MediationMixedResult result = MediationAnalysis.fitPedigree(
                input.outcome(), input.treatment(), input.mediator(),
                input.covariates(), List.of(pedigree),
                input.randomEffects(), RemlOptions.defaults(),
                options.confidence, options.backend);
            return Fit.mixed(result, "pedigree-reml",
                loaded.individuals().size(), matchedObservations,
                aliasesResolved, singletons.size(), singletonObservations);
        }
        if (!input.randomEffects().isEmpty()) {
            MediationMixedResult result = MediationAnalysis.fitMixed(
                input.outcome(), input.treatment(), input.mediator(),
                input.covariates(), input.randomEffects(),
                RemlOptions.defaults(), options.confidence, options.backend);
            return Fit.mixed(result, "grouped-reml", 0, 0);
        }
        OlsOptions ols = new OlsOptions(RankDeficiencyStrategy.ERROR,
            options.confidence, MissingDataPolicy.ERROR);
        MediationResult result = MediationAnalysis.fit(
            input.outcome(), input.treatment(), input.mediator(),
            input.covariates(), ols, options.backend);
        return Fit.ols(result);
    }

    private static String effects(Fit fit, char delimiter) {
        StringBuilder result = new StringBuilder();
        row(result, delimiter, "effect", "estimate", "standard_error",
            "statistic", "p_value", "ci_lower", "ci_upper",
            "degrees_of_freedom");
        for (MediationEffect effect : fit.effects()) {
            row(result, delimiter, effect.name(),
                Double.toString(effect.estimate()),
                Double.toString(effect.standardError()),
                Double.toString(effect.statistic()),
                Double.toString(effect.pValue()),
                Double.toString(effect.confidenceLower()),
                Double.toString(effect.confidenceUpper()),
                Double.toString(effect.degreesOfFreedom()));
        }
        return result.toString();
    }

    private static String runLog(Input input, Fit fit, Options options) {
        return "command=mediation\n"
            + "input=" + options.input.toAbsolutePath() + "\n"
            + "model=" + fit.model() + "\n"
            + "outcome=" + options.outcome + "\n"
            + "treatment=" + options.treatment + "\n"
            + "mediator=" + options.mediator + "\n"
            + "covariates=" + String.join(",", options.covariates) + "\n"
            + "groups=" + String.join(",", options.groups) + "\n"
            + "input_rows=" + input.originalRows() + "\n"
            + "analysis_rows=" + input.rows() + "\n"
            + "missing_rows_omitted="
            + (input.originalRows() - input.rows()) + "\n"
            + "pedigree_file_members=" + fit.pedigreeFileMembers() + "\n"
            + "pedigree_file_observations_matched="
            + fit.pedigreeFileObservationsMatched() + "\n"
            + "pedigree_unqualified_aliases_resolved="
            + fit.pedigreeAliasesResolved() + "\n"
            + "pedigree_singletons_added=" + fit.singletonFamilies() + "\n"
            + "pedigree_singleton_observations="
            + fit.singletonObservations() + "\n"
            + "confidence_level=" + options.confidence + "\n"
            + "backend=" + options.backend + "\n"
            + "converged=" + fit.converged() + "\n"
            + "output=" + options.output.toAbsolutePath() + "\n";
    }

    private record Input(
            double[] outcome, double[] treatment, double[] mediator,
            double[][] covariates, List<RandomEffectTerm> randomEffects,
            List<String> individualIds, int rows, int originalRows) {
        static Input read(Options options) throws IOException {
            DelimitedData table = DelimitedData.read(options.input);
            int outcomeColumn = column(table, options.outcome);
            int treatmentColumn = column(table, options.treatment);
            int mediatorColumn = column(table, options.mediator);
            int[] covariateColumns = options.covariates.stream()
                .mapToInt(name -> column(table, name)).toArray();
            int[] groupColumns = options.groups.stream()
                .mapToInt(name -> column(table, name)).toArray();
            int individualColumn = options.individualId == null ? -1
                : column(table, options.individualId);
            List<Integer> retained = new ArrayList<>();
            for (int row = 0; row < table.rows().size(); row++) {
                String[] fields = table.rows().get(row);
                boolean complete = present(fields[outcomeColumn])
                    && present(fields[treatmentColumn])
                    && present(fields[mediatorColumn]);
                for (int index : covariateColumns) {
                    complete &= present(fields[index]);
                }
                for (int index : groupColumns) {
                    complete &= present(fields[index]);
                }
                if (individualColumn >= 0) {
                    complete &= present(fields[individualColumn]);
                }
                if (complete) retained.add(row);
            }
            if (retained.isEmpty()) {
                throw new IllegalArgumentException(
                    "no complete mediation observations remain");
            }
            int rows = retained.size();
            double[] outcome = new double[rows];
            double[] treatment = new double[rows];
            double[] mediator = new double[rows];
            double[][] covariates = covariateColumns.length == 0
                ? null : new double[rows][covariateColumns.length];
            List<List<String>> groupValues = new ArrayList<>();
            for (int ignored : groupColumns) groupValues.add(new ArrayList<>());
            List<String> individualIds = new ArrayList<>();
            for (int target = 0; target < rows; target++) {
                int source = retained.get(target);
                String[] fields = table.rows().get(source);
                outcome[target] = numeric(fields[outcomeColumn],
                    options.outcome, source + 2);
                treatment[target] = numeric(fields[treatmentColumn],
                    options.treatment, source + 2);
                mediator[target] = numeric(fields[mediatorColumn],
                    options.mediator, source + 2);
                for (int index = 0; index < covariateColumns.length; index++) {
                    covariates[target][index] = numeric(
                        fields[covariateColumns[index]],
                        options.covariates.get(index), source + 2);
                }
                for (int index = 0; index < groupColumns.length; index++) {
                    groupValues.get(index).add(
                        fields[groupColumns[index]].trim());
                }
                if (individualColumn >= 0) {
                    individualIds.add(fields[individualColumn].trim());
                }
            }
            List<RandomEffectTerm> effects = new ArrayList<>();
            for (int index = 0; index < options.groups.size(); index++) {
                effects.add(RandomEffectTerm.randomIntercept(
                    options.groups.get(index), groupValues.get(index)));
            }
            return new Input(outcome, treatment, mediator, covariates,
                List.copyOf(effects), List.copyOf(individualIds), rows,
                table.rows().size());
        }
    }

    private record Fit(
            String model, List<MediationEffect> effects, boolean converged,
            int pedigreeFileMembers, int pedigreeFileObservationsMatched,
            int pedigreeAliasesResolved, int singletonFamilies,
            int singletonObservations) {
        static Fit ols(MediationResult result) {
            return new Fit("ols", effects(result.aPath(), result.bPath(),
                result.indirectEffect(), result.directEffect(),
                result.totalEffect()), true, 0, 0, 0, 0, 0);
        }
        static Fit mixed(MediationMixedResult result, String model,
                int singletonFamilies, int singletonObservations) {
            return mixed(result, model, 0, 0, 0, singletonFamilies,
                singletonObservations);
        }
        static Fit mixed(MediationMixedResult result, String model,
                int pedigreeFileMembers,
                int pedigreeFileObservationsMatched,
                int pedigreeAliasesResolved, int singletonFamilies,
                int singletonObservations) {
            return new Fit(model, effects(result.aPath(), result.bPath(),
                result.indirectEffect(), result.directEffect(),
                result.totalEffect()), result.converged(),
                pedigreeFileMembers, pedigreeFileObservationsMatched,
                pedigreeAliasesResolved, singletonFamilies,
                singletonObservations);
        }
        private static List<MediationEffect> effects(
                MediationEffect... values) {
            return List.of(values);
        }
    }

    private static int column(DelimitedData table, String name) {
        for (int index = 0; index < table.header().size(); index++) {
            if (table.header().get(index).equalsIgnoreCase(name)) return index;
        }
        throw new IllegalArgumentException("column is absent: " + name);
    }

    private static boolean present(String value) {
        if (value == null) return false;
        String normalized = value.trim().toLowerCase(Locale.ROOT);
        return !normalized.isEmpty() && !normalized.equals("na")
            && !normalized.equals("n/a") && !normalized.equals("nan")
            && !normalized.equals("null") && !normalized.equals(".");
    }

    private static double numeric(String value, String column, int row) {
        try {
            double result = Double.parseDouble(value.trim());
            if (!Double.isFinite(result)) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException("non-numeric value in column "
                + column + " at row " + row + ": " + value, exception);
        }
    }

    private static final class Options {
        Path input, output, log, pedigree;
        String outcome, treatment, mediator, individualId;
        String pedigreeId, sireId, damId, pedigreeFamilyId;
        final List<String> covariates = new ArrayList<>();
        final List<String> groups = new ArrayList<>();
        double confidence = 0.95;
        BackendPolicy backend = BackendPolicy.PREFERRED;
        boolean overwrite, help;

        static Options parse(String[] arguments) {
            Options result = new Options();
            for (int index = 0; index < arguments.length; index++) {
                String option = arguments[index];
                switch (option) {
                    case "--input", "--in" -> result.input = Path.of(
                        value(arguments, ++index, option));
                    case "--out", "--output" -> result.output = Path.of(
                        value(arguments, ++index, option));
                    case "--log" -> result.log = Path.of(
                        value(arguments, ++index, option));
                    case "--outcome" -> result.outcome =
                        value(arguments, ++index, option);
                    case "--treatment", "--exposure" -> result.treatment =
                        value(arguments, ++index, option);
                    case "--mediator" -> result.mediator =
                        value(arguments, ++index, option);
                    case "--covariates" -> addColumns(result.covariates,
                        value(arguments, ++index, option));
                    case "--group" -> addColumns(result.groups,
                        value(arguments, ++index, option));
                    case "--individual-id" -> result.individualId =
                        value(arguments, ++index, option);
                    case "--pedigree" -> result.pedigree = Path.of(
                        value(arguments, ++index, option));
                    case "--pedigree-id" -> result.pedigreeId =
                        value(arguments, ++index, option);
                    case "--sire-id" -> result.sireId =
                        value(arguments, ++index, option);
                    case "--dam-id" -> result.damId =
                        value(arguments, ++index, option);
                    case "--pedigree-family-id" -> result.pedigreeFamilyId =
                        value(arguments, ++index, option);
                    case "--confidence" -> result.confidence =
                        probability(value(arguments, ++index, option), option);
                    case "--backend" -> result.backend =
                        backend(value(arguments, ++index, option));
                    case "--overwrite" -> result.overwrite = true;
                    case "--help", "-h" -> result.help = true;
                    default -> throw new IllegalArgumentException(
                        "unknown mediation option: " + option);
                }
            }
            if (!result.help && (result.input == null || result.output == null
                    || result.outcome == null || result.treatment == null
                    || result.mediator == null)) {
                throw new IllegalArgumentException(
                    "--input, --outcome, --treatment, --mediator, and --out "
                        + "are required");
            }
            if (result.pedigree != null && (result.individualId == null
                    || result.pedigreeId == null || result.sireId == null
                    || result.damId == null)) {
                throw new IllegalArgumentException(
                    "--pedigree requires --individual-id, --pedigree-id, "
                        + "--sire-id, and --dam-id");
            }
            if (result.pedigree == null && (result.individualId != null
                    || result.pedigreeId != null || result.sireId != null
                    || result.damId != null
                    || result.pedigreeFamilyId != null)) {
                throw new IllegalArgumentException(
                    "pedigree ID options require --pedigree");
            }
            Set<String> variables = new HashSet<>();
            variables.add(result.outcome);
            variables.add(result.treatment);
            variables.add(result.mediator);
            for (String covariate : result.covariates) {
                if (!variables.add(covariate)) {
                    throw new IllegalArgumentException(
                        "mediation variable columns must be distinct: "
                            + covariate);
                }
            }
            return result;
        }
    }

    private static void prepareOutputs(Options options, Path log)
            throws IOException {
        Path output = options.output.toAbsolutePath().normalize();
        Path logPath = log.toAbsolutePath().normalize();
        if (output.equals(logPath)) {
            throw new IllegalArgumentException(
                "output and log paths must be distinct");
        }
        for (Path input : List.of(options.input)) {
            if (output.equals(input.toAbsolutePath().normalize())
                    || logPath.equals(input.toAbsolutePath().normalize())) {
                throw new IllegalArgumentException(
                    "output or log path aliases the input");
            }
        }
        if (options.pedigree != null
                && (output.equals(options.pedigree.toAbsolutePath().normalize())
                    || logPath.equals(
                        options.pedigree.toAbsolutePath().normalize()))) {
            throw new IllegalArgumentException(
                "output or log path aliases the pedigree");
        }
        if (!options.overwrite) {
            PipelinePaths.requireFreshOutputs(options.output, log);
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

    private static char delimiter(Path path) {
        return path.toString().toLowerCase(Locale.ROOT).endsWith(".csv")
            ? ',' : '\t';
    }

    private static void row(StringBuilder output, char delimiter,
            String... values) {
        for (int index = 0; index < values.length; index++) {
            if (index > 0) output.append(delimiter);
            String value = values[index];
            if (value.indexOf(delimiter) >= 0 || value.indexOf('"') >= 0
                    || value.indexOf('\n') >= 0 || value.indexOf('\r') >= 0) {
                output.append("\"").append(value.replace("\"", "\"\""))
                    .append("\"");
            } else {
                output.append(value);
            }
        }
        output.append('\n');
    }

    private static void addColumns(List<String> target, String value) {
        for (String column : value.split(",", -1)) {
            String name = column.trim();
            if (!name.isEmpty() && !target.contains(name)) target.add(name);
        }
        if (target.isEmpty()) {
            throw new IllegalArgumentException(
                "column list must not be empty");
        }
    }

    private static double probability(String value, String option) {
        try {
            double result = Double.parseDouble(value);
            if (!(result > 0.0 && result < 1.0)
                    || !Double.isFinite(result)) {
                throw new NumberFormatException();
            }
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                option + " must be in (0,1)", exception);
        }
    }

    private static BackendPolicy backend(String value) {
        try {
            return BackendPolicy.valueOf(value.toUpperCase(Locale.ROOT));
        } catch (IllegalArgumentException exception) {
            throw new IllegalArgumentException(
                "unknown backend policy: " + value, exception);
        }
    }

    private static String value(String[] values, int index, String option) {
        if (index >= values.length || values[index].isBlank()) {
            throw new IllegalArgumentException(option + " requires a value");
        }
        return values[index];
    }

    private static String help() {
        return """
            Usage: mediation --input FILE --outcome Y --treatment X
              --mediator M [--covariates C1,C2] [--group COLUMN]
              [--confidence 0.95] --out FILE

            Add one or more --group columns for shared grouped-REML random
            intercepts. For pedigree REML add:
              --individual-id COLUMN --pedigree FILE --pedigree-id COLUMN
              --sire-id COLUMN --dam-id COLUMN
              [--pedigree-family-id COLUMN]

            Input is observation-by-variable CSV/TSV. Outcome, treatment,
            mediator, and covariates must be numeric. One common complete-case
            sample is used. The output reports a, b, indirect (Sobel),
            direct, and total effects. This is Gaussian linear mediation;
            it does not establish causal identification or bootstrap the
            indirect effect.
            """;
    }
}
