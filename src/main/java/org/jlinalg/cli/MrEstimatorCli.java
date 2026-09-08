/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.jlinalg.mr.HarmonizedInstrument;
import org.jlinalg.mr.MendelianRandomization;
import org.jlinalg.mr.MrEstimate;
import org.jlinalg.mr.MrOptions;
import org.jlinalg.mr.MrPlot;
import org.jlinalg.mr.CorrelatedMendelianRandomization;
import org.jlinalg.mr.OverlapAwareMendelianRandomization;
import org.jlinalg.mr.MrEggerResult;
import org.jlinalg.genetics.ConditionalAssociation;
import org.jlinalg.genetics.ConditionalAssociationModel;
import org.jlinalg.genetics.ConditionalAssociationResult;
import org.jlinalg.genetics.ConditionalSignalSelection;
import org.jlinalg.genetics.SecondarySignalClumper;

/** General summary-data MR estimator command. */
final class MrEstimatorCli {
    private MrEstimatorCli() { }

    static int run(String[] arguments, java.io.PrintStream output,
            java.io.PrintStream error) {
        try {
            Options options = Options.parse(arguments);
            if (options.help) { output.println(help()); return 0; }
            if (options.method.equals("conditional") || options.method.equals("secondary-signals"))
                return conditional(options, output);
            if (options.joint || options.conditionOn != null)
                throw new IllegalArgumentException("--joint and --condition-on require --method conditional; causal MR inputs are not automatically conditioned");
            List<HarmonizedInstrument> instruments = read(options.input);
            List<MrEstimate> estimates = new ArrayList<>();
            MrEggerResult plotEgger = null;
            switch (options.method) {
                case "ivw-fixed" -> estimates.add(MendelianRandomization.ivw(instruments, false, options.confidence));
                case "ivw-random" -> estimates.add(MendelianRandomization.ivw(instruments, true, options.confidence));
                case "egger" -> { plotEgger = MendelianRandomization.egger(instruments, options.confidence); estimates.add(plotEgger.slope()); }
                case "weighted-median" -> estimates.add(MendelianRandomization.weightedMedian(instruments,
                    new MrOptions(options.confidence, options.bootstrap, options.seed)));
                case "ivw-generalized-fixed" -> estimates.add(CorrelatedMendelianRandomization.ivw(instruments, requireLd(options), false, options.confidence, options.backend).estimate());
                case "ivw-generalized-random" -> estimates.add(CorrelatedMendelianRandomization.ivw(instruments, requireLd(options), true, options.confidence, options.backend).estimate());
                case "egger-generalized" -> { plotEgger = CorrelatedMendelianRandomization.egger(instruments, requireLd(options), options.confidence, options.backend).estimate(); estimates.add(plotEgger.slope()); }
                case "overlap-aware" -> estimates.add(overlap(instruments, options));
                case "all" -> { estimates.add(MendelianRandomization.ivw(instruments, false, options.confidence)); estimates.add(MendelianRandomization.ivw(instruments, true, options.confidence)); estimates.add(MendelianRandomization.egger(instruments, options.confidence).slope()); estimates.add(MendelianRandomization.weightedMedian(instruments, new MrOptions(options.confidence, options.bootstrap, options.seed))); if (options.ld != null) estimates.add(CorrelatedMendelianRandomization.ivw(instruments, options.ld, false, options.confidence, options.backend).estimate()); if (options.samplingCovariance != null) estimates.add(overlap(instruments, options)); }
                default -> throw new IllegalArgumentException("unknown MR estimator: " + options.method);
            }
            StringBuilder result = new StringBuilder("method\testimate\tstandard_error\tcausal_p_value\tci_lower\tci_upper\tinstruments\n");
            for (MrEstimate estimate : estimates) result.append(estimate.method()).append('\t').append(estimate.estimate()).append('\t').append(estimate.standardError()).append('\t').append(estimate.pValue()).append('\t').append(estimate.confidenceLower()).append('\t').append(estimate.confidenceUpper()).append('\t').append(estimate.instrumentCount()).append('\n');
            if (options.output == null) output.print(result); else Files.writeString(options.output, result.toString());
            if (options.plot != null) {
                if (plotEgger != null) MrPlot.writeSvg(options.plot, instruments, plotEgger);
                else MrPlot.writeSvg(options.plot, instruments, estimates.get(0));
            }
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage()); return 2;
        }
    }

    private static MrEstimate overlap(List<HarmonizedInstrument> instruments, Options options) throws IOException {
        var fit = OverlapAwareMendelianRandomization.ivw(instruments, requireSamplingCovariance(options), options.confidence);
        if (!fit.converged()) throw new IllegalArgumentException("overlap-aware IVW did not converge in " + fit.iterations() + " iterations");
        return fit.estimate();
    }

    private static int conditional(Options options, java.io.PrintStream output) throws IOException {
        if (options.plot != null) throw new IllegalArgumentException("--plot is a causal MR scatter; conditional SNP tests have no causal slope");
        if (options.samplingCovariance != null) throw new IllegalArgumentException("conditional single-trait inference does not use exposure/outcome sampling covariance");
        List<String> lines = Files.readAllLines(options.input);
        if (lines.size() < 2) throw new IllegalArgumentException("conditional input needs header and associations");
        String separator = options.input.toString().toLowerCase(Locale.ROOT).endsWith(".csv") ? "," : "\t";
        String[] header = lines.get(0).split(separator, -1);
        int id = column(header, "variant_id", "SNP", "rsid");
        int beta = column(header, "beta_" + options.association, "beta." + options.association, "beta");
        int se = column(header, "se_" + options.association, "se." + options.association, "se");
        List<ConditionalAssociation> values = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            if (lines.get(i).isBlank()) continue;
            String[] row = lines.get(i).split(separator, -1);
            if (row.length != header.length) throw new IllegalArgumentException("inconsistent conditional input row " + (i + 1));
            values.add(new ConditionalAssociation(row[id], number(row[beta]), number(row[se]), "locus"));
        }
        List<ConditionalAssociationResult> tests;
        List<String> selected = List.of();
        if (options.method.equals("secondary-signals")) {
            if (options.joint || options.conditionOn != null) throw new IllegalArgumentException("secondary-signals selects its own conditioning set");
            ConditionalSignalSelection selection = SecondarySignalClumper.select(values, requireLd(options),
                options.pThreshold, options.maxSignals, options.confidence);
            if (!selection.converged()) throw new IllegalArgumentException("conditional signal selection cycled or exhausted its iteration limit");
            tests = selection.associations(); selected = selection.selectedVariantIds();
        } else {
            if (options.joint == (options.conditionOn != null))
                throw new IllegalArgumentException("conditional requires exactly one of --joint or --condition-on SNP1,SNP2");
            int[] indices;
            if (options.joint) indices = java.util.stream.IntStream.range(0, values.size()).toArray();
            else {
                String[] ids = options.conditionOn.split(",", -1); indices = new int[ids.length];
                for (int i = 0; i < ids.length; i++) {
                    indices[i] = -1;
                    for (int j = 0; j < values.size(); j++) if (values.get(j).variantId().equals(ids[i].trim())) indices[i] = j;
                    if (indices[i] < 0) throw new IllegalArgumentException("unknown conditioning SNP: " + ids[i]);
                }
            }
            tests = new ConditionalAssociationModel(values, requireLd(options)).condition(indices, options.confidence);
        }
        StringBuilder result = new StringBuilder("variant_id\tassociation\tconditional_beta\tconditional_se\tconditional_instrument_p_value\tci_lower\tci_upper\tconditioned_on\tselected\n");
        for (ConditionalAssociationResult test : tests) result.append(test.variantId()).append('\t').append(options.association)
            .append('\t').append(test.effect()).append('\t').append(test.standardError()).append('\t').append(test.pValue())
            .append('\t').append(test.confidenceLower()).append('\t').append(test.confidenceUpper())
            .append('\t').append(String.join(",", test.conditionedOn())).append('\t').append(selected.contains(test.variantId())).append('\n');
        if (options.output == null) output.print(result); else Files.writeString(options.output, result.toString());
        return 0;
    }

    private static List<HarmonizedInstrument> read(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path); if (lines.size() < 2) throw new IOException("MR input needs a header and at least one row");
        char delimiter = path.toString().toLowerCase(Locale.ROOT).endsWith(".csv") ? ',' : '\t';
        String[] header = lines.get(0).split(java.util.regex.Pattern.quote(String.valueOf(delimiter)), -1);
        int id = column(header, "variant_id", "SNP", "rsid");
        int bx = column(header, "beta_exposure", "beta.exposure", "exposure_effect");
        int sx = column(header, "se_exposure", "se.exposure", "exposure_se");
        int by = column(header, "beta_outcome", "beta.outcome", "outcome_effect");
        int sy = column(header, "se_outcome", "se.outcome", "outcome_se");
        List<HarmonizedInstrument> result = new ArrayList<>();
        for (int line = 1; line < lines.size(); line++) { if (lines.get(line).isBlank()) continue; String[] v = lines.get(line).split(java.util.regex.Pattern.quote(String.valueOf(delimiter)), -1); if (v.length != header.length) throw new IOException("inconsistent MR input row " + (line + 1)); result.add(new HarmonizedInstrument(v[id], "", "", number(v[bx]), number(v[sx]), number(v[by]), number(v[sy]), Double.NaN, Double.NaN, false, false)); }
        return List.copyOf(result);
    }

    private static int column(String[] header, String... names) { for (int i = 0; i < header.length; i++) for (String name : names) if (header[i].trim().equalsIgnoreCase(name)) return i; throw new IllegalArgumentException("MR input is missing one of " + Arrays.toString(names)); }
    private static double number(String value) { try { return Double.parseDouble(value.trim()); } catch (NumberFormatException e) { throw new IllegalArgumentException("invalid numeric MR value: " + value, e); } }
    private static double[][] requireLd(Options options) throws IOException { if (options.ld == null) throw new IllegalArgumentException("--ld is required for generalized MR or conditional association analysis"); return options.ld; }
    private static double[] requireSamplingCovariance(Options options) throws IOException { if (options.samplingCovariance == null) throw new IllegalArgumentException("--sampling-covariance is required for overlap-aware MR"); return options.samplingCovariance; }
    private static double[][] readMatrix(Path path) throws IOException { List<String> lines = Files.readAllLines(path); if (lines.isEmpty()) throw new IOException("LD matrix is empty"); char delimiter = path.toString().toLowerCase(Locale.ROOT).endsWith(".csv") ? ',' : '\t'; int size = lines.size(); double[][] result = new double[size][]; for (int row = 0; row < size; row++) { String[] values = lines.get(row).split(java.util.regex.Pattern.quote(String.valueOf(delimiter)), -1); if (values.length != size) throw new IOException("LD matrix must be square"); result[row] = new double[size]; for (int column = 0; column < size; column++) result[row][column] = number(values[column]); } return result; }
    private static double[] readVector(Path path) throws IOException { List<String> lines = Files.readAllLines(path); List<Double> values = new ArrayList<>(); for (String line : lines) for (String value : line.split("[,\\t]")) if (!value.isBlank()) values.add(number(value)); double[] result = new double[values.size()]; for (int index = 0; index < result.length; index++) result[index] = values.get(index); return result; }
    private static String help() { return "Usage: mr-estimate --input FILE [--method ivw-fixed|ivw-random|egger|weighted-median|ivw-generalized-fixed|ivw-generalized-random|egger-generalized|overlap-aware|conditional|secondary-signals|all] [--ld MATRIX] [--sampling-covariance FILE] [--confidence 0.95] [--backend preferred|cpu|...] [--output FILE] [--plot FILE]\nInput columns: variant_id, beta_exposure, se_exposure, beta_outcome, se_outcome (CSV or TSV). LD is signed, allele-aligned and in exact input order. Conditional SNP association tests require --joint or --condition-on SNP1,SNP2, with --association exposure|outcome (default exposure; beta/se columns also accepted). secondary-signals uses --p-threshold 5e-8 and --max-signals 10. Conditional p-values test SNP association, not causal MR. --sampling-covariance supplies one exposure/outcome sampling covariance per instrument."; }

    private static final class Options {
        Path input, output, plot; double[][] ld; double[] samplingCovariance; String method = "all"; double confidence = 0.95; int bootstrap = 1000; long seed = 20260831L; boolean help;
        org.jlinalg.compute.BackendPolicy backend = org.jlinalg.compute.BackendPolicy.PREFERRED;
        String conditionOn, association = "exposure"; boolean joint; double pThreshold = 5e-8; int maxSignals = 10;
        static Options parse(String[] args) throws IOException { Options o = new Options(); for (int i = 0; i < args.length; i++) switch (args[i]) { case "--backend" -> o.backend = org.jlinalg.compute.BackendPolicy.valueOf(value(args, ++i, "--backend").toUpperCase(Locale.ROOT)); case "--joint" -> o.joint = true; case "--condition-on" -> o.conditionOn = value(args, ++i, "--condition-on"); case "--association" -> o.association = value(args, ++i, "--association").toLowerCase(Locale.ROOT); case "--p-threshold" -> o.pThreshold = Double.parseDouble(value(args, ++i, "--p-threshold")); case "--max-signals" -> o.maxSignals = Integer.parseInt(value(args, ++i, "--max-signals")); case "--input", "--in" -> o.input = Path.of(value(args, ++i, "--input")); case "--output", "--out" -> o.output = Path.of(value(args, ++i, "--output")); case "--plot" -> o.plot = Path.of(value(args, ++i, "--plot")); case "--ld", "--ld-correlation" -> o.ld = readMatrix(Path.of(value(args, ++i, "--ld"))); case "--sampling-covariance", "--overlap-covariance" -> o.samplingCovariance = readVector(Path.of(value(args, ++i, "--sampling-covariance"))); case "--method" -> o.method = value(args, ++i, "--method").toLowerCase(Locale.ROOT); case "--confidence" -> o.confidence = Double.parseDouble(value(args, ++i, "--confidence")); case "--bootstrap" -> o.bootstrap = Integer.parseInt(value(args, ++i, "--bootstrap")); case "--seed" -> o.seed = Long.parseLong(value(args, ++i, "--seed")); case "--help", "-h" -> o.help = true; default -> throw new IllegalArgumentException("unknown mr-estimate option: " + args[i]); } if (!o.help && o.input == null) throw new IllegalArgumentException("--input is required"); if (!(o.confidence > 0.0 && o.confidence < 1.0) || !Double.isFinite(o.confidence)) throw new IllegalArgumentException("--confidence must be in (0,1)"); if (o.bootstrap < 2) throw new IllegalArgumentException("--bootstrap must be at least 2"); if (!o.association.equals("exposure") && !o.association.equals("outcome")) throw new IllegalArgumentException("--association must be exposure or outcome"); if (!(o.pThreshold > 0.0 && o.pThreshold < 1.0) || o.maxSignals < 1) throw new IllegalArgumentException("invalid selection threshold or maximum signals"); return o; }
        static String value(String[] a, int i, String option) { if (i >= a.length || a[i].isBlank()) throw new IllegalArgumentException(option + " requires a value"); return a[i]; }
    }
}
