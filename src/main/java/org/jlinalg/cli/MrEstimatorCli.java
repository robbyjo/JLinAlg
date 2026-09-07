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

/** General summary-data MR estimator command. */
final class MrEstimatorCli {
    private MrEstimatorCli() { }

    static int run(String[] arguments, java.io.PrintStream output,
            java.io.PrintStream error) {
        try {
            Options options = Options.parse(arguments);
            if (options.help) { output.println(help()); return 0; }
            List<HarmonizedInstrument> instruments = read(options.input);
            List<MrEstimate> estimates = new ArrayList<>();
            switch (options.method) {
                case "ivw-fixed" -> estimates.add(MendelianRandomization.ivw(instruments, false, options.confidence));
                case "ivw-random" -> estimates.add(MendelianRandomization.ivw(instruments, true, options.confidence));
                case "egger" -> estimates.add(MendelianRandomization.egger(instruments, options.confidence).slope());
                case "weighted-median" -> estimates.add(MendelianRandomization.weightedMedian(instruments,
                    new MrOptions(options.confidence, options.bootstrap, options.seed)));
                case "all" -> { estimates.add(MendelianRandomization.ivw(instruments, false, options.confidence)); estimates.add(MendelianRandomization.ivw(instruments, true, options.confidence)); estimates.add(MendelianRandomization.egger(instruments, options.confidence).slope()); estimates.add(MendelianRandomization.weightedMedian(instruments, new MrOptions(options.confidence, options.bootstrap, options.seed))); }
                default -> throw new IllegalArgumentException("unknown MR estimator: " + options.method);
            }
            StringBuilder result = new StringBuilder("method\testimate\tstandard_error\tp_value\tci_lower\tci_upper\tinstruments\n");
            for (MrEstimate estimate : estimates) result.append(estimate.method()).append('\t').append(estimate.estimate()).append('\t').append(estimate.standardError()).append('\t').append(estimate.pValue()).append('\t').append(estimate.confidenceLower()).append('\t').append(estimate.confidenceUpper()).append('\t').append(estimate.instrumentCount()).append('\n');
            if (options.output == null) output.print(result); else Files.writeString(options.output, result.toString());
            if (options.plot != null) MrPlot.writeSvg(options.plot, instruments, estimates.get(0));
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage()); return 2;
        }
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
    private static String help() { return "Usage: mr-estimate --input FILE [--method ivw-fixed|ivw-random|egger|weighted-median|all] [--output FILE] [--plot FILE]\nInput columns: variant_id, beta_exposure, se_exposure, beta_outcome, se_outcome (CSV or TSV)."; }

    private static final class Options {
        Path input, output, plot; String method = "all"; double confidence = 0.95; int bootstrap = 1000; long seed = 20260831L; boolean help;
        static Options parse(String[] args) { Options o = new Options(); for (int i = 0; i < args.length; i++) switch (args[i]) { case "--input", "--in" -> o.input = Path.of(value(args, ++i, "--input")); case "--output", "--out" -> o.output = Path.of(value(args, ++i, "--output")); case "--plot" -> o.plot = Path.of(value(args, ++i, "--plot")); case "--method" -> o.method = value(args, ++i, "--method").toLowerCase(Locale.ROOT); case "--confidence" -> o.confidence = Double.parseDouble(value(args, ++i, "--confidence")); case "--bootstrap" -> o.bootstrap = Integer.parseInt(value(args, ++i, "--bootstrap")); case "--seed" -> o.seed = Long.parseLong(value(args, ++i, "--seed")); case "--help", "-h" -> o.help = true; default -> throw new IllegalArgumentException("unknown mr-estimate option: " + args[i]); } if (!o.help && o.input == null) throw new IllegalArgumentException("--input is required"); if (!(o.confidence > 0.0 && o.confidence < 1.0) || !Double.isFinite(o.confidence)) throw new IllegalArgumentException("--confidence must be in (0,1)"); if (o.bootstrap < 2) throw new IllegalArgumentException("--bootstrap must be at least 2"); return o; }
        static String value(String[] a, int i, String option) { if (i >= a.length || a[i].isBlank()) throw new IllegalArgumentException(option + " requires a value"); return a[i]; }
    }
}
