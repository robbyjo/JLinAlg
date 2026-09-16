/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import jdistlib.ChiSquare;
import jdistlib.Normal;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mr.MultivariableInstrument;
import org.jlinalg.mr.MultivariableMendelianRandomization;
import org.jlinalg.mr.MultivariableMrResult;
import org.jlinalg.mr.MultivariateInstrument;

/** Multiple-exposure, single-outcome multivariable MR command. */
final class MultivariableMrCli {
    private MultivariableMrCli() { }

    static int run(String[] arguments, PrintStream output, PrintStream error) {
        try {
            Options options = Options.parse(arguments);
            if (options.help) { output.println(help()); return 0; }
            PipelinePaths.requireFreshOutputs(options.output);
            List<MultivariateInstrument> wide = MrWideTable.read(options.input,
                options.exposures, List.of(options.outcome));
            List<MultivariableInstrument> instruments = new ArrayList<>();
            for (MultivariateInstrument value : wide)
                instruments.add(new MultivariableInstrument(value.variantId(),
                    value.exposureEffects(), value.exposureStandardErrors(),
                    value.outcomeEffects()[0], value.outcomeStandardErrors()[0]));
            boolean egger = options.method.equals("egger");
            MultivariableMrResult result =
                MultivariableMendelianRandomization.fit(instruments,
                    options.exposures, egger, options.backend);
            double critical = Normal.quantile(.5 + options.confidence / 2,
                0, 1, true, false);
            double[] beta = result.beta(), se = result.standardErrors();
            double[] p = result.pValues(), strength = result.marginalFStatistics();
            double qP = ChiSquare.cumulative(result.cochranQ(),
                result.heterogeneityDegreesOfFreedom(), false, false);
            StringBuilder text = new StringBuilder(
                "method\texposure\toutcome\ttest\ttest_scope\ttest_df\ttest_p_value\t"
                + "estimate\tstandard_error\tcausal_p_value\tci_lower\tci_upper\t"
                + "marginal_f\tcochran_q\theterogeneity_df\theterogeneity_p_value\t"
                + "intercept\tintercept_se\tinstruments\n");
            for (int index = 0; index < beta.length; index++)
                text.append(egger ? "MVMR_EGGER" : "MVMR_IVW").append('\t')
                    .append(options.exposures.get(index)).append('\t')
                    .append(options.outcome).append("\tcoefficient\texposure\t1\t")
                    .append(p[index]).append('\t').append(beta[index]).append('\t')
                    .append(se[index]).append('\t').append(p[index]).append('\t')
                    .append(beta[index] - critical * se[index]).append('\t')
                    .append(beta[index] + critical * se[index]).append('\t')
                    .append(strength[index]).append('\t')
                    .append(result.cochranQ()).append('\t')
                    .append(result.heterogeneityDegreesOfFreedom()).append('\t')
                    .append(qP).append('\t').append(result.intercept()).append('\t')
                    .append(result.interceptStandardError()).append('\t')
                    .append(instruments.size()).append('\n');
            write(options.output, text.toString(), output);
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static void write(Path path, String text, PrintStream output)
            throws IOException {
        if (path == null) output.print(text);
        else Files.writeString(path, text,
            java.nio.file.StandardOpenOption.CREATE_NEW);
    }

    private static String help() {
        return "Usage: mr-mvmr --input FILE --exposures NAME1,NAME2 "
            + "[--outcome NAME] [--method ivw|egger] [--confidence 0.95] "
            + "[--backend preferred|cpu|...] [--output FILE]\n"
            + "This is multivariable MR: multiple exposures and exactly one outcome. "
            + "Wide input uses variant_id plus beta_exposure_NAME/se_exposure_NAME "
            + "and beta_outcome/se_outcome (or outcome-suffixed columns). "
            + "Use mr-multivariate for multiple correlated outcomes.";
    }

    private static final class Options {
        Path input, output; List<String> exposures; String outcome = "outcome";
        String method = "ivw"; double confidence = .95; boolean help;
        BackendPolicy backend = BackendPolicy.PREFERRED;
        static Options parse(String[] args) {
            Options result = new Options();
            String exposureNames = null;
            for (int index = 0; index < args.length; index++) switch (args[index]) {
                case "--input", "--in" -> result.input = Path.of(value(args, ++index, "--input"));
                case "--output", "--out" -> result.output = Path.of(value(args, ++index, "--output"));
                case "--exposures" -> exposureNames = value(args, ++index, "--exposures");
                case "--outcome" -> result.outcome = value(args, ++index, "--outcome");
                case "--method" -> result.method = value(args, ++index, "--method").toLowerCase(Locale.ROOT);
                case "--confidence" -> result.confidence = Double.parseDouble(value(args, ++index, "--confidence"));
                case "--backend" -> result.backend = BackendPolicy.valueOf(value(args, ++index, "--backend").toUpperCase(Locale.ROOT));
                case "--help", "-h" -> result.help = true;
                default -> throw new IllegalArgumentException("unknown mr-mvmr option: " + args[index]);
            }
            if (result.help) return result;
            if (result.input == null) throw new IllegalArgumentException("--input is required");
            result.exposures = MrWideTable.names(exposureNames, "--exposures");
            if (result.exposures.size() < 2)
                throw new IllegalArgumentException("mr-mvmr requires at least two exposure names");
            if (!result.method.equals("ivw") && !result.method.equals("egger"))
                throw new IllegalArgumentException("--method must be ivw or egger");
            if (!(result.confidence > 0 && result.confidence < 1))
                throw new IllegalArgumentException("--confidence must be in (0,1)");
            return result;
        }
        private static String value(String[] args, int index, String option) {
            if (index >= args.length || args[index].isBlank())
                throw new IllegalArgumentException(option + " requires a value");
            return args[index];
        }
    }
}
