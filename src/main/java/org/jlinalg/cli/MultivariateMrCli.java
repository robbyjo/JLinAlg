/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Locale;
import jdistlib.Normal;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mr.MultivariateInstrument;
import org.jlinalg.mr.MultivariateMendelianRandomization;
import org.jlinalg.mr.MultivariateMrJointTest;
import org.jlinalg.mr.MultivariateMrOutlier;
import org.jlinalg.mr.MultivariateMrPresso;
import org.jlinalg.mr.MultivariateMrPressoResult;
import org.jlinalg.mr.MultivariateMrPlot;
import org.jlinalg.mr.MultivariateMrResult;

/** Joint multiple-correlated-outcome summary-data MR command. */
final class MultivariateMrCli {
    private MultivariateMrCli() { }

    static int run(String[] arguments, PrintStream output, PrintStream error) {
        try {
            Options options = Options.parse(arguments);
            if (options.help) { output.println(help()); return 0; }
            PipelinePaths.requireFreshOutputs(options.output, options.diagnostics,
                options.plot);
            List<MultivariateInstrument> instruments = MrWideTable.read(
                options.input, options.exposures, options.outcomes);
            double[][] correlation = MrWideTable.readMatrix(
                options.outcomeCorrelation);
            String result;
            MultivariateMrResult plotFit;
            if (options.method.equals("presso")) {
                MultivariateMrPressoResult fit = MultivariateMrPresso.analyze(
                    instruments, options.exposures, options.outcomes,
                    correlation, options.outlierThreshold, options.simulations,
                    options.seed, options.backend);
                result = effects("MULTIVARIATE_PRESSO_RAW", fit.rawEstimate(),
                    options.confidence, instruments.size());
                if (fit.outlierCorrectedEstimate() != null)
                    result += effects("MULTIVARIATE_PRESSO_CORRECTED",
                        fit.outlierCorrectedEstimate(), options.confidence,
                        instruments.size() - fit.outlierVariants().size())
                        .substring(header().length());
                plotFit = fit.outlierCorrectedEstimate() == null
                    ? fit.rawEstimate() : fit.outlierCorrectedEstimate();
                Files.writeString(options.diagnostics, diagnostics(fit),
                    java.nio.file.StandardOpenOption.CREATE_NEW);
            } else {
                boolean random = options.method.equals("ivw-random");
                MultivariateMrResult fit =
                    MultivariateMendelianRandomization.fit(instruments,
                        options.exposures, options.outcomes, correlation,
                        random, options.backend);
                result = effects(random ? "MULTIVARIATE_IVW_RANDOM"
                    : "MULTIVARIATE_IVW_FIXED", fit, options.confidence,
                    instruments.size());
                plotFit = fit;
            }
            if (options.output == null) output.print(result);
            else Files.writeString(options.output, result,
                java.nio.file.StandardOpenOption.CREATE_NEW);
            if (options.plot != null)
                MultivariateMrPlot.writeForestSvg(options.plot, plotFit,
                    options.confidence);
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static String effects(String method, MultivariateMrResult result,
            double confidence, int instruments) {
        double critical = Normal.quantile(.5 + confidence / 2,
            0, 1, true, false);
        double[] beta = result.beta(), se = result.standardErrors();
        double[] p = result.pValues();
        StringBuilder text = new StringBuilder(header());
        MultivariateMrJointTest overall = result.overallTest();
        for (int outcome = 0; outcome < result.outcomeNames().size(); outcome++) {
            MultivariateMrJointTest outcomeTest =
                result.outcomeJointTests().get(outcome);
            for (int exposure = 0;
                    exposure < result.exposureNames().size(); exposure++) {
                int index = result.coefficientIndex(exposure, outcome);
                MultivariateMrJointTest exposureTest =
                    result.exposureJointTests().get(exposure);
                text.append(method).append('\t')
                    .append(result.exposureNames().get(exposure)).append('\t')
                    .append(result.outcomeNames().get(outcome))
                    .append("\tcoefficient\texposure_outcome\t1\t")
                    .append(p[index]).append('\t')
                    .append(beta[index]).append('\t').append(se[index]).append('\t')
                    .append(p[index]).append('\t')
                    .append(beta[index] - critical * se[index]).append('\t')
                    .append(beta[index] + critical * se[index]).append('\t')
                    .append(exposureTest.chiSquare()).append('\t')
                    .append(exposureTest.degreesOfFreedom()).append('\t')
                    .append(exposureTest.pValue()).append('\t')
                    .append(outcomeTest.chiSquare()).append('\t')
                    .append(outcomeTest.degreesOfFreedom()).append('\t')
                    .append(outcomeTest.pValue()).append('\t')
                    .append(overall.chiSquare()).append('\t')
                    .append(overall.degreesOfFreedom()).append('\t')
                    .append(overall.pValue()).append('\t')
                    .append(result.cochranQ()).append('\t')
                    .append(result.heterogeneityDegreesOfFreedom()).append('\t')
                    .append(result.heterogeneityPValue()).append('\t')
                    .append(result.dispersion()).append('\t')
                    .append(result.marginalFStatistics()[exposure]).append('\t')
                    .append(instruments).append('\n');
            }
        }
        return text.toString();
    }

    private static String header() {
        return "method\texposure\toutcome\ttest\ttest_scope\ttest_df\ttest_p_value\t"
            + "estimate\tstandard_error\tcausal_p_value\tci_lower\tci_upper\t"
            + "exposure_joint_chi_square\texposure_joint_df\texposure_joint_p_value\t"
            + "outcome_joint_chi_square\toutcome_joint_df\toutcome_joint_p_value\t"
            + "overall_chi_square\toverall_df\toverall_p_value\tcochran_q\t"
            + "heterogeneity_df\theterogeneity_p_value\tdispersion\tmarginal_f\t"
            + "instruments\n";
    }

    private static String diagnostics(MultivariateMrPressoResult result) {
        StringBuilder text = new StringBuilder(
            "variant_id\tmahalanobis_distance\tempirical_p_value\t"
            + "bonferroni_p_value\toutlier\tglobal_statistic\tglobal_p_value\t"
            + "simulations\tseed\n");
        for (MultivariateMrOutlier test : result.instrumentTests())
            text.append(test.variantId()).append('\t')
                .append(test.mahalanobisDistance()).append('\t')
                .append(test.empiricalPValue()).append('\t')
                .append(test.bonferroniPValue()).append('\t')
                .append(test.outlier()).append('\t')
                .append(result.globalStatistic()).append('\t')
                .append(result.globalPValue()).append('\t')
                .append(result.simulations()).append('\t')
                .append(result.seed()).append('\n');
        return text.toString();
    }

    private static String help() {
        return "Usage: mr-multivariate --input FILE --exposures NAME[,NAME...] "
            + "--outcomes NAME1,NAME2 --outcome-correlation MATRIX "
            + "[--method ivw-fixed|ivw-random|presso] [--output FILE] "
            + "[--diagnostics FILE] [--plot FOREST.svg]\n"
            + "This is multivariate (multi-response) MR: at least two correlated "
            + "outcomes, with one or more exposures. Wide input uses variant_id, "
            + "beta_exposure_NAME/se_exposure_NAME, and "
            + "beta_outcome_NAME/se_outcome_NAME. The correlation matrix is "
            + "headerless, symmetric, positive definite, and ordered exactly as "
            + "--outcomes. PRESSO additionally requires --diagnostics and supports "
            + "--simulations 1000 --seed 20260916 --outlier-threshold 0.05. "
            + "Use mr-mvmr for multiple exposures with one outcome.";
    }

    private static final class Options {
        Path input, output, diagnostics, plot, outcomeCorrelation;
        List<String> exposures, outcomes; String method = "ivw-fixed";
        double confidence = .95, outlierThreshold = .05;
        int simulations = 1000; long seed = 20260916L; boolean help;
        BackendPolicy backend = BackendPolicy.PREFERRED;
        static Options parse(String[] args) {
            Options result = new Options();
            String exposureNames = null, outcomeNames = null;
            for (int index = 0; index < args.length; index++) switch (args[index]) {
                case "--input", "--in" -> result.input = Path.of(value(args, ++index, "--input"));
                case "--output", "--out" -> result.output = Path.of(value(args, ++index, "--output"));
                case "--diagnostics" -> result.diagnostics = Path.of(value(args, ++index, "--diagnostics"));
                case "--plot" -> result.plot = Path.of(value(args, ++index, "--plot"));
                case "--outcome-correlation" -> result.outcomeCorrelation = Path.of(value(args, ++index, "--outcome-correlation"));
                case "--exposures" -> exposureNames = value(args, ++index, "--exposures");
                case "--outcomes" -> outcomeNames = value(args, ++index, "--outcomes");
                case "--method" -> result.method = value(args, ++index, "--method").toLowerCase(Locale.ROOT);
                case "--confidence" -> result.confidence = Double.parseDouble(value(args, ++index, "--confidence"));
                case "--outlier-threshold" -> result.outlierThreshold = Double.parseDouble(value(args, ++index, "--outlier-threshold"));
                case "--simulations" -> result.simulations = Integer.parseInt(value(args, ++index, "--simulations"));
                case "--seed" -> result.seed = Long.parseLong(value(args, ++index, "--seed"));
                case "--backend" -> result.backend = BackendPolicy.valueOf(value(args, ++index, "--backend").toUpperCase(Locale.ROOT));
                case "--help", "-h" -> result.help = true;
                default -> throw new IllegalArgumentException("unknown mr-multivariate option: " + args[index]);
            }
            if (result.help) return result;
            if (result.input == null || result.outcomeCorrelation == null)
                throw new IllegalArgumentException(
                    "--input and --outcome-correlation are required");
            result.exposures = MrWideTable.names(exposureNames, "--exposures");
            result.outcomes = MrWideTable.names(outcomeNames, "--outcomes");
            if (result.outcomes.size() < 2)
                throw new IllegalArgumentException(
                    "mr-multivariate requires at least two outcome names");
            if (!List.of("ivw-fixed", "ivw-random", "presso").contains(result.method))
                throw new IllegalArgumentException(
                    "--method must be ivw-fixed, ivw-random, or presso");
            if (result.method.equals("presso") && result.diagnostics == null)
                throw new IllegalArgumentException(
                    "--diagnostics is required with --method presso");
            if (!result.method.equals("presso") && result.diagnostics != null)
                throw new IllegalArgumentException(
                    "--diagnostics is only valid with --method presso");
            if (!(result.confidence > 0 && result.confidence < 1)
                    || !(result.outlierThreshold > 0 && result.outlierThreshold < 1)
                    || result.simulations < 20)
                throw new IllegalArgumentException(
                    "invalid confidence, outlier threshold, or simulation count");
            return result;
        }
        private static String value(String[] args, int index, String option) {
            if (index >= args.length || args[index].isBlank())
                throw new IllegalArgumentException(option + " requires a value");
            return args[index];
        }
    }
}
