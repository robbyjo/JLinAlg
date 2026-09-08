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
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.BetaRegression;
import org.jlinalg.distributional.BetaRegressionOptions;
import org.jlinalg.distributional.BetaRegressionResult;
import org.jlinalg.penalized.ElasticNetOptions;
import org.jlinalg.penalized.PenalizedCrossValidationResult;
import org.jlinalg.penalized.PenalizedRegression;
import org.jlinalg.penalized.PenalizedRegressionCrossValidation;
import org.jlinalg.penalized.PenalizedRegressionResult;

/** File-oriented beta and Gaussian penalized regression commands. */
final class RegressionCli {
    private RegressionCli() { }

    static int run(String command, String[] arguments,
            java.io.PrintStream output, java.io.PrintStream error) {
        try {
            Options options = Options.parse(command, arguments);
            if (options.help) { output.println(help(command)); return 0; }
            PipelinePaths.requireFreshOutputs(options.output);
            Table table = Table.read(options.input);
            String result = command.equals("beta-regression")
                ? beta(table, options) : penalized(table, options);
            if (options.output == null) output.print(result);
            else Files.writeString(options.output, result, java.nio.file.StandardOpenOption.CREATE_NEW);
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static String beta(Table table, Options options) {
        double[] response = table.column(options.response);
        double[][] mean = table.design(options.meanColumns, options.intercept);
        List<String> meanNames = table.designNames(options.meanColumns, options.intercept);
        List<String> precisionColumns = options.precisionColumns.isEmpty()
            ? List.of() : options.precisionColumns;
        double[][] precision = table.design(precisionColumns, options.intercept);
        BetaRegressionOptions betaOptions = precisionColumns.isEmpty()
            ? BetaRegressionOptions.constantPrecisionDefaults()
            : BetaRegressionOptions.variablePrecisionDefaults();
        BetaRegressionResult fit = BetaRegression.fit(response, mean, precision,
            betaOptions, BackendPolicy.CPU);
        StringBuilder result = new StringBuilder(
            "model\tblock\tterm\testimate\tstandard_error\tstatistic\tp_value\tconverged\n");
        double[] errors = fit.standardErrors(), statistics = fit.statistics(), pValues = fit.pValues();
        double[] meanBeta = fit.meanCoefficients(), precisionBeta = fit.precisionCoefficients();
        for (int i = 0; i < meanBeta.length; i++) row(result, "beta", "mean",
            meanNames.get(i), meanBeta[i], errors[i], statistics[i], pValues[i], fit.converged());
        List<String> precisionNames = table.designNames(precisionColumns, options.intercept);
        for (int i = 0; i < precisionBeta.length; i++) row(result, "beta", "precision",
            precisionNames.get(i), precisionBeta[i], errors[meanBeta.length + i],
            statistics[meanBeta.length + i], pValues[meanBeta.length + i], fit.converged());
        return result.toString();
    }

    private static String penalized(Table table, Options options) {
        double[] response = table.column(options.response);
        double[][] predictors = table.matrix(options.predictorColumns);
        ElasticNetOptions fitOptions = ElasticNetOptions.builder()
            .alpha(options.alpha).fitIntercept(options.intercept)
            .standardize(options.standardize).build();
        PenalizedRegressionResult fit;
        double selectedLambda = options.lambda;
        double selectedAlpha = options.alpha;
        if (options.cvFolds > 1) {
            PenalizedCrossValidationResult cv = PenalizedRegressionCrossValidation.fit(
                response, predictors, options.lambdas, options.cvFolds,
                options.seed, fitOptions);
            fit = cv.minimumErrorFit();
            selectedLambda = cv.lambdaMinimum();
        } else {
            fit = options.model.equals("ridge")
                ? PenalizedRegression.ridge(response, predictors, options.lambda, fitOptions)
                : options.model.equals("lasso")
                    ? PenalizedRegression.lasso(response, predictors, options.lambda, fitOptions)
                    : PenalizedRegression.fit(response, predictors, options.lambda, fitOptions);
        }
        StringBuilder result = new StringBuilder(
            "model\tterm\testimate\tstandard_error\tstatistic\tp_value\tlambda\talpha\tactive\tconverged\n");
        row(result, options.model, "(Intercept)", fit.intercept(), selectedLambda,
            selectedAlpha, fit.activeCoefficientCount(), fit.converged());
        double[] coefficients = fit.coefficients();
        for (int i = 0; i < coefficients.length; i++)
            row(result, options.model, options.predictorColumns.get(i), coefficients[i],
                selectedLambda, selectedAlpha, fit.activeCoefficientCount(), fit.converged());
        return result.toString();
    }

    private static void row(StringBuilder output, String model, String block,
            String term, double estimate, double error, double statistic,
            double pValue, boolean converged) {
        output.append(model).append('\t').append(block).append('\t').append(term)
            .append('\t').append(number(estimate)).append('\t').append(number(error))
            .append('\t').append(number(statistic)).append('\t').append(number(pValue))
            .append('\t').append(converged).append('\n');
    }

    private static void row(StringBuilder output, String model, String term,
            double estimate, double lambda, double alpha, int active,
            boolean converged) {
        output.append(model).append('\t').append(term).append('\t')
            .append(number(estimate)).append("\t\t\t\t")
            .append(number(lambda)).append('\t').append(number(alpha)).append('\t')
            .append(active).append('\t').append(converged).append('\n');
    }

    private static String number(double value) {
        return Double.isFinite(value) ? Double.toString(value) : "";
    }

    private static String help(String command) {
        if (command.equals("beta-regression")) return """
            Usage: beta-regression --input FILE --response COLUMN
              [--mean COLUMNS] [--precision COLUMNS] [--out FILE]

            Numeric CSV/TSV input. An intercept is included by default in both
            mean and precision designs; use --no-intercept to disable it.
            Responses must lie strictly between zero and one.
            """;
        return """
            Usage: penalized-regression --input FILE --response COLUMN
              --predictors COLUMNS --model ridge|lasso|elastic-net
              [--lambda VALUE | --lambda-grid a,b,c --cv-folds K]
              [--alpha VALUE] [--out FILE]

            Numeric CSV/TSV input. Intercept and standardization are enabled by
            default; use --no-intercept or --no-standardize to disable them.
            """;
    }

    private static final class Options {
        Path input, output; String response, model = "ridge";
        List<String> meanColumns = List.of(), precisionColumns = List.of(), predictorColumns = List.of();
        double lambda = 1.0, alpha = 0.5; double[] lambdas = {1.0};
        int cvFolds = 0; long seed = 20260831L; boolean intercept = true, standardize = true, help;

        static Options parse(String command, String[] arguments) {
            Options result = new Options();
            for (int i = 0; i < arguments.length; i++) {
                String option = arguments[i];
                switch (option) {
                    case "--input", "--in" -> result.input = Path.of(value(arguments, ++i, option));
                    case "--out", "--output" -> result.output = Path.of(value(arguments, ++i, option));
                    case "--response" -> result.response = value(arguments, ++i, option);
                    case "--mean" -> result.meanColumns = columns(value(arguments, ++i, option));
                    case "--precision" -> result.precisionColumns = columns(value(arguments, ++i, option));
                    case "--predictors" -> result.predictorColumns = columns(value(arguments, ++i, option));
                    case "--model" -> result.model = value(arguments, ++i, option).toLowerCase(Locale.ROOT);
                    case "--lambda" -> result.lambda = nonnegative(value(arguments, ++i, option), option);
                    case "--lambda-grid" -> result.lambdas = nonnegativeList(value(arguments, ++i, option), option);
                    case "--alpha" -> result.alpha = probability(value(arguments, ++i, option), option);
                    case "--cv-folds" -> result.cvFolds = positiveInteger(value(arguments, ++i, option), option);
                    case "--seed" -> result.seed = Long.parseLong(value(arguments, ++i, option));
                    case "--no-intercept" -> result.intercept = false;
                    case "--no-standardize" -> result.standardize = false;
                    case "--help", "-h" -> result.help = true;
                    default -> throw new IllegalArgumentException("unknown " + command + " option: " + option);
                }
            }
            if (!result.help && (result.input == null || result.response == null))
                throw new IllegalArgumentException("--input and --response are required");
            if (result.cvFolds == 1)
                throw new IllegalArgumentException("--cv-folds must be at least two");
            if (!result.help && command.equals("beta-regression")) {
                if (!result.meanColumns.isEmpty() && !result.precisionColumns.isEmpty()) { /* accepted */ }
            } else if (!result.help) {
                if (result.predictorColumns.isEmpty()) throw new IllegalArgumentException("--predictors is required");
                if (!List.of("ridge", "lasso", "elastic-net").contains(result.model))
                    throw new IllegalArgumentException("--model must be ridge, lasso, or elastic-net");
                if (result.model.equals("ridge")) result.alpha = 0.0;
                if (result.model.equals("lasso")) result.alpha = 1.0;
                if (result.cvFolds > 1 && result.lambdas.length < 2)
                    throw new IllegalArgumentException("--lambda-grid needs at least two values with --cv-folds");
            }
            return result;
        }
        private static String value(String[] values, int index, String option) { if (index >= values.length || values[index].isBlank()) throw new IllegalArgumentException(option + " requires a value"); return values[index]; }
        private static List<String> columns(String value) { List<String> result = Arrays.stream(value.split(",", -1)).map(String::trim).filter(v -> !v.isBlank()).toList(); if (result.isEmpty()) throw new IllegalArgumentException("column list must not be empty"); return result; }
        private static double nonnegative(String value, String option) { double result = Double.parseDouble(value); if (!(result >= 0) || !Double.isFinite(result)) throw new IllegalArgumentException(option + " must be finite and nonnegative"); return result; }
        private static double probability(String value, String option) { double result = Double.parseDouble(value); if (!(result >= 0 && result <= 1) || !Double.isFinite(result)) throw new IllegalArgumentException(option + " must be in [0,1]"); return result; }
        private static double[] nonnegativeList(String value, String option) { double[] result = Arrays.stream(value.split(",", -1)).mapToDouble(v -> nonnegative(v.trim(), option)).toArray(); if (result.length == 0) throw new IllegalArgumentException(option + " must not be empty"); return result; }
        private static int positiveInteger(String value, String option) { int result = Integer.parseInt(value); if (result < 1) throw new IllegalArgumentException(option + " must be positive"); return result; }
    }

    private static final class Table {
        final String[] header; final double[][] values; final char delimiter;
        private Table(String[] header, double[][] values, char delimiter) { this.header = header; this.values = values; this.delimiter = delimiter; }
        static Table read(Path path) throws IOException {
            List<String> lines = Files.readAllLines(path); if (lines.size() < 2) throw new IOException("input needs a header and at least one row");
            char delimiter = path.toString().toLowerCase(Locale.ROOT).endsWith(".csv") ? ',' : '\t';
            String[] header = split(lines.get(0), delimiter); if (Arrays.stream(header).anyMatch(String::isBlank)) throw new IOException("column names must not be blank");
            List<double[]> rows = new ArrayList<>();
            for (int line = 1; line < lines.size(); line++) { if (lines.get(line).isBlank()) continue; String[] fields = split(lines.get(line), delimiter); if (fields.length != header.length) throw new IOException("inconsistent field count at line " + (line + 1)); double[] row = new double[fields.length]; for (int i = 0; i < row.length; i++) { try { row[i] = Double.parseDouble(fields[i].trim()); } catch (NumberFormatException exception) { throw new IOException("non-numeric value at line " + (line + 1) + ", column " + header[i], exception); } } rows.add(row); }
            if (rows.isEmpty()) throw new IOException("input contains no rows");
            return new Table(header, rows.toArray(double[][]::new), delimiter);
        }
        double[] column(String name) { int index = index(name); double[] result = new double[values.length]; for (int i = 0; i < result.length; i++) result[i] = values[i][index]; return result; }
        double[][] matrix(List<String> names) { double[][] result = new double[values.length][names.size()]; for (int row = 0; row < result.length; row++) for (int column = 0; column < names.size(); column++) result[row][column] = values[row][index(names.get(column))]; return result; }
        double[][] design(List<String> names, boolean intercept) { double[][] result = new double[values.length][names.size() + (intercept ? 1 : 0)]; for (int row = 0; row < result.length; row++) { int offset = 0; if (intercept) result[row][offset++] = 1.0; for (String name : names) result[row][offset++] = values[row][index(name)]; } return result; }
        List<String> designNames(List<String> names, boolean intercept) { List<String> result = new ArrayList<>(); if (intercept) result.add("(Intercept)"); result.addAll(names); return result; }
        private int index(String name) { for (int i = 0; i < header.length; i++) if (header[i].trim().equalsIgnoreCase(name)) return i; throw new IllegalArgumentException("column is absent: " + name); }
        private static String[] split(String line, char delimiter) { List<String> fields = new ArrayList<>(); StringBuilder current = new StringBuilder(); boolean quoted = false; for (int i = 0; i < line.length(); i++) { char value = line.charAt(i); if (value == '"') quoted = !quoted; else if (value == delimiter && !quoted) { fields.add(current.toString()); current.setLength(0); } else current.append(value); } if (quoted) throw new IllegalArgumentException("unterminated quoted field"); fields.add(current.toString()); return fields.toArray(String[]::new); }
    }
}
