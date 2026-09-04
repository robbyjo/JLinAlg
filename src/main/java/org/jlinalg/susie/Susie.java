/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.susie;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Sum of Single Effects regression using the IBSS algorithm. */
public final class Susie {
    private Susie() { }

    /** Fits individual-level data after one-time centering and scaling. */
    public static SusieResult fit(
            double[] response, double[][] design, List<String> variableNames) {
        return fit(response, design, variableNames,
            SusieOptions.defaults(), BackendPolicy.CPU);
    }

    public static SusieResult fit(
            double[] response, double[][] design, List<String> variableNames,
            SusieOptions options, BackendPolicy backendPolicy) {
        if (response == null || design == null || response.length != design.length) {
            throw new IllegalArgumentException("response and design dimensions are invalid");
        }
        int rows = response.length;
        if (rows < 2 || design[0] == null || design[0].length < 1) {
            throw new IllegalArgumentException("response and design dimensions are invalid");
        }
        int columns = design[0].length;
        for (double[] row : design) {
            if (row == null || row.length != columns) {
                throw new IllegalArgumentException("design must be rectangular");
            }
        }
        List<String> names = names(variableNames, columns);
        double responseMean = Arrays.stream(response).average().orElseThrow();
        if (!Double.isFinite(responseMean)) throw new IllegalArgumentException("response must be finite");
        double[] centeredResponse = new double[rows];
        double[] matrix = new double[rows * columns];
        double[] means = new double[columns];
        double[] scales = new double[columns];
        for (int row = 0; row < rows; row++) centeredResponse[row] = response[row] - responseMean;
        for (int column = 0; column < columns; column++) {
            int offset = column * rows;
            for (int row = 0; row < rows; row++) {
                double value = design[row][column];
                if (!Double.isFinite(value)) throw new IllegalArgumentException("design must be finite");
                means[column] += value;
            }
            means[column] /= rows;
            for (int row = 0; row < rows; row++) {
                double centered = design[row][column] - means[column];
                matrix[offset + row] = centered;
                scales[column] += centered * centered;
            }
            scales[column] = Math.sqrt(scales[column] / (rows - 1.0));
            if (!(scales[column] > 0.0)) {
                throw new IllegalArgumentException("constant design column: " + names.get(column));
            }
            for (int row = 0; row < rows; row++) {
                matrix[offset + row] /= scales[column];
            }
        }
        double yty = 0.0;
        for (double value : centeredResponse) yty += value * value;
        double[] xtx = new double[columns * columns];
        double[] xty = new double[columns];
        IntStream columnStream = IntStream.range(0, columns);
        if ((long) rows * columns * columns >= 10_000_000L) {
            columnStream = columnStream.parallel();
        }
        columnStream.forEach(column -> {
            int offset = column * rows;
            xty[column] = dot(matrix, offset, centeredResponse, 0, rows);
            for (int other = 0; other <= column; other++) {
                double product = dot(matrix, offset, matrix, other * rows, rows);
                xtx[column * columns + other] = product;
                xtx[other * columns + column] = product;
            }
        });
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            Core core = fitCore(xtx, xty, yty, rows, names, options, backend);
            double[] beta = core.posteriorMean().clone();
            double[] effectMean = core.mu().clone();
            for (int column = 0; column < columns; column++) beta[column] /= scales[column];
            for (int effect = 0; effect < core.effects(); effect++) {
                for (int column = 0; column < columns; column++)
                    effectMean[effect * columns + column] /= scales[column];
            }
            double intercept = responseMean;
            for (int column = 0; column < columns; column++) intercept -= beta[column] * means[column];
            return result(core, beta, effectMean, intercept,
                names, xtx, options, context);
        }
    }

    /** Fits standardized summary statistics, with z and allele-aligned LD order matching. */
    public static SusieResult fitSummary(
            double[] zScores, double[][] ldCorrelation, double sampleSize,
            List<String> variableNames, SusieOptions options,
            BackendPolicy backendPolicy) {
        if (zScores == null || ldCorrelation == null
                || ldCorrelation.length != zScores.length
                || !(sampleSize > 1.0)) {
            throw new IllegalArgumentException("summary dimensions or sample size are invalid");
        }
        int columns = zScores.length;
        double[] ld = MatrixOps.rowMajor(ldCorrelation, columns);
        validateLd(ld, columns);
        double[] xtx = ld.clone();
        double degreesOfFreedom = sampleSize - 1.0;
        for (int index = 0; index < xtx.length; index++) xtx[index] *= degreesOfFreedom;
        double[] xty = new double[columns];
        for (int column = 0; column < columns; column++) {
            if (!Double.isFinite(zScores[column])) throw new IllegalArgumentException("z scores must be finite");
            double adjustment = degreesOfFreedom
                / (zScores[column] * zScores[column] + sampleSize - 2.0);
            xty[column] = Math.sqrt(degreesOfFreedom * adjustment) * zScores[column];
        }
        List<String> names = names(variableNames, columns);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            Core core = fitCore(xtx, xty, degreesOfFreedom, (int) Math.round(sampleSize),
                names, options, context.backend());
            return result(core, core.posteriorMean(), core.mu(), 0.0,
                names, xtx, options, context);
        }
    }

    /** Fits caller-supplied sufficient statistics X'X, X'y, and y'y. */
    public static SusieResult fitSufficientStatistics(
            double[] xtx, double[] xty, double yty, int observations,
            List<String> variableNames, SusieOptions options,
            BackendPolicy backendPolicy) {
        int columns = xty == null ? 0 : xty.length;
        if (xtx == null || xtx.length != columns * columns
                || !(yty > 0.0) || observations < 2) {
            throw new IllegalArgumentException("sufficient statistics are invalid");
        }
        List<String> names = names(variableNames, columns);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            Core core = fitCore(xtx.clone(), xty.clone(), yty, observations,
                names, options, context.backend());
            return result(core, core.posteriorMean(), core.mu(), 0.0,
                names, xtx, options, context);
        }
    }

    private static Core fitCore(
            double[] xtx, double[] xty, double yty, int observations,
            List<String> names, SusieOptions options, ComputeBackend backend) {
        if (options == null) throw new IllegalArgumentException("options are required");
        int variables = xty.length;
        int effects = Math.min(options.effects(), variables);
        double[] alpha = new double[effects * variables];
        double[] mu = new double[effects * variables];
        double[] second = new double[effects * variables];
        double[] logBayesFactors = new double[effects * variables];
        double[] total = new double[variables];
        double[] fittedCross = new double[variables];
        double[] effectCross = new double[effects * variables];
        double[] logBf = new double[variables];
        double[] conditionalMean = new double[variables];
        double[] conditionalVariance = new double[variables];
        double[] kl = new double[effects];
        double residualVariance = Math.max(1e-8, yty / (observations - 1.0));
        boolean converged = false;
        int iterations = 0;
        double objective = Double.NEGATIVE_INFINITY;
        for (int iteration = 1; iteration <= options.maximumIterations(); iteration++) {
            iterations = iteration;
            double previousObjective = objective;
            for (int effect = 0; effect < effects; effect++) {
                int offset = effect * variables;
                for (int variable = 0; variable < variables; variable++) {
                    total[variable] -= mu[offset + variable];
                    fittedCross[variable] -= effectCross[offset + variable];
                }
                double maximumLogBf = Double.NEGATIVE_INFINITY;
                for (int variable = 0; variable < variables; variable++) {
                    double diagonal = xtx[variable * variables + variable];
                    if (!(diagonal > 0.0)) throw new IllegalArgumentException("X'X has nonpositive diagonal");
                    double bhat = (xty[variable] - fittedCross[variable]) / diagonal;
                    double se2 = residualVariance / diagonal;
                    conditionalVariance[variable] = options.priorVariance() * se2
                        / (options.priorVariance() + se2);
                    conditionalMean[variable] = options.priorVariance()
                        / (options.priorVariance() + se2) * bhat;
                    logBf[variable] = 0.5 * (Math.log(se2
                        / (se2 + options.priorVariance()))
                        + bhat * bhat * options.priorVariance()
                            / (se2 * (se2 + options.priorVariance())));
                    maximumLogBf = Math.max(maximumLogBf, logBf[variable]);
                }
                double sum = 0.0;
                for (int variable = 0; variable < variables; variable++)
                    sum += Math.exp(logBf[variable] - maximumLogBf);
                double modelLogBayesFactor = maximumLogBf
                    + Math.log(sum * (1.0 / variables
                        + Math.sqrt(Math.ulp(1.0))));
                for (int variable = 0; variable < variables; variable++) {
                    double probability = Math.exp(logBf[variable] - maximumLogBf) / sum;
                    alpha[offset + variable] = probability;
                    logBayesFactors[offset + variable] = logBf[variable];
                    mu[offset + variable] = probability * conditionalMean[variable];
                    second[offset + variable] = probability
                        * (conditionalVariance[variable]
                            + conditionalMean[variable] * conditionalMean[variable]);
                    total[variable] += mu[offset + variable];
                }
                if ("cpu".equals(backend.selectedBackend())) {
                    multiply(xtx, variables, mu, offset, effectCross, offset);
                } else {
                    backend.dgemv(jdistlib.accelerator.MatrixTranspose.NONE,
                        variables, variables, 1.0, xtx, 0, variables,
                        mu, offset, 1, 0.0, effectCross, offset, 1);
                }
                double posteriorExpectedLogLikelihood = 0.0;
                for (int variable = 0; variable < variables; variable++) {
                    double residualCross = xty[variable] - fittedCross[variable];
                    posteriorExpectedLogLikelihood +=
                        -2.0 * mu[offset + variable] * residualCross
                        + xtx[variable * variables + variable]
                            * second[offset + variable];
                    fittedCross[variable] += effectCross[offset + variable];
                }
                posteriorExpectedLogLikelihood *= -0.5 / residualVariance;
                kl[effect] = -modelLogBayesFactor
                    + posteriorExpectedLogLikelihood;
            }
            double expectedResidualSumSquares = yty;
            for (int variable = 0; variable < variables; variable++) {
                expectedResidualSumSquares += total[variable]
                    * (fittedCross[variable] - 2.0 * xty[variable]);
            }
            for (int effect = 0; effect < effects; effect++) {
                int offset = effect * variables;
                for (int variable = 0; variable < variables; variable++) {
                    expectedResidualSumSquares +=
                        xtx[variable * variables + variable]
                            * second[offset + variable]
                        - mu[offset + variable]
                            * effectCross[offset + variable];
                }
            }
            objective = -0.5 * observations
                * Math.log(2.0 * Math.PI * residualVariance)
                - 0.5 * expectedResidualSumSquares / residualVariance;
            for (double divergence : kl) objective -= divergence;
            if (objective - previousObjective < options.convergenceTolerance()) {
                converged = true;
                break;
            }
            if (options.estimateResidualVariance()) {
                residualVariance = Math.max(1e-8,
                    expectedResidualSumSquares / observations);
            }
        }
        return new Core(alpha, mu, logBayesFactors, total, residualVariance,
            effects, iterations, converged, objective);
    }

    private static SusieResult result(
            Core core, double[] beta, double[] effectMean,
            double intercept, List<String> names,
            double[] xtx, SusieOptions options, BackendContext context) {
        int variables = names.size();
        double[] pip = new double[variables];
        Arrays.fill(pip, 1.0);
        for (int variable = 0; variable < variables; variable++) {
            for (int effect = 0; effect < core.effects(); effect++)
                pip[variable] *= 1.0 - core.alpha()[effect * variables + variable];
            pip[variable] = 1.0 - pip[variable];
        }
        List<CredibleSet> sets = credibleSets(
            core.alpha(), core.effects(), names, xtx, options);
        return new SusieResult(names, pip, beta, core.alpha(), effectMean,
            core.logBayesFactors(), sets,
            intercept, core.residualVariance(), core.effects(), core.iterations(),
            core.converged(), core.objective(), context.provenance());
    }

    private static double dot(
            double[] first, int firstOffset, double[] second,
            int secondOffset, int length) {
        double result = 0.0;
        for (int index = 0; index < length; index++) {
            result = Math.fma(first[firstOffset + index],
                second[secondOffset + index], result);
        }
        return result;
    }

    private static void multiply(
            double[] matrix, int size, double[] vector, int vectorOffset,
            double[] result, int resultOffset) {
        IntStream rows = IntStream.range(0, size);
        if (size >= 256) rows = rows.parallel();
        rows.forEach(row -> {
            double sum = 0.0;
            int offset = row * size;
            for (int column = 0; column < size; column++) {
                sum = Math.fma(matrix[offset + column],
                    vector[vectorOffset + column], sum);
            }
            result[resultOffset + row] = sum;
        });
    }

    private static List<CredibleSet> credibleSets(
            double[] alpha, int effects, List<String> names,
            double[] xtx, SusieOptions options) {
        int variables = names.size();
        List<CredibleSet> result = new ArrayList<>();
        for (int effect = 0; effect < effects; effect++) {
            int offset = effect * variables;
            Integer[] order = new Integer[variables];
            for (int index = 0; index < variables; index++) order[index] = index;
            Arrays.sort(order, Comparator.comparingDouble(
                (Integer index) -> alpha[offset + index]).reversed());
            List<Integer> selected = new ArrayList<>();
            double coverage = 0.0;
            for (int index : order) {
                selected.add(index);
                coverage += alpha[offset + index];
                if (coverage >= options.credibleSetCoverage()) break;
            }
            double purity = 1.0;
            for (int first : selected) {
                for (int second : selected) {
                    double correlation = xtx[first * variables + second]
                        / Math.sqrt(xtx[first * variables + first]
                            * xtx[second * variables + second]);
                    purity = Math.min(purity, Math.abs(correlation));
                }
            }
            if (purity >= options.minimumCredibleSetPurity()) {
                result.add(new CredibleSet(effect,
                    selected.stream().map(names::get).toList(), coverage, purity));
            }
        }
        return result;
    }

    private static void validateLd(double[] ld, int size) {
        for (int row = 0; row < size; row++) {
            if (Math.abs(ld[row * size + row] - 1.0) > 1e-8)
                throw new IllegalArgumentException("LD diagonal must equal one");
            for (int column = 0; column < row; column++) {
                if (Math.abs(ld[row * size + column] - ld[column * size + row]) > 1e-10
                        || Math.abs(ld[row * size + column]) > 1.0 + 1e-10)
                    throw new IllegalArgumentException("LD must be a symmetric correlation matrix");
            }
        }
    }

    private static List<String> names(List<String> names, int columns) {
        if (names == null) {
            List<String> result = new ArrayList<>(columns);
            for (int index = 0; index < columns; index++) result.add("variable" + (index + 1));
            return result;
        }
        if (names.size() != columns || names.stream().anyMatch(value -> value == null || value.isBlank()))
            throw new IllegalArgumentException("one nonblank name is required per variable");
        return List.copyOf(names);
    }

    private record Core(double[] alpha, double[] mu, double[] logBayesFactors,
                        double[] posteriorMean,
                        double residualVariance, int effects, int iterations,
                        boolean converged, double objective) { }
}
