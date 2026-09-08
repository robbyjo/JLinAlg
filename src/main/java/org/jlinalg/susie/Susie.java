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
import org.jlinalg.genetics.GeneticCovarianceValidation;

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
            if (!(scales[column] > 0.0) || !Double.isFinite(scales[column])) {
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
                || zScores.length == 0 || !(sampleSize > 2.0)
                || !Double.isFinite(sampleSize) || sampleSize > Integer.MAX_VALUE
                || sampleSize != Math.rint(sampleSize)) {
            throw new IllegalArgumentException("summary dimensions or sample size are invalid");
        }
        int columns = zScores.length;
        for (double[] row : ldCorrelation)
            if (row == null || row.length != columns)
                throw new IllegalArgumentException("LD must be square");
        double[] ld = MatrixOps.rowMajor(ldCorrelation, columns);
        validateLd(ld, columns);
        double[] xtx = ld.clone();
        double degreesOfFreedom = sampleSize - 1.0;
        for (int index = 0; index < xtx.length; index++) xtx[index] *= degreesOfFreedom;
        double[] xty = new double[columns];
        for (int column = 0; column < columns; column++) {
            if (!Double.isFinite(zScores[column])) throw new IllegalArgumentException("z scores must be finite");
            xty[column] = degreesOfFreedom
                * (zScores[column] / Math.hypot(zScores[column], Math.sqrt(sampleSize - 2.0)));
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
        if (columns == 0 || xtx == null || (long) columns * columns != xtx.length
                || !(yty > 0.0) || !Double.isFinite(yty) || observations < 2) {
            throw new IllegalArgumentException("sufficient statistics are invalid");
        }
        List<String> names = names(variableNames, columns);
        for (double value : xty)
            if (!Double.isFinite(value))
                throw new IllegalArgumentException("X'y must be finite");
        // These are claimed to be actual cross-products, unlike approximate
        // external-LD summary inputs: the joint [X y]'[X y] must also be PSD.
        int jointSize = columns + 1;
        double[] joint = new double[jointSize * jointSize];
        for (int row = 0; row < columns; row++) {
            System.arraycopy(xtx, row * columns, joint, row * jointSize, columns);
            joint[row * jointSize + columns] = joint[columns * jointSize + row] = xty[row];
        }
        joint[joint.length - 1] = yty;
        GeneticCovarianceValidation.requirePositiveSemidefinite(joint, jointSize);
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
        if (!(yty > 0) || !Double.isFinite(yty))
            throw new IllegalArgumentException("response sum of squares must be finite and positive");
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
                    double residualCross = xty[variable] - fittedCross[variable];
                    double se2 = residualVariance / diagonal;
                    double ratio = options.priorVariance() / se2;
                    double shrink = ratio <= 1 ? ratio / (1 + ratio)
                        : 1 / (1 + se2 / options.priorVariance());
                    conditionalVariance[variable] = ratio <= 1
                        ? options.priorVariance() / (1 + ratio) : se2 * shrink;
                    conditionalMean[variable] = ratio <= 1
                        ? multiplyDivide(conditionalVariance[variable], residualCross, residualVariance)
                        : multiplyDivide(shrink, residualCross, diagonal);
                    double logPenalty = ratio <= 1 ? Math.log1p(ratio)
                        : Math.log(options.priorVariance()) + Math.log(diagonal) - Math.log(residualVariance)
                            + Math.log1p(se2 / options.priorVariance());
                    logBf[variable] = 0.5 * (-logPenalty
                        + multiplyDivide(residualCross, conditionalMean[variable], residualVariance));
                    if (!Double.isFinite(logBf[variable]))
                        throw new IllegalArgumentException("single-effect posterior exceeds numerical range");
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
            if (!Double.isFinite(objective) || !Double.isFinite(expectedResidualSumSquares)
                    || expectedResidualSumSquares < 0)
                throw new IllegalArgumentException("inconsistent or numerically unrepresentable sufficient statistics");
            if (Math.abs(objective - previousObjective) < options.convergenceTolerance()) {
                converged = true;
                break;
            }
            if (options.estimateResidualVariance() && iteration < options.maximumIterations()) {
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
        for (int variable = 0; variable < variables; variable++) {
            for (int effect = 0; effect < core.effects(); effect++)
                pip[variable] += Math.log1p(-core.alpha()[effect * variables + variable]);
            pip[variable] = -Math.expm1(pip[variable]);
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

    /** Preserve representable products when either direct operation would overflow/underflow. */
    private static double multiplyDivide(double first, double second, double divisor) {
        int a = Math.getExponent(first), b = Math.getExponent(second), c = Math.getExponent(divisor);
        return Math.scalb((Math.scalb(first, -a) * Math.scalb(second, -b))
            / Math.scalb(divisor, -c), a + b - c);
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
                        / Math.sqrt(xtx[first * variables + first])
                        / Math.sqrt(xtx[second * variables + second]);
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
        GeneticCovarianceValidation.requirePositiveSemidefinite(ld, size);
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
        if (names.size() != columns || names.stream().anyMatch(value -> value == null || value.isBlank())
                || new java.util.HashSet<>(names).size() != columns)
            throw new IllegalArgumentException("one unique nonblank name is required per variable");
        return List.copyOf(names);
    }

    private record Core(double[] alpha, double[] mu, double[] logBayesFactors,
                        double[] posteriorMean,
                        double residualVariance, int effects, int iterations,
                        boolean converged, double objective) { }
}
