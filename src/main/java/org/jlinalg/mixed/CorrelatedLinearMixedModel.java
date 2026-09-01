/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceEstimation;

/** Cholesky-parameterized Gaussian LMM with correlated random blocks. */
public final class CorrelatedLinearMixedModel {
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    private static final double INVALID = Double.MAX_VALUE / 16.0;

    private CorrelatedLinearMixedModel() { }

    public static CorrelatedLinearMixedModelResult fit(
            double[] response, double[][] fixedEffects,
            List<CorrelatedRandomEffectBlock> blocks,
            RemlOptions options, BackendPolicy backendPolicy) {
        if (response == null)
            throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, fixed, response.length, fixedEffects[0].length,
            blocks, options, backendPolicy);
    }

    public static CorrelatedLinearMixedModelResult fit(
            double[] response, double[] fixed, int rows, int columns,
            List<CorrelatedRandomEffectBlock> blocks,
            RemlOptions options, BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixed, rows, columns);
        if (blocks == null || blocks.isEmpty() || options == null
                || backendPolicy == null)
            throw new IllegalArgumentException("blocks and options are required");
        if (options.degreesOfFreedomMethod()
                != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION)
            throw new IllegalArgumentException(
                "correlated-block fitting currently supports residual-approximation DF");
        int parameters = 1;
        for (CorrelatedRandomEffectBlock block : blocks) {
            if (block == null || block.observations() != rows)
                throw new IllegalArgumentException(
                    "correlated block rows must equal observations");
            parameters += block.effectCount() * (block.effectCount() + 1) / 2;
        }
        double sampleVariance = sampleVariance(response);
        double[] initial = new double[parameters];
        double[] lower = new double[parameters];
        double[] upper = new double[parameters];
        int parameter = 0;
        double initialSd = Math.sqrt(sampleVariance / (blocks.size() + 1));
        for (CorrelatedRandomEffectBlock block : blocks)
            for (int row = 0; row < block.effectCount(); row++)
                for (int column = 0; column <= row; column++) {
                    boolean diagonal = row == column;
                    initial[parameter] = diagonal ? Math.log(initialSd) : 0.0;
                    lower[parameter] = diagonal
                        ? 0.5 * Math.log(options.minimumVariance())
                        : -Math.sqrt(options.maximumVariance());
                    upper[parameter] = diagonal
                        ? 0.5 * Math.log(options.maximumVariance())
                        : Math.sqrt(options.maximumVariance());
                    parameter++;
                }
        initial[parameter] = Math.log(initialSd);
        lower[parameter] = 0.5 * Math.log(options.minimumVariance());
        upper[parameter] = 0.5 * Math.log(options.maximumVariance());
        int maximumEvaluations = Math.max(200,
            options.maximumIterations() * 20);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            Objective objective = new Objective(response, fixed, rows, columns,
                blocks, options.varianceEstimation(), backend);
            Optimum optimized = optimize(objective, initial, lower, upper,
                Math.min(0.5, options.maximumLogVarianceStep()),
                Math.max(1e-7, options.relativeTolerance()),
                maximumEvaluations);
            Evaluation fitted = objective.evaluate(optimized.parameters());
            double[] standardErrors = new double[columns];
            for (int column = 0; column < columns; column++)
                standardErrors[column] = Math.sqrt(Math.max(0.0,
                    fitted.fixedCovariance()[column * columns + column]));
            double degrees = rows - columns - 1.0;
            AssociationStatistics association = AssociationStatistics.studentT(
                fitted.beta(), standardErrors, degrees,
                DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
            return new CorrelatedLinearMixedModelResult(association,
                fitted.fixedCovariance(), fitted.estimates(),
                fitted.residualVariance(), fitted.fitted(), fitted.residuals(),
                fitted.logLikelihood(), options.varianceEstimation(),
                optimized.evaluations(), optimized.converged(),
                context.provenance());
        }
    }

    private static Optimum optimize(
            Objective objective, double[] initial, double[] lower,
            double[] upper, double initialStep, double tolerance,
            int maximumEvaluations) {
        double[] current = initial.clone();
        double currentValue = objective.value(current);
        int evaluations = 1;
        double step = initialStep;
        while (evaluations < maximumEvaluations && step > tolerance) {
            boolean improved = false;
            for (int parameter = 0; parameter < current.length
                    && evaluations < maximumEvaluations; parameter++) {
                double bestValue = currentValue;
                double bestParameter = current[parameter];
                for (int direction : new int[] {-1, 1}) {
                    double candidateValue = Math.max(lower[parameter],
                        Math.min(upper[parameter],
                            current[parameter] + direction * step));
                    if (candidateValue == current[parameter]) continue;
                    double[] candidate = current.clone();
                    candidate[parameter] = candidateValue;
                    double value = objective.value(candidate);
                    evaluations++;
                    if (value < bestValue) {
                        bestValue = value;
                        bestParameter = candidateValue;
                    }
                }
                if (bestValue < currentValue) {
                    current[parameter] = bestParameter;
                    currentValue = bestValue;
                    improved = true;
                }
            }
            if (!improved) step *= 0.5;
        }
        return new Optimum(current, evaluations, step <= tolerance);
    }

    private static final class Objective {
        private final double[] response;
        private final double[] fixed;
        private final int rows;
        private final int columns;
        private final List<CorrelatedRandomEffectBlock> blocks;
        private final VarianceEstimation estimation;
        private final ComputeBackend backend;

        Objective(double[] response, double[] fixed, int rows, int columns,
                List<CorrelatedRandomEffectBlock> blocks,
                VarianceEstimation estimation, ComputeBackend backend) {
            this.response = response;
            this.fixed = fixed;
            this.rows = rows;
            this.columns = columns;
            this.blocks = List.copyOf(blocks);
            this.estimation = estimation;
            this.backend = backend;
        }

        double value(double[] parameters) {
            try {
                double likelihood = evaluate(parameters).logLikelihood();
                return Double.isFinite(likelihood) ? -likelihood : INVALID;
            } catch (IllegalArgumentException | IllegalStateException exception) {
                return INVALID;
            }
        }

        Evaluation evaluate(double[] parameters) {
            Decoded decoded = decode(parameters, blocks);
            double[] covariance = new double[rows * rows];
            for (int row = 0; row < rows; row++)
                covariance[row * rows + row] = decoded.residualVariance();
            for (int blockIndex = 0; blockIndex < blocks.size(); blockIndex++) {
                CorrelatedRandomEffectBlock block = blocks.get(blockIndex);
                int[] groups = block.groupIndices();
                double[] design = block.effectDesign();
                double[] randomCovariance = decoded.covariances().get(blockIndex);
                int effects = block.effectCount();
                for (int row = 0; row < rows; row++)
                    for (int column = 0; column <= row; column++) {
                        if (groups[row] != groups[column]) continue;
                        double value = bilinear(design, row, column,
                            effects, randomCovariance);
                        covariance[row * rows + column] += value;
                        if (row != column)
                            covariance[column * rows + row] += value;
                    }
            }
            CholeskyFactor factor = backend.dpotrf(covariance, rows);
            double[] inverseX = factor.solve(fixed, columns);
            double[] information = MatrixOps.transposeMultiply(backend,
                fixed, rows, columns, inverseX, columns);
            CholeskyFactor fixedFactor = backend.dpotrf(information, columns);
            double[] inverseY = factor.solve(response);
            double[] right = new double[columns];
            backend.dgemv(MatrixTranspose.TRANSPOSE, rows, columns,
                1.0, fixed, inverseY, 0.0, right);
            double[] beta = fixedFactor.solve(right);
            double[] residual = MatrixOps.subtract(response,
                MatrixOps.multiply(backend, fixed, rows, columns, beta));
            double[] projected = factor.solve(residual);
            double quadratic = backend.ddot(rows, residual, 0, 1,
                projected, 0, 1);
            boolean restricted = estimation == VarianceEstimation.REML;
            double likelihood = -0.5 * (
                (restricted ? rows - columns : rows) * LOG_TWO_PI
                    + factor.logDeterminant()
                    + (restricted ? fixedFactor.logDeterminant() : 0.0)
                    + quadratic);
            double[] fixedCovariance = fixedFactor.solve(
                MatrixOps.identity(columns), columns);
            List<CorrelatedRandomEffectEstimates> estimates = new ArrayList<>();
            double[] fitted = MatrixOps.multiply(
                backend, fixed, rows, columns, beta);
            for (int blockIndex = 0; blockIndex < blocks.size(); blockIndex++) {
                CorrelatedRandomEffectBlock block = blocks.get(blockIndex);
                int effects = block.effectCount();
                int[] groups = block.groupIndices();
                double[] design = block.effectDesign();
                double[] randomCovariance = decoded.covariances().get(blockIndex);
                double[] modes = new double[block.groupCount() * effects];
                for (int row = 0; row < rows; row++)
                    for (int left = 0; left < effects; left++) {
                        double contribution = 0.0;
                        for (int rightEffect = 0; rightEffect < effects; rightEffect++)
                            contribution += randomCovariance[left * effects + rightEffect]
                                * design[row * effects + rightEffect];
                        modes[groups[row] * effects + left] +=
                            contribution * projected[row];
                    }
                for (int row = 0; row < rows; row++)
                    for (int effect = 0; effect < effects; effect++)
                        fitted[row] += design[row * effects + effect]
                            * modes[groups[row] * effects + effect];
                estimates.add(new CorrelatedRandomEffectEstimates(block.name(),
                    block.groupNames(), block.effectNames(), randomCovariance,
                    modes));
            }
            return new Evaluation(beta, fixedCovariance, estimates,
                decoded.residualVariance(), fitted,
                MatrixOps.subtract(response, fitted), likelihood);
        }
    }

    private static Decoded decode(
            double[] parameters, List<CorrelatedRandomEffectBlock> blocks) {
        List<double[]> covariance = new ArrayList<>(blocks.size());
        int parameter = 0;
        for (CorrelatedRandomEffectBlock block : blocks) {
            int effects = block.effectCount();
            double[] lower = new double[effects * effects];
            for (int row = 0; row < effects; row++)
                for (int column = 0; column <= row; column++)
                    lower[row * effects + column] = row == column
                        ? Math.exp(parameters[parameter++])
                        : parameters[parameter++];
            double[] value = new double[effects * effects];
            for (int row = 0; row < effects; row++)
                for (int column = 0; column <= row; column++) {
                    double sum = 0.0;
                    for (int index = 0; index <= Math.min(row, column); index++)
                        sum += lower[row * effects + index]
                            * lower[column * effects + index];
                    value[row * effects + column] = sum;
                    value[column * effects + row] = sum;
                }
            covariance.add(value);
        }
        double residual = Math.exp(2.0 * parameters[parameter]);
        return new Decoded(List.copyOf(covariance), residual);
    }

    private static double bilinear(
            double[] design, int first, int second, int effects,
            double[] covariance) {
        double value = 0.0;
        for (int left = 0; left < effects; left++)
            for (int right = 0; right < effects; right++)
                value += design[first * effects + left]
                    * covariance[left * effects + right]
                    * design[second * effects + right];
        return value;
    }

    private static double sampleVariance(double[] values) {
        double mean = Arrays.stream(values).average().orElse(0.0);
        double result = 0.0;
        for (double value : values) result += (value - mean) * (value - mean);
        return Math.max(1e-6, result / Math.max(1, values.length - 1));
    }

    private record Decoded(List<double[]> covariances,
                           double residualVariance) { }
    private record Evaluation(double[] beta, double[] fixedCovariance,
                              List<CorrelatedRandomEffectEstimates> estimates,
                              double residualVariance, double[] fitted,
                              double[] residuals, double logLikelihood) { }
    private record Optimum(double[] parameters, int evaluations,
                           boolean converged) { }
}
