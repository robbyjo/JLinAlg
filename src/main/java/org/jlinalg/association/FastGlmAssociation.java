/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;

/** Prepared-null, block-score GLM association scanner. */
public final class FastGlmAssociation {
    private static final double MINIMUM_WORKING_WEIGHT = 1e-12;
    private static final double MAXIMUM_WORKING_WEIGHT = 1e12;

    private final int observations;
    private final int covariateCount;
    private final double[] weightedCovariates;
    private final double[] informationInverse;
    private final double[] workingResidual;
    private final double[] squareRootWorkingWeights;
    private final double dispersion;
    private final GlmResult nullModel;

    private FastGlmAssociation(
            int observations,
            int covariateCount,
            double[] weightedCovariates,
            double[] informationInverse,
            double[] workingResidual,
            double[] squareRootWorkingWeights,
            double dispersion,
            GlmResult nullModel) {
        this.observations = observations;
        this.covariateCount = covariateCount;
        this.weightedCovariates = weightedCovariates;
        this.informationInverse = informationInverse;
        this.workingResidual = workingResidual;
        this.squareRootWorkingWeights = squareRootWorkingWeights;
        this.dispersion = dispersion;
        this.nullModel = nullModel;
    }

    /** Fits the shared-covariate null model and caches its efficient-score metric. */
    public static FastGlmAssociation prepare(
            double[] response,
            double[][] covariates,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            GlmOptions glmOptions,
            AssociationEngineOptions engineOptions) {
        if (response == null)
            throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(covariates, response.length);
        return prepare(response, fixed, response.length, covariates[0].length,
            family, priorWeights, offset, glmOptions, engineOptions);
    }

    /** Contiguous row-major null-model preparation. */
    public static FastGlmAssociation prepare(
            double[] response,
            double[] covariates,
            int rows,
            int columns,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            GlmOptions glmOptions,
            AssociationEngineOptions engineOptions) {
        MatrixOps.validateModelData(response, covariates, rows, columns);
        if (family == null || glmOptions == null || engineOptions == null)
            throw new IllegalArgumentException("family and options are required");
        double[] weights = positiveWeights(priorWeights, rows);
        double[] offsets = offsets(offset, rows);
        GlmResult nullModel = Glm.fit(response, covariates, rows, columns,
            family, weights, offsets, glmOptions, engineOptions.backendPolicy());
        if (!nullModel.converged())
            throw new IllegalArgumentException(
                "GLM null model did not converge: "
                    + nullModel.convergenceMessage());
        double[] means = nullModel.fittedMeans();
        double[] eta = nullModel.linearPredictor();
        double[] squareRootWorkingWeights = new double[rows];
        double[] workingResidual = new double[rows];
        for (int row = 0; row < rows; row++) {
            double derivative = family.meanDerivative(eta[row]);
            double variance = family.variance(means[row]);
            if (!Double.isFinite(derivative) || derivative == 0.0
                    || !Double.isFinite(variance) || !(variance > 0.0))
                throw new IllegalArgumentException(
                    "family produced invalid null-model score weights");
            double workingWeight = clamp(
                weights[row] * derivative * derivative / variance,
                MINIMUM_WORKING_WEIGHT, MAXIMUM_WORKING_WEIGHT);
            squareRootWorkingWeights[row] = Math.sqrt(workingWeight);
            workingResidual[row] = squareRootWorkingWeights[row]
                * (response[row] - means[row]) / derivative;
        }
        double[] weightedCovariates = weight(
            covariates, rows, columns, squareRootWorkingWeights);
        double[] informationInverse = nullModel.covariance();
        double dispersion = nullModel.dispersion();
        if (!(dispersion > 0.0) || !Double.isFinite(dispersion))
            throw new IllegalArgumentException(
                "null-model dispersion must be finite and positive");
        for (int index = 0; index < informationInverse.length; index++)
            informationInverse[index] /= dispersion;
        return new FastGlmAssociation(rows, columns, weightedCovariates,
            informationInverse, workingResidual, squareRootWorkingWeights,
            dispersion, nullModel);
    }

    /** Scans an observation-by-predictor matrix in prepared score blocks. */
    public AssociationBatchResult scan(
            double[][] predictors,
            List<String> predictorNames,
            AssociationEngineOptions options) {
        if (predictors == null || predictors.length != observations)
            throw new IllegalArgumentException(
                "predictor row count must equal null-model observations");
        double[] values = MatrixOps.rowMajorUnchecked(predictors, observations);
        return scan(values, predictors[0].length, predictorNames, options);
    }

    /** Scans contiguous row-major predictors using efficient one-step estimates. */
    public AssociationBatchResult scan(
            double[] predictors,
            int predictorCount,
            List<String> predictorNames,
            AssociationEngineOptions options) {
        if (predictors == null || predictorCount < 1
                || predictors.length != observations * predictorCount
                || options == null)
            throw new IllegalArgumentException(
                "predictor dimensions and options are required");
        List<String> names = names(predictorNames, predictorCount);
        double[] beta = nan(predictorCount);
        double[] standardErrors = nan(predictorCount);
        ConcurrentLinkedQueue<AssociationFailure> failures =
            new ConcurrentLinkedQueue<>();
        long started = System.nanoTime();
        int blocks = (predictorCount + options.chunkSize() - 1)
            / options.chunkSize();
        execute(blocks, options.parallelism(), blockIndex -> {
            int first = blockIndex * options.chunkSize();
            int count = Math.min(options.chunkSize(), predictorCount - first);
            PreparedBlock prepared = prepareBlock(predictors, predictorCount,
                first, count, names, options, failures);
            try (BackendContext context = BackendContext.select(
                    options.backendPolicy())) {
                ComputeBackend backend = context.backend();
                double[] fixedCross = MatrixOps.transposeMultiply(backend,
                    weightedCovariates, observations, covariateCount,
                    prepared.values(), count);
                double[] projectionCoefficients = MatrixOps.multiply(backend,
                    informationInverse, covariateCount, covariateCount,
                    fixedCross, count);
                double[] fixedProjection = MatrixOps.multiply(backend,
                    weightedCovariates, observations, covariateCount,
                    projectionCoefficients, count);
                double[] block = prepared.values();
                for (int index = 0; index < block.length; index++)
                    block[index] -= fixedProjection[index];
                for (int variable = 0; variable < count; variable++) {
                    if (!prepared.valid()[variable]) continue;
                    int destination = first + variable;
                    double score = 0.0;
                    double information = 0.0;
                    for (int row = 0; row < observations; row++) {
                        double value = block[row * count + variable];
                        score += value * workingResidual[row];
                        information += value * value;
                    }
                    if (!(information > 1e-14)
                            || !Double.isFinite(information)) {
                        fail(destination, names.get(destination),
                            "predictor is constant or collinear with null covariates",
                            options.failurePolicy(), failures);
                        continue;
                    }
                    beta[destination] = score / information;
                    standardErrors[destination] = Math.sqrt(
                        dispersion / information);
                }
            }
        });
        double degreesOfFreedom = observations - covariateCount - 1.0;
        if (!(degreesOfFreedom > 0.0))
            throw new IllegalArgumentException(
                "GLM association scan requires positive approximate DFE");
        AssociationStatistics statistics = AssociationStatistics.studentT(
            beta, standardErrors, degreesOfFreedom,
            DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
        List<AssociationFailure> ordered = failures.stream()
            .sorted(Comparator.comparingInt(AssociationFailure::index)).toList();
        return new AssociationBatchResult(names, statistics.beta(),
            statistics.standardErrors(), statistics.statistics(),
            statistics.pValues(), statistics.degreesOfFreedom(),
            statistics.negativeLog10PValues(), ordered, covariateCount,
            Math.min(blocks, options.parallelism()),
            System.nanoTime() - started);
    }

    /** The fitted covariate-only GLM used by every marker score test. */
    public GlmResult nullModel() { return nullModel; }

    private PreparedBlock prepareBlock(
            double[] predictors, int total, int first, int count,
            List<String> names, AssociationEngineOptions options,
            ConcurrentLinkedQueue<AssociationFailure> failures) {
        double[] result = new double[observations * count];
        boolean[] valid = new boolean[count];
        Arrays.fill(valid, true);
        for (int variable = 0; variable < count; variable++) {
            int source = first + variable;
            double sum = 0.0;
            int finite = 0;
            int firstMissing = -1;
            for (int row = 0; row < observations; row++) {
                double value = predictors[row * total + source];
                if (Double.isFinite(value)) {
                    sum += value;
                    finite++;
                } else if (firstMissing < 0) firstMissing = row;
            }
            String message = null;
            if (finite == 0) message = "predictor has no finite values";
            else if (firstMissing >= 0
                    && options.predictorMissingPolicy()
                        == VariableMissingPolicy.ERROR)
                message = "non-finite predictor at row " + firstMissing;
            if (message != null) {
                fail(source, names.get(source), message,
                    options.failurePolicy(), failures);
                valid[variable] = false;
                continue;
            }
            double mean = sum / finite;
            for (int row = 0; row < observations; row++) {
                double value = predictors[row * total + source];
                result[row * count + variable] =
                    squareRootWorkingWeights[row]
                        * (Double.isFinite(value) ? value : mean);
            }
        }
        return new PreparedBlock(result, valid);
    }

    private static void fail(
            int index, String name, String message,
            AssociationFailurePolicy policy,
            ConcurrentLinkedQueue<AssociationFailure> failures) {
        if (policy == AssociationFailurePolicy.FAIL_FAST)
            throw new IllegalArgumentException(name + ": " + message);
        failures.add(new AssociationFailure(index, name,
            "NonEstimablePredictor", message));
    }

    private static double[] weight(
            double[] matrix, int rows, int columns, double[] squareRoots) {
        double[] result = new double[matrix.length];
        for (int row = 0; row < rows; row++)
            for (int column = 0; column < columns; column++)
                result[row * columns + column] =
                    squareRoots[row] * matrix[row * columns + column];
        return result;
    }

    private static double[] positiveWeights(double[] supplied, int rows) {
        double[] result = new double[rows];
        if (supplied == null) Arrays.fill(result, 1.0);
        else {
            if (supplied.length != rows)
                throw new IllegalArgumentException(
                    "prior weight length must equal rows");
            for (int row = 0; row < rows; row++) {
                if (!(supplied[row] > 0.0) || !Double.isFinite(supplied[row]))
                    throw new IllegalArgumentException(
                        "prior weights must be finite and positive");
                result[row] = supplied[row];
            }
        }
        return result;
    }

    private static double[] offsets(double[] supplied, int rows) {
        if (supplied == null) return new double[rows];
        if (supplied.length != rows)
            throw new IllegalArgumentException("offset length must equal rows");
        return MatrixOps.finiteCopy(supplied, "offset");
    }

    private static List<String> names(List<String> supplied, int count) {
        if (supplied == null) {
            List<String> result = new ArrayList<>(count);
            for (int index = 0; index < count; index++)
                result.add("variable" + (index + 1));
            return List.copyOf(result);
        }
        if (supplied.size() != count
                || supplied.stream().anyMatch(
                    value -> value == null || value.isBlank()))
            throw new IllegalArgumentException(
                "one nonblank name is required per predictor");
        return List.copyOf(supplied);
    }

    private static void execute(
            int blocks, int parallelism, BlockOperation operation) {
        ForkJoinPool pool = new ForkJoinPool(Math.min(blocks, parallelism));
        try {
            pool.submit(() -> IntStream.range(0, blocks).parallel()
                .forEach(operation::run)).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException(
                "GLM association scan was interrupted", exception);
        } catch (ExecutionException exception) {
            Throwable cause = exception.getCause();
            if (cause instanceof RuntimeException runtime) throw runtime;
            throw new IllegalStateException("GLM association scan failed", cause);
        } finally {
            pool.shutdownNow();
        }
    }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    private static double[] nan(int size) {
        double[] result = new double[size];
        Arrays.fill(result, Double.NaN);
        return result;
    }

    @FunctionalInterface
    private interface BlockOperation { void run(int block); }
    private record PreparedBlock(double[] values, boolean[] valid) { }
}
