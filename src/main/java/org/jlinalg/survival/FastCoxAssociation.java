/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import org.jlinalg.association.AssociationBatchResult;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.association.AssociationFailure;
import org.jlinalg.association.AssociationFailurePolicy;
import org.jlinalg.association.VariableMissingPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.MatrixOps;

/** Prepared-null, bounded-block score scanner for Cox association tests. */
public final class FastCoxAssociation {
    private final CoxSurvivalData survival;
    private final int observations;
    private final int covariateCount;
    private final double[] covariates;
    private final double[] offset;
    private final CoxOptions coxOptions;
    private final CoxResult nullModel;
    private final CoxRiskSetPlan rightCensoredPlan;
    private final CoxCountingProcessPlan countingProcessPlan;
    private final List<String> clusterIds;
    private final double[] scoreCorrelation;

    private FastCoxAssociation(CoxSurvivalData survival, int observations,
            int covariateCount, double[] covariates, double[] offset,
            CoxOptions coxOptions, CoxResult nullModel,
            List<String> clusterIds, double[] scoreCorrelation) {
        this.survival = survival;
        this.observations = observations;
        this.covariateCount = covariateCount;
        this.covariates = covariates;
        this.offset = offset;
        this.coxOptions = coxOptions;
        this.nullModel = nullModel;
        this.rightCensoredPlan = CoxRiskSetPlan.prepare(survival);
        this.countingProcessPlan = rightCensoredPlan == null
            ? CoxCountingProcessPlan.prepare(survival) : null;
        this.clusterIds = clusterIds;
        this.scoreCorrelation = scoreCorrelation;
    }

    public static FastCoxAssociation prepare(
            CoxSurvivalData survival, double[][] covariates,
            double[] offset, CoxOptions coxOptions,
            AssociationEngineOptions engineOptions) {
        return prepare(survival, covariates, offset, coxOptions,
            engineOptions, null, null);
    }

    /**
     * Prepares one null fit. Cluster IDs enable sandwich score variance;
     * an optional row-major observation correlation enables relatedness-aware
     * score variance.
     */
    public static FastCoxAssociation prepare(
            CoxSurvivalData survival, double[][] covariates,
            double[] offset, CoxOptions coxOptions,
            AssociationEngineOptions engineOptions,
            List<String> clusterIds, double[] scoreCorrelation) {
        if (survival == null || covariates == null || coxOptions == null
                || engineOptions == null
                || covariates.length != survival.observations())
            throw new IllegalArgumentException(
                "Cox score-scan inputs are required and row-aligned");
        int rows = survival.observations();
        double[] fixed = MatrixOps.rowMajor(covariates, rows);
        double[] offsets = offset == null ? new double[rows] : offset.clone();
        if (offsets.length != rows)
            throw new IllegalArgumentException("one Cox offset is required per row");
        if (clusterIds != null && (clusterIds.size() != rows
                || clusterIds.stream().anyMatch(value -> value == null)))
            throw new IllegalArgumentException(
                "one nonnull cluster id is required per row");
        if (scoreCorrelation != null
                && scoreCorrelation.length != Math.multiplyExact(rows, rows))
            throw new IllegalArgumentException(
                "score correlation must be an observation-square matrix");
        CoxResult nullModel = CoxRegression.fit(survival, covariates,
            offsets, coxOptions, engineOptions.backendPolicy());
        if (!nullModel.converged())
            throw new IllegalArgumentException("Cox null model did not converge: "
                + nullModel.convergenceMessage());
        return new FastCoxAssociation(survival, rows, covariates[0].length,
            fixed, offsets, coxOptions, nullModel,
            clusterIds == null ? null : List.copyOf(clusterIds),
            scoreCorrelation == null ? null : scoreCorrelation.clone());
    }

    public AssociationBatchResult scan(double[][] predictors,
            List<String> predictorNames, AssociationEngineOptions options,
            CoxScoreVariance variance) {
        if (predictors == null || predictors.length != observations)
            throw new IllegalArgumentException(
                "predictor rows must equal null-model observations");
        return scan(MatrixOps.rowMajorUnchecked(predictors, observations),
            predictors[0].length, predictorNames, options, variance);
    }

    public AssociationBatchResult scan(double[] predictors, int predictorCount,
            List<String> predictorNames, AssociationEngineOptions options,
            CoxScoreVariance variance) {
        if (predictors == null || predictorCount < 1
                || predictors.length != observations * predictorCount
                || options == null || variance == null)
            throw new IllegalArgumentException(
                "predictor dimensions and score controls are required");
        if (variance == CoxScoreVariance.CLUSTER_ROBUST && clusterIds == null)
            throw new IllegalArgumentException(
                "cluster-robust score variance requires cluster IDs");
        if (variance == CoxScoreVariance.RELATEDNESS
                && scoreCorrelation == null)
            throw new IllegalArgumentException(
                "relatedness score variance requires a score correlation");
        List<String> names = names(predictorNames, predictorCount);
        double[] beta = nan(predictorCount);
        double[] standardErrors = nan(predictorCount);
        ConcurrentLinkedQueue<AssociationFailure> failures =
            new ConcurrentLinkedQueue<>();
        long started = System.nanoTime();
        int blocks = (predictorCount + options.chunkSize() - 1)
            / options.chunkSize();
        execute(blocks, options.parallelism(), block -> {
            int first = block * options.chunkSize();
            int count = Math.min(options.chunkSize(), predictorCount - first);
            try {
                evaluateBlock(predictors, predictorCount, first, count,
                    names, options, variance, beta, standardErrors, failures);
            } catch (RuntimeException exception) {
                if (options.failurePolicy() == AssociationFailurePolicy.FAIL_FAST)
                    throw exception;
                for (int index = first; index < first + count; index++)
                    failures.add(new AssociationFailure(index, names.get(index),
                        exception.getClass().getSimpleName(), message(exception)));
            }
        });
        AssociationStatistics statistics = AssociationStatistics.normal(
            beta, standardErrors);
        List<AssociationFailure> ordered = failures.stream()
            .sorted(Comparator.comparingInt(AssociationFailure::index)).toList();
        return new AssociationBatchResult(names, statistics.beta(),
            statistics.standardErrors(), statistics.statistics(),
            statistics.pValues(), statistics.degreesOfFreedom(),
            statistics.negativeLog10PValues(), ordered, covariateCount,
            Math.min(blocks, options.parallelism()),
            System.nanoTime() - started);
    }

    public CoxResult nullModel() { return nullModel; }
    public BackendProvenance backend() { return nullModel.backend(); }
    public int observations() { return observations; }

    private void evaluateBlock(double[] predictors, int total, int first,
            int count, List<String> names, AssociationEngineOptions options,
            CoxScoreVariance variance, double[] beta, double[] standardErrors,
            ConcurrentLinkedQueue<AssociationFailure> failures) {
        int columns = covariateCount + count;
        double[] design = new double[observations * columns];
        boolean[] valid = new boolean[count];
        Arrays.fill(valid, true);
        for (int row = 0; row < observations; row++)
            System.arraycopy(covariates, row * covariateCount, design,
                row * columns, covariateCount);
        for (int variable = 0; variable < count; variable++) {
            int source = first + variable;
            double sum = 0.0;
            int finite = 0;
            int missing = -1;
            for (int row = 0; row < observations; row++) {
                double value = predictors[row * total + source];
                if (Double.isFinite(value)) { sum += value; finite++; }
                else if (missing < 0) missing = row;
            }
            String failure = finite == 0 ? "predictor has no finite values"
                : missing >= 0 && options.predictorMissingPolicy()
                    == VariableMissingPolicy.ERROR
                    ? "non-finite predictor at row " + missing : null;
            if (failure != null) {
                fail(source, names.get(source), failure, options, failures);
                valid[variable] = false;
                continue;
            }
            double mean = sum / finite;
            for (int row = 0; row < observations; row++) {
                double value = predictors[row * total + source];
                design[row * columns + covariateCount + variable] =
                    Double.isFinite(value) ? value : mean;
            }
        }
        double[] coefficients = new double[columns];
        System.arraycopy(nullModel.beta(), 0, coefficients, 0, covariateCount);
        CoxPartialLikelihood.Evaluation evaluation =
            CoxPartialLikelihood.evaluate(survival, design, columns,
                coefficients, offset, coxOptions.ties(), rightCensoredPlan,
                countingProcessPlan);
        double[] nullCovariance = nullModel.covariance();
        double[] residuals = variance == CoxScoreVariance.MODEL_BASED ? null
            : CoxDiagnostics.scoreResiduals(survival, design, columns,
                coefficients, offset, coxOptions.ties());
        for (int variable = 0; variable < count; variable++) {
            if (!valid[variable]) continue;
            int destination = first + variable;
            int predictor = covariateCount + variable;
            double[] projection = new double[covariateCount];
            for (int left = 0; left < covariateCount; left++)
                for (int right = 0; right < covariateCount; right++)
                    projection[left] += evaluation.information()[
                        predictor * columns + right]
                        * nullCovariance[right * covariateCount + left];
            double information = evaluation.information()[
                predictor * columns + predictor];
            for (int column = 0; column < covariateCount; column++)
                information -= projection[column]
                    * evaluation.information()[column * columns + predictor];
            if (!(information > 1e-14) || !Double.isFinite(information)) {
                fail(destination, names.get(destination),
                    "predictor is constant or collinear with null covariates",
                    options, failures);
                continue;
            }
            double score = evaluation.score()[predictor];
            beta[destination] = score / information;
            double varianceValue = 1.0 / information;
            if (residuals != null) {
                double[] efficient = new double[observations];
                for (int row = 0; row < observations; row++) {
                    efficient[row] = residuals[row * columns + predictor];
                    for (int column = 0; column < covariateCount; column++)
                        efficient[row] -= projection[column]
                            * residuals[row * columns + column];
                }
                double meat = variance == CoxScoreVariance.CLUSTER_ROBUST
                    ? clusterMeat(efficient) : relatednessMeat(efficient);
                varianceValue = meat / (information * information);
            }
            standardErrors[destination] = Math.sqrt(Math.max(0.0,
                varianceValue));
        }
    }

    private double clusterMeat(double[] score) {
        Map<String, Double> sums = new LinkedHashMap<>();
        for (int row = 0; row < observations; row++)
            sums.merge(clusterIds.get(row), score[row], Double::sum);
        double result = sums.values().stream()
            .mapToDouble(value -> value * value).sum();
        return result;
    }

    private double relatednessMeat(double[] score) {
        double result = 0.0;
        for (int row = 0; row < observations; row++)
            for (int column = 0; column < observations; column++)
                result += score[row]
                    * scoreCorrelation[row * observations + column]
                    * score[column];
        return result;
    }

    private static void fail(int index, String name, String message,
            AssociationEngineOptions options,
            ConcurrentLinkedQueue<AssociationFailure> failures) {
        if (options.failurePolicy() == AssociationFailurePolicy.FAIL_FAST)
            throw new IllegalArgumentException(name + ": " + message);
        failures.add(new AssociationFailure(index, name,
            "NonEstimablePredictor", message));
    }

    private static List<String> names(List<String> supplied, int count) {
        if (supplied == null) {
            List<String> result = new ArrayList<>(count);
            for (int index = 0; index < count; index++)
                result.add("variable" + (index + 1));
            return List.copyOf(result);
        }
        if (supplied.size() != count || supplied.stream()
                .anyMatch(value -> value == null || value.isBlank()))
            throw new IllegalArgumentException(
                "one nonblank name is required per predictor");
        return List.copyOf(supplied);
    }

    private static void execute(int blocks, int parallelism,
            java.util.function.IntConsumer operation) {
        ForkJoinPool pool = new ForkJoinPool(Math.min(blocks, parallelism));
        try {
            pool.submit(() -> IntStream.range(0, blocks).parallel()
                .forEach(operation)).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("Cox score scan was interrupted",
                exception);
        } catch (ExecutionException exception) {
            Throwable cause = exception.getCause();
            if (cause instanceof RuntimeException runtime) throw runtime;
            throw new IllegalStateException("Cox score scan failed", cause);
        } finally {
            pool.shutdownNow();
        }
    }

    private static String message(RuntimeException exception) {
        return exception.getMessage() == null ? "" : exception.getMessage();
    }

    private static double[] nan(int size) {
        double[] result = new double[size];
        Arrays.fill(result, Double.NaN);
        return result;
    }
}
