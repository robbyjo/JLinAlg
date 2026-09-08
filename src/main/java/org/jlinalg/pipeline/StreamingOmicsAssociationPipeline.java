/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.jlinalg.association.AssociationBatchResult;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.association.AssociationEstimate;
import org.jlinalg.association.AssociationFailure;
import org.jlinalg.association.FastGlmAssociation;
import org.jlinalg.association.FastOlsAssociation;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.survival.CoxOptions;
import org.jlinalg.survival.CoxScoreVariance;
import org.jlinalg.survival.CoxSurvivalData;
import org.jlinalg.survival.FastCoxAssociation;

/** Block-streamed fast OLS scans for transformed omics feature matrices. */
public final class StreamingOmicsAssociationPipeline {
    private StreamingOmicsAssociationPipeline() { }

    /**
     * Treats every matrix row as a changing response and returns the selected
     * fixed-design coefficient. This is the natural high-throughput QTL shape.
     */
    public static OmicsAssociationResult scanResponses(
            NumericMatrixSource source,
            List<String> analysisSampleIds,
            double[][] design,
            int testedCoefficient,
            OmicsTransform transform,
            OmicsMissingPolicy missingPolicy,
            int featureBlockSize,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions) throws IOException {
        validate(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy);
        if (design == null || design.length != analysisSampleIds.size())
            throw new IllegalArgumentException(
                "design rows must match analysis samples");
        return scan(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy, (matrix, names) ->
                FastOlsAssociation.scanResponses(matrix, design,
                    testedCoefficient, names, weights, offset,
                    olsOptions, engineOptions));
    }

    /**
     * Treats every matrix row as a changing predictor of one response. This
     * supports predictor-oriented TWAS, EWAS, PWAS, and metabolome scans.
     */
    public static OmicsAssociationResult scanPredictors(
            NumericMatrixSource source,
            List<String> analysisSampleIds,
            double[] response,
            double[][] covariates,
            OmicsTransform transform,
            OmicsMissingPolicy missingPolicy,
            int featureBlockSize,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions) throws IOException {
        validate(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy);
        if (response == null || covariates == null
                || response.length != analysisSampleIds.size()
                || covariates.length != analysisSampleIds.size())
            throw new IllegalArgumentException(
                "response and covariate rows must match analysis samples");
        return scan(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy, (matrix, names) ->
                FastOlsAssociation.scanPredictors(response, covariates,
                    matrix, names, weights, offset, olsOptions, engineOptions));
    }

    /** Streams predictor-oriented OLS results without retaining all features. */
    public static OmicsAssociationSummary scanPredictorsTo(
            NumericMatrixSource source,
            List<String> analysisSampleIds,
            double[] response,
            double[][] covariates,
            OmicsTransform transform,
            OmicsMissingPolicy missingPolicy,
            int featureBlockSize,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions,
            OmicsAssociationSink sink) throws IOException {
        validate(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy);
        if (sink == null || response == null || covariates == null
                || response.length != analysisSampleIds.size()
                || covariates.length != analysisSampleIds.size())
            throw new IllegalArgumentException(
                "response, covariates, and omics sink must match analysis samples");
        return scanTo(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy, (matrix, names) ->
                FastOlsAssociation.scanPredictors(response, covariates,
                    matrix, names, weights, offset, olsOptions, engineOptions),
            sink);
    }

    /** Streams prepared-null GLM score tests for generic numeric omics rows. */
    public static OmicsAssociationSummary scanPredictorsGlmTo(
            NumericMatrixSource source,
            List<String> analysisSampleIds,
            double[] response,
            double[][] covariates,
            GlmFamily family,
            OmicsTransform transform,
            OmicsMissingPolicy missingPolicy,
            int featureBlockSize,
            double[] weights,
            double[] offset,
            GlmOptions glmOptions,
            AssociationEngineOptions engineOptions,
            OmicsAssociationSink sink) throws IOException {
        validate(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy);
        if (sink == null || response == null || covariates == null
                || response.length != analysisSampleIds.size()
                || covariates.length != analysisSampleIds.size())
            throw new IllegalArgumentException(
                "response, covariates, and omics sink must match analysis samples");
        FastGlmAssociation prepared = FastGlmAssociation.prepare(response,
            covariates, family, weights, offset, glmOptions, engineOptions);
        return scanTo(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy,
            (matrix, names) -> prepared.scan(matrix, names, engineOptions),
            sink);
    }

    /** Streams prepared-null Cox score tests for file-backed omics rows. */
    public static OmicsAssociationSummary scanPredictorsCoxTo(
            NumericMatrixSource source,
            List<String> analysisSampleIds,
            CoxSurvivalData survival,
            double[][] covariates,
            CoxScoreVariance scoreVariance,
            List<String> clusterIds,
            double[] scoreCorrelation,
            OmicsTransform transform,
            OmicsMissingPolicy missingPolicy,
            int featureBlockSize,
            double[] offset,
            CoxOptions coxOptions,
            AssociationEngineOptions engineOptions,
            OmicsAssociationSink sink) throws IOException {
        validate(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy);
        if (sink == null || survival == null || covariates == null
                || survival.observations() != analysisSampleIds.size()
                || covariates.length != analysisSampleIds.size())
            throw new IllegalArgumentException(
                "survival, covariates, and sink must match analysis samples");
        FastCoxAssociation prepared = FastCoxAssociation.prepare(survival,
            covariates, offset, coxOptions, engineOptions,
            clusterIds, scoreCorrelation);
        return scanTo(source, analysisSampleIds, featureBlockSize,
            transform, missingPolicy, (matrix, names) -> prepared.scan(
                matrix, names, engineOptions, scoreVariance), sink);
    }

    private static OmicsAssociationResult scan(
            NumericMatrixSource source, List<String> analysisSampleIds,
            int blockSize, OmicsTransform transform,
            OmicsMissingPolicy missingPolicy, BlockScanner scanner)
            throws IOException {
        int[] sampleOrder = SampleAlignment.requireOrder(
            source.metadata().sampleIds(), analysisSampleIds);
        List<OmicsAssociationEstimate> estimates = new ArrayList<>();
        List<AssociationPipelineFailure> failures = new ArrayList<>();
        long sourceFeatures = 0;
        try (NumericBlockReader reader = source.open(sampleOrder)) {
            for (NumericBlock block; (block = reader.read(blockSize)) != null;) {
                sourceFeatures += block.rows().size();
                List<String> names = block.rows().stream()
                    .map(NumericRow::id).toList();
                double[][] matrix = new double[analysisSampleIds.size()]
                    [block.rows().size()];
                for (int feature = 0; feature < block.rows().size(); feature++) {
                    double[] values = transformed(block.rows().get(feature),
                        transform, analysisSampleIds.size());
                    prepareMissing(values, missingPolicy,
                        block.rows().get(feature).id());
                    for (int sample = 0; sample < values.length; sample++)
                        matrix[sample][feature] = values[sample];
                }
                AssociationBatchResult result = scanner.scan(matrix, names);
                AssociationFailure[] failedRows = failureRows(result);
                for (int index = 0; index < result.size(); index++) {
                    AssociationEstimate estimate = result.estimate(index);
                    AssociationPipelineFailure failure = failure(estimate, failedRows[index]);
                    if (failure != null) { failures.add(failure); continue; }
                    estimates.add(new OmicsAssociationEstimate(estimate.name(),
                        estimate.beta(), estimate.standardError(),
                        estimate.statistic(), estimate.degreesOfFreedom(),
                        estimate.pValue(), estimate.log10PValue(),
                        estimate.negativeLog10PValue()));
                }
            }
        }
        return new OmicsAssociationResult(
            sourceFeatures, estimates, failures);
    }

    private static OmicsAssociationSummary scanTo(
            NumericMatrixSource source, List<String> analysisSampleIds,
            int blockSize, OmicsTransform transform,
            OmicsMissingPolicy missingPolicy, BlockScanner scanner,
            OmicsAssociationSink sink) throws IOException {
        int[] sampleOrder = SampleAlignment.requireOrder(
            source.metadata().sampleIds(), analysisSampleIds);
        long sourceFeatures = 0;
        long tested = 0;
        long failed = 0;
        try (NumericBlockReader reader = source.open(sampleOrder)) {
            for (NumericBlock block; (block = reader.read(blockSize)) != null;) {
                sourceFeatures += block.rows().size();
                List<String> names = block.rows().stream()
                    .map(NumericRow::id).toList();
                double[][] matrix = new double[analysisSampleIds.size()]
                    [block.rows().size()];
                for (int feature = 0; feature < block.rows().size(); feature++) {
                    double[] values = transformed(block.rows().get(feature),
                        transform, analysisSampleIds.size());
                    prepareMissing(values, missingPolicy,
                        block.rows().get(feature).id());
                    for (int sample = 0; sample < values.length; sample++)
                        matrix[sample][feature] = values[sample];
                }
                AssociationBatchResult result = scanner.scan(matrix, names);
                AssociationFailure[] failedRows = failureRows(result);
                for (int index = 0; index < result.size(); index++) {
                    AssociationEstimate estimate = result.estimate(index);
                    AssociationPipelineFailure failure = failure(estimate, failedRows[index]);
                    if (failure != null) {
                        sink.acceptFailure(failure);
                        failed++;
                        continue;
                    }
                    sink.acceptEstimate(new OmicsAssociationEstimate(
                        estimate.name(), estimate.beta(),
                        estimate.standardError(), estimate.statistic(),
                        estimate.degreesOfFreedom(), estimate.pValue(),
                        estimate.log10PValue(),
                        estimate.negativeLog10PValue()));
                    tested++;
                }
            }
        }
        return new OmicsAssociationSummary(sourceFeatures, tested, failed);
    }

    private static void prepareMissing(
            double[] values, OmicsMissingPolicy policy, String id) {
        double sum = 0;
        int finite = 0;
        for (double value : values)
            if (Double.isFinite(value)) {
                sum += value;
                finite++;
            }
        if (finite == values.length) return;
        if (finite == 0)
            throw new IllegalArgumentException(
                "omics feature has no finite values: " + id);
        if (policy == OmicsMissingPolicy.ERROR)
            throw new IllegalArgumentException(
                "omics feature contains missing values: " + id);
        double mean = sum / finite;
        for (int index = 0; index < values.length; index++)
            if (!Double.isFinite(values[index])) values[index] = mean;
    }

    private static double[] transformed(NumericRow row, OmicsTransform transform, int samples) {
        if (row.valuesView().length != samples)
            throw new IllegalArgumentException("feature length does not match aligned samples: " + row.id());
        double[] result = transform.apply(row.values());
        if (result == null || result.length != samples)
            throw new IllegalArgumentException("transform must preserve sample count: " + row.id());
        return result;
    }

    private static AssociationFailure[] failureRows(AssociationBatchResult result) {
        AssociationFailure[] failures = new AssociationFailure[result.size()];
        for (AssociationFailure failure : result.failures()) failures[failure.index()] = failure;
        return failures;
    }

    private static AssociationPipelineFailure failure(AssociationEstimate estimate, AssociationFailure failure) {
        if (failure != null) return new AssociationPipelineFailure(estimate.name(),
            failure.exceptionType(), failure.message());
        if (!Double.isFinite(estimate.beta()) || !Double.isFinite(estimate.standardError())
                || !Double.isFinite(estimate.pValue()))
            return new AssociationPipelineFailure(estimate.name(), "NonEstimablePredictor",
                "feature is not estimable after adjustment");
        return null;
    }

    private static void validate(
            NumericMatrixSource source, List<String> analysisSampleIds,
            int blockSize, OmicsTransform transform,
            OmicsMissingPolicy missingPolicy) {
        if (source == null || analysisSampleIds == null
                || analysisSampleIds.isEmpty() || blockSize < 1
                || transform == null || missingPolicy == null)
            throw new IllegalArgumentException(
                "omics pipeline inputs and positive block size are required");
    }

    @FunctionalInterface
    private interface BlockScanner {
        AssociationBatchResult scan(double[][] matrix, List<String> names);
    }
}
