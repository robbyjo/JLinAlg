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
import org.jlinalg.gwas.AssociationScanOptions;
import org.jlinalg.gwas.AssociationScanResult;
import org.jlinalg.gwas.RemlAssociationScanner;
import org.jlinalg.ols.OlsOptions;

/**
 * Bounded-memory single-trait GWAS scanner over any {@link VariantSource}.
 *
 * <p>Filtering is calculated after sample alignment and before mean imputation.
 * Results remain in source order. The fast OLS path reuses covariates within
 * each block; the REML path reuses one P3D/EMMAX null model across all
 * blocks.</p>
 */
public final class StreamingAssociationPipeline {
    private StreamingAssociationPipeline() { }

    public static AssociationPipelineResult fastOls(
            VariantSource source,
            List<String> analysisSampleIds,
            double[] response,
            double[][] covariates,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions,
            AssociationPipelineOptions pipelineOptions) throws IOException {
        validate(source, analysisSampleIds, response, covariates,
            pipelineOptions);
        return collect(source, analysisSampleIds, pipelineOptions,
            (variants, names) -> BlockResult.of(
                FastOlsAssociation.scanPredictors(
                    response, covariates, variants, names, weights, offset,
                    olsOptions, engineOptions)));
    }

    public static AssociationPipelineSummary fastOlsTo(
            VariantSource source,
            List<String> analysisSampleIds,
            double[] response,
            double[][] covariates,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions,
            AssociationPipelineOptions pipelineOptions,
            AssociationPipelineSink sink) throws IOException {
        validate(source, analysisSampleIds, response, covariates,
            pipelineOptions);
        return scanTo(source, analysisSampleIds, pipelineOptions,
            (variants, names) -> BlockResult.of(
                FastOlsAssociation.scanPredictors(
                    response, covariates, variants, names, weights, offset,
                    olsOptions, engineOptions)), sink);
    }

    /**
     * Fits one covariate-only GLM and streams efficient-score tests for all
     * variants. The prepared null model is reused across every input block.
     */
    public static AssociationPipelineSummary fastGlmTo(
            VariantSource source,
            List<String> analysisSampleIds,
            double[] response,
            double[][] covariates,
            GlmFamily family,
            double[] weights,
            double[] offset,
            GlmOptions glmOptions,
            AssociationEngineOptions engineOptions,
            AssociationPipelineOptions pipelineOptions,
            AssociationPipelineSink sink) throws IOException {
        validate(source, analysisSampleIds, response, covariates,
            pipelineOptions);
        FastGlmAssociation prepared = FastGlmAssociation.prepare(response,
            covariates, family, weights, offset, glmOptions, engineOptions);
        return scanTo(source, analysisSampleIds, pipelineOptions,
            (variants, names) -> BlockResult.of(
                prepared.scan(variants, names, engineOptions)), sink);
    }

    public static AssociationPipelineResult remlP3d(
            VariantSource source,
            List<String> analysisSampleIds,
            RemlAssociationScanner scanner,
            AssociationScanOptions scanOptions,
            AssociationPipelineOptions pipelineOptions) throws IOException {
        if (scanner == null || scanOptions == null)
            throw new IllegalArgumentException("REML scanner and options are required");
        if (source == null || analysisSampleIds == null
                || analysisSampleIds.isEmpty() || pipelineOptions == null)
            throw new IllegalArgumentException(
                "source, analysis samples, and pipeline options are required");
        return collect(source, analysisSampleIds, pipelineOptions,
            (variants, names) -> {
                AssociationScanResult result =
                    scanner.scan(variants, names, scanOptions);
                return BlockResult.of(result);
            });
    }

    public static AssociationPipelineSummary remlP3dTo(
            VariantSource source,
            List<String> analysisSampleIds,
            RemlAssociationScanner scanner,
            AssociationScanOptions scanOptions,
            AssociationPipelineOptions pipelineOptions,
            AssociationPipelineSink sink) throws IOException {
        if (scanner == null || scanOptions == null)
            throw new IllegalArgumentException("REML scanner and options are required");
        if (source == null || analysisSampleIds == null
                || analysisSampleIds.isEmpty() || pipelineOptions == null)
            throw new IllegalArgumentException(
                "source, analysis samples, and pipeline options are required");
        return scanTo(source, analysisSampleIds, pipelineOptions,
            (variants, names) -> BlockResult.of(
                scanner.scan(variants, names, scanOptions)), sink);
    }

    private static AssociationPipelineResult collect(
            VariantSource source,
            List<String> analysisSampleIds,
            AssociationPipelineOptions options,
            BlockScanner scanner) throws IOException {
        List<AssociationPipelineEstimate> estimates = new ArrayList<>();
        List<VariantFilterResult> excluded = new ArrayList<>();
        List<AssociationPipelineFailure> failures = new ArrayList<>();
        AssociationPipelineSummary summary = scanTo(source, analysisSampleIds,
            options, scanner, new AssociationPipelineSink() {
                @Override public void acceptEstimate(
                        AssociationPipelineEstimate estimate) {
                    estimates.add(estimate);
                }
                @Override public void acceptExcluded(
                        VariantFilterResult result) {
                    excluded.add(result);
                }
                @Override public void acceptFailure(
                        AssociationPipelineFailure failure) {
                    failures.add(failure);
                }
            });
        return new AssociationPipelineResult(summary.sourceVariants(), estimates,
            excluded, failures);
    }

    private static AssociationPipelineSummary scanTo(
            VariantSource source,
            List<String> analysisSampleIds,
            AssociationPipelineOptions options,
            BlockScanner scanner,
            AssociationPipelineSink sink) throws IOException {
        if (sink == null)
            throw new IllegalArgumentException("association sink is required");
        int[] sampleOrder = SampleAlignment.requireOrder(
            source.metadata().sampleIds(), analysisSampleIds);
        long sourceVariants = 0;
        long tested = 0;
        long excludedCount = 0;
        long failureCount = 0;
        try (VariantBlockReader reader = source.open(sampleOrder)) {
            for (VariantBlock block;
                    (block = reader.read(options.variantBlockSize())) != null;) {
                sourceVariants += block.variants().size();
                List<VariantFilterResult> included = new ArrayList<>();
                for (VariantRecord variant : block.variants()) {
                    VariantFilterResult result = VariantFilters.evaluate(
                        variant, options.variantFilter());
                    if (result.included()) included.add(result);
                    else {
                        sink.acceptExcluded(result);
                        excludedCount++;
                    }
                }
                if (included.isEmpty()) continue;
                double[][] matrix = transpose(included, analysisSampleIds.size());
                List<String> names = included.stream()
                    .map(result -> result.variant().id()).toList();
                BlockResult blockResult = scanner.scan(matrix, names);
                for (int index = 0; index < included.size(); index++) {
                    int estimateIndex = index;
                    VariantFilterResult filtered = included.get(index);
                    AssociationEstimate estimate = blockResult.estimate(index);
                    AssociationPipelineEstimate output =
                        new AssociationPipelineEstimate(
                        filtered.variant(), filtered.statistics(),
                        estimate.beta(), estimate.standardError(),
                        estimate.statistic(), estimate.degreesOfFreedom(),
                        estimate.pValue(), estimate.log10PValue(),
                        estimate.negativeLog10PValue());
                    sink.acceptEstimate(output);
                    tested++;
                    if (!Double.isFinite(estimate.beta())
                            && blockResult.failures().stream()
                                .noneMatch(value ->
                                    value.index() == estimateIndex)) {
                        sink.acceptFailure(new AssociationPipelineFailure(
                            filtered.variant().id(), "NonEstimablePredictor",
                            "variant is constant or collinear after adjustment"));
                        failureCount++;
                    }
                }
                for (AssociationFailure failure : blockResult.failures()) {
                    sink.acceptFailure(new AssociationPipelineFailure(
                        included.get(failure.index()).variant().id(),
                        failure.exceptionType(), failure.message()));
                    failureCount++;
                }
            }
        }
        return new AssociationPipelineSummary(sourceVariants, tested,
            excludedCount, failureCount);
    }

    private static double[][] transpose(
            List<VariantFilterResult> variants, int samples) {
        double[][] result = new double[samples][variants.size()];
        for (int variant = 0; variant < variants.size(); variant++) {
            double[] dosages = variants.get(variant).variant().dosagesView();
            if (dosages.length != samples)
                throw new IllegalArgumentException(
                    "variant dosage count does not match aligned samples");
            for (int sample = 0; sample < samples; sample++)
                result[sample][variant] = dosages[sample];
        }
        return result;
    }

    private static void validate(
            VariantSource source, List<String> analysisSampleIds,
            double[] response, double[][] covariates,
            AssociationPipelineOptions options) {
        if (source == null || analysisSampleIds == null
                || analysisSampleIds.isEmpty() || response == null
                || covariates == null || options == null)
            throw new IllegalArgumentException("pipeline inputs are required");
        if (response.length != analysisSampleIds.size()
                || covariates.length != analysisSampleIds.size())
            throw new IllegalArgumentException(
                "response and covariate rows must match analysis samples");
    }

    @FunctionalInterface
    private interface BlockScanner {
        BlockResult scan(double[][] variants, List<String> names);
    }

    private record BlockResult(
            List<AssociationEstimate> estimates,
            List<AssociationFailure> failures) {
        private AssociationEstimate estimate(int index) {
            return estimates.get(index);
        }
        private static BlockResult of(AssociationBatchResult result) {
            List<AssociationEstimate> estimates =
                new ArrayList<>(result.names().size());
            for (int index = 0; index < result.names().size(); index++)
                estimates.add(result.estimate(index));
            return new BlockResult(estimates, result.failures());
        }
        private static BlockResult of(AssociationScanResult result) {
            List<AssociationEstimate> estimates =
                new ArrayList<>(result.markerNames().size());
            double[] beta = result.beta();
            double[] standardErrors = result.standardErrors();
            double[] statistics = result.tStatistics();
            double[] p = result.pValues();
            double[] logP = result.log10PValues();
            double[] minusLogP = result.negativeLog10PValues();
            double[] degrees = result.statistics().degreesOfFreedom();
            for (int index = 0; index < beta.length; index++) {
                estimates.add(new AssociationEstimate(
                    result.markerNames().get(index), beta[index],
                    standardErrors[index], statistics[index], p[index],
                    logP[index], minusLogP[index], degrees[index]));
            }
            return new BlockResult(estimates, List.of());
        }
    }
}
