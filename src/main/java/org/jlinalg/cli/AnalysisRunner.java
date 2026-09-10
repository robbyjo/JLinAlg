/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.jlinalg.association.AssociationFitter;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.association.AssociationFailurePolicy;
import org.jlinalg.association.VariableMissingPolicy;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.CompiledFormula;
import org.jlinalg.formula.CompiledMixedFormula;
import org.jlinalg.formula.Formula;
import org.jlinalg.formula.MixedFormula;
import org.jlinalg.formula.MixedFormulaOptions;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.glmm.GlmmLaplace;
import org.jlinalg.glmm.GlmmLaplaceOptions;
import org.jlinalg.glmm.GlmmLaplaceResult;
import org.jlinalg.glmm.SparseGlmmLaplace;
import org.jlinalg.genetics.GenomicRelationshipMatrix;
import org.jlinalg.gwas.AssociationScanOptions;
import org.jlinalg.gwas.GenotypeMissingPolicy;
import org.jlinalg.gwas.RemlAssociationScanner;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.OlsResult;
import org.jlinalg.pipeline.AssociationPipelineOptions;
import org.jlinalg.pipeline.AssociationPipelineSummary;
import org.jlinalg.pipeline.DataFormat;
import org.jlinalg.pipeline.DelimitedMatrixSource;
import org.jlinalg.pipeline.NumericMatrixSource;
import org.jlinalg.pipeline.OmicsAssociationSummary;
import org.jlinalg.pipeline.OmicsMissingPolicy;
import org.jlinalg.pipeline.OmicsTransform;
import org.jlinalg.pipeline.StreamingAssociationPipeline;
import org.jlinalg.pipeline.StreamingOmicsAssociationPipeline;
import org.jlinalg.pipeline.VariantFilterOptions;
import org.jlinalg.pipeline.VariantSource;
import org.jlinalg.pipeline.VariantSources;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.reml.Reml;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.RemlResult;
import org.jlinalg.reml.VarianceComponent;
import org.jlinalg.survival.CoxKinshipFrailty;
import org.jlinalg.survival.CoxMixedOptions;
import org.jlinalg.survival.CoxMixedResult;
import org.jlinalg.survival.CoxOptions;
import org.jlinalg.survival.CoxRegression;
import org.jlinalg.survival.CoxResult;
import org.jlinalg.survival.CoxSurvivalData;
import org.jlinalg.survival.CoxTies;

/** Resolves CLI configuration into library model and streaming pipeline calls. */
final class AnalysisRunner {
    private final CliOptions options;
    private final FormulaPlan plan;
    private final RunLog log;
    private final PrintStream output;
    private String resolvedVarianceComponents = "not-applicable";
    private String resolvedMixedFit = "not-applicable";
    private String resolvedStatisticType = "not-applicable";
    private String resolvedDfMethod = "not-applicable";
    private String resolvedPartialR2Method = "not-applicable";

    AnalysisRunner(CliOptions options, FormulaPlan plan, RunLog log,
            PrintStream output) {
        this.options = options;
        this.plan = plan;
        this.log = log;
        this.output = output;
    }

    int execute() throws IOException {
        return options.omics == null ? phenotypeOnly() : omics();
    }

    private int omics() throws IOException {
        OmicsTypeDetector.Detection detection = OmicsTypeDetector.detect(
            options.omics, options.omicsType);
        info("omics_type=" + detection.type());
        info("omics_type_source=" + detection.source());
        info("omics_type_confidence=" + detection.confidence());
        if (options.annotation != null && detection.type().equals("gwas"))
            throw new IllegalArgumentException(
                "--annot is intended for non-WGS/GWAS omics");
        AnnotationLookup annotation = options.annotation == null
            ? AnnotationLookup.empty()
            : AnnotationLookup.read(options.annotation, options.annotationId,
                options.annotationColumns);
        info("annotation_rows=" + annotation.size());

        boolean genotype = detection.type().equals("gwas");
        List<String> sourceIds;
        VariantSource variantSource = null;
        NumericMatrixSource numericSource = null;
        if (genotype) {
            DataFormat format = DataFormat.infer(options.omics);
            List<String> external = format == DataFormat.BGEN
                    && options.bgenSamples != null
                ? SampleFiles.read(options.bgenSamples) : null;
            variantSource = VariantSources.open(options.omics, format, external);
            sourceIds = variantSource.metadata().sampleIds();
        } else {
            numericSource = DelimitedMatrixSource.open(options.omics);
            sourceIds = numericSource.metadata().sampleIds();
        }

        boolean binomial = options.family.equals("binomial");
        PhenotypeData phenotype = PhenotypeData.read(
            options.phenotype, options.idColumn);
        PhenotypeData.Prepared prepared = phenotype.prepare(sourceIds,
            plan.response(), binomial, options.caseValue, options.controlValue);
        int idAlignedSamples = prepared.ids().size();
        prepared = completeCases(
            phenotype, prepared, plan.response(), binomial);
        reportSampleAlignment(sourceIds.size(), phenotype.originalIds().size(),
            idAlignedSamples, prepared.ids().size());
        logBinary(prepared);
        GrmContext grm = grm(phenotype, prepared);
        PedigreeContext pedigree = pedigree(phenotype, prepared);
        String model = resolveModel();
        CompiledFormula fixed;
        CompiledMixedFormula mixed = null;
        if (plan.hasRandomEffects()) {
            mixed = compileMixed(
                plan.withoutOmics(), prepared.modelTable(), pedigree);
            fixed = mixed.fixed();
        } else {
            fixed = Formula.compile(
                plan.withoutOmics(), prepared.modelTable());
        }
        double[][] covariates = matrix(
            fixed.design(), fixed.rows(), fixed.columns());
        int blockSize = AdaptiveBlockSizer.choose(prepared.ids().size(),
            options.blockSize);
        resolvedVarianceComponents =
            resolveVarianceComponents(model, genotype);
        resolvedMixedFit = resolveMixedFit(
            model, genotype, resolvedVarianceComponents);
        resolveOutputMetadata(model, genotype);
        info("resolved_model=" + model);
        info("variance_components=" + resolvedVarianceComponents);
        info("mixed_fit=" + resolvedMixedFit);
        logOutputMetadata();
        info("block_size=" + blockSize);
        info("threads=" + options.threads);
        info("backend=" + options.backend);
        if (options.explain || options.dryRun) {
            output.println("omics type: " + detection.type()
                + " (" + detection.source() + ")");
            output.println("model: " + model);
            output.println("variance components: "
                + resolvedVarianceComponents);
            output.println("mixed fit: " + resolvedMixedFit);
            output.println("output format: " + outputFormat());
            output.println("adaptive block size: " + blockSize);
            output.println("backend: " + options.backend);
            output.println("output: " + options.output);
            if (options.dryRun) return 0;
        }
        if (options.resume && Files.exists(Path.of(options.output + ".partial")))
            throw new IllegalArgumentException(
                "partial output lacks resumable block metadata; "
                    + "use --overwrite to restart");
        OmicsTransform transform = TransformParser.parse(
            options.transforms, options.transformPlugins);
        if (genotype && !options.transforms.isEmpty())
            throw new IllegalArgumentException(
                "genotype transforms are not supported in the variant pipeline");
        AssociationEngineOptions engine = new AssociationEngineOptions(
            options.threads, Math.max(1, Math.min(256, blockSize)),
            options.backend, AssociationFailurePolicy.RECORD_NAN,
            VariableMissingPolicy.MEAN_IMPUTE);
        VariantFilterOptions filters = VariantFilterOptions.builder()
            .minimumMaf(options.minimumMaf)
            .minimumMac(options.minimumMac)
            .maximumMissingRate(options.maximumMissingRate)
            .minimumImputationQuality(options.minimumInfo)
            .build();
        AssociationPipelineOptions pipeline =
            new AssociationPipelineOptions(blockSize, filters);
        Counts counts;
        try (CliResultSink sink = new CliResultSink(options.output,
                options.overwrite, genotype, resolvedStatisticType,
                annotation, genotype ? prepared.caseControlGroups() : null)) {
            counts = switch (model) {
                case "ols" -> scanOls(variantSource, numericSource, genotype,
                    prepared.ids(), fixed, covariates, transform, blockSize,
                    engine, pipeline, sink);
                case "glm" -> scanGlm(variantSource, numericSource, genotype,
                    prepared.ids(), fixed, covariates, transform, blockSize,
                    engine, pipeline, sink);
                case "lmm" -> scanLmm(variantSource, numericSource, genotype,
                    prepared.ids(), fixed, mixed, grm, pedigree, covariates,
                    transform, blockSize, engine, pipeline, sink);
                case "glmm" -> scanGlmm(variantSource, numericSource, genotype,
                    prepared.ids(), fixed, mixed, grm, pedigree, covariates,
                    transform, blockSize, engine, sink);
                default -> throw new IllegalArgumentException(
                    "omics model is not yet supported: " + model);
            };
            sink.finish();
            info("fdr_tests=" + sink.adjustedTests());
        }
        info("source_features=" + counts.source());
        info("tested_features=" + counts.tested());
        info("failed_features=" + counts.failed());
        manifest(model, detection.type(), blockSize, counts);
        return 0;
    }

    private Counts scanOls(
            VariantSource variant, NumericMatrixSource numeric,
            boolean genotype, List<String> ids, CompiledFormula fixed,
            double[][] covariates, OmicsTransform transform, int blockSize,
            AssociationEngineOptions engine,
            AssociationPipelineOptions pipeline, CliResultSink sink)
            throws IOException {
        if (genotype) {
            AssociationPipelineSummary summary =
                StreamingAssociationPipeline.fastOlsTo(
                    variant, ids, fixed.response(), covariates,
                    fixed.weights(), fixed.offset(), OlsOptions.defaults(),
                    engine, pipeline, sink);
            return Counts.of(summary);
        }
        OmicsAssociationSummary summary =
            StreamingOmicsAssociationPipeline.scanPredictorsTo(
                numeric, ids, fixed.response(), covariates, transform,
                OmicsMissingPolicy.MEAN_IMPUTE, blockSize,
                fixed.weights(), fixed.offset(), OlsOptions.defaults(),
                engine, sink);
        return Counts.of(summary);
    }

    private Counts scanGlm(
            VariantSource variant, NumericMatrixSource numeric,
            boolean genotype, List<String> ids, CompiledFormula fixed,
            double[][] covariates, OmicsTransform transform, int blockSize,
            AssociationEngineOptions engine,
            AssociationPipelineOptions pipeline, CliResultSink sink)
            throws IOException {
        GlmFamily family = family(options.family);
        if (genotype) {
            AssociationPipelineSummary summary =
                StreamingAssociationPipeline.fastGlmTo(
                    variant, ids, fixed.response(), covariates, family,
                    fixed.weights(), fixed.offset(), GlmOptions.defaults(),
                    engine, pipeline, sink);
            return Counts.of(summary);
        }
        OmicsAssociationSummary summary =
            StreamingOmicsAssociationPipeline.scanPredictorsGlmTo(
                numeric, ids, fixed.response(), covariates, family,
                transform, OmicsMissingPolicy.MEAN_IMPUTE, blockSize,
                fixed.weights(), fixed.offset(), GlmOptions.defaults(),
                engine, sink);
        return Counts.of(summary);
    }

    private Counts scanLmm(
            VariantSource variant, NumericMatrixSource numeric,
            boolean genotype, List<String> ids,
            CompiledFormula fixed, CompiledMixedFormula mixed,
            GrmContext grm, PedigreeContext pedigree,
            double[][] covariates, OmicsTransform transform, int blockSize,
            AssociationEngineOptions engine,
            AssociationPipelineOptions pipeline, CliResultSink sink)
            throws IOException {
        if ((mixed == null && grm == null)
                || mixed != null
                    && !mixed.correlatedRandomEffects().isEmpty())
            throw new IllegalArgumentException(
                "omics LMM requires a GRM or independent random-effect terms");
        if (!genotype) {
            requireRefit("numeric LMM");
            return scanExactReml(numeric, ids, fixed, mixed, grm, pedigree,
                covariates, transform, blockSize, engine, sink);
        }
        if (!resolvedVarianceComponents.equals("null-model"))
            throw new IllegalArgumentException(
                "genotype LMM per-marker refitting is not yet available; "
                    + "use --variance-components null-model");
        if (pedigree != null)
            throw new IllegalArgumentException(
                "pedigree genotype scans require per-marker refitting");
        if (!options.degreesOfFreedom.equals("auto"))
            throw new IllegalArgumentException(
                "Satterthwaite/KR marker tests require refit mode; "
                    + "null-model scans use residual approximation");
        RemlAssociationScanner scanner = RemlAssociationScanner.prepare(
            fixed.response(), covariates,
            components(mixed == null ? List.of() : mixed.randomEffects(),
                fixed.rows(), grm == null ? null : grm.component()),
            RemlOptions.defaults(), options.backend);
        AssociationPipelineSummary summary =
            StreamingAssociationPipeline.remlP3dTo(
                variant, ids, scanner,
                new AssociationScanOptions(
                    Math.max(1, Math.min(256, blockSize)),
                    GenotypeMissingPolicy.MEAN_IMPUTE, options.threads),
                pipeline, sink);
        return Counts.of(summary);
    }

    private Counts scanExactReml(
            NumericMatrixSource numeric, List<String> ids,
            CompiledFormula fixed, CompiledMixedFormula mixed,
            GrmContext grm, PedigreeContext pedigree,
            double[][] covariates, OmicsTransform transform, int blockSize,
            AssociationEngineOptions engine, CliResultSink sink)
            throws IOException {
        if (grm != null && pedigree != null)
            throw new IllegalArgumentException(
                "a pedigree and --grm cannot be combined in one CLI model");
        if (fixed.weights() != null)
            throw new IllegalArgumentException(
                "weighted streamed LMM refits are not yet supported");
        double[] response = adjustedResponse(fixed);
        RemlOptions reml = RemlOptions.builder()
            .degreesOfFreedomMethod(
                dfMethod(options.degreesOfFreedom)).build();
        if (grm != null) {
            AssociationFitter fitter = (values, design, rows, columns, backend) -> {
                RemlResult fit = Reml.fit(values, design, rows, columns,
                    components(mixed == null ? List.of() : mixed.randomEffects(),
                        rows, grm.component()),
                    reml, backend);
                requireConverged(fit.converged(), "REML");
                return fit.associationStatistics();
            };
            return exactOmicsScan(numeric, ids, response, covariates,
                transform, blockSize, fitter, engine, sink);
        }
        if (mixed == null)
            throw new IllegalArgumentException(
                "numeric LMM requires a random-effect formula term");
        SparseRandomStructure structure =
            sparseStructure(mixed, pedigree);
        try (SparseLinearMixedModel.Prepared prepared =
                SparseLinearMixedModel.prepareWithPrecision(
                    fixed.rows(), structure.terms(), structure.precisions(),
                    reml, options.backend)) {
            AssociationFitter fitter = (values, design, rows, columns, backend) -> {
                SparseLinearMixedModelResult fit =
                    prepared.fit(values, design, columns);
                requireConverged(fit.converged(), "REML");
                return fit.associationStatistics();
            };
            return exactOmicsScan(numeric, ids, response, covariates,
                transform, blockSize, fitter, engine, sink);
        }
    }

    private Counts scanGlmm(
            VariantSource variant, NumericMatrixSource numeric,
            boolean genotype, List<String> ids,
            CompiledFormula fixed, CompiledMixedFormula mixed,
            GrmContext grm, PedigreeContext pedigree,
            double[][] covariates, OmicsTransform transform, int blockSize,
            AssociationEngineOptions engine, CliResultSink sink)
            throws IOException {
        if (genotype)
            throw new IllegalArgumentException(
                "genotype Laplace GLMM scans are not yet available");
        requireRefit("numeric GLMM");
        if ((mixed == null && grm == null)
                || mixed != null
                    && !mixed.correlatedRandomEffects().isEmpty())
            throw new IllegalArgumentException(
                "omics GLMM requires a GRM or independent random-effect terms");
        if (grm != null && pedigree != null)
            throw new IllegalArgumentException(
                "a pedigree and --grm cannot be combined in one CLI model");
        GlmFamily family = family(options.family);
        GlmmLaplaceOptions laplace = GlmmLaplaceOptions.defaults();
        AssociationFitter fitter;
        if (grm != null) {
            fitter = (values, design, rows, columns, backend) -> {
                GlmmLaplaceResult fit = GlmmLaplace.fit(
                    values, design, rows, columns, family,
                    randomComponents(
                        mixed == null ? List.of() : mixed.randomEffects(),
                        rows, grm.component()),
                    fixed.weights(), fixed.offset(), laplace, backend);
                requireConverged(fit.converged(), "Laplace GLMM");
                return fit.associationStatistics();
            };
            return exactOmicsScan(numeric, ids, fixed.response(), covariates,
                transform, blockSize, fitter, engine, sink);
        }
        SparseRandomStructure structure =
            sparseStructure(mixed, pedigree);
        try (SparseGlmmLaplace.Prepared prepared =
                SparseGlmmLaplace.prepareWithPrecision(
                    fixed.rows(), family, structure.terms(),
                    structure.precisions(), laplace, options.backend)) {
            fitter = (values, design, rows, columns, backend) -> {
                GlmmLaplaceResult fit = prepared.fit(
                    values, design, columns, fixed.weights(), fixed.offset());
                requireConverged(fit.converged(), "Laplace GLMM");
                return fit.associationStatistics();
            };
            return exactOmicsScan(numeric, ids, fixed.response(), covariates,
                transform, blockSize, fitter, engine, sink);
        }
    }

    private Counts exactOmicsScan(
            NumericMatrixSource numeric, List<String> ids,
            double[] response, double[][] covariates,
            OmicsTransform transform, int blockSize,
            AssociationFitter fitter, AssociationEngineOptions engine,
            CliResultSink sink) throws IOException {
        OmicsAssociationSummary summary =
            StreamingOmicsAssociationPipeline.scanPredictorsRefitTo(
                numeric, ids, response, covariates, transform,
                OmicsMissingPolicy.MEAN_IMPUTE, blockSize,
                fitter, engine, sink);
        return Counts.of(summary);
    }

    private int phenotypeOnly() throws IOException {
        String model = resolveModel();
        String preparedResponse = plan.isCox()
            ? plan.survival().event() : plan.response();
        PhenotypeData phenotype = PhenotypeData.read(
            options.phenotype, options.idColumn);
        PhenotypeData.Prepared prepared = phenotype.prepare(null,
            preparedResponse, options.family.equals("binomial") || plan.isCox(),
            options.caseValue, options.controlValue);
        int phenotypeSamples = prepared.ids().size();
        prepared = completeCases(phenotype, prepared, preparedResponse,
            options.family.equals("binomial") || plan.isCox());
        reportPhenotypeSamples(phenotypeSamples, prepared.ids().size());
        GrmContext grm = grm(phenotype, prepared);
        PedigreeContext pedigree = pedigree(phenotype, prepared);
        resolvedVarianceComponents =
            resolveVarianceComponents(model, false);
        resolvedMixedFit = resolveMixedFit(
            model, false, resolvedVarianceComponents);
        resolveOutputMetadata(model, false);
        info("omics_type=none");
        info("resolved_model=" + model);
        info("variance_components=" + resolvedVarianceComponents);
        info("mixed_fit=" + resolvedMixedFit);
        logOutputMetadata();
        long tests = switch (model) {
            case "ols" -> phenotypeOls(prepared);
            case "glm" -> phenotypeGlm(prepared);
            case "lmm" -> phenotypeLmm(prepared, grm, pedigree);
            case "glmm" -> phenotypeGlmm(prepared, grm, pedigree);
            case "cox" -> phenotypeCox(prepared, grm);
            default -> throw new IllegalArgumentException(
                "unsupported model: " + model);
        };
        info("fdr_tests=" + tests);
        manifest(model, "none", 0,
            new Counts(prepared.ids().size(), tests, 0));
        return 0;
    }

    private long phenotypeOls(PhenotypeData.Prepared prepared)
            throws IOException {
        CompiledFormula compiled = Formula.compile(
            plan.withoutOmics(), prepared.modelTable());
        OlsResult fit = compiled.fitOls(
            OlsOptions.defaults(), options.backend);
        return CoefficientOutput.write(options.output, options.overwrite,
            compiled.coefficientNames(), fit.beta(), fit.standardErrors(),
            fit.tStatistics(), filled(compiled.columns(),
                fit.residualDegreesOfFreedom()), fit.pValues(),
            "t", null, "transformed_effect");
    }

    private long phenotypeGlm(PhenotypeData.Prepared prepared)
            throws IOException {
        CompiledFormula compiled = Formula.compile(
            plan.withoutOmics(), prepared.modelTable());
        GlmResult fit = compiled.fitGlm(family(options.family),
            GlmOptions.defaults(), options.backend);
        double df = fit.observations() - fit.rank();
        AssociationStatistics statistics = AssociationStatistics.studentT(
            fit.beta(), fit.standardErrors(), df,
            DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
        return CoefficientOutput.write(options.output, options.overwrite,
            compiled.coefficientNames(), statistics.beta(),
            statistics.standardErrors(), statistics.statistics(),
            statistics.degreesOfFreedom(), statistics.pValues(),
            "t_approx", null, "transformed_effect");
    }

    private long phenotypeLmm(
            PhenotypeData.Prepared prepared, GrmContext grm,
            PedigreeContext pedigree)
            throws IOException {
        CompiledFormula fixed;
        CompiledMixedFormula mixed = null;
        if (plan.hasRandomEffects()) {
            mixed = compileMixed(
                plan.withoutOmics(), prepared.modelTable(), pedigree);
            fixed = mixed.fixed();
        } else {
            fixed = Formula.compile(
                plan.withoutOmics(), prepared.modelTable());
        }
        if (grm == null && mixed == null)
            throw new IllegalArgumentException(
                "LMM requires --grm or a random-effect formula term");
        requireRefit("LMM");
        RemlOptions reml = RemlOptions.builder()
            .degreesOfFreedomMethod(dfMethod(options.degreesOfFreedom)).build();
        AssociationStatistics statistics;
        if (grm == null) {
            if (pedigree == null) {
                statistics = mixed.fitSparse(reml, options.backend)
                    .associationStatistics();
            } else {
                if (!mixed.correlatedRandomEffects().isEmpty()
                        || fixed.weights() != null)
                    throw new IllegalArgumentException(
                        "pedigree REML currently requires independent random "
                            + "terms and unweighted observations");
                SparseRandomStructure structure =
                    sparseStructure(mixed, pedigree);
                SparseLinearMixedModelResult fit =
                    SparseLinearMixedModel.fitWithPrecision(
                        adjustedResponse(fixed), fixed.design(),
                        fixed.rows(), fixed.columns(), structure.terms(),
                        structure.precisions(), reml, options.backend);
                requireConverged(fit.converged(), "REML");
                statistics = fit.associationStatistics();
            }
        } else {
            if (!mixedCompatibleWithDenseKinship(mixed, pedigree))
                throw new IllegalArgumentException(
                    "--grm cannot currently be combined with pedigree or "
                        + "correlated random-slope terms");
            if (fixed.weights() != null)
                throw new IllegalArgumentException(
                    "weighted GRM LMM is not supported by this CLI path");
            RemlResult fit = Reml.fit(adjustedResponse(fixed), fixed.design(),
                fixed.rows(), fixed.columns(),
                components(mixed == null ? List.of() : mixed.randomEffects(),
                    fixed.rows(), grm.component()),
                reml, options.backend);
            requireConverged(fit.converged(), "REML");
            statistics = fit.associationStatistics();
        }
        return CoefficientOutput.write(options.output, options.overwrite,
            fixed.coefficientNames(), statistics.beta(),
            statistics.standardErrors(), statistics.statistics(),
            statistics.degreesOfFreedom(), statistics.pValues(), "t",
            null, "transformed_effect");
    }

    private long phenotypeGlmm(
            PhenotypeData.Prepared prepared, GrmContext grm,
            PedigreeContext pedigree)
            throws IOException {
        CompiledFormula fixed;
        CompiledMixedFormula mixed = null;
        if (plan.hasRandomEffects()) {
            mixed = compileMixed(
                plan.withoutOmics(), prepared.modelTable(), pedigree);
            if (!mixed.correlatedRandomEffects().isEmpty())
                throw new IllegalArgumentException(
                    "correlated random slopes are not yet supported by "
                        + "the Laplace GLMM CLI");
            fixed = mixed.fixed();
        } else {
            fixed = Formula.compile(
                plan.withoutOmics(), prepared.modelTable());
        }
        if (grm == null && mixed == null)
            throw new IllegalArgumentException(
                "GLMM requires --grm or a random-effect formula term");
        requireRefit("GLMM");
        GlmmLaplaceResult fit;
        if (grm == null) {
            SparseRandomStructure structure =
                sparseStructure(mixed, pedigree);
            fit = SparseGlmmLaplace.fitWithPrecision(
                fixed.response(), fixed.design(), fixed.rows(), fixed.columns(),
                family(options.family), structure.terms(),
                structure.precisions(), fixed.weights(), fixed.offset(),
                GlmmLaplaceOptions.defaults(), options.backend);
        } else {
            if (!mixedCompatibleWithDenseKinship(mixed, pedigree))
                throw new IllegalArgumentException(
                    "--grm cannot currently be combined with pedigree or "
                        + "correlated random-slope terms");
            fit = GlmmLaplace.fit(
                fixed.response(), fixed.design(), fixed.rows(), fixed.columns(),
                family(options.family),
                randomComponents(
                    mixed == null ? List.of() : mixed.randomEffects(),
                    fixed.rows(), grm.component()),
                fixed.weights(), fixed.offset(),
                GlmmLaplaceOptions.defaults(), options.backend);
        }
        requireConverged(fit.converged(), "Laplace GLMM");
        return CoefficientOutput.write(options.output, options.overwrite,
            fixed.coefficientNames(), fit.beta(),
            fit.standardErrors(), fit.statistics(),
            filled(fixed.columns(), Double.POSITIVE_INFINITY),
            fit.pValues(), "z", null, "transformed_effect");
    }

    private long phenotypeCox(
            PhenotypeData.Prepared prepared, GrmContext grm)
            throws IOException {
        if (plan.hasRandomEffects())
            throw new IllegalArgumentException(
                "Cox frailty terms are not yet exposed by the CLI");
        String rhs = plan.withoutOmics().substring(
            plan.withoutOmics().indexOf('~') + 1);
        FormulaPlan.Survival survival = plan.survival();
        CompiledFormula design = Formula.compile(
            survival.stop() + "~0+" + rhs, prepared.modelTable());
        double[] stop = design.response();
        double[] start = survival.start() == null ? new double[stop.length]
            : Formula.compile(survival.start() + "~1",
                prepared.modelTable()).response();
        double[] eventValues = Formula.compile(
            survival.event() + "~1", prepared.modelTable()).response();
        boolean[] event = new boolean[eventValues.length];
        for (int index = 0; index < event.length; index++)
            event[index] = eventValues[index] == 1.0;
        CoxTies ties = CoxTies.valueOf(
            options.ties.toUpperCase(Locale.ROOT));
        CoxSurvivalData survivalData =
            new CoxSurvivalData(start, stop, event, null);
        double[] beta;
        double[] standardErrors;
        double[] statistics;
        double[] pValues;
        double[] hazardRatios;
        boolean converged;
        String convergenceMessage;
        if (grm == null) {
            CoxResult fit = CoxRegression.fit(
                survivalData,
                matrix(design.design(), design.rows(), design.columns()),
                design.offset(), CoxOptions.defaults().withTies(ties),
                options.backend);
            beta = fit.beta();
            standardErrors = fit.standardErrors();
            statistics = fit.zStatistics();
            pValues = fit.pValues();
            hazardRatios = fit.hazardRatios();
            converged = fit.converged();
            convergenceMessage = fit.convergenceMessage();
        } else {
            CoxMixedOptions defaults = CoxMixedOptions.defaults();
            CoxMixedOptions mixedOptions = new CoxMixedOptions(
                CoxOptions.defaults().withTies(ties),
                defaults.initialVariances(),
                defaults.maximumVarianceIterations(),
                defaults.logVarianceTolerance(), defaults.minimumVariance(),
                defaults.maximumVariance());
            CoxMixedResult fit = CoxKinshipFrailty.fit(
                survivalData,
                matrix(design.design(), design.rows(), design.columns()),
                grm.observationIds(), grm.matrix(), design.offset(),
                mixedOptions, 1e-8, options.backend);
            beta = fit.beta();
            standardErrors = fit.standardErrors();
            statistics = fit.zStatistics();
            pValues = fit.pValues();
            hazardRatios = fit.hazardRatios();
            converged = fit.converged();
            convergenceMessage = fit.convergenceMessage();
        }
        if (!converged)
            warning("cox_convergence=" + convergenceMessage);
        return CoefficientOutput.write(options.output, options.overwrite,
            design.coefficientNames(), beta, standardErrors,
            statistics, filled(design.columns(), Double.POSITIVE_INFINITY),
            pValues, "z", hazardRatios, "hazard_ratio");
    }

    private String resolveModel() {
        String resolved;
        if (!options.model.equals("auto")) resolved = options.model;
        else if (plan.isCox()) resolved = "cox";
        else if (plan.hasRandomEffects() || options.grm != null
                || options.pedigree != null)
            resolved = options.family.equals("gaussian") ? "lmm" : "glmm";
        else resolved = options.family.equals("gaussian") ? "ols" : "glm";
        if (options.grm != null && !resolved.equals("lmm")
                && !resolved.equals("glmm") && !resolved.equals("cox"))
            throw new IllegalArgumentException(
                "--grm is applicable to lmm, glmm, and Cox models");
        return resolved;
    }

    private String resolveVarianceComponents(
            String model, boolean genotype) {
        if (!model.equals("lmm") && !model.equals("glmm"))
            return "not-applicable";
        if (!options.varianceComponents.equals("auto"))
            return options.varianceComponents;
        return genotype && model.equals("lmm") ? "null-model" : "refit";
    }

    private void requireRefit(String model) {
        if (!resolvedVarianceComponents.equals("refit"))
            throw new IllegalArgumentException(
                model + " uses per-feature variance-component refitting; "
                    + "use --variance-components refit or auto");
    }

    private static String resolveMixedFit(
            String model, boolean genotype, String varianceComponents) {
        if (model.equals("lmm"))
            return genotype && varianceComponents.equals("null-model")
                ? "p3d-null-model" : "exact-reml-refit";
        if (model.equals("glmm")) return "laplace-marginal-refit";
        return "not-applicable";
    }

    private void resolveOutputMetadata(
            String model, boolean genotype) {
        resolvedStatisticType = switch (model) {
            case "glmm", "cox" -> "z";
            case "glm" -> "t_approx";
            default -> "t";
        };
        resolvedDfMethod = switch (model) {
            case "glmm", "cox" -> "asymptotic";
            case "glm" -> "residual-approximation";
            case "lmm" -> genotype
                ? "residual-approximation"
                : dfMethod(options.degreesOfFreedom).name()
                    .toLowerCase(Locale.ROOT).replace('_', '-');
            default -> "residual";
        };
        resolvedPartialR2Method = resolvedStatisticType.equals("t")
            ? "test-statistic" : "not-applicable";
    }

    private void logOutputMetadata() throws IOException {
        info("statistic_type=" + resolvedStatisticType);
        info("df_method=" + resolvedDfMethod);
        info("partial_r2_method=" + resolvedPartialR2Method);
        info("output_format=" + outputFormat());
    }

    private String outputFormat() {
        return options.output.getFileName().toString()
            .toLowerCase(Locale.ROOT).endsWith(".csv") ? "csv" : "tsv";
    }

    private static GlmFamily family(String name) {
        return switch (name) {
            case "gaussian" -> GlmFamilies.gaussian();
            case "binomial" -> GlmFamilies.binomial();
            case "poisson" -> GlmFamilies.poisson();
            case "gamma" -> GlmFamilies.gamma();
            case "inverse-gaussian" -> GlmFamilies.inverseGaussian();
            case "quasi-binomial" -> GlmFamilies.quasiBinomial();
            case "quasi-poisson" -> GlmFamilies.quasiPoisson();
            default -> throw new IllegalArgumentException(
                "unsupported family: " + name);
        };
    }

    private static DegreesOfFreedomMethod dfMethod(String name) {
        return switch (name) {
            case "auto", "residual", "residual-approximation" ->
                DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION;
            case "satterth", "satterthwaite" ->
                DegreesOfFreedomMethod.SATTERTHWAITE;
            case "kr", "kenward-roger" ->
                DegreesOfFreedomMethod.KENWARD_ROGER;
            default -> throw new IllegalArgumentException(
                "unknown denominator-DF method: " + name);
        };
    }

    private static List<VarianceComponent> components(
            List<RandomEffectTerm> terms, int rows,
            VarianceComponent grm) {
        List<VarianceComponent> result =
            new ArrayList<>(randomComponents(terms, rows, grm));
        result.add(VarianceComponent.identity("residual", rows));
        return result;
    }

    private static List<VarianceComponent> randomComponents(
            List<RandomEffectTerm> terms, int rows,
            VarianceComponent grm) {
        List<VarianceComponent> result = new ArrayList<>();
        for (RandomEffectTerm term : terms) {
            double[] design = term.design();
            int columns = term.coefficients();
            double[] covariance = new double[rows * rows];
            for (int first = 0; first < rows; first++)
                for (int second = 0; second < rows; second++) {
                    double sum = 0.0;
                    for (int column = 0; column < columns; column++)
                        sum += design[first * columns + column]
                            * design[second * columns + column];
                    covariance[first * rows + second] = sum;
                }
            result.add(new VarianceComponent(
                term.name(), rows, covariance));
        }
        if (grm != null) result.add(grm);
        return result;
    }

    private static SparseRandomStructure sparseStructure(
            CompiledMixedFormula mixed, PedigreeContext pedigree) {
        if (mixed == null)
            throw new IllegalArgumentException(
                "a random-effect formula term is required");
        List<RandomEffectTerm> terms = new ArrayList<>();
        List<SparsePrecisionMatrix> precisions = new ArrayList<>();
        String mappedName = pedigree == null ? null
            : "1|" + pedigree.matchingColumn();
        if (pedigree != null) {
            PedigreeRandomEffectTerm value =
                PedigreeRandomEffectTerm.ofSparse(
                    mappedName, pedigree.observationIds(),
                    pedigree.loaded().individuals(),
                    pedigree.loaded().inbreedingCoefficients());
            terms.add(value.randomEffect());
            precisions.add(value.precision());
        }
        for (PedigreeRandomEffectTerm value :
                mixed.pedigreeRandomEffects()) {
            terms.add(value.randomEffect());
            precisions.add(value.precision());
        }
        for (RandomEffectTerm value : mixed.randomEffects()) {
            if (value.name().equals(mappedName)) continue;
            terms.add(value);
            precisions.add(
                SparsePrecisionMatrix.identity(value.coefficients()));
        }
        return new SparseRandomStructure(terms, precisions);
    }

    private static boolean mixedCompatibleWithDenseKinship(
            CompiledMixedFormula mixed, PedigreeContext pedigree) {
        return pedigree == null && (mixed == null
            || mixed.pedigreeRandomEffects().isEmpty()
                && mixed.correlatedRandomEffects().isEmpty());
    }

    private static double[] adjustedResponse(CompiledFormula fixed) {
        double[] response = fixed.response();
        double[] offset = fixed.offset();
        if (offset != null)
            for (int row = 0; row < response.length; row++)
                response[row] -= offset[row];
        return response;
    }

    private static void requireConverged(
            boolean converged, String model) {
        if (!converged)
            throw new IllegalArgumentException(
                model + " association model did not converge");
    }

    private CompiledMixedFormula compileMixed(
            String formula, org.jlinalg.formula.ModelTable table,
            PedigreeContext pedigree) {
        CompiledMixedFormula result = MixedFormula.compile(formula, table);
        if (pedigree != null) {
            String expected = "1|" + pedigree.matchingColumn();
            long matches = result.randomEffects().stream()
                .filter(term -> term.name().equals(expected)).count();
            if (matches != 1)
                throw new IllegalArgumentException(
                    "--pedigree requires the formula term (1|"
                        + pedigree.matchingColumn() + ")");
        }
        return result;
    }

    private GrmContext grm(
            PhenotypeData phenotype, PhenotypeData.Prepared prepared)
            throws IOException {
        if (options.grm == null) return null;
        GrmReader.Loaded loaded = GrmReader.read(options.grm);
        String matchingColumn = options.individualId == null
            ? options.idColumn : options.individualId;
        List<String> observationIds = matchingColumn.equals(options.idColumn)
            ? prepared.ids()
            : phenotype.alignedValues(prepared.ids(), matchingColumn);
        VarianceComponent component = loaded.matrix()
            .varianceComponent("grm", observationIds);
        info("grm=" + options.grm.toAbsolutePath());
        info("grm_format=" + loaded.format());
        info("grm_samples=" + loaded.matrix().samples());
        info("grm_match_column=" + matchingColumn);
        info("grm_dense_bytes=" + Math.multiplyExact(
            (long) loaded.matrix().samples() * loaded.matrix().samples(),
            Double.BYTES));
        return new GrmContext(
            loaded.matrix(), observationIds, component, matchingColumn);
    }

    private PedigreeContext pedigree(
            PhenotypeData phenotype, PhenotypeData.Prepared prepared)
            throws IOException {
        if (options.pedigree == null) return null;
        String matchingColumn = options.individualId == null
            ? options.idColumn : options.individualId;
        List<String> observationIds = matchingColumn.equals(options.idColumn)
            ? prepared.ids()
            : phenotype.alignedValues(prepared.ids(), matchingColumn);
        PedigreeReader.Loaded value = PedigreeReader.read(
            options.pedigree, options.pedigreeId, options.sireId,
            options.damId, options.pedigreeFamilyId);
        java.util.Set<String> pedigreeIds = value.individuals().stream()
            .map(org.jlinalg.pedigree.PedigreeIndividual::id)
            .collect(java.util.stream.Collectors.toUnmodifiableSet());
        for (String id : observationIds) {
            if (!pedigreeIds.contains(id))
                throw new IllegalArgumentException(
                    "phenotype pedigree ID is absent from the pedigree: "
                        + id + " (matching column " + matchingColumn + ")");
        }
        info("pedigree=" + options.pedigree.toAbsolutePath());
        info("pedigree_members=" + value.individuals().size());
        info("pedigree_match_column=" + matchingColumn);
        info("pedigree_precision=sparse_additive_relationship_inverse");
        return new PedigreeContext(value, observationIds, matchingColumn);
    }

    private static double[][] matrix(double[] values, int rows, int columns) {
        double[][] result = new double[rows][columns];
        for (int row = 0; row < rows; row++)
            System.arraycopy(values, row * columns,
                result[row], 0, columns);
        return result;
    }

    private static double[] filled(int size, double value) {
        double[] result = new double[size];
        Arrays.fill(result, value);
        return result;
    }

    private void logBinary(PhenotypeData.Prepared prepared)
            throws IOException {
        if (prepared.binaryMapping() == null) return;
        info("case_value=" + prepared.binaryMapping().caseValue());
        info("control_value=" + prepared.binaryMapping().controlValue());
    }

    private void reportSampleAlignment(
            int omicsSamples, int phenotypeSamples, int alignedSamples,
            int analysisSamples)
            throws IOException {
        int omicsOnly = omicsSamples - alignedSamples;
        int phenotypeOnly = phenotypeSamples - alignedSamples;
        int phenotypeMissing = alignedSamples - analysisSamples;
        info("omics_samples=" + omicsSamples);
        info("phenotype_samples=" + phenotypeSamples);
        info("aligned_samples=" + alignedSamples);
        info("omics_only_samples=" + omicsOnly);
        info("phenotype_only_samples=" + phenotypeOnly);
        info("phenotype_missing_samples_omitted=" + phenotypeMissing);
        info("analysis_samples=" + analysisSamples);
        output.println("Aligned samples: " + alignedSamples
            + " (omics=" + omicsSamples + ", phenotype=" + phenotypeSamples
            + ", omics-only=" + omicsOnly
            + ", phenotype-only=" + phenotypeOnly + ")");
        output.println("Analysis samples: " + analysisSamples
            + " (phenotype-missing omitted=" + phenotypeMissing + ")");
    }

    private void reportPhenotypeSamples(
            int phenotypeSamples, int analysisSamples) throws IOException {
        int phenotypeMissing = phenotypeSamples - analysisSamples;
        info("phenotype_samples=" + phenotypeSamples);
        info("phenotype_missing_samples_omitted=" + phenotypeMissing);
        info("analysis_samples=" + analysisSamples);
        output.println("Analysis samples: " + analysisSamples
            + " (phenotype=" + phenotypeSamples
            + ", phenotype-missing omitted=" + phenotypeMissing + ")");
    }

    private PhenotypeData.Prepared completeCases(
            PhenotypeData phenotype, PhenotypeData.Prepared prepared,
            String response, boolean encodeBinomial) {
        int[] retained = completeCaseRows(prepared.modelTable());
        if (retained.length == prepared.ids().size()) return prepared;
        List<String> ids = new ArrayList<>(retained.length);
        for (int row : retained) ids.add(prepared.ids().get(row));
        return phenotype.prepare(ids, response, encodeBinomial,
            options.caseValue, options.controlValue);
    }

    private int[] completeCaseRows(org.jlinalg.formula.ModelTable table) {
        if (plan.isCox()) {
            String rhs = plan.withoutOmics().substring(
                plan.withoutOmics().indexOf('~') + 1);
            FormulaPlan.Survival survival = plan.survival();
            int[] retained = Formula.compile(
                survival.stop() + "~0+" + rhs, table,
                MissingDataPolicy.OMIT).retainedRows();
            retained = intersectRows(retained, Formula.compile(
                survival.event() + "~1", table,
                MissingDataPolicy.OMIT).retainedRows());
            if (survival.start() != null) {
                retained = intersectRows(retained, Formula.compile(
                    survival.start() + "~1", table,
                    MissingDataPolicy.OMIT).retainedRows());
            }
            return retained;
        }
        if (plan.hasRandomEffects()) {
            return MixedFormula.compile(plan.withoutOmics(), table,
                MixedFormulaOptions.builder()
                    .missingDataPolicy(MissingDataPolicy.OMIT)
                    .build()).retainedRows();
        }
        return Formula.compile(plan.withoutOmics(), table,
            MissingDataPolicy.OMIT).retainedRows();
    }

    private static int[] intersectRows(int[] left, int[] right) {
        int[] result = new int[Math.min(left.length, right.length)];
        int leftIndex = 0;
        int rightIndex = 0;
        int count = 0;
        while (leftIndex < left.length && rightIndex < right.length) {
            if (left[leftIndex] == right[rightIndex]) {
                result[count++] = left[leftIndex];
                leftIndex++;
                rightIndex++;
            } else if (left[leftIndex] < right[rightIndex]) {
                leftIndex++;
            } else {
                rightIndex++;
            }
        }
        if (count == 0) {
            throw new IllegalArgumentException(
                "no complete observations remain");
        }
        return Arrays.copyOf(result, count);
    }

    private void manifest(
            String model, String omicsType, int blockSize, Counts counts)
            throws IOException {
        new ManifestWriter()
            .put("run_id", log == null ? null : log.runId())
            .put("version", JLinAlgCli.version())
            .put("phenotype", options.phenotype.toAbsolutePath())
            .put("omics", options.omics == null ? null
                : options.omics.toAbsolutePath())
            .put("grm", options.grm == null ? null
                : options.grm.toAbsolutePath())
            .put("grm_match_column", options.grm == null ? null
                : options.individualId == null
                    ? options.idColumn : options.individualId)
            .put("pedigree", options.pedigree == null ? null
                : options.pedigree.toAbsolutePath())
            .put("pedigree_match_column", options.pedigree == null ? null
                : options.individualId == null
                    ? options.idColumn : options.individualId)
            .put("pedigree_precision", options.pedigree == null ? null
                : "sparse_additive_relationship_inverse")
            .put("formula", options.formula)
            .put("model", model)
            .put("family", options.family)
            .put("omics_type", omicsType)
            .put("variance_components", resolvedVarianceComponents)
            .put("mixed_fit", resolvedMixedFit)
            .put("df", options.degreesOfFreedom)
            .put("statistic_type", resolvedStatisticType)
            .put("df_method", resolvedDfMethod)
            .put("partial_r2_method", resolvedPartialR2Method)
            .put("output_format", outputFormat())
            .put("block_size", blockSize)
            .put("threads", options.threads)
            .put("transform_plugins", options.transformPlugins)
            .put("source_features", counts.source())
            .put("tested_features", counts.tested())
            .put("failed_features", counts.failed())
            .put("output", options.output.toAbsolutePath())
            .write(options.manifestPath());
    }

    private void info(String message) throws IOException {
        if (log != null) log.info(message);
    }
    private void warning(String message) throws IOException {
        if (log != null) log.warning(message);
    }

    private record Counts(long source, long tested, long failed) {
        private static Counts of(AssociationPipelineSummary value) {
            return new Counts(value.sourceVariants(), value.testedVariants(),
                value.failures());
        }
        private static Counts of(OmicsAssociationSummary value) {
            return new Counts(value.sourceFeatures(), value.testedFeatures(),
                value.failedFeatures());
        }
    }

    private record GrmContext(
        GenomicRelationshipMatrix matrix, List<String> observationIds,
        VarianceComponent component, String matchingColumn) {
        private GrmContext {
            observationIds = List.copyOf(observationIds);
        }
    }

    private record PedigreeContext(
        PedigreeReader.Loaded loaded, List<String> observationIds,
        String matchingColumn) {
        private PedigreeContext {
            observationIds = List.copyOf(observationIds);
        }
    }

    private record SparseRandomStructure(
        List<RandomEffectTerm> terms,
        List<SparsePrecisionMatrix> precisions) {
        private SparseRandomStructure {
            terms = List.copyOf(terms);
            precisions = List.copyOf(precisions);
        }
    }
}
