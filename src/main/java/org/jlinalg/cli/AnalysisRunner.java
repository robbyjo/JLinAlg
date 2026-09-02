/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.association.AssociationFailurePolicy;
import org.jlinalg.association.VariableMissingPolicy;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.CompiledFormula;
import org.jlinalg.formula.CompiledMixedFormula;
import org.jlinalg.formula.Formula;
import org.jlinalg.formula.MixedFormula;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.glmm.GlmmPql;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.glmm.GlmmPqlResult;
import org.jlinalg.genetics.GenomicRelationshipMatrix;
import org.jlinalg.gwas.AssociationScanOptions;
import org.jlinalg.gwas.GenotypeMissingPolicy;
import org.jlinalg.gwas.RemlAssociationScanner;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.mixed.RandomEffectTerm;
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

    AnalysisRunner(CliOptions options, FormulaPlan plan, RunLog log) {
        this.options = options;
        this.plan = plan;
        this.log = log;
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
        logBinary(prepared);
        GrmContext grm = grm(phenotype, prepared);
        String model = resolveModel();
        CompiledFormula fixed;
        CompiledMixedFormula mixed = null;
        if (plan.hasRandomEffects()) {
            mixed = MixedFormula.compile(
                plan.withoutOmics(), prepared.modelTable());
            fixed = mixed.fixed();
        } else {
            fixed = Formula.compile(
                plan.withoutOmics(), prepared.modelTable());
        }
        double[][] covariates = matrix(
            fixed.design(), fixed.rows(), fixed.columns());
        int blockSize = AdaptiveBlockSizer.choose(sourceIds.size(),
            options.blockSize);
        info("resolved_model=" + model);
        info("variance_components=" + options.varianceComponents);
        info("block_size=" + blockSize);
        info("threads=" + options.threads);
        info("backend=" + options.backend);
        if (options.explain || options.dryRun) {
            System.out.println("omics type: " + detection.type()
                + " (" + detection.source() + ")");
            System.out.println("model: " + model);
            System.out.println("samples: " + sourceIds.size());
            System.out.println("adaptive block size: " + blockSize);
            System.out.println("backend: " + options.backend);
            System.out.println("output: " + options.output);
            if (options.dryRun) return 0;
        }
        if (options.resume && Files.exists(Path.of(options.output + ".partial")))
            throw new IllegalArgumentException(
                "partial output lacks resumable block metadata; "
                    + "use --overwrite to restart");
        if (options.varianceComponents.equals("refit"))
            throw new IllegalArgumentException(
                "per-feature variance-component refitting is not yet wired "
                    + "to the streaming CLI");
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
        String statisticType = model.equals("glm") ? "t_approx" : "t";
        String dfMethod = model.equals("glm") || model.equals("lmm")
            ? "residual-approximation" : "residual";
        Counts counts;
        try (CliResultSink sink = new CliResultSink(options.output,
                options.overwrite, detection.type(), statisticType, dfMethod,
                annotation, genotype ? prepared.caseControlGroups() : null)) {
            counts = switch (model) {
                case "ols" -> scanOls(variantSource, numericSource, genotype,
                    sourceIds, fixed, covariates, transform, blockSize,
                    engine, pipeline, sink);
                case "glm" -> scanGlm(variantSource, numericSource, genotype,
                    sourceIds, fixed, covariates, transform, blockSize,
                    engine, pipeline, sink);
                case "lmm" -> scanLmm(variantSource, genotype, sourceIds,
                    fixed, mixed, grm, covariates, blockSize, pipeline, sink);
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
            VariantSource variant, boolean genotype, List<String> ids,
            CompiledFormula fixed, CompiledMixedFormula mixed,
            GrmContext grm, double[][] covariates, int blockSize,
            AssociationPipelineOptions pipeline, CliResultSink sink)
            throws IOException {
        if (!genotype)
            throw new IllegalArgumentException(
                "streamed non-genotype LMM scans are not yet available");
        if ((mixed == null && grm == null)
                || mixed != null
                    && !mixed.correlatedRandomEffects().isEmpty())
            throw new IllegalArgumentException(
                "omics LMM requires a GRM or independent random-effect terms");
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

    private int phenotypeOnly() throws IOException {
        String model = resolveModel();
        String preparedResponse = plan.isCox()
            ? plan.survival().event() : plan.response();
        PhenotypeData phenotype = PhenotypeData.read(
            options.phenotype, options.idColumn);
        PhenotypeData.Prepared prepared = phenotype.prepare(null,
            preparedResponse, options.family.equals("binomial") || plan.isCox(),
            options.caseValue, options.controlValue);
        GrmContext grm = grm(phenotype, prepared);
        info("resolved_model=" + model);
        long tests = switch (model) {
            case "ols" -> phenotypeOls(prepared);
            case "glm" -> phenotypeGlm(prepared);
            case "lmm" -> phenotypeLmm(prepared, grm);
            case "glmm" -> phenotypeGlmm(prepared, grm);
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
            "t", "residual", null, "transformed_effect");
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
            "t_approx", "residual-approximation",
            null, "transformed_effect");
    }

    private long phenotypeLmm(
            PhenotypeData.Prepared prepared, GrmContext grm)
            throws IOException {
        CompiledFormula fixed;
        List<RandomEffectTerm> randomEffects;
        if (plan.hasRandomEffects()) {
            CompiledMixedFormula compiled = MixedFormula.compile(
                plan.withoutOmics(), prepared.modelTable());
            if (!compiled.correlatedRandomEffects().isEmpty())
                throw new IllegalArgumentException(
                    "correlated random slopes are not yet exposed by the CLI");
            fixed = compiled.fixed();
            randomEffects = compiled.randomEffects();
        } else {
            fixed = Formula.compile(
                plan.withoutOmics(), prepared.modelTable());
            randomEffects = List.of();
        }
        if (grm == null && randomEffects.isEmpty())
            throw new IllegalArgumentException(
                "LMM requires --grm or a random-effect formula term");
        if (fixed.weights() != null || fixed.offset() != null)
            throw new IllegalArgumentException(
                "weighted or offset LMM is not supported by this CLI path");
        RemlOptions reml = RemlOptions.builder()
            .degreesOfFreedomMethod(dfMethod(options.degreesOfFreedom)).build();
        RemlResult fit = Reml.fit(fixed.response(), fixed.design(),
            fixed.rows(), fixed.columns(),
            components(randomEffects, fixed.rows(),
                grm == null ? null : grm.component()),
            reml, options.backend);
        AssociationStatistics statistics = fit.associationStatistics();
        return CoefficientOutput.write(options.output, options.overwrite,
            fixed.coefficientNames(), statistics.beta(),
            statistics.standardErrors(), statistics.statistics(),
            statistics.degreesOfFreedom(), statistics.pValues(), "t",
            statistics.degreesOfFreedomMethod().name()
                .toLowerCase(Locale.ROOT),
            null, "transformed_effect");
    }

    private long phenotypeGlmm(
            PhenotypeData.Prepared prepared, GrmContext grm)
            throws IOException {
        CompiledFormula fixed;
        List<RandomEffectTerm> randomEffects;
        if (plan.hasRandomEffects()) {
            CompiledMixedFormula compiled = MixedFormula.compile(
                plan.withoutOmics(), prepared.modelTable());
            if (!compiled.correlatedRandomEffects().isEmpty())
                throw new IllegalArgumentException(
                    "correlated random slopes are not supported by PQL CLI");
            fixed = compiled.fixed();
            randomEffects = compiled.randomEffects();
        } else {
            fixed = Formula.compile(
                plan.withoutOmics(), prepared.modelTable());
            randomEffects = List.of();
        }
        if (grm == null && randomEffects.isEmpty())
            throw new IllegalArgumentException(
                "GLMM requires --grm or a random-effect formula term");
        GlmmPqlOptions pql = GlmmPqlOptions.builder()
            .remlOptions(RemlOptions.builder()
                .degreesOfFreedomMethod(
                    DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION)
                .build())
            .build();
        GlmmPqlResult fit = GlmmPql.fit(
            fixed.response(),
            matrix(fixed.design(), fixed.rows(), fixed.columns()),
            family(options.family),
            components(randomEffects, fixed.rows(),
                grm == null ? null : grm.component()),
            fixed.weights(), fixed.offset(),
            pql, options.backend);
        return CoefficientOutput.write(options.output, options.overwrite,
            fixed.coefficientNames(), fit.beta(),
            fit.standardErrors(), fit.statistics(), fit.degreesOfFreedom(),
            fit.pValues(), "t_approx", "residual-approximation",
            null, "transformed_effect");
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
            pValues, "z", "asymptotic", hazardRatios, "hazard_ratio");
    }

    private String resolveModel() {
        String resolved;
        if (!options.model.equals("auto")) resolved = options.model;
        else if (plan.isCox()) resolved = "cox";
        else if (plan.hasRandomEffects() || options.grm != null)
            resolved = options.family.equals("gaussian") ? "lmm" : "glmm";
        else resolved = options.family.equals("gaussian") ? "ols" : "glm";
        if (options.grm != null && !resolved.equals("lmm")
                && !resolved.equals("glmm") && !resolved.equals("cox"))
            throw new IllegalArgumentException(
                "--grm is applicable to lmm, glmm, and Cox models");
        return resolved;
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
        result.add(VarianceComponent.identity("residual", rows));
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
            .put("formula", options.formula)
            .put("model", model)
            .put("family", options.family)
            .put("omics_type", omicsType)
            .put("variance_components", options.varianceComponents)
            .put("df", options.degreesOfFreedom)
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
}
