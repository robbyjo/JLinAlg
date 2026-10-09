/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.plot.DiagramPlots;
import org.jlinalg.plot.LegendPosition;
import org.jlinalg.plot.PlotOptions;
import org.jlinalg.plot.Plots;
import org.jlinalg.study.StudyOperationExecutor;
import org.jlinalg.study.StudyNeed;
import org.yaml.snakeyaml.Yaml;

/** Concrete adapters from study operations to existing JLinAlg CLI workflows. */
final class StudyCliExecutor implements StudyOperationExecutor {
    private static final Set<String> BULK_ASSOCIATION_TYPES = Set.of(
        "generic", "transcriptomics", "expression", "mirna", "lincrna",
        "lncrna", "extracellular_rna", "epigenomics", "ewas", "gwas",
        "genomics", "proteomics", "metabolomics", "phenotype");
    private static final Set<String> CUSTOM_COMMANDS = Set.of(
        "enrichment", "variant-score", "susie", "coloc", "mr-estimate", "mr-xwas",
        "clump", "mr-instruments", "ldsc", "twas", "pwas", "genomic-factor", "network",
        "differential", "ewas-regions", "multiple-test",
        "multiple-impute", "mi-pool");

    @Override
    public Result execute(Context context) throws Exception {
        Files.createDirectories(context.outputDirectory());
        String type = text(context.operation().get("type"));
        return switch (type) {
            case "validate_harmonize" -> harmonize(context);
            case "association" -> association(context);
            case "effect_comparison" -> effectComparison(context);
            case "discovery_replication" -> replication(context);
            case "meta_analysis" -> meta(context);
            case "diagnostic_review" -> diagnostics(context);
            case "report" -> report(context);
            case "participant_characteristics" -> new StudyParticipantCharacteristics().execute(context);
            case "irqtl_candidate_review" -> new StudyIrQtlReview().execute(context);
            case "candidate_prioritization" -> prioritization(context);
            case "source_discovery", "bioinformatics_investigation", "literature_search" -> new StudyBioinformatics().execute(context);
            case "literature_support", "advisor_review", "related_phenotype_advice", "manuscript_draft" -> new StudyAssistance().execute(context);
            case "sensitivity_analysis" -> sensitivity(context);
            case "gene_ontology" -> automaticOrConfiguredEnrichment(
                context, "GO", "gene-ontology.tsv");
            case "pathway_enrichment" -> automaticOrConfiguredEnrichment(
                context, "Reactome", "pathway-enrichment.tsv");
            case "variant_scoring", "fine_mapping", "colocalization",
                    "mendelian_randomization" -> genetic(context);
            default -> throw new IllegalArgumentException(
                "Unsupported study operation: " + type);
        };
    }

    private Result harmonize(Context context) throws IOException {
        Map<String, Map<String, Object>> datasets = datasets(context.study());
        List<String> selected = strings(context.operation().get("datasets"));
        List<String> outputs = new ArrayList<>();
        Map<String, Path> standardized = new LinkedHashMap<>();
        Map<String, Standardization> summaries = new LinkedHashMap<>();
        Map<String, Object> hypothesis = hypothesis(context.study(),
            text(context.operation().get("hypothesis")));
        Map<String, Object> analysis = map(
            hypothesis.get("analysis"), "hypothesis.analysis");
        Map<String, Object> harmonization = optionalMap(
            analysis.get("harmonization"), "analysis.harmonization");
        String invalidRows = text(harmonization.getOrDefault(
            "invalid_rows", "error")).toLowerCase(Locale.ROOT);
        if (!Set.of("error", "exclude").contains(invalidRows))
            throw new IllegalArgumentException(
                "analysis.harmonization.invalid_rows must be error or exclude");
        String commonUniverse = text(harmonization.getOrDefault(
            "common_universe", "none")).toLowerCase(Locale.ROOT);
        if (!Set.of("none", "intersection").contains(commonUniverse))
            throw new IllegalArgumentException(
                "analysis.harmonization.common_universe must be none or intersection");
        Set<String> exactCommon = null;
        Set<String> validCommon = null;
        StringBuilder table = new StringBuilder(
            "dataset_id\tfile_role\tpath\tsha256\n");
        for (String id : selected) {
            Map<String, Object> dataset = required(datasets, id, "dataset");
            Map<String, Object> files = map(dataset.get("files"), "dataset.files");
            for (Map.Entry<String, Object> file : files.entrySet()) {
                Path path = context.project().resolve(text(file.getValue())).normalize();
                if (!Files.isRegularFile(path))
                    throw new IOException("Dataset file is absent: " + path);
                table.append(id).append('\t').append(file.getKey()).append('\t')
                    .append(portable(context.project(), path)).append('\t')
                    .append(FollowupSupport.hash(path)).append('\n');
            }
            String stage = text(dataset.get("stage"));
            if (stage.endsWith("_results")) {
                Path result = context.project().resolve(
                    text(files.get("results"))).normalize();
                Path standard = context.outputDirectory()
                    .resolve("standard-results").resolve(id + ".tsv");
                Files.createDirectories(standard.getParent());
                Standardization summary = standardize(result,
                    optionalMap(dataset.get("columns"), "dataset.columns"),
                    standard, invalidRows.equals("exclude"));
                summaries.put(id, summary);
                standardized.put(id, standard);
                if (exactCommon == null) {
                    exactCommon = new LinkedHashSet<>(summary.allFeatures());
                    validCommon = new LinkedHashSet<>(summary.eligibleFeatures());
                } else {
                    exactCommon.retainAll(summary.allFeatures());
                    validCommon.retainAll(summary.eligibleFeatures());
                }
                outputs.add("standard-results/" + id + ".tsv");
            }
        }
        if (commonUniverse.equals("intersection")) {
            if (validCommon == null || validCommon.isEmpty())
                throw new IOException("Strict common feature universe is empty");
            for (Path standard : standardized.values())
                retainFeatures(standard, validCommon);
        }
        Files.writeString(context.outputDirectory().resolve("inputs.tsv"),
            table, StandardCharsets.UTF_8);
        outputs.add("inputs.tsv");
        if (!summaries.isEmpty()) {
            StringBuilder summary = new StringBuilder(
                "dataset_id\tinput_rows\teligible_rows\texcluded_invalid"
                    + "\texcluded_status\tretained_common\n");
            int retained = commonUniverse.equals("intersection")
                ? validCommon.size() : -1;
            for (Map.Entry<String, Standardization> entry : summaries.entrySet()) {
                Standardization value = entry.getValue();
                summary.append(entry.getKey()).append('\t')
                    .append(value.inputRows()).append('\t')
                    .append(value.eligibleFeatures().size()).append('\t')
                    .append(value.excludedInvalid()).append('\t')
                    .append(value.excludedStatus()).append('\t')
                    .append(retained < 0 ? value.eligibleFeatures().size() : retained)
                    .append('\n');
            }
            summary.append("# exact_common_features\t")
                .append(exactCommon == null ? 0 : exactCommon.size()).append('\n')
                .append("# strict_common_features\t")
                .append(validCommon == null ? 0 : validCommon.size()).append('\n')
                .append("# common_universe_policy\t").append(commonUniverse)
                .append('\n');
            Files.writeString(context.outputDirectory()
                .resolve("harmonization-summary.tsv"), summary,
                StandardCharsets.UTF_8);
            outputs.add("harmonization-summary.tsv");
        }
        outputs.addAll(StudyFeatureAnnotations.write(context, standardized));
        boolean retryable;
        try (var reader = Files.newBufferedReader(context.outputDirectory().resolve("annotation-summary.yaml"))) {
            retryable = list(map(new Yaml().load(reader), "annotation summary").get("input_requests"), "input_requests")
                .stream().anyMatch(item -> Boolean.TRUE.equals(map(item, "request").get("retryable")));
        }
        return new Result(State.COMPLETE, "Validated and hashed " + selected.size()
            + " dataset inputs", outputs, Map.of("retryable", retryable));
    }

    private Result association(Context context) throws IOException {
        String datasetId = text(context.operation().get("dataset_id"));
        Map<String, Object> dataset = required(
            datasets(context.study()), datasetId, "dataset");
        Map<String, Object> hypothesis = hypothesis(context.study(),
            text(context.operation().get("hypothesis")));
        String datasetType = text(dataset.get("type")).toLowerCase(Locale.ROOT);
        if (!BULK_ASSOCIATION_TYPES.contains(datasetType))
            return waiting(context, "No automatic association adapter is registered for "
                + datasetType + "; configure a modality-specific operation");
        Map<String, Object> files = map(dataset.get("files"), "dataset.files");
        Map<String, Object> columns = optionalMap(
            dataset.get("columns"), "dataset.columns");
        Map<String, Object> analysis = map(
            hypothesis.get("analysis"), "hypothesis.analysis");
        Map<String, Object> options = optionalMap(
            analysis.get("association"), "analysis.association");
        String stage = text(dataset.get("stage")).toLowerCase(Locale.ROOT);
        if (!stage.equals("processed_matrix") && !options.containsKey("model"))
            return waiting(context, "Association from " + stage
                + " requires an explicit analysis.association model and preprocessing decision");
        if (!files.containsKey("phenotype"))
            return waiting(context,
                "Matrix association requires files.phenotype");
        if (hypothesis.get("outcome") == null)
            return waiting(context,
                "Matrix association requires hypothesis.outcome");
        Path matrix = context.project().resolve(text(files.get("matrix"))).normalize();
        Path phenotype = context.project().resolve(
            text(files.get("phenotype"))).normalize();
        if (!Files.isRegularFile(matrix) || !Files.isRegularFile(phenotype))
            throw new IOException("Association requires existing matrix and phenotype files");

        String outcomeKey = text(hypothesis.get("outcome"));
        Map<String, Object> phenotypes = optionalMap(
            columns.get("phenotypes"), "columns.phenotypes");
        String outcome = phenotypes.containsKey(outcomeKey)
            ? text(phenotypes.get(outcomeKey)) : outcomeKey;
        String id = text(columns.getOrDefault("sample_id", "id"));
        Map<String, Object> covariateColumns = optionalMap(
            columns.get("covariates"), "columns.covariates");
        List<String> covariates = new ArrayList<>();
        for (String logical : strings(hypothesis.get("covariates")))
            covariates.add(text(covariateColumns.getOrDefault(logical, logical)));
        String formula = outcome + " ~ "
            + (covariates.isEmpty() ? "" : String.join(" + ", covariates) + " + ")
            + "<omics>";
        String model = text(options.getOrDefault("model", "ols"));
        Path raw = context.outputDirectory().resolve("association.tsv");
        List<String> arguments = new ArrayList<>(List.of(
            "--pheno", phenotype.toString(), "--omics", matrix.toString(),
            "--id", id, "--formula", formula, "--model", model,
            "--omics-type", omicsType(datasetType),
            "--out", raw.toString(), "--threads",
            text(options.getOrDefault("threads", "1")), "--no-log",
            "--no-config"));
        if (options.get("family") != null)
            arguments.addAll(List.of("--family", text(options.get("family"))));
        if (options.get("link") != null)
            arguments.addAll(List.of("--link", text(options.get("link"))));
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();
        int status = JLinAlgCli.run(arguments.toArray(String[]::new),
            new PrintStream(stdout), new PrintStream(stderr));
        Files.writeString(context.outputDirectory().resolve("console.txt"),
            stdout.toString(StandardCharsets.UTF_8)
                + stderr.toString(StandardCharsets.UTF_8));
        if (status != 0)
            throw new IOException("Association command failed: "
                + stderr.toString(StandardCharsets.UTF_8).trim());
        Path standard = context.outputDirectory()
            .resolve("standard-results").resolve(datasetId + ".tsv");
        Files.createDirectories(standard.getParent());
        standardize(raw, Map.of("feature_id", "id", "effect", "beta",
            "standard_error", "standard_error", "p_value", "p_value",
            "status", "status"), standard, false);
        List<String> outputs = new ArrayList<>(List.of("association.tsv",
            "standard-results/" + datasetId + ".tsv", "console.txt"));
        outputs.addAll(StudyFeatureAnnotations.write(context, Map.of(datasetId, standard)));
        var method = StudyMethods.record(model, Map.of("dataset", datasetId, "formula", formula, "model", model,
            "covariates", covariates, "command_arguments", arguments, "threads", options.getOrDefault("threads", "1")),
            StudyPublicEvidence.list(options.get("references")), options.get("rationale"),
            model.equals("ols") ? List.of("linear regression") : List.of(model, "regression"));
        return new Result(State.COMPLETE, "Fitted association model for " + datasetId, outputs, Map.of("method_report", method));
    }

    private Result replication(Context context) throws IOException {
        return replication(context, results(context));
    }

    private Result replication(Context context, Map<String, Path> results)
            throws IOException {
        Map<String, Object> roles = map(
            context.operation().get("roles"), "operation.roles");
        List<String> discovery = role(roles, "discovery");
        List<String> replication = role(roles, "replication");
        if (discovery.isEmpty() || replication.isEmpty())
            return waiting(context,
                "Discovery-replication comparison requires both discovery and replication roles");
        Map<String, Object> hypothesis = hypothesis(context.study(),
            text(context.operation().get("hypothesis")));
        Map<String, Object> analysis = map(
            hypothesis.get("analysis"), "hypothesis.analysis");
        Map<String, Object> testing = optionalMap(
            analysis.get("multiple_testing"), "analysis.multiple_testing");
        String method = text(testing.getOrDefault("method", "BH"))
            .toUpperCase(Locale.ROOT);
        if (!Set.of("BH", "NONE").contains(method))
            throw new IllegalArgumentException(
                "Discovery-replication currently supports multiple_testing method BH or NONE");
        double familyAlpha = number(testing.getOrDefault("alpha", 0.05),
            "multiple_testing.alpha");
        double discoveryAlpha = number(
            analysis.getOrDefault("discovery_alpha", familyAlpha),
            "discovery_alpha");
        double replicationAlpha = number(
            analysis.getOrDefault("replication_alpha", familyAlpha),
            "replication_alpha");
        Map<String, Object> replicationOptions = optionalMap(
            analysis.get("replication"), "analysis.replication");
        String adjustmentScope = text(replicationOptions.getOrDefault(
            "adjustment_scope", "all_results")).toLowerCase(Locale.ROOT);
        if (!Set.of("all_results", "discovery_subset").contains(adjustmentScope))
            throw new IllegalArgumentException(
                "analysis.replication.adjustment_scope must be all_results or discovery_subset");
        boolean requireSameSign = booleanValue(replicationOptions.getOrDefault(
            "require_same_sign", Boolean.TRUE),
            "analysis.replication.require_same_sign");
        String thresholdComparison = text(replicationOptions.getOrDefault(
            "threshold_comparison", "less_than_or_equal"))
            .toLowerCase(Locale.ROOT);
        if (!Set.of("less_than", "less_than_or_equal")
                .contains(thresholdComparison))
            throw new IllegalArgumentException(
                "analysis.replication.threshold_comparison must be less_than or less_than_or_equal");
        boolean strictThreshold = thresholdComparison.equals("less_than");
        probability(discoveryAlpha, "discovery_alpha");
        probability(replicationAlpha, "replication_alpha");
        Path output = context.outputDirectory().resolve("replication.tsv");
        Path vennDirectory = context.outputDirectory().resolve("replication-venn");
        Files.createDirectories(vennDirectory);
        List<String> outputs = new ArrayList<>(List.of("replication.tsv"));
        StringBuilder vennSummary = new StringBuilder(
            "discovery_dataset\treplication_dataset\tdiscovery_only"
                + "\treplicated_overlap\treplication_criterion_only"
                + "\tdiscovery_set_size\treplication_criterion_set_size"
                + "\tadjustment_scope\trequire_same_sign\n");
        int vennIndex = 0;
        try (BufferedWriter writer = Files.newBufferedWriter(output)) {
            writer.write("feature_id\tdiscovery_dataset\treplication_dataset"
                + "\tdiscovery_beta\tdiscovery_p\tdiscovery_adjusted_p"
                + "\treplication_beta\treplication_p\treplication_adjusted_p"
                + "\tdiscovery_selected\tsame_sign\treplicated\n");
            for (String discoveryId : discovery) {
                Map<String, Estimate> left = estimates(required(
                    results, discoveryId, "harmonized result"));
                Map<String, Double> leftAdjusted = adjusted(left, method);
                Set<String> selected = new LinkedHashSet<>();
                for (Map.Entry<String, Double> value : leftAdjusted.entrySet())
                    if (significant(value.getValue(), discoveryAlpha,
                            strictThreshold))
                        selected.add(value.getKey());
                for (String replicationId : replication) {
                    Map<String, Estimate> right = estimates(required(
                        results, replicationId, "harmonized result"));
                    Map<String, Estimate> adjustmentFamily = right;
                    if (adjustmentScope.equals("discovery_subset")) {
                        adjustmentFamily = new LinkedHashMap<>();
                        for (String feature : selected) {
                            Estimate estimate = right.get(feature);
                            if (estimate != null)
                                adjustmentFamily.put(feature, estimate);
                        }
                    }
                    Map<String, Double> rightAdjusted = adjusted(
                        adjustmentFamily, method);
                    Set<String> replicationCriterion = new LinkedHashSet<>();
                    Set<String> replicatedFeatures = new LinkedHashSet<>();
                    for (Map.Entry<String, Estimate> row : left.entrySet()) {
                        Estimate other = right.get(row.getKey());
                        if (other == null) continue;
                        Estimate first = row.getValue();
                        double discoveryQ = leftAdjusted.get(row.getKey());
                        double replicationQ = rightAdjusted.getOrDefault(
                            row.getKey(), Double.NaN);
                        boolean discoverySelected = significant(discoveryQ,
                            discoveryAlpha, strictThreshold);
                        boolean same = Math.signum(first.beta)
                            == Math.signum(other.beta);
                        boolean replicationSelected = significant(replicationQ,
                            replicationAlpha, strictThreshold);
                        if (replicationSelected && (!requireSameSign || same))
                            replicationCriterion.add(row.getKey());
                        boolean replicated = discoverySelected
                            && replicationSelected && (!requireSameSign || same);
                        if (replicated) replicatedFeatures.add(row.getKey());
                        writer.write(row.getKey() + "\t" + discoveryId + "\t"
                            + replicationId + "\t" + first.beta + "\t"
                            + value(first.p) + "\t" + value(discoveryQ) + "\t"
                            + other.beta + "\t" + value(other.p) + "\t"
                            + value(replicationQ) + "\t" + discoverySelected + "\t"
                            + same + "\t" + replicated + "\n");
                    }
                    Set<String> discoveryOnly = new LinkedHashSet<>(selected);
                    discoveryOnly.removeAll(replicationCriterion);
                    Set<String> replicationOnly = new LinkedHashSet<>(replicationCriterion);
                    replicationOnly.removeAll(selected);
                    vennSummary.append(discoveryId).append('\t')
                        .append(replicationId).append('\t')
                        .append(discoveryOnly.size()).append('\t')
                        .append(replicatedFeatures.size()).append('\t')
                        .append(replicationOnly.size()).append('\t')
                        .append(selected.size()).append('\t')
                        .append(replicationCriterion.size()).append('\t')
                        .append(adjustmentScope).append('\t')
                        .append(requireSameSign).append('\n');
                    if (!selected.isEmpty() || !replicationCriterion.isEmpty()) {
                        String stem = String.format(Locale.ROOT, "%02d-%s-vs-%s",
                            ++vennIndex, fileToken(discoveryId),
                            fileToken(replicationId));
                        PlotOptions plotOptions = PlotOptions.builder().size(850, 620)
                            .title("Discovery and replication overlap")
                            .build();
                        var venn = DiagramPlots.venn(new DiagramPlots.VennData(
                            discoveryId + " selected", replicationId
                                + " replication criterion",
                            discoveryOnly.size(), replicatedFeatures.size(),
                            replicationOnly.size()), plotOptions);
                        Path png = vennDirectory.resolve(stem + ".png");
                        Path pdf = vennDirectory.resolve(stem + ".pdf");
                        venn.save(png);
                        venn.save(pdf);
                        outputs.add(context.outputDirectory().relativize(png)
                            .toString().replace('\\', '/'));
                        outputs.add(context.outputDirectory().relativize(pdf)
                            .toString().replace('\\', '/'));
                    }
                }
            }
        }
        Files.writeString(context.outputDirectory().resolve("replication-venn.tsv"),
            vennSummary, StandardCharsets.UTF_8);
        outputs.add(1, "replication-venn.tsv");
        return Result.complete("Compared discovery and replication effects and generated overlap diagrams",
            outputs);
    }

    private Result effectComparison(Context context) throws IOException {
        Map<String, Path> resultPaths = results(context);
        if (resultPaths.size() < 2)
            return new Result(State.NOT_APPLICABLE,
                "Effect comparison requires at least two result tables",
                List.of(), Map.of("available_result_tables", resultPaths.size()));
        Map<String, Map<String, Object>> described = datasets(context.study());
        List<String> ids = new ArrayList<>(resultPaths.keySet());
        List<String> outputs = new ArrayList<>();
        List<Object> comparisons = new ArrayList<>();
        Path plotDirectory = context.outputDirectory().resolve("effect-scatterplots");
        Files.createDirectories(plotDirectory);
        StringBuilder summary = new StringBuilder(
            "left_dataset\tright_dataset\tleft_cohort\tright_cohort"
                + "\tcomparison_type\tcommon_features\tpearson_r"
                + "\tpearson_r_display\n");
        int plotIndex = 0;
        for (int leftIndex = 0; leftIndex < ids.size(); leftIndex++) {
            String leftId = ids.get(leftIndex);
            Map<String, Estimate> left = estimates(resultPaths.get(leftId));
            for (int rightIndex = leftIndex + 1; rightIndex < ids.size(); rightIndex++) {
                String rightId = ids.get(rightIndex);
                Map<String, Estimate> right = estimates(resultPaths.get(rightId));
                List<Double> xValues = new ArrayList<>(), yValues = new ArrayList<>();
                Moment moment = new Moment();
                for (Map.Entry<String, Estimate> entry : left.entrySet()) {
                    Estimate other = right.get(entry.getKey());
                    if (other == null) continue;
                    xValues.add(entry.getValue().beta);
                    yValues.add(other.beta);
                    moment.add(entry.getValue().beta, other.beta, false, false);
                }
                String leftCohort = cohort(described.get(leftId));
                String rightCohort = cohort(described.get(rightId));
                String comparisonType = leftCohort.equals(rightCohort)
                    ? "within_cohort" : "cross_cohort";
                double correlation = moment.correlation();
                summary.append(leftId).append('\t').append(rightId).append('\t')
                    .append(leftCohort).append('\t').append(rightCohort).append('\t')
                    .append(comparisonType).append('\t').append(moment.n).append('\t')
                    .append(value(correlation)).append('\t')
                    .append(Double.isFinite(correlation)
                        ? String.format(Locale.ROOT, "%.4f", correlation) : "NA")
                    .append('\n');
                Map<String, Object> comparison = new LinkedHashMap<>();
                comparison.put("left_dataset", leftId);
                comparison.put("right_dataset", rightId);
                comparison.put("left_cohort", leftCohort);
                comparison.put("right_cohort", rightCohort);
                comparison.put("comparison_type", comparisonType);
                comparison.put("common_features", moment.n);
                comparison.put("pearson_r", Double.isFinite(correlation)
                    ? correlation : "NA");
                if (!xValues.isEmpty()) {
                    String stem = String.format(Locale.ROOT, "%02d-%s-vs-%s",
                        ++plotIndex, fileToken(leftId), fileToken(rightId));
                    Path png = plotDirectory.resolve(stem + ".png");
                    Path pdf = plotDirectory.resolve(stem + ".pdf");
                    double[] x = xValues.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] y = yValues.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] bounds = commonBounds(x, y);
                    String legend = "Pearson r = " + (Double.isFinite(correlation)
                        ? String.format(Locale.ROOT, "%.4f", correlation) : "NA")
                        + " (n = " + x.length + ")";
                    PlotOptions options = PlotOptions.builder().size(850, 700)
                        .title((comparisonType.equals("within_cohort")
                            ? "Within-cohort" : "Cross-cohort")
                            + " effect comparison")
                        .xLabel(leftId + " effect").yLabel(rightId + " effect")
                        .legend(LegendPosition.RIGHT).markerSize(7)
                        .xRange(bounds[0], bounds[1]).yRange(bounds[0], bounds[1])
                        .build();
                    var plot = Plots.xy().options(options)
                        .points(legend, x, y, new Color(0x00, 0x72, 0xB2))
                        .line("Identity (y = x)",
                            new double[] {bounds[0], bounds[1]},
                            new double[] {bounds[0], bounds[1]}, Color.DARK_GRAY)
                        .build();
                    plot.save(png);
                    plot.save(pdf);
                    String pngRelative = context.outputDirectory().relativize(png)
                        .toString().replace('\\', '/');
                    String pdfRelative = context.outputDirectory().relativize(pdf)
                        .toString().replace('\\', '/');
                    outputs.add(pngRelative);
                    outputs.add(pdfRelative);
                    comparison.put("plot_png", pngRelative);
                    comparison.put("plot_pdf", pdfRelative);
                }
                comparisons.add(comparison);
            }
        }
        Path table = context.outputDirectory().resolve("effect-comparisons.tsv");
        Files.writeString(table, summary, StandardCharsets.UTF_8);
        outputs.add(0, "effect-comparisons.tsv");
        return new Result(State.COMPLETE,
            "Generated " + comparisons.size() + " pairwise effect comparison(s)",
            outputs, Map.of("effect_comparisons", comparisons,
                "correlation_display", "Pearson correlation is shown to four decimal places in each plot legend"));
    }

    private Result meta(Context context) throws IOException {
        return meta(context, results(context));
    }

    private Result meta(Context context, Map<String, Path> results) throws IOException {
        if (results.isEmpty())
            return waiting(context, "Meta-analysis requires association results");
        Map<String, Object> configuration = optionalMap(
            context.operation().get("configuration"), "operation.configuration");
        String model = text(configuration.getOrDefault("model", "fixed"))
            .toLowerCase(Locale.ROOT);
        if (!Set.of("fixed", "random").contains(model))
            throw new IllegalArgumentException(
                "meta_analysis model must be fixed or random");
        int defaultMinimum = results.size() > 1 ? 2 : 1;
        int minimum = positiveInteger(
            configuration.getOrDefault("min_cohorts", defaultMinimum),
            "meta_analysis.min_cohorts");
        List<String> arguments = new ArrayList<>();
        for (Map.Entry<String, Path> result : results.entrySet())
            arguments.addAll(List.of("--cohort",
                result.getKey() + "=" + result.getValue()));
        Path output = context.outputDirectory().resolve("meta-analysis.tsv");
        arguments.addAll(List.of("--out", output.toString(),
            "--id-column", "feature_id", "--effect-column", "beta",
            "--se-column", "se", "--min-cohorts", Integer.toString(minimum),
            "--model", model));
        ByteArrayOutputStream console = new ByteArrayOutputStream();
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        int status = MetaCli.run("meta-analysis",
            arguments.toArray(String[]::new), new PrintStream(console),
            new PrintStream(error));
        Files.writeString(context.outputDirectory().resolve("console.txt"),
            console.toString(StandardCharsets.UTF_8)
                + error.toString(StandardCharsets.UTF_8));
        if (status != 0)
            throw new IOException("Meta-analysis failed: "
                + error.toString(StandardCharsets.UTF_8).trim());
        var defaults = org.jlinalg.meta.MetaAnalysisOptions.builder().build();
        Map<String, Object> parameters = new LinkedHashMap<>();
        parameters.put("model", model); parameters.put("weighting", "inverse variance"); parameters.put("min_cohorts", minimum);
        parameters.put("inference", defaults.inferenceMethod().toString()); parameters.put("confidence_level", defaults.confidenceLevel());
        parameters.put("tau_squared_estimator", model.equals("random") ? defaults.tauSquaredEstimator().toString() : "not applicable");
        parameters.put("maximum_iterations", defaults.maximumIterations()); parameters.put("tolerance", defaults.tolerance());
        parameters.put("command_arguments", arguments); parameters.put("effective_defaults", "Backend defaults captured at execution, not reconstructed for older runs.");
        parameters.put("overlap_handling", "No covariance-based correction; confirm cohort independence/overlap declarations.");
        return new Result(State.COMPLETE, "Pooled " + results.size() + " cohort result tables",
            List.of("meta-analysis.tsv", "meta-analysis.tsv.log", "console.txt"), Map.of("method_report",
                StudyMethods.record(model + "-effect inverse-variance meta-analysis", parameters,
                    StudyPublicEvidence.list(configuration.get("references")), configuration.get("rationale"), List.of("meta-analysis", "inverse variance"))));
    }

    private Result sensitivity(Context context) throws IOException {
        Map<String, Path> results = results(context);
        Map<String, Object> assistance = optionalMap(
            context.study().get("assistance"), "assistance");
        Map<String, Object> configuration = optionalMap(
            context.operation().get("configuration"), "operation.configuration");
        Set<String> actions = new LinkedHashSet<>(
            strings(assistance.get("allow_analyses")));
        if (Boolean.TRUE.equals(configuration.get(
                "automatic_leave_one_cohort_out")))
            actions.add("leave_one_cohort_out");
        List<Object> dispositions = new ArrayList<>();
        List<String> outputs = new ArrayList<>();
        boolean incomplete = false;
        // These are supplementary influence diagnostics, never a silent
        // replacement for the all-cohort primary analysis.
        boolean primaryReady = false;
        Path primaryMeta = null;
        Path primaryReplication = null;
        for (Path dependency : context.dependencyDirectories().values()) {
            if (Files.isRegularFile(dependency.resolve("replication.tsv")))
                primaryReplication = dependency.resolve("replication.tsv");
            if (Files.isRegularFile(dependency.resolve("meta-analysis.tsv"))) {
                primaryMeta = dependency.resolve("meta-analysis.tsv");
                try (var reader = Files.newBufferedReader(
                        dependency.resolve("manifest.yaml"))) {
                    primaryReady |= "complete".equals(map(new Yaml().load(reader),
                        "manifest").get("state"));
                }
            }
        }
        ReplicationMetrics primaryReplicationMetrics = replicationMetrics(
            primaryReplication);
        MetaMetrics primaryMetaMetrics = metaMetrics(primaryMeta);
        Map<String, Map<String, Object>> described = datasets(context.study());
        Map<String, List<String>> cohorts = new LinkedHashMap<>();
        boolean unambiguous = true;
        for (String id : results.keySet()) {
            Map<String, Object> dataset = described.get(id);
            Object cohort = dataset == null ? null : dataset.get("cohort");
            if (!(cohort instanceof String name) || name.isBlank()) { unambiguous = false; continue; }
            cohorts.computeIfAbsent((String) cohort, ignored -> new ArrayList<>()).add(id);
        }
        StringBuilder comparisonTable = new StringBuilder(
            "action\tomitted_cohort\tincluded_datasets"
                + "\tprimary_replication_fraction\tsensitivity_replication_fraction"
                + "\treplication_fraction_delta\tprimary_mean_i_squared"
                + "\tsensitivity_mean_i_squared\tmean_i_squared_delta"
                + "\tprimary_high_i_squared_fraction"
                + "\tsensitivity_high_i_squared_fraction"
                + "\todd_one_candidate\tmanuscript_disposition\n");
        List<Object> improvementCandidates = new ArrayList<>();
        for (String action : actions) {
            if (!primaryReady || results.size() < 2 || action.equals("leave_one_cohort_out") && (!unambiguous || cohorts.size() < 3)) {
                incomplete = true;
                dispositions.add(Map.of("action", action, "status", "waiting_for_input", "required_input",
                    action.equals("leave_one_cohort_out") ? "Completed primary meta-analysis and at least three explicitly identified cohorts"
                        : "Completed primary meta-analysis with at least two comparable cohort result tables"));
                continue;
            }
            List<String> variants = action.equals("fixed_random_comparison") ? List.of("fixed", "random") : new ArrayList<>(cohorts.keySet());
            for (int i = 0; i < variants.size(); i++) {
                String variant = variants.get(i);
                Path folder = context.outputDirectory().resolve(action + "/analysis-" + (i + 1));
                Files.createDirectories(folder);
                Map<String, Path> selected = new LinkedHashMap<>(results);
                if (action.equals("leave_one_cohort_out")) cohorts.get(variant).forEach(selected::remove);
                Map<String, Object> operation = new LinkedHashMap<>(context.operation());
                Map<String, Object> primary = optionalMap(map(hypothesis(context.study(), text(operation.get("hypothesis")))
                    .get("analysis"), "analysis").get("operations"), "operations");
                Map<String, Object> settings = new LinkedHashMap<>(optionalMap(primary.get("meta_analysis"), "meta_analysis"));
                if (positiveInteger(settings.getOrDefault("min_cohorts", 2), "min_cohorts") > selected.size()) {
                    incomplete = true;
                    dispositions.add(Map.of("action", action, "variant", variant, "status", "waiting_for_input",
                        "required_input", "The remaining cohorts do not meet the approved min_cohorts setting; changing that setting requires an amendment"));
                    continue;
                }
                if (action.equals("fixed_random_comparison")) settings.put("model", variant);
                operation.put("configuration", settings);
                Result result = meta(new Context(context.project(), folder, operation, context.study(), context.plan(), context.dependencyDirectories()), selected);
                for (String output : result.outputs()) outputs.add(context.outputDirectory().relativize(folder.resolve(output)).toString().replace('\\', '/'));
                Path replicationOutput = null;
                if (action.equals("leave_one_cohort_out")) {
                    Map<String, Object> sourceHypothesis = hypothesis(context.study(),
                        text(operation.get("hypothesis")));
                    Map<String, Object> sourceRoles = optionalMap(
                        sourceHypothesis.get("inputs"), "hypothesis.inputs");
                    Map<String, Object> filteredRoles = new LinkedHashMap<>();
                    for (Map.Entry<String, Object> role : sourceRoles.entrySet()) {
                        List<String> retained = strings(role.getValue()).stream()
                            .filter(selected::containsKey).toList();
                        if (!retained.isEmpty()) filteredRoles.put(role.getKey(), retained);
                    }
                    if (!role(filteredRoles, "discovery").isEmpty()
                            && !role(filteredRoles, "replication").isEmpty()) {
                        Map<String, Object> replicationOperation = new LinkedHashMap<>(operation);
                        replicationOperation.put("type", "discovery_replication");
                        replicationOperation.put("roles", filteredRoles);
                        Result replicationResult = replication(new Context(
                            context.project(), folder, replicationOperation,
                            context.study(), context.plan(),
                            context.dependencyDirectories()), selected);
                        for (String output : replicationResult.outputs()) {
                            String relative = context.outputDirectory().relativize(
                                folder.resolve(output)).toString().replace('\\', '/');
                            if (!outputs.contains(relative)) outputs.add(relative);
                        }
                        replicationOutput = folder.resolve("replication.tsv");
                    }
                }
                MetaMetrics sensitivityMeta = metaMetrics(
                    folder.resolve("meta-analysis.tsv"));
                ReplicationMetrics sensitivityReplication = replicationMetrics(
                    replicationOutput);
                double replicationDelta = difference(
                    sensitivityReplication.fraction(),
                    primaryReplicationMetrics.fraction());
                double heterogeneityDelta = difference(
                    sensitivityMeta.meanISquared(),
                    primaryMetaMetrics.meanISquared());
                boolean replicationImproved = Double.isFinite(replicationDelta)
                    && replicationDelta > 1e-12;
                boolean heterogeneityImproved = Double.isFinite(heterogeneityDelta)
                    && heterogeneityDelta < -1e-12;
                boolean oddOneCandidate = action.equals("leave_one_cohort_out")
                    && replicationImproved && heterogeneityImproved;
                String manuscriptDisposition = oddOneCandidate
                    ? "feature_in_sensitivity_section_primary_unchanged"
                    : "retain_in_supplementary_sensitivity_results";
                Map<String, Object> disposition = new LinkedHashMap<>();
                disposition.put("action", action);
                disposition.put("variant", variant);
                disposition.put("status", "complete");
                disposition.put("result", context.outputDirectory().relativize(
                    folder.resolve("meta-analysis.tsv")).toString().replace('\\', '/'));
                disposition.put("included_datasets", new ArrayList<>(selected.keySet()));
                disposition.put("primary_metrics", metrics(
                    primaryReplicationMetrics, primaryMetaMetrics));
                disposition.put("sensitivity_metrics", metrics(
                    sensitivityReplication, sensitivityMeta));
                disposition.put("replication_fraction_delta", numeric(replicationDelta));
                disposition.put("mean_i_squared_delta", numeric(heterogeneityDelta));
                disposition.put("replication_improved", replicationImproved);
                disposition.put("heterogeneity_improved", heterogeneityImproved);
                disposition.put("odd_one_candidate", oddOneCandidate);
                disposition.put("manuscript_disposition", manuscriptDisposition);
                dispositions.add(disposition);
                if (oddOneCandidate) improvementCandidates.add(Map.of(
                    "omitted_cohort", variant,
                    "replication_fraction_delta", replicationDelta,
                    "mean_i_squared_delta", heterogeneityDelta,
                    "manuscript_disposition", manuscriptDisposition));
                comparisonTable.append(action).append('\t')
                    .append(action.equals("leave_one_cohort_out") ? variant : "")
                    .append('\t').append(String.join(",", selected.keySet()))
                    .append('\t').append(value(primaryReplicationMetrics.fraction()))
                    .append('\t').append(value(sensitivityReplication.fraction()))
                    .append('\t').append(value(replicationDelta))
                    .append('\t').append(value(primaryMetaMetrics.meanISquared()))
                    .append('\t').append(value(sensitivityMeta.meanISquared()))
                    .append('\t').append(value(heterogeneityDelta))
                    .append('\t').append(value(primaryMetaMetrics.highISquaredFraction()))
                    .append('\t').append(value(sensitivityMeta.highISquaredFraction()))
                    .append('\t').append(oddOneCandidate)
                    .append('\t').append(manuscriptDisposition).append('\n');
            }
        }
        Files.writeString(context.outputDirectory().resolve(
            "sensitivity-comparison.tsv"), comparisonTable, StandardCharsets.UTF_8);
        outputs.add(0, "sensitivity-comparison.tsv");
        Map<String, Object> decision = new LinkedHashMap<>();
        decision.put("analyses", dispositions);
        decision.put("primary_metrics", metrics(primaryReplicationMetrics,
            primaryMetaMetrics));
        decision.put("odd_one_candidates", improvementCandidates);
        decision.put("automatic_rule",
            "Flag an omission only when its replication fraction increases and its mean I-squared decreases versus the unchanged all-cohort primary analysis");
        decision.put("authorization",
            "Automatic leave-one-cohort-out influence diagnostics are part of the approved deterministic plan; other analyses require assistance.allow_analyses");
        decision.put("manuscript_guidance", improvementCandidates.isEmpty()
            ? "No omission improved both registered criteria; retain all-cohort primary results and report the complete sensitivity table in supplementary material."
            : "Feature flagged omissions in a labeled sensitivity subsection with exact deltas and full-cohort primary results alongside them. Promotion to a primary analysis requires a prespecified rationale or a reviewed study amendment; outcome-selected exclusion is not automatic.");
        decision.put("interpretation",
            "Exploratory influence diagnostics; primary results unchanged. Every cohort omission is retained, and no effect-scale correction, selective exclusion, or new hypothesis is inferred.");
        Files.writeString(context.outputDirectory().resolve("decision.yaml"), new Yaml().dump(decision));
        outputs.add("decision.yaml");
        return new Result(incomplete ? State.WAITING_FOR_INPUT : State.COMPLETE,
            incomplete ? "Some preapproved sensitivity analyses need additional inputs" : "Completed preapproved sensitivity analyses", outputs, decision);
    }

    private Result automaticOrConfiguredEnrichment(Context context,
            String collection, String outputName) throws Exception {
        Map<String, Object> configuration = optionalMap(
            context.operation().get("configuration"), "operation.configuration");
        if (configuration.containsKey("command")) return configured(context);
        Path source = dependencyFile(context, "meta-analysis.tsv");
        if (source == null)
            return waiting(context,
                collection + " enrichment needs a completed association or meta-analysis table",
                Map.of("required_input", "feature-level association results"));

        DelimitedData table = DelimitedData.read(source);
        int feature = table.column("feature_id");
        int pValue = table.column("p_value");
        int status = table.header().contains("status")
            ? table.column("status") : -1;
        int term = table.header().contains("term") ? table.column("term") : -1;
        Map<String, Double> pValues = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            if (status >= 0 && !Set.of("ok", "single_cohort", "tested")
                    .contains(row[status].toLowerCase(Locale.ROOT))) continue;
            if (term >= 0 && !Set.of("pooled", "na", "")
                    .contains(row[term].toLowerCase(Locale.ROOT))) continue;
            if (missing(row[pValue])) continue;
            pValues.putIfAbsent(row[feature], probability(row[pValue], source));
        }
        if (pValues.isEmpty())
            return waiting(context,
                "No eligible feature-level p-values were available for enrichment",
                Map.of("source", source.toString()));
        Map<String, Double> fdr = bh(pValues);
        Map<String, Object> hypothesis = hypothesis(context.study(),
            text(context.operation().get("hypothesis")));
        Map<String, Object> analysis = map(
            hypothesis.get("analysis"), "hypothesis.analysis");
        Map<String, Object> testing = optionalMap(
            analysis.get("multiple_testing"), "analysis.multiple_testing");
        double alpha = number(testing.getOrDefault("alpha", 0.05),
            "multiple_testing.alpha");
        Set<String> selected = new LinkedHashSet<>();
        Path selection = context.outputDirectory().resolve(
            "feature-selection.tsv");
        try (BufferedWriter writer = Files.newBufferedWriter(selection)) {
            writer.write("feature_id\tp_value\tFDR\tselected\n");
            for (Map.Entry<String, Double> entry : pValues.entrySet()) {
                boolean keep = fdr.get(entry.getKey()) <= alpha;
                if (keep) selected.add(entry.getKey());
                writer.write(entry.getKey() + "\t" + entry.getValue()
                    + "\t" + fdr.get(entry.getKey()) + "\t" + keep + "\n");
            }
        }

        Map<String, GeneMapping> mappings = geneMappings(context,
            pValues.keySet());
        boolean symbols = collection.equals("GO");
        Set<String> backgroundGenes = new java.util.TreeSet<>();
        Set<String> selectedGenes = new java.util.TreeSet<>();
        Path mappingOutput = context.outputDirectory().resolve(
            "feature-to-gene.tsv");
        try (BufferedWriter writer = Files.newBufferedWriter(mappingOutput)) {
            writer.write("feature_id\tgene_id\tgene_symbol\tmapped\n");
            for (String id : pValues.keySet()) {
                GeneMapping mapping = mappings.get(id);
                String gene = mapping == null ? null
                    : symbols ? mapping.symbol() : mapping.ensembl();
                boolean mapped = gene != null && !gene.isBlank();
                writer.write(id + "\t" + (mapping == null ? "" : mapping.ensembl())
                    + "\t" + (mapping == null ? "" : mapping.symbol())
                    + "\t" + mapped + "\n");
                if (mapped) {
                    backgroundGenes.add(gene);
                    if (selected.contains(id)) selectedGenes.add(gene);
                }
            }
        }
        if (backgroundGenes.isEmpty())
            return waiting(context,
                "Feature identifiers could not be mapped to the gene namespace required by "
                    + collection,
                Map.of("required_input",
                    "dataset annotation_release (for GENCODE IDs) or files.annotation",
                    "unmapped_features", pValues.size()));

        Path selectedFile = context.outputDirectory().resolve(
            "selected-genes.txt");
        Path backgroundFile = context.outputDirectory().resolve(
            "background-genes.txt");
        Files.writeString(selectedFile, String.join("\n", selectedGenes)
            + (selectedGenes.isEmpty() ? "" : "\n"));
        Files.writeString(backgroundFile, String.join("\n", backgroundGenes)
            + "\n");
        if (selectedGenes.isEmpty()) {
            Map<String, Object> summary = new LinkedHashMap<>();
            summary.put("state", "complete");
            summary.put("collection", collection);
            summary.put("tested_features", pValues.size());
            summary.put("selected_features", selected.size());
            summary.put("selected_genes", 0);
            summary.put("interpretation",
                "No features passed the prespecified FDR threshold; enrichment was not tested");
            Files.writeString(context.outputDirectory().resolve(
                "enrichment-status.yaml"), new Yaml().dump(summary));
            return new Result(State.COMPLETE, "No genes met the prespecified enrichment selection threshold",
                List.of("feature-selection.tsv", "feature-to-gene.tsv",
                    "selected-genes.txt", "background-genes.txt",
                    "enrichment-status.yaml"), Map.of("method_report", Map.of("performed", false, "reason", "No selected genes; enrichment was not tested", "selection_alpha", alpha)));
        }

        String directoryName = collection.equals("GO")
            ? "go-human" : "reactome-human";
        Path database = StudyResourceManager.resourceRoot(
            context.project(), context.study()).resolve("enrichment")
            .resolve(directoryName);
        ByteArrayOutputStream download = new ByteArrayOutputStream();
        StudyPublicEvidence.write(context.outputDirectory().resolve("resource-acquisition-intent.yaml"), Map.of(
            "hypothesis", context.operation().get("hypothesis"), "collection", collection,
            "scientific_rationale", "Test whether the hypothesis-selected genes are overrepresented in " + collection
                + " using the mapped tested genes as background. This is functional context, not independent replication or causal evidence.",
            "selection_source_sha256", FollowupSupport.hash(source), "selected_genes", selectedGenes.size(),
            "background_genes", backgroundGenes.size(), "resource_directory", database.toString(),
            "request_needed", !Files.isRegularFile(database.resolve("manifest.properties"))));
        if (!Files.isRegularFile(database.resolve("manifest.properties"))) {
            if (!StudyResourceManager.automatic(context.study()))
                return waiting(context,
                    collection + " database is missing and automatic downloads are disabled",
                    Map.of("expected_resource", database.toString()));
            Files.createDirectories(database.getParent());
            List<String> downloadArguments = new ArrayList<>(List.of(
                "--enrichment", collection, "--download", database.toString(),
                "--species", "human", "--gene-id-type",
                symbols ? "symbol" : "ensembl"));
            int downloadStatus = EnrichmentCli.run(
                downloadArguments.toArray(String[]::new),
                new PrintStream(download), new PrintStream(download));
            if (downloadStatus != 0)
                throw new IOException(collection + " resource download failed: "
                    + download.toString(StandardCharsets.UTF_8).trim());
        }
        Files.writeString(context.outputDirectory().resolve(
            "resource-download.log"),
            download.toString(StandardCharsets.UTF_8));
        Files.copy(database.resolve("manifest.properties"),
            context.outputDirectory().resolve("reference-manifest.properties"));

        Path output = context.outputDirectory().resolve(outputName);
        ByteArrayOutputStream console = new ByteArrayOutputStream();
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        int result = EnrichmentCli.run(new String[] {
            "--enrichment", collection,
            "--enrichment-db", database.toString(),
            "--input", selectedFile.toString(),
            "--background", backgroundFile.toString(),
            "--out", output.toString(),
            "--gene-id-type", symbols ? "symbol" : "ensembl"
        }, new PrintStream(console), new PrintStream(error));
        Files.writeString(context.outputDirectory().resolve("console.txt"),
            console.toString(StandardCharsets.UTF_8)
                + error.toString(StandardCharsets.UTF_8));
        if (result != 0)
            throw new IOException(collection + " enrichment failed: "
                + error.toString(StandardCharsets.UTF_8).trim());
        return new Result(State.COMPLETE, "Completed local " + collection
            + " enrichment for " + selectedGenes.size() + " selected genes",
            List.of("feature-selection.tsv", "feature-to-gene.tsv",
                "selected-genes.txt", "background-genes.txt", outputName,
                outputName + ".log", outputName + ".features.tsv",
                outputName + ".manifest.properties", "console.txt",
                "resource-download.log", "reference-manifest.properties"), Map.of("method_report",
                    StudyMethods.record(collection + " gene-set overrepresentation", Map.of("collection", collection, "selection_alpha", alpha,
                        "feature_selection", "Benjamini-Hochberg across all eligible tested features", "background", "All mapped tested genes",
                        "feature_mapping", "feature-to-gene.tsv", "selected_gene_count", selectedGenes.size(), "background_gene_count", backgroundGenes.size(),
                        "database_release_and_hashes", "reference-manifest.properties", "test_settings", outputName + ".manifest.properties", "command_log", outputName + ".log"),
                        StudyPublicEvidence.list(configuration.get("references")), configuration.get("rationale"), List.of(collection.equals("GO") ? "Gene Ontology" : "Reactome", "enrichment"))));
    }

    private Result genetic(Context context) throws IOException {
        Map<String, Object> configuration = optionalMap(
            context.operation().get("configuration"), "operation.configuration");
        if (configuration.containsKey("command")) return configured(context);
        String operation = text(context.operation().get("type"));
        if (!operation.equals("variant_scoring") && (!StudyPublicEvidence.list(StudyBioinformatics.config(context).get("genetic_pairs")).isEmpty()
                || operation.equals("mendelian_randomization") && !StudyPublicEvidence.map(StudyBioinformatics.config(context).get("mr")).isEmpty()
                || operation.equals("colocalization") && !StudyPublicEvidence.map(StudyBioinformatics.config(context).get("colocalization")).isEmpty()))
            return new StudyGeneticEvidence().execute(context);
        boolean variantInput = selectedDatasets(context).stream().anyMatch(dataset -> {
            String type = String.valueOf(dataset.get("type"))
                .toLowerCase(Locale.ROOT);
            Map<String, Object> columns = optionalMap(dataset.get("columns"),
                "dataset.columns");
            return Set.of("gwas", "genomics").contains(type)
                || columns.containsKey("variant_id")
                || columns.containsKey("rsid");
        });
        Map<String, Object> details = new LinkedHashMap<>();
        if (operation.equals("variant_scoring") && !variantInput)
            return new Result(State.NOT_APPLICABLE,
                "Variant scoring is not applicable to the declared non-variant omics results",
                List.of(), Map.of("reason", "No variant-level association input"));
        details.put("eligible_variant_summary_input", variantInput);
        details.put("annotation_requirements", annotationRequests(context));
        if (!variantInput) {
            details.put("molecular_followup", List.of(
                "Resolve probe/transcript coordinates, annotation release, and genome build",
                "Identify same-phenotype GWAS studies and nearby loci; nearby lead variants alone are insufficient",
                "Obtain tissue-appropriate molecular-QTL effects for the actual measured feature and full regional outcome-GWAS statistics",
                "Check allele/build/ancestry alignment, instrument validity, and sample overlap before fine mapping, colocalization, or MR"));
            details.put("catalog_search_status", context.dependencyDirectories().values().stream()
                .anyMatch(path -> Files.isRegularFile(path.resolve("source-discovery.yaml"))) ? "see_source_discovery" : "not_run");
            details.put("catalog_search_reason", "Inspect source_discovery evidence and citations. For supported public MR sources use concise analysis.bioinformatics.mr; use genetic_pairs for custom reviewed regional sources. Candidate catalog hits alone cannot identify a valid causal analysis.");
            details.put("ai_advisor_status", "not_configured");
        }
        details.put("default_ld_database", "1000g-phase3 (GRCh37); automatic public GRCh38 MR can select 1000g-30x-grch38 explicitly");
        details.put("ld_download_policy", "automatic_when_variant_inputs_are_eligible");
        details.put("privacy_note",
            "Reference downloads are inbound; no study result is uploaded");
        details.put("required_inputs", switch (operation) {
            case "fine_mapping" -> List.of(
                "regional variant-level summary statistics with beta/se or z",
                "sample size and genome build", "ancestry or LD population",
                "variant-matched LD matrix or extractable LD panel");
            case "colocalization" -> List.of(
                "two traits with regional variant-level summary statistics",
                "harmonized build, alleles, and overlapping variants",
                "fine-mapping effects for both traits");
            case "mendelian_randomization" -> List.of(
                "exposure and outcome variant-level summary statistics",
                "effect alleles and allele frequencies", "ancestry",
                "valid instrument definition and sample-overlap declaration");
            default -> List.of(
                "variant-level association results and annotation policy");
        });
        if (variantInput && !operation.equals("variant_scoring")) {
            ByteArrayOutputStream resourceLog = new ByteArrayOutputStream();
            try {
                StudyResourceManager.LdReference reference =
                    StudyResourceManager.ld(context.project(), context.study(),
                        new PrintStream(resourceLog));
                details.put("ld_reference", reference == null ? "not_installed"
                    : reference.directory().toString());
            } catch (InterruptedException interrupted) {
                Thread.currentThread().interrupt();
                throw new IOException("LD reference download was interrupted",
                    interrupted);
            }
            Files.writeString(context.outputDirectory().resolve(
                "ld-reference-download.log"),
                resourceLog.toString(StandardCharsets.UTF_8));
            details.put("download_log", "ld-reference-download.log");
        }
        return waiting(context,
            operation + " cannot be inferred from the currently declared inputs",
            details, geneticNeeds(context, operation));
    }

    private static List<StudyNeed> geneticNeeds(Context context, String operation) throws IOException {
        List<StudyNeed> needs = new ArrayList<>();
        for (Object raw : annotationRequests(context)) {
            Map<String, Object> request = optionalMap(raw, "annotation request");
            String id = text(request.get("dataset"));
            if (id.isBlank()) {
                needs.add(StudyNeed.review("annotation-review:pooled:" + context.operation().get("hypothesis"),
                    "How should conflicting feature annotations across datasets be reconciled?",
                    "Inspect annotation-summary.yaml and resolve pooled identity, build and interval conflicts before genomic inference.", true));
                continue;
            }
            Map<String, Object> dataset = selectedDatasets(context).stream()
                .filter(item -> id.equals(item.get("id"))).findFirst().orElseThrow();
            Map<String, Object> files = optionalMap(dataset.get("files"), "files");
            if (!files.containsKey("annotation"))
                needs.add(StudyNeed.dataset(id, "files.annotation",
                    "Which local feature annotation file should be used for dataset " + id + "?",
                    "Genomic follow-up needs unambiguous feature coordinates, genome build, and coordinate convention. Supply the mapping source and release; do not infer them from nearby GWAS hits.",
                    "local_file", List.of()));
            else needs.add(StudyNeed.review("annotation-review:" + id,
                "How should missing or ambiguous annotation for dataset " + id + " be resolved?",
                "The supplied annotation leaves unresolved mappings. Inspect annotation-summary.yaml and document corrected mappings, source/release, build, and coordinate convention in an amendment.", true));
        }
        String hypothesis = text(context.operation().get("hypothesis"));
        needs.add(StudyNeed.review("genetic-inputs:" + hypothesis + ":" + operation,
            "Which scientifically relevant local GWAS/QTL inputs and analysis specification should support " + operation + " for " + hypothesis + "?",
            "Document phenotype and molecular feature matching, tissue relevance, source citations, ancestry, genome build, alleles, and sample overlap. Lead SNP proximity alone is insufficient; review request.yaml for operation-specific requirements. Acquisition and adapter configuration remain separate steps.", true));
        return needs;
    }

    private Result prioritization(Context context) throws IOException {
        Path meta = dependencyFile(context, "meta-analysis.tsv");
        Path output = context.outputDirectory().resolve(
            "candidate-prioritization.tsv");
        Map<String, String> evidence = dependencyStates(context);
        Map<String, List<String>> featureEvidence = new LinkedHashMap<>();
        for (Path directory : context.dependencyDirectories().values()) {
            for (String name : List.of("bioinformatics-investigation.yaml", "literature-search.yaml", "genetic-evidence.yaml")) {
                Path path = directory.resolve(name);
                if (!Files.isRegularFile(path)) continue;
                Map<String, Object> summary;
                try (var reader = Files.newBufferedReader(path)) { summary = StudyPublicEvidence.map(new Yaml().load(reader)); }
                for (Object raw : StudyPublicEvidence.list(summary.getOrDefault("evidence", summary.get("pairs")))) {
                    Map<String, Object> row = StudyPublicEvidence.map(raw);
                    String feature = row.get("feature_id") == null
                        ? "" : String.valueOf(row.get("feature_id"));
                    if (feature.isBlank()) continue;
                    String kind = row.get("kind") == null
                        ? "" : String.valueOf(row.get("kind"));
                    if (kind.isBlank()) kind = directory.getFileName().toString();
                    String status = row.get("status") == null
                        ? "unknown" : String.valueOf(row.get("status"));
                    if (status.isBlank()) status = "unknown";
                    featureEvidence.computeIfAbsent(feature, ignored -> new ArrayList<>()).add(kind + ":" + status);
                }
            }
        }
        try (BufferedWriter writer = Files.newBufferedWriter(output)) {
            writer.write("rank\tfeature_id\tbeta\tp_value\tFDR\ti_squared"
                + "\tevidence_layers_complete\tevidence_layers_blocked"
                + "\tfeature_evidence\tinterpretation\n");
            if (meta != null) {
                DelimitedData table = DelimitedData.read(meta);
                int id = table.column("feature_id"), beta = table.column("beta");
                int p = table.column("p_value");
                int i2 = table.header().contains("i_squared")
                    ? table.column("i_squared") : -1;
                int status = table.header().contains("status")
                    ? table.column("status") : -1;
                Map<String, Candidate> candidates = new LinkedHashMap<>();
                Map<String, Double> probabilities = new LinkedHashMap<>();
                for (String[] row : table.rows()) {
                    if (status >= 0 && !Set.of("ok", "single_cohort")
                            .contains(row[status].toLowerCase(Locale.ROOT))) continue;
                    if (missing(row[p])) continue;
                    double probability = probability(row[p], meta);
                    candidates.putIfAbsent(row[id], new Candidate(row[id],
                        row[beta], probability, i2 < 0 ? "NA" : row[i2]));
                    probabilities.putIfAbsent(row[id], probability);
                }
                Map<String, Double> adjusted = bh(probabilities);
                List<Candidate> ranked = candidates.values().stream()
                    .sorted(Comparator.comparingDouble(Candidate::pValue)
                        .thenComparing(Candidate::feature)).toList();
                String complete = evidence.entrySet().stream()
                    .filter(entry -> entry.getValue().equals("complete"))
                    .map(Map.Entry::getKey).sorted().collect(
                        java.util.stream.Collectors.joining(";"));
                String blocked = evidence.entrySet().stream()
                    .filter(entry -> !Set.of("complete", "not_applicable").contains(entry.getValue()))
                    .map(Map.Entry::getKey).sorted().collect(
                        java.util.stream.Collectors.joining(";"));
                int rank = 0;
                for (Candidate candidate : ranked)
                    writer.write(++rank + "\t" + candidate.feature() + "\t"
                        + candidate.beta() + "\t" + candidate.pValue() + "\t"
                        + adjusted.get(candidate.feature()) + "\t"
                        + candidate.iSquared() + "\t" + complete + "\t"
                        + blocked + "\t" + String.join(";", featureEvidence.getOrDefault(candidate.feature(), List.of())) + "\t"
                        + "association-ranked; operation completion is not feature-specific biological support; review feature_evidence and source records"
                        + "\n");
            } else {
                int rank = 0;
                for (StudyBioinformatics.Feature feature : StudyBioinformatics.features(context))
                    writer.write(++rank + "\t" + feature.id() + "\tNA\t" + feature.p() + "\t" + feature.fdr()
                        + "\tNA\t\t\t" + String.join(";", featureEvidence.getOrDefault(feature.id(), List.of()))
                        + "\tprimary-association-ranked exploratory candidate; no meta-analysis was available\n");
            }
        }
        String methodology = "# Candidate prioritization method\n\n"
            + "Candidates are ordered by meta-analysis p-value and retain FDR, "
            + "effect, heterogeneity, and evidence-layer availability as separate "
            + "columns. No opaque composite score is used. Association ranking "
            + "and feature-specific contextual evidence are kept separate from operation-level availability. A protein structure or keyword literature match is not target validation. Without meta-analysis, bounded FDR-selected annotated primary findings are listed with their original full-family FDR. "
            + "alone is not causal evidence, target validation, safety evidence, "
            + "or a clinical recommendation.\n";
        Files.writeString(context.outputDirectory().resolve(
            "prioritization-method.md"), methodology);
        return Result.complete("Generated a transparent candidate evidence table",
            List.of("candidate-prioritization.tsv", "prioritization-method.md"));
    }

    private Result configured(Context context) throws IOException {
        Map<String, Object> configuration = optionalMap(
            context.operation().get("configuration"), "operation.configuration");
        if (!configuration.containsKey("command"))
            return waiting(context, text(context.operation().get("type"))
                + " requires analysis.operations configuration and reference resources");
        String command = text(configuration.get("command"));
        if (!CUSTOM_COMMANDS.contains(command) || command.equals("study"))
            throw new IllegalArgumentException(
                "Command is not allowed for configured study execution: " + command);
        List<String> configuredArguments = strings(configuration.get("arguments"));
        if ("offline".equals(optionalMap(context.study().get("resources"), "resources").get("network"))
                && (!Set.of("susie", "coloc", "mr-estimate", "differential", "ewas-regions", "multiple-test", "multiple-impute", "mi-pool").contains(command)
                    || configuredArguments.stream().anyMatch(argument -> argument.contains("://")
                        || argument.startsWith("--download") || argument.startsWith("--config"))))
            return waiting(context, "Offline policy prevents a potentially network-enabled configured adapter",
                Map.of("connection_status", "offline_policy", "command", command,
                    "required_input", "Use an explicitly local supported adapter or repeat in an environment with permitted network access"));
        boolean needsLd = Boolean.TRUE.equals(configuration.get("requires_ld"))
            || configuredArguments.stream().anyMatch(argument ->
                argument.contains("@ld_database@")
                    || argument.contains("@ld_population@"));
        StudyResourceManager.LdReference ld = null;
        if (needsLd) {
            ByteArrayOutputStream resourceLog = new ByteArrayOutputStream();
            try {
                ld = StudyResourceManager.ld(context.project(), context.study(),
                    new PrintStream(resourceLog));
            } catch (InterruptedException interrupted) {
                Thread.currentThread().interrupt();
                throw new IOException("LD reference download was interrupted",
                    interrupted);
            }
            Files.writeString(context.outputDirectory().resolve(
                "ld-reference-download.log"),
                resourceLog.toString(StandardCharsets.UTF_8));
            if (ld == null)
                return waiting(context,
                    "An LD reference is required and automatic downloads are disabled",
                    Map.of("default_database", "1000g-phase3",
                        "resource_policy", "manual_or_denied"));
            if ("auto".equalsIgnoreCase(ld.population())
                    && configuredArguments.stream().anyMatch(argument ->
                        argument.contains("@ld_population@")))
                return waiting(context,
                    "The 1000 Genomes reference is installed, but ancestry must be declared before LD analysis",
                    Map.of("installed_database", ld.database(),
                        "available_populations", List.of("AFR", "AMR", "EAS", "EUR", "SAS"),
                        "required_field", "resources.ld.population"));
        }
        List<String> arguments = new ArrayList<>();
        if (!command.equals("enrichment")) arguments.add(command);
        arguments.add("--no-config");
        for (String argument : configuredArguments) {
            String expanded = argument.replace("@project@",
                context.project().toString()).replace("@output@",
                context.outputDirectory().toString());
            if (expanded.contains("@upstream@"))
                expanded = expanded.replace("@upstream@", upstream(context));
            if (ld != null) expanded = expanded
                .replace("@ld_database@", ld.directory().toString())
                .replace("@ld_population@", ld.population());
            arguments.add(expanded);
        }
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();
        int status = JLinAlgCli.run(arguments.toArray(String[]::new),
            new PrintStream(stdout), new PrintStream(stderr));
        Files.writeString(context.outputDirectory().resolve("console.txt"),
            stdout.toString(StandardCharsets.UTF_8)
                + stderr.toString(StandardCharsets.UTF_8));
        if (status != 0)
            throw new IOException(command + " failed: "
                + stderr.toString(StandardCharsets.UTF_8).trim());
        List<String> outputs = strings(configuration.get("outputs"));
        outputs.add("console.txt");
        if (needsLd) outputs.add("ld-reference-download.log");
        return Result.complete("Executed configured " + command + " operation",
            outputs);
    }

    private static String upstream(Context context) {
        if (context.dependencyDirectories().size() != 1)
            throw new IllegalArgumentException(
                "@upstream@ requires exactly one immediate dependency; found "
                    + context.dependencyDirectories().keySet());
        return context.dependencyDirectories().values().iterator().next()
            .toString();
    }

    private Result diagnostics(Context context) throws IOException {
        Map<String, Object> document = new LinkedHashMap<>();
        document.put("schema_version", 1);
        document.put("hypothesis", context.operation().get("hypothesis"));
        document.put("status", "standard_diagnostics_complete");
        List<Object> observations = new ArrayList<>();
        List<String> outputs = new ArrayList<>();
        Path replication = dependencyFile(context, "replication.tsv");
        if (replication != null) {
            DelimitedData table = DelimitedData.read(replication);
            int feature = table.column("feature_id");
            int discoveryDataset = table.column("discovery_dataset");
            int replicationDataset = table.column("replication_dataset");
            int discoveryBeta = table.column("discovery_beta");
            int replicationBeta = table.column("replication_beta");
            int selected = table.column("discovery_selected");
            int replicated = table.column("replicated");
            Map<String, Moment> moments = new LinkedHashMap<>();
            Set<String> selectedFeatures = new LinkedHashSet<>();
            Set<String> replicatedFeatures = new LinkedHashSet<>();
            for (String[] row : table.rows()) {
                String pair = row[discoveryDataset] + " -> "
                    + row[replicationDataset];
                Moment moment = moments.computeIfAbsent(pair,
                    ignored -> new Moment());
                boolean selectedValue = Boolean.parseBoolean(row[selected]);
                boolean replicatedValue = Boolean.parseBoolean(row[replicated]);
                moment.add(finite(row[discoveryBeta], "discovery_beta", replication),
                    finite(row[replicationBeta], "replication_beta", replication),
                    selectedValue, replicatedValue);
                if (selectedValue) selectedFeatures.add(row[feature]);
                if (replicatedValue) replicatedFeatures.add(row[feature]);
            }
            Path correlations = context.outputDirectory().resolve(
                "effect-correlations.tsv");
            List<Object> comparisons = new ArrayList<>();
            try (BufferedWriter writer = Files.newBufferedWriter(correlations)) {
                writer.write("comparison\tn\tpearson_r\tols_intercept"
                    + "\tols_slope\tsd_ratio\tsame_sign_fraction"
                    + "\tdiscovery_selected\treplicated"
                    + "\treplication_fraction\n");
                for (Map.Entry<String, Moment> entry : moments.entrySet()) {
                    Moment moment = entry.getValue();
                    double correlation = entry.getValue().correlation();
                    double slope = moment.slope();
                    double replicationFraction = moment.replicationFraction();
                    Map<String, Object> summary = new LinkedHashMap<>();
                    summary.put("comparison", entry.getKey());
                    summary.put("n", moment.n);
                    summary.put("pearson_r", correlation);
                    summary.put("ols_intercept", moment.intercept());
                    summary.put("ols_slope", slope);
                    summary.put("sd_ratio", moment.standardDeviationRatio());
                    summary.put("same_sign_fraction", moment.sameSignFraction());
                    summary.put("discovery_selected", moment.selected);
                    summary.put("replicated", moment.replicated);
                    summary.put("replication_fraction", replicationFraction);
                    comparisons.add(summary);
                    writer.write(entry.getKey() + "\t" + moment.n
                        + "\t" + value(correlation) + "\t"
                        + value(moment.intercept()) + "\t"
                        + value(slope) + "\t"
                        + value(moment.standardDeviationRatio()) + "\t"
                        + value(moment.sameSignFraction()) + "\t"
                        + moment.selected + "\t" + moment.replicated + "\t"
                        + value(replicationFraction) + "\n");
                    if (moment.n >= 100
                            && Double.isFinite(correlation)
                            && Math.abs(correlation) < 0.2) {
                        Map<String, Object> finding = new LinkedHashMap<>();
                        finding.put("type", "low_cross_cohort_effect_correlation");
                        finding.put("comparison", entry.getKey());
                        finding.put("pearson_r", correlation);
                        finding.put("interpretation",
                            "May reflect weak shared signal, phenotype or tissue heterogeneity, or processing artifacts; inspect plots and cohort metadata");
                        observations.add(finding);
                    }
                    if (moment.selected >= 100
                            && Double.isFinite(replicationFraction)
                            && replicationFraction < 0.01) {
                        Map<String, Object> finding = new LinkedHashMap<>();
                        finding.put("type", "near_zero_replication_rate");
                        finding.put("comparison", entry.getKey());
                        finding.put("discovery_selected", moment.selected);
                        finding.put("replicated", moment.replicated);
                        finding.put("replication_fraction",
                            replicationFraction);
                        finding.put("pearson_r", correlation);
                        finding.put("interpretation",
                            Math.abs(correlation) >= 0.8
                                ? "Effects are strongly correlated but almost no discovery-selected features replicate; review uncertainty, sample size, phenotype/tissue definitions, scale, and multiple-testing inputs"
                                : "Almost no discovery-selected features replicate; review power, phenotype/tissue definitions, harmonization, and processing");
                        observations.add(finding);
                    }
                    double absoluteSlope = Math.abs(slope);
                    if (moment.n >= 100 && Math.abs(correlation) >= 0.5
                            && Double.isFinite(absoluteSlope)
                            && (absoluteSlope > 4.0 || absoluteSlope < 0.25)) {
                        observations.add(Map.of(
                            "type", "cross_cohort_effect_scale_mismatch",
                            "comparison", entry.getKey(),
                            "pearson_r", correlation,
                            "ols_slope", slope,
                            "interpretation",
                            "Effects track one another but differ substantially in scale; verify phenotype scaling, units, transformations, and model parameterization before pooling"));
                    }
                }
            }
            document.put("discovery_replication_summary", Map.of(
                "discovery_selected_features", selectedFeatures.size(),
                "replicated_in_at_least_one", replicatedFeatures.size(),
                "replication_fraction_any", selectedFeatures.isEmpty()
                    ? "NA" : replicatedFeatures.size()
                        / (double) selectedFeatures.size()));
            document.put("discovery_replication", comparisons);
            outputs.add("effect-correlations.tsv");
        }
        Path meta = dependencyFile(context, "meta-analysis.tsv");
        if (meta != null) {
            DelimitedData table = DelimitedData.read(meta);
            int p = table.column("p_value");
            int feature = table.column("feature_id");
            int i2 = table.header().contains("i_squared")
                ? table.column("i_squared") : -1;
            Map<String, Double> probabilities = new LinkedHashMap<>();
            int highHeterogeneity = 0;
            for (String[] row : table.rows()) {
                if (!missing(row[p])) probabilities.putIfAbsent(row[feature],
                    probability(row[p], meta));
                if (i2 >= 0 && !missing(row[i2])
                        && finite(row[i2], "i_squared", meta) >= 50.0)
                    highHeterogeneity++;
            }
            long significant = bh(probabilities).values().stream()
                .filter(value -> value <= 0.05).count();
            document.put("meta_analysis", Map.of(
                "tested_features", probabilities.size(),
                "BH_FDR_0_05_features", significant,
                "I2_at_least_50_rows", highHeterogeneity));
            if (!probabilities.isEmpty()
                    && highHeterogeneity > probabilities.size() / 10) {
                observations.add(Map.of(
                    "type", "widespread_meta_analysis_heterogeneity",
                    "I2_at_least_50_rows", highHeterogeneity,
                    "interpretation",
                    "Review cohort, tissue, sex, ancestry, phenotype, and processing differences and consider random-effects sensitivity analysis"));
            }
        }
        document.put("unusual_observations", observations);
        List<Object> nextSteps = new ArrayList<>();
        {
            List<?> requests = annotationRequests(context);
            if (!requests.isEmpty()) nextSteps.add(Map.of(
                "id", "resolve-feature-annotation", "status", "input_needed_for_genomic_followup",
                "trigger", "Missing or ambiguous feature annotation",
                "required_inputs", requests, "source", "annotation-summary.yaml"));
        }
        for (Object observation : observations) {
            Map<String, Object> finding = optionalMap(observation, "observation");
            String type = text(finding.get("type"));
            String action = switch (type) {
                case "cross_cohort_effect_scale_mismatch" -> "Verify phenotype units, transformations, outcome coding, and model parameterization; define a common estimand before pooling. Do not rescale using the fitted between-cohort slope.";
                case "near_zero_replication_rate", "low_cross_cohort_effect_correlation" -> "Audit standard errors, sample sizes, selection thresholds, tissue, and processing. Review cohort effect plots and assess replication power.";
                case "widespread_meta_analysis_heterogeneity" -> "After confirming comparable effect definitions, compare fixed and random effects and perform leave-one-cohort-out analysis; report changes in direction, uncertainty, FDR, and heterogeneity.";
                default -> "Review the diagnostic observation and propose a documented sensitivity analysis.";
            };
            Map<String, Object> next = new LinkedHashMap<>();
            next.put("id", "review-" + (nextSteps.size() + 1));
            next.put("status", "proposed");
            next.put("trigger", finding);
            next.put("suggested_action", action);
            next.put("source", "diagnostics.yaml");
            next.put("required_inputs", List.of("cohort-level effect estimates and standard errors",
                "phenotype transformations, units, covariates, and cohort relationship metadata"));
            next.put("execution", "Record an analysis_requests amendment with rationale and multiplicity plan; approve a new plan before execution");
            nextSteps.add(next);
        }
        document.put("next_steps", nextSteps);
        Files.writeString(context.outputDirectory().resolve("next-steps.yaml"),
            new Yaml().dump(Map.of("schema_version", 1, "hypothesis", context.operation().get("hypothesis"),
                "producer", "jlinalg_deterministic_rules", "proposals", nextSteps)));
        outputs.add("next-steps.yaml");
        document.put("visual_review_boundary",
            "Aggregate screening cannot reliably classify streaks, clusters, or other plot morphology; record them with study observe and a reviewed amendment");
        document.put("note",
            "Automated flags are prompts for scientific review, not determinations of artifact or biology");
        String text = new Yaml().dump(document);
        Files.writeString(context.outputDirectory().resolve("diagnostics.yaml"),
            text, StandardCharsets.UTF_8);
        outputs.add("diagnostics.yaml");
        List<StudyNeed> needs = new ArrayList<>();
        for (Object raw : observations) {
            Map<String, Object> finding = optionalMap(raw, "observation");
            String kind = text(finding.get("type"));
            String hypothesis = text(context.operation().get("hypothesis"));
            String comparison = String.valueOf(finding.getOrDefault("comparison", "all cohorts"));
            needs.add(StudyNeed.review("diagnostic:" + hypothesis + ":" + kind + ":" + comparison,
                "What documented metadata or sensitivity analysis should be reviewed for " + kind + " in " + hypothesis + " (" + comparison + ")?",
                "See diagnostics.yaml and next-steps.yaml. The flag does not establish artifact or biology; examine units, transformations, model definitions, and study differences before proposing a new analysis. Record exploratory status and multiplicity considerations.", false));
        }
        return new Result(State.COMPLETE, "Completed registered deterministic diagnostics", outputs,
            needs.isEmpty() ? Map.of() : Map.of("needs", needs.stream().map(StudyNeed::document).toList()));
    }

    private Result report(Context context) throws IOException {
        String hypothesis = text(context.operation().get("hypothesis"));
        StringBuilder report = new StringBuilder("# Hypothesis ")
            .append(hypothesis).append(" evidence report\n\n")
            .append("This report separates completed evidence from unavailable "
                + "evidence. Every statement resolves to an operation manifest "
                + "with input and output hashes.\n\n")
            .append("## Evidence ledger\n\n");
        int complete = 0, blocked = 0, failed = 0, skipped = 0;
        for (Map.Entry<String, Path> dependency
                : context.dependencyDirectories().entrySet()) {
            Path manifestPath = dependency.getValue().resolve("manifest.yaml");
            Map<String, Object> manifest = Files.isRegularFile(manifestPath)
                ? loadYaml(manifestPath) : Map.of();
            String state = String.valueOf(manifest.getOrDefault(
                "state", "unknown"));
            if (state.equals("complete")) complete++;
            else if (state.equals("not_applicable")) skipped++;
            else if (state.equals("failed")) failed++;
            else blocked++;
            report.append("- `").append(manifest.getOrDefault(
                    "operation_type", dependency.getKey()))
                .append("` — **").append(state).append("**: ")
                .append(manifest.getOrDefault("summary", "No summary"))
                .append(". Result ID: `")
                .append(manifest.getOrDefault("run_id", "unknown")).append('/')
                .append(dependency.getKey()).append("`; manifest: `")
                .append(portable(context.project(), manifestPath)).append("`\n");
        }
        report.append("\n## Readiness\n\n")
            .append("- Completed evidence operations: ").append(complete).append("\n")
            .append("- Blocked or waiting evidence operations: ").append(blocked).append("\n")
            .append("- Failed evidence operations: ").append(failed).append("\n\n")
            .append("- Inapplicable operations: ").append(skipped).append("\n\n")
            .append("A blocked MR, fine-mapping, colocalization, or literature "
                + "branch is not negative evidence. It means the required input "
                + "or approved source was unavailable.\n\n")
            .append("## Reuse\n\nUse `index/LLM_CONTEXT.md` for a compact, "
                + "machine-readable drafting handoff and `index/index.json` for "
                + "artifact-level retrieval. Verify the project before relying "
                + "on either.\n");
        Files.writeString(context.outputDirectory().resolve("report.md"),
            report, StandardCharsets.UTF_8);
        String publication = "# Publication evidence package for " + hypothesis
            + "\n\n" + report.substring(report.indexOf("## Evidence ledger"))
            + "\nThis is an evidence inventory, not manuscript prose. Effect "
            + "definitions, cohort comparability, model assumptions, and any "
            + "post hoc amendments must be reviewed by the scientific team.\n";
        Files.writeString(context.outputDirectory().resolve(
            "publication-evidence.md"), publication, StandardCharsets.UTF_8);
        List<String> outputs = new ArrayList<>(List.of("report.md",
            "publication-evidence.md"));
        outputs.addAll(StudyPublication.write(context));
        if (objectives(context.study()).stream().anyMatch(
                value -> value.contains("drug"))) {
            String drug = "# Drug-discovery evidence package for " + hypothesis
                + "\n\nThe candidate table, when present, ranks association "
                + "evidence transparently. Genetic causality, colocalization, "
                + "tractability, safety, tissue relevance, and literature are "
                + "separate evidence layers; missing layers are shown as blockers "
                + "and must not be interpreted as support.\n\n"
                + report.substring(report.indexOf("## Evidence ledger"));
            Files.writeString(context.outputDirectory().resolve(
                "drug-discovery-evidence.md"), drug, StandardCharsets.UTF_8);
            outputs.add("drug-discovery-evidence.md");
        }
        return Result.complete("Generated traceable hypothesis report",
            outputs);
    }

    private Map<String, GeneMapping> geneMappings(Context context,
            Set<String> features) throws IOException, InterruptedException {
        Map<String, GeneMapping> resolved = new LinkedHashMap<>();
        Set<String> excluded = new LinkedHashSet<>();
        for (Path assessed : dependencyFiles(context, "feature-annotation.tsv")) {
            DelimitedData table = DelimitedData.read(assessed);
            int id = table.column("feature_id"), gene = table.column("gene_id");
            int symbol = table.column("gene_symbol"), status = table.column("mapping_status");
            for (String[] row : table.rows()) {
                if (!features.contains(row[id])) continue;
                if (!Set.of("located", "missing_location").contains(row[status])) {
                    excluded.add(row[id]);
                    continue;
                }
                GeneMapping mapping = new GeneMapping(row[gene].replaceFirst("\\.\\d+$", ""), row[symbol]);
                GeneMapping previous = resolved.putIfAbsent(row[id], mapping);
                if (previous != null && !previous.equals(mapping)) excluded.add(row[id]);
            }
        }
        excluded.forEach(resolved::remove);
        return resolved;
    }

    private static List<Object> annotationRequests(Context context) throws IOException {
        List<Object> requests = new ArrayList<>();
        for (Path path : dependencyFiles(context, "annotation-summary.yaml")) {
            Object values = loadYaml(path).get("input_requests");
            if (values instanceof List<?> list) requests.addAll(list);
        }
        return requests;
    }

    private static List<Path> dependencyFiles(Context context, String name) throws IOException {
        List<Path> files = new ArrayList<>();
        for (Path directory : context.dependencyDirectories().values())
            try (var paths = Files.walk(directory)) {
                files.addAll(paths.filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().equals(name)).sorted().toList());
            }
        return files;
    }

    private static Map<String, Double> bh(Map<String, Double> pValues) {
        List<Map.Entry<String, Double>> sorted = pValues.entrySet().stream()
            .sorted(Map.Entry.comparingByValue()).toList();
        Map<String, Double> result = new LinkedHashMap<>();
        double previous = 1.0;
        for (int index = sorted.size() - 1; index >= 0; index--) {
            Map.Entry<String, Double> entry = sorted.get(index);
            previous = Math.min(previous, Math.min(1.0,
                entry.getValue() * sorted.size() / (index + 1.0)));
            result.put(entry.getKey(), previous);
        }
        return result;
    }

    private static Path dependencyFile(Context context, String name)
            throws IOException {
        String producer = switch (name) {
            case "replication.tsv" -> "discovery_replication";
            case "meta-analysis.tsv" -> "meta_analysis";
            default -> throw new IllegalArgumentException("No registered dependency artifact role for " + name);
        };
        List<Path> candidates = new ArrayList<>();
        for (Path directory : context.dependencyDirectories().values()) {
            Path manifestPath = directory.resolve("manifest.yaml");
            if (!Files.isRegularFile(manifestPath)) continue;
            Map<String, Object> manifest = loadYaml(manifestPath);
            if (!producer.equals(manifest.get("operation_type")) || !"complete".equals(manifest.get("state"))) continue;
            Map<String, Object> outputs = optionalMap(manifest.get("output_sha256"), "output_sha256");
            if (outputs.containsKey(name) && Files.isRegularFile(directory.resolve(name))) {
                candidates.add(directory.resolve(name));
                continue;
            }
            for (String output : outputs.keySet()) {
                Path artifact = directory.resolve(output).normalize();
                if (!artifact.startsWith(directory)) throw new IOException("Unsafe dependency artifact: " + output);
                if (artifact.getFileName().toString().equals(name) && Files.isRegularFile(artifact)) candidates.add(artifact);
            }
        }
        if (candidates.size() > 1) throw new IOException("Ambiguous " + producer + " dependency artifacts: " + candidates);
        return candidates.isEmpty() ? null : candidates.get(0);
    }

    private static List<Map<String, Object>> selectedDatasets(Context context) {
        Map<String, Map<String, Object>> all = datasets(context.study());
        List<Map<String, Object>> result = new ArrayList<>();
        for (String id : strings(context.operation().get("datasets"))) {
            Map<String, Object> dataset = all.get(id);
            if (dataset != null) result.add(dataset);
        }
        return result;
    }

    private static Map<String, String> dependencyStates(Context context)
            throws IOException {
        Map<String, String> result = new LinkedHashMap<>();
        for (Map.Entry<String, Path> dependency
                : context.dependencyDirectories().entrySet()) {
            Path manifest = dependency.getValue().resolve("manifest.yaml");
            if (!Files.isRegularFile(manifest)) continue;
            Map<String, Object> document = loadYaml(manifest);
            String name = String.valueOf(document.getOrDefault(
                "operation_type", dependency.getKey()));
            if (result.containsKey(name)) name += ":" + dependency.getKey();
            result.put(name, String.valueOf(document.getOrDefault(
                "state", "unknown")));
        }
        return result;
    }

    private static Set<String> objectives(Map<String, Object> study) {
        Map<String, Object> project = map(study.get("project"), "project");
        Object value = project.get("objective");
        Set<String> result = new LinkedHashSet<>();
        if (value instanceof Collection<?> values)
            values.forEach(item -> result.add(String.valueOf(item)
                .toLowerCase(Locale.ROOT).replace('-', '_').replace(' ', '_')));
        else if (value != null) result.add(String.valueOf(value)
            .toLowerCase(Locale.ROOT).replace('-', '_').replace(' ', '_'));
        return result;
    }

    @SuppressWarnings("unchecked")
    private static Map<String, Object> loadYaml(Path path) throws IOException {
        Object value = new Yaml().load(Files.readString(path));
        if (!(value instanceof Map<?, ?> raw))
            throw new IOException("Expected YAML mapping: " + path);
        Map<String, Object> result = new LinkedHashMap<>();
        raw.forEach((key, item) -> result.put(String.valueOf(key), item));
        return result;
    }

    private Result waiting(Context context, String message) throws IOException {
        return waiting(context, message, Map.of());
    }

    private Result waiting(Context context, String message,
            Map<String, Object> details) throws IOException {
        return waiting(context, message, details, List.of());
    }

    private Result waiting(Context context, String message,
            Map<String, Object> details, List<StudyNeed> needs) throws IOException {
        Map<String, Object> request = new LinkedHashMap<>();
        request.put("schema_version", 1);
        request.put("operation", context.operation().get("type"));
        request.put("hypothesis", context.operation().get("hypothesis"));
        request.put("state", "waiting_for_input");
        request.put("reason", message);
        request.putAll(details);
        if (!needs.isEmpty()) request.put("needs", needs.stream().map(StudyNeed::document).toList());
        String content = new Yaml().dump(request);
        Files.writeString(context.outputDirectory().resolve("request.yaml"),
            content, StandardCharsets.UTF_8);
        Map<String, Object> resultDetails = new LinkedHashMap<>(Map.of("request", "request.yaml", "requirements", details));
        if (!needs.isEmpty()) resultDetails.put("needs", needs.stream().map(StudyNeed::document).toList());
        return new Result(State.WAITING_FOR_INPUT, message, List.of("request.yaml"), resultDetails);
    }

    private static Standardization standardize(Path input,
            Map<String, Object> columns, Path output, boolean excludeInvalid)
            throws IOException {
        DelimitedData table = DelimitedData.readResult(input);
        String idName = configured(columns, "feature_id",
            table.header(), "feature_id", "id", "probe_id", "variant",
            "__row_name__");
        String betaName = configured(columns, "effect",
            table.header(), "beta", "effect", "estimate");
        String seName = configured(columns, "standard_error",
            table.header(), "se", "standard_error", "stderr");
        String pName = configured(columns, "p_value",
            table.header(), "p_value", "p", "pvalue");
        String statusName = columns.get("status") == null ? null
            : text(columns.get("status"));
        int id = table.column(idName), beta = table.column(betaName);
        int se = table.column(seName), p = table.column(pName);
        int status = statusName == null || !table.header().contains(statusName)
            ? -1 : table.column(statusName);
        Set<String> seen = new LinkedHashSet<>();
        Set<String> eligible = new LinkedHashSet<>();
        Map<String, Object> validation = optionalMap(
            columns.get("validation"), "columns.validation");
        List<Integer> finiteColumns = columnIndexes(table,
            validation.get("finite"), "validation.finite");
        List<Integer> positiveColumns = columnIndexes(table,
            validation.get("positive"), "validation.positive");
        List<Integer> nonnegativeColumns = columnIndexes(table,
            validation.get("nonnegative"), "validation.nonnegative");
        List<Integer> probabilityColumns = columnIndexes(table,
            validation.get("probability"), "validation.probability");
        int inputRows = 0, excludedInvalid = 0, excludedStatus = 0;
        try (BufferedWriter writer = Files.newBufferedWriter(output)) {
            writer.write("feature_id\tbeta\tse\tp_value\n");
            for (String[] row : table.rows()) {
                inputRows++;
                String feature = row[id].trim();
                if (feature.isEmpty() || !seen.add(feature))
                    throw new IOException("Blank or duplicate feature ID in " + input
                        + ": " + feature);
                if (status >= 0 && !Set.of("ok", "tested", "complete")
                        .contains(row[status].toLowerCase(Locale.ROOT))) {
                    excludedStatus++;
                    continue;
                }
                try {
                    validate(row, finiteColumns, positiveColumns,
                        nonnegativeColumns, probabilityColumns, table, input);
                    double effect = finite(row[beta], "effect", input);
                    double error = finite(row[se], "standard error", input);
                    if (!(error > 0.0))
                        throw new IOException(
                            "Standard errors must be positive in " + input);
                    String probability = missing(row[p]) ? "NA"
                        : Double.toString(probability(row[p], input));
                    writer.write(feature + "\t" + effect + "\t" + error
                        + "\t" + probability + "\n");
                    eligible.add(feature);
                } catch (IOException invalid) {
                    if (!excludeInvalid) throw invalid;
                    excludedInvalid++;
                }
            }
        }
        if (eligible.isEmpty())
            throw new IOException("No analyzable association rows in " + input);
        return new Standardization(inputRows, excludedInvalid, excludedStatus,
            Set.copyOf(seen), Set.copyOf(eligible));
    }

    private static List<Integer> columnIndexes(DelimitedData table,
            Object names, String label) {
        List<Integer> result = new ArrayList<>();
        for (String name : strings(names)) {
            if (!table.header().contains(name))
                throw new IllegalArgumentException(label
                    + " column is absent: " + name);
            result.add(table.column(name));
        }
        return result;
    }

    private static void validate(String[] row, List<Integer> finiteColumns,
            List<Integer> positiveColumns, List<Integer> nonnegativeColumns,
            List<Integer> probabilityColumns, DelimitedData table, Path input)
            throws IOException {
        for (int column : finiteColumns)
            finite(row[column], table.header().get(column), input);
        for (int column : positiveColumns)
            if (!(finite(row[column], table.header().get(column), input) > 0.0))
                throw new IOException(table.header().get(column)
                    + " must be positive in " + input);
        for (int column : nonnegativeColumns)
            if (finite(row[column], table.header().get(column), input) < 0.0)
                throw new IOException(table.header().get(column)
                    + " must be nonnegative in " + input);
        for (int column : probabilityColumns)
            probability(row[column], input);
    }

    private static void retainFeatures(Path path, Set<String> retained)
            throws IOException {
        DelimitedData table = DelimitedData.read(path);
        int id = table.column("feature_id");
        try (BufferedWriter writer = Files.newBufferedWriter(path)) {
            writer.write(String.join("\t", table.header()));
            writer.write('\n');
            for (String[] row : table.rows())
                if (retained.contains(row[id])) {
                    writer.write(String.join("\t", row));
                    writer.write('\n');
                }
        }
    }

    private static Map<String, Estimate> estimates(Path path) throws IOException {
        DelimitedData table = DelimitedData.read(path);
        int id = table.column("feature_id"), beta = table.column("beta");
        int se = table.column("se"), p = table.column("p_value");
        Map<String, Estimate> result = new LinkedHashMap<>();
        for (String[] row : table.rows())
            result.put(row[id], new Estimate(
                finite(row[beta], "beta", path),
                finite(row[se], "se", path),
                missing(row[p]) ? Double.NaN : probability(row[p], path)));
        return result;
    }

    private static Map<String, Double> adjusted(
            Map<String, Estimate> estimates, String method) {
        Map<String, Double> result = new LinkedHashMap<>();
        estimates.forEach((id, estimate) -> result.put(id, estimate.p));
        if (method.equals("NONE")) return result;
        List<Map.Entry<String, Estimate>> finite = estimates.entrySet().stream()
            .filter(entry -> Double.isFinite(entry.getValue().p))
            .sorted(java.util.Comparator.comparingDouble(
                entry -> entry.getValue().p)).toList();
        double previous = 1.0;
        for (int index = finite.size() - 1; index >= 0; index--) {
            Map.Entry<String, Estimate> entry = finite.get(index);
            double value = Math.min(previous,
                entry.getValue().p * finite.size() / (index + 1.0));
            previous = Math.min(1.0, value);
            result.put(entry.getKey(), previous);
        }
        return result;
    }

    private static Map<String, Path> results(Context context) throws IOException {
        Map<String, Path> result = new LinkedHashMap<>();
        for (Path directory : context.dependencyDirectories().values()) {
            Path standard = directory.resolve("standard-results");
            if (!Files.isDirectory(standard)) continue;
            try (var paths = Files.list(standard)) {
                for (Path path : paths.filter(Files::isRegularFile)
                        .sorted().toList()) {
                    String name = path.getFileName().toString();
                    if (name.endsWith(".tsv"))
                        result.put(name.substring(0, name.length() - 4), path);
                }
            }
        }
        return result;
    }

    private static Map<String, Map<String, Object>> datasets(
            Map<String, Object> study) {
        Map<String, Map<String, Object>> result = new LinkedHashMap<>();
        for (Object value : list(study.get("datasets"), "datasets")) {
            Map<String, Object> dataset = map(value, "dataset");
            result.put(text(dataset.get("id")), dataset);
        }
        return result;
    }

    private static String cohort(Map<String, Object> dataset) {
        if (dataset == null || dataset.get("cohort") == null
                || String.valueOf(dataset.get("cohort")).isBlank())
            return "unspecified";
        return String.valueOf(dataset.get("cohort"));
    }

    private static String fileToken(String value) {
        String token = value.toLowerCase(Locale.ROOT)
            .replaceAll("[^a-z0-9._-]+", "-")
            .replaceAll("^-+|-+$", "");
        return token.isEmpty() ? "dataset" : token;
    }

    private static double[] commonBounds(double[] left, double[] right) {
        double minimum = Double.POSITIVE_INFINITY;
        double maximum = Double.NEGATIVE_INFINITY;
        for (double value : left) {
            minimum = Math.min(minimum, value);
            maximum = Math.max(maximum, value);
        }
        for (double value : right) {
            minimum = Math.min(minimum, value);
            maximum = Math.max(maximum, value);
        }
        minimum = Math.min(minimum, 0.0);
        maximum = Math.max(maximum, 0.0);
        double spread = maximum - minimum;
        if (!(spread > 0.0)) spread = Math.max(1.0, Math.abs(maximum) * 0.2);
        double padding = spread * 0.06;
        return new double[] {minimum - padding, maximum + padding};
    }

    private static ReplicationMetrics replicationMetrics(Path path)
            throws IOException {
        if (path == null || !Files.isRegularFile(path))
            return new ReplicationMetrics(0, 0, Double.NaN);
        DelimitedData table = DelimitedData.read(path);
        int selectedColumn = table.column("discovery_selected");
        int replicatedColumn = table.column("replicated");
        int selected = 0, replicated = 0;
        for (String[] row : table.rows()) {
            if (Boolean.parseBoolean(row[selectedColumn])) selected++;
            if (Boolean.parseBoolean(row[replicatedColumn])) replicated++;
        }
        return new ReplicationMetrics(selected, replicated,
            selected == 0 ? Double.NaN : replicated / (double) selected);
    }

    private static MetaMetrics metaMetrics(Path path) throws IOException {
        if (path == null || !Files.isRegularFile(path))
            return new MetaMetrics(0, 0, Double.NaN, Double.NaN);
        DelimitedData table = DelimitedData.read(path);
        int feature = table.column("feature_id");
        int p = table.column("p_value");
        int iSquared = table.header().contains("i_squared")
            ? table.column("i_squared") : -1;
        Map<String, Double> pValues = new LinkedHashMap<>();
        double iSquaredSum = 0.0;
        int iSquaredCount = 0, highISquared = 0;
        for (String[] row : table.rows()) {
            if (!missing(row[p])) pValues.putIfAbsent(row[feature],
                probability(row[p], path));
            if (iSquared >= 0 && !missing(row[iSquared])) {
                double value = finite(row[iSquared], "i_squared", path);
                iSquaredSum += value;
                iSquaredCount++;
                if (value >= 50.0) highISquared++;
            }
        }
        int significant = (int) bh(pValues).values().stream()
            .filter(value -> value <= 0.05).count();
        return new MetaMetrics(pValues.size(), significant,
            iSquaredCount == 0 ? Double.NaN : iSquaredSum / iSquaredCount,
            iSquaredCount == 0 ? Double.NaN
                : highISquared / (double) iSquaredCount);
    }

    private static Map<String, Object> metrics(ReplicationMetrics replication,
            MetaMetrics meta) {
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("discovery_selected_pairwise_tests", replication.selected());
        result.put("replicated_pairwise_tests", replication.replicated());
        result.put("replication_fraction", numeric(replication.fraction()));
        result.put("meta_tested_features", meta.tested());
        result.put("meta_BH_FDR_0_05_features", meta.significant());
        result.put("mean_i_squared", numeric(meta.meanISquared()));
        result.put("high_i_squared_fraction", numeric(
            meta.highISquaredFraction()));
        return result;
    }

    private static double difference(double value, double reference) {
        return Double.isFinite(value) && Double.isFinite(reference)
            ? value - reference : Double.NaN;
    }

    private static Object numeric(double value) {
        return Double.isFinite(value) ? value : "NA";
    }

    private static Map<String, Object> hypothesis(Map<String, Object> study,
            String id) {
        for (Object value : list(study.get("hypotheses"), "hypotheses")) {
            Map<String, Object> hypothesis = map(value, "hypothesis");
            if (id.equals(hypothesis.get("id"))) return hypothesis;
        }
        throw new IllegalArgumentException("Unknown hypothesis: " + id);
    }

    private static List<String> role(Map<String, Object> roles, String name) {
        return roles.get(name) == null ? List.of() : strings(roles.get(name));
    }

    private static String configured(Map<String, Object> columns, String key,
            List<String> header, String... defaults) {
        if (columns.get(key) != null) return text(columns.get(key));
        for (String candidate : defaults)
            if (header.contains(candidate)) return candidate;
        throw new IllegalArgumentException("No column mapping for " + key
            + "; available columns: " + header);
    }

    private static String omicsType(String type) {
        return switch (type.toLowerCase(Locale.ROOT)) {
            case "transcriptomics", "expression", "mirna", "lincrna",
                    "lncrna", "extracellular_rna", "single_cell_rna",
                    "spatial_transcriptomics" -> "expression";
            case "epigenomics", "ewas" -> "ewas";
            case "gwas", "genomics" -> "gwas";
            case "proteomics", "single_cell_protein",
                    "spatial_proteomics" -> "proteomics";
            default -> "generic";
        };
    }

    private static double finite(String value, String name, Path path)
            throws IOException {
        try {
            double result = Double.parseDouble(value);
            if (!Double.isFinite(result)) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException invalid) {
            throw new IOException("Invalid " + name + " in " + path + ": " + value);
        }
    }

    private static double probability(String value, Path path)
            throws IOException {
        double result = finite(value, "p-value", path);
        if (result < 0.0 || result > 1.0)
            throw new IOException("P-values must be in [0,1] in " + path);
        return result;
    }

    private static boolean missing(String value) {
        return Set.of("", ".", "na", "n/a", "null", "nan")
            .contains(value.trim().toLowerCase(Locale.ROOT));
    }

    private static boolean significant(
            double adjusted, double alpha, boolean strict) {
        return Double.isFinite(adjusted)
            && (strict ? adjusted < alpha : adjusted <= alpha);
    }

    private static boolean booleanValue(Object value, String name) {
        if (value instanceof Boolean result) return result;
        if (value instanceof String text) {
            if (text.equalsIgnoreCase("true")) return true;
            if (text.equalsIgnoreCase("false")) return false;
        }
        throw new IllegalArgumentException(name + " must be true or false");
    }

    private static String value(double value) {
        return Double.isFinite(value) ? Double.toString(value) : "NA";
    }

    private static double number(Object value, String name) {
        double result = Double.parseDouble(String.valueOf(value));
        if (!Double.isFinite(result))
            throw new IllegalArgumentException(name + " must be finite");
        return result;
    }

    private static void probability(double value, String name) {
        if (value < 0.0 || value > 1.0)
            throw new IllegalArgumentException(name + " must be in [0,1]");
    }

    private static int positiveInteger(Object value, String name) {
        try {
            int result = Integer.parseInt(String.valueOf(value));
            if (result < 1) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException invalid) {
            throw new IllegalArgumentException(name + " must be a positive integer");
        }
    }

    private static Map<String, Object> map(Object value, String name) {
        if (!(value instanceof Map<?, ?> raw))
            throw new IllegalArgumentException(name + " must be a mapping");
        Map<String, Object> result = new LinkedHashMap<>();
        raw.forEach((key, item) -> result.put(String.valueOf(key), item));
        return result;
    }

    private static Map<String, Object> optionalMap(Object value, String name) {
        return value == null ? new LinkedHashMap<>() : map(value, name);
    }

    private static List<Object> list(Object value, String name) {
        if (!(value instanceof Collection<?> raw))
            throw new IllegalArgumentException(name + " must be a list");
        return new ArrayList<>(raw);
    }

    private static List<String> strings(Object value) {
        if (value == null) return new ArrayList<>();
        if (value instanceof String text) return new ArrayList<>(List.of(text));
        return list(value, "value").stream().map(String::valueOf)
            .collect(java.util.stream.Collectors.toCollection(ArrayList::new));
    }

    private static String text(Object value) {
        if (value == null || String.valueOf(value).isBlank())
            throw new IllegalArgumentException("Required text is absent");
        return String.valueOf(value);
    }

    private static <T> T required(Map<String, T> values, String key,
            String name) {
        T result = values.get(key);
        if (result == null)
            throw new IllegalArgumentException("Unknown " + name + ": " + key);
        return result;
    }

    private static String portable(Path root, Path path) {
        try {
            return root.relativize(path.toAbsolutePath().normalize())
                .toString().replace('\\', '/');
        } catch (IllegalArgumentException differentRoot) {
            return path.toAbsolutePath().normalize().toString();
        }
    }

    private record Estimate(double beta, double se, double p) { }

    private record ReplicationMetrics(int selected, int replicated,
            double fraction) { }

    private record MetaMetrics(int tested, int significant,
            double meanISquared, double highISquaredFraction) { }

    private record GeneMapping(String ensembl, String symbol) { }

    private record Candidate(String feature, String beta, double pValue,
            String iSquared) { }

    private static final class Moment {
        private long n, selected, replicated, sameSign;
        private double x, y, xx, yy, xy;

        void add(double left, double right, boolean selectedValue,
                boolean replicatedValue) {
            n++; x += left; y += right; xx += left * left;
            yy += right * right; xy += left * right;
            if (Math.signum(left) == Math.signum(right)) sameSign++;
            if (selectedValue) selected++;
            if (replicatedValue) replicated++;
        }

        double correlation() {
            double numerator = n * xy - x * y;
            double denominator = Math.sqrt((n * xx - x * x)
                * (n * yy - y * y));
            return denominator > 0 ? numerator / denominator : Double.NaN;
        }

        double replicationFraction() {
            return selected == 0 ? Double.NaN
                : replicated / (double) selected;
        }

        double slope() {
            double denominator = n * xx - x * x;
            return denominator > 0 ? (n * xy - x * y) / denominator
                : Double.NaN;
        }

        double intercept() {
            double slopeValue = slope();
            return n > 0 && Double.isFinite(slopeValue)
                ? (y - slopeValue * x) / n : Double.NaN;
        }

        double standardDeviationRatio() {
            double left = n * xx - x * x;
            double right = n * yy - y * y;
            return left > 0 && right >= 0 ? Math.sqrt(right / left)
                : Double.NaN;
        }

        double sameSignFraction() {
            return n == 0 ? Double.NaN : sameSign / (double) n;
        }
    }

    private record Standardization(int inputRows, int excludedInvalid,
            int excludedStatus, Set<String> allFeatures,
            Set<String> eligibleFeatures) { }
}
