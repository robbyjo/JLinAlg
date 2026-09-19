/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.differential.DifferentialFit;
import org.jlinalg.differential.DifferentialResult;
import org.jlinalg.differential.EmpiricalBayesDifferential;
import org.jlinalg.ewas.EwasProbe;
import org.jlinalg.ewas.EwasRegion;
import org.jlinalg.ewas.RegionLevelEwas;
import org.jlinalg.imputation.ImputationDiagnostic;
import org.jlinalg.imputation.MiceImputer;
import org.jlinalg.imputation.MiceOptions;
import org.jlinalg.imputation.MiceResult;
import org.jlinalg.imputation.RubinPooling;
import org.jlinalg.imputation.VariableType;
import org.jlinalg.multipletesting.HierarchicalTesting;
import org.jlinalg.multipletesting.IndependentHypothesisWeighting;
import org.jlinalg.multipletesting.WeightedBenjaminiHochberg;
import org.jlinalg.pipeline.DelimitedMatrixSource;
import org.jlinalg.pipeline.NumericBlock;
import org.jlinalg.pipeline.NumericBlockReader;

/** Differential, regional EWAS, multiple-testing, and imputation workflows. */
final class XwasInferenceCli {
    private XwasInferenceCli() { }

    static int run(String command, String[] arguments, PrintStream output,
            PrintStream error) {
        try {
            if (Arrays.asList(arguments).contains("--help")) {
                output.println(help(command));
                return 0;
            }
            switch (command) {
                case "differential" -> differential(arguments, output);
                case "ewas-regions" -> regions(arguments, output);
                case "multiple-test" -> multipleTesting(arguments, output);
                case "multiple-impute" -> impute(arguments, output);
                case "mi-pool" -> pool(arguments, output);
                default -> throw new IllegalArgumentException(
                    "unknown xWAS inference command: " + command);
            }
            return 0;
        } catch (IOException | RuntimeException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static void differential(String[] arguments, PrintStream output)
            throws IOException {
        Map<String, String> options = XwasFiles.options(arguments,
            "--omics", "--pheno", "--id", "--group", "--covariates",
            "--method", "--out");
        String method = XwasFiles.required(options, "--method").toLowerCase(Locale.ROOT);
        DelimitedData phenotype = DelimitedData.read(XwasFiles.path(options, "--pheno"));
        String idColumn = XwasFiles.required(options, "--id");
        String groupColumn = XwasFiles.required(options, "--group");
        OmicsInput omics = omics(XwasFiles.path(options, "--omics"), phenotype, idColumn);
        List<String> covariates = options.containsKey("--covariates")
            ? XwasFiles.names(options.get("--covariates")) : List.of();
        Design design = design(omics.phenotypeRows, phenotype, groupColumn, covariates);
        DifferentialFit fit = switch (method) {
            case "limma" -> EmpiricalBayesDifferential.fitContinuous(
                omics.values, design.matrix, design.contrast);
            case "voom" -> EmpiricalBayesDifferential.fitVoom(
                omics.values, design.matrix, design.contrast);
            case "negative-binomial", "nb" ->
                EmpiricalBayesDifferential.fitNegativeBinomial(
                    omics.values, design.matrix, design.contrast);
            default -> throw new IllegalArgumentException(
                "--method must be limma, voom, or negative-binomial");
        };
        // Retain the full tested family conservatively; unavailable tests enter
        // adjustment as one, but remain NaN in the published inference columns.
        double[] adjusted = WeightedBenjaminiHochberg.adjust(fit.results().stream()
            .mapToDouble(row -> Double.isFinite(row.pValue()) ? row.pValue() : 1.0).toArray());
        for (int i = 0; i < adjusted.length; i++)
            if (!Double.isFinite(fit.results().get(i).pValue())) adjusted[i] = Double.NaN;
        StringBuilder results = new StringBuilder(
            "feature_id\teffect\tlog2_fold_change\tse\tstatistic\tdf\tp\tfdr_bh"
            + "\tbase_mean\traw_variance\tmoderated_variance\tdispersion"
            + "\tmean_precision_weight\tconverged\titerations\n");
        for (int index = 0; index < fit.results().size(); index++) {
            DifferentialResult row = fit.results().get(index);
            results.append(omics.featureIds.get(index)).append('\t')
                .append(row.effect()).append('\t').append(row.log2FoldChange()).append('\t')
                .append(row.standardError()).append('\t').append(row.statistic()).append('\t')
                .append(row.degreesOfFreedom()).append('\t').append(row.pValue()).append('\t')
                .append(adjusted[index]).append('\t').append(row.baseMean()).append('\t')
                .append(row.rawVariance()).append('\t').append(row.moderatedVariance()).append('\t')
                .append(row.dispersion()).append('\t').append(row.meanPrecisionWeight()).append('\t')
                .append(row.converged()).append('\t').append(row.iterations()).append('\n');
        }
        Path destination = XwasFiles.path(options, "--out");
        StringBuilder factors = new StringBuilder("sample_id\tsize_factor\n");
        double[] sizeFactors = fit.sizeFactors();
        for (int sample = 0; sample < sizeFactors.length; sample++)
            factors.append(omics.sampleIds.get(sample)).append('\t')
                .append(sizeFactors[sample]).append('\n');
        String metadata = "key\tvalue\nmethod\t" + fit.method()
            + "\nreference_group\t" + design.referenceGroup
            + "\ntested_group\t" + design.testedGroup
            + "\naligned_samples\t" + omics.sampleIds.size()
            + "\nfeatures\t" + omics.featureIds.size()
            + "\nprior_df\t" + fit.priorDegreesOfFreedom()
            + "\nprior_variance_or_dispersion\t" + fit.priorVariance()
            + "\nfdr_family\tall fitted features, including nonsignificant rows\n";
        Map<Path, String> files = new LinkedHashMap<>();
        files.put(destination, results.toString());
        files.put(suffix(destination, ".size-factors.tsv"), factors.toString());
        files.put(suffix(destination, ".metadata.tsv"), metadata);
        XwasFiles.publish(files);
        output.printf(Locale.ROOT,
            "method=%s aligned_samples=%d features=%d reference=%s tested=%s%n",
            fit.method(), omics.sampleIds.size(), omics.featureIds.size(),
            design.referenceGroup, design.testedGroup);
    }

    private static void regions(String[] arguments, PrintStream output)
            throws IOException {
        Map<String, String> options = XwasFiles.options(arguments,
            "--input", "--probe", "--chromosome", "--position", "--effect",
            "--p", "--genome-build", "--max-gap", "--min-probes",
            "--correlation-length", "--out");
        DelimitedData input = DelimitedData.read(XwasFiles.path(options, "--input"));
        String probeColumn = options.getOrDefault("--probe", "probe_id");
        String chromosomeColumn = options.getOrDefault("--chromosome", "chromosome");
        String positionColumn = options.getOrDefault("--position", "position");
        String effectColumn = options.getOrDefault("--effect", "effect");
        String pColumn = options.getOrDefault("--p", "p");
        List<EwasProbe> probes = new ArrayList<>();
        for (String[] row : input.rows()) probes.add(new EwasProbe(
            XwasFiles.id(row[input.column(probeColumn)]),
            XwasFiles.id(row[input.column(chromosomeColumn)]),
            positiveLong(row[input.column(positionColumn)], positionColumn),
            XwasFiles.number(row[input.column(effectColumn)]),
            probability(row[input.column(pColumn)], pColumn)));
        RegionLevelEwas.Options regionOptions = new RegionLevelEwas.Options(
            longOption(options, "--max-gap", 500L),
            integerOption(options, "--min-probes", 3),
            doubleOption(options, "--correlation-length", 200.0));
        List<EwasRegion> regions = RegionLevelEwas.scan(probes,
            XwasFiles.required(options, "--genome-build"), regionOptions);
        StringBuilder text = new StringBuilder(
            "region_id\tgenome_build\tchromosome\tstart\tend\tprobe_count\tprobes"
            + "\tmean_effect\tspatial_z\tp\tfdr_bh\tprobes_per_kb"
            + "\tmaximum_observed_gap\tcorrelation_length\n");
        for (int index = 0; index < regions.size(); index++) {
            EwasRegion region = regions.get(index);
            text.append("region_").append(index + 1).append('\t')
                .append(region.genomeBuild()).append('\t').append(region.chromosome()).append('\t')
                .append(region.start()).append('\t').append(region.end()).append('\t')
                .append(region.probes().size()).append('\t')
                .append(String.join(";", region.probes())).append('\t')
                .append(region.meanEffect()).append('\t').append(region.statistic()).append('\t')
                .append(region.pValue()).append('\t').append(region.adjustedPValue()).append('\t')
                .append(region.probesPerKilobase()).append('\t')
                .append(region.maximumObservedGap()).append('\t')
                .append(region.correlationLength()).append('\n');
        }
        Path destination = XwasFiles.path(options, "--out");
        Map<Path, String> files = new LinkedHashMap<>();
        files.put(destination, text.toString());
        files.put(suffix(destination, ".metadata.tsv"),
            "key\tvalue\ngenome_build\t" + XwasFiles.required(options, "--genome-build")
            + "\nregion_definition\tnonoverlapping same-chromosome runs"
            + "\nmaximum_gap\t" + regionOptions.maximumGap()
            + "\nminimum_probes\t" + regionOptions.minimumProbes()
            + "\nspatial_correlation\texp(-distance/" + regionOptions.correlationLength() + ")"
            + "\nregion_fdr_family\tall retained regions\n");
        XwasFiles.publish(files);
        output.printf(Locale.ROOT, "probes=%d regions=%d genome_build=%s%n",
            probes.size(), regions.size(), XwasFiles.required(options, "--genome-build"));
    }

    private static void multipleTesting(String[] arguments, PrintStream output)
            throws IOException {
        Map<String, String> options = XwasFiles.options(arguments,
            "--input", "--method", "--id", "--p", "--status", "--covariate",
            "--fold", "--bins", "--lambda", "--maximum-weight",
            "--null-independent", "--parent", "--alpha", "--out");
        String method = XwasFiles.required(options, "--method").toLowerCase(Locale.ROOT);
        DelimitedData input = DelimitedData.read(XwasFiles.path(options, "--input"));
        String idColumn = options.getOrDefault("--id", "id");
        String pColumn = options.getOrDefault("--p", "p");
        int idIndex = input.column(idColumn);
        int pIndex = input.column(pColumn);
        int statusIndex = options.containsKey("--status")
            ? input.column(options.get("--status")) : -1;
        String[] ids = new String[input.rows().size()];
        String[] status = new String[ids.length];
        double[] pValues = new double[ids.length];
        for (int row = 0; row < ids.length; row++) {
            String[] fields = input.rows().get(row);
            ids[row] = XwasFiles.id(fields[idIndex]);
            status[row] = statusIndex < 0 ? "ok" : fields[statusIndex].trim();
            pValues[row] = status[row].equalsIgnoreCase("ok")
                ? probability(fields[pIndex], pColumn) : 1.0;
        }
        Path destination = XwasFiles.path(options, "--out");
        StringBuilder text;
        if (method.equals("ihw")) {
            if (!Boolean.parseBoolean(XwasFiles.required(options, "--null-independent")))
                throw new IllegalArgumentException(
                    "IHW requires --null-independent true after design-based verification");
            String covariateColumn = XwasFiles.required(options, "--covariate");
            String foldColumn = XwasFiles.required(options, "--fold");
            double[] covariate = new double[ids.length];
            int[] folds = new int[ids.length];
            Map<String, Integer> foldCodes = new LinkedHashMap<>();
            for (int row = 0; row < ids.length; row++) {
                String[] fields = input.rows().get(row);
                covariate[row] = XwasFiles.number(fields[input.column(covariateColumn)]);
                String fold = XwasFiles.id(fields[input.column(foldColumn)]);
                folds[row] = foldCodes.computeIfAbsent(fold, ignored -> foldCodes.size());
            }
            IndependentHypothesisWeighting.Options ihwOptions =
                new IndependentHypothesisWeighting.Options(
                    integerOption(options, "--bins", 5),
                    doubleOption(options, "--lambda", 0.5),
                    doubleOption(options, "--maximum-weight", 5.0));
            IndependentHypothesisWeighting.Result result =
                IndependentHypothesisWeighting.adjust(pValues, covariate, folds, ihwOptions);
            text = new StringBuilder(
                "id\tstatus\tp_for_family\tcovariate\tfold\tbin\tweight\tfdr_ihw\n");
            double[] weights = result.weights();
            int[] bins = result.bins();
            double[] adjusted = result.adjustedPValues();
            for (int row = 0; row < ids.length; row++)
                text.append(ids[row]).append('\t').append(status[row]).append('\t')
                    .append(pValues[row]).append('\t').append(covariate[row]).append('\t')
                    .append(folds[row]).append('\t').append(bins[row]).append('\t')
                    .append(weights[row]).append('\t').append(adjusted[row]).append('\n');
        } else if (method.equals("hierarchy")) {
            String parentColumn = XwasFiles.required(options, "--parent");
            String[] parents = new String[ids.length];
            int parentIndex = input.column(parentColumn);
            for (int row = 0; row < ids.length; row++)
                parents[row] = input.rows().get(row)[parentIndex].trim();
            double alpha = doubleOption(options, "--alpha", 0.05);
            HierarchicalTesting.Result result = HierarchicalTesting.adjust(
                ids, parents, pValues, alpha);
            double[] weights = result.weights();
            double[] adjusted = result.adjustedPValues();
            boolean[] rejected = result.rejected();
            text = new StringBuilder(
                "id\tparent\tstatus\tp_for_family\tweight\thierarchical_fwer\trejected\n");
            for (int row = 0; row < ids.length; row++)
                text.append(ids[row]).append('\t').append(parents[row]).append('\t')
                    .append(status[row]).append('\t').append(pValues[row]).append('\t')
                    .append(weights[row]).append('\t').append(adjusted[row]).append('\t')
                    .append(rejected[row]).append('\n');
        } else {
            throw new IllegalArgumentException("--method must be ihw or hierarchy");
        }
        long failed = Arrays.stream(status)
            .filter(value -> !value.equalsIgnoreCase("ok")).count();
        Map<Path, String> files = new LinkedHashMap<>();
        files.put(destination, text.toString());
        files.put(suffix(destination, ".metadata.tsv"),
            "key\tvalue\nmethod\t" + method + "\nhypotheses\t" + ids.length
            + "\nfailed_rows\t" + failed
            + "\nfailure_policy\tretained in family with p=1\n");
        XwasFiles.publish(files);
        output.printf(Locale.ROOT, "method=%s hypotheses=%d failed_rows=%d%n",
            method, ids.length, failed);
    }

    private static void impute(String[] arguments, PrintStream output)
            throws IOException {
        Map<String, String> options = XwasFiles.options(arguments,
            "--input", "--id", "--types", "--imputations", "--iterations",
            "--donors", "--seed", "--ridge", "--out");
        DelimitedData input = DelimitedData.read(XwasFiles.path(options, "--input"));
        LinkedHashMap<String, VariableType> requested = types(
            XwasFiles.required(options, "--types"));
        int rows = input.rows().size();
        double[][] data = new double[rows][requested.size()];
        VariableType[] variableTypes = requested.values().toArray(VariableType[]::new);
        List<String> columns = new ArrayList<>(requested.keySet());
        List<Map<String, Integer>> encoders = new ArrayList<>();
        for (int column = 0; column < columns.size(); column++) {
            VariableType type = variableTypes[column];
            int source = input.column(columns.get(column));
            Map<String, Integer> encoder = new LinkedHashMap<>();
            encoders.add(encoder);
            for (int row = 0; row < rows; row++) {
                String value = input.rows().get(row)[source].trim();
                if (missing(value)) data[row][column] = Double.NaN;
                else if (type == VariableType.CATEGORICAL)
                    data[row][column] = encoder.computeIfAbsent(value,
                        ignored -> encoder.size());
                else data[row][column] = XwasFiles.number(value);
            }
        }
        MiceOptions miceOptions = new MiceOptions(
            integerOption(options, "--imputations", 5),
            integerOption(options, "--iterations", 10),
            integerOption(options, "--donors", 5),
            longOption(options, "--seed", 1L),
            doubleOption(options, "--ridge", 1e-6));
        MiceResult result = MiceImputer.impute(data, variableTypes, miceOptions);
        Path prefix = Path.of(XwasFiles.required(options, "--out")).toAbsolutePath().normalize();
        int idIndex = options.containsKey("--id") ? input.column(options.get("--id")) : -1;
        Map<Path, String> files = new LinkedHashMap<>();
        List<double[][]> datasets = result.datasets();
        for (int imputation = 0; imputation < datasets.size(); imputation++) {
            StringBuilder text = new StringBuilder();
            if (idIndex >= 0) text.append(options.get("--id")).append('\t');
            text.append(String.join("\t", columns)).append('\n');
            double[][] completed = datasets.get(imputation);
            for (int row = 0; row < rows; row++) {
                if (idIndex >= 0) text.append(input.rows().get(row)[idIndex]).append('\t');
                for (int column = 0; column < columns.size(); column++) {
                    if (column > 0) text.append('\t');
                    if (variableTypes[column] == VariableType.CATEGORICAL)
                        text.append(label(encoders.get(column),
                            (int) Math.rint(completed[row][column])));
                    else text.append(completed[row][column]);
                }
                text.append('\n');
            }
            files.put(Path.of(prefix + ".imp" + (imputation + 1) + ".tsv"),
                text.toString());
        }
        StringBuilder diagnostics = new StringBuilder(
            "variable\ttype\tobserved\tmissing\tbetween_chain_variance"
            + "\taverage_within_chain_variance\tchain_means\tchain_variances\n");
        for (ImputationDiagnostic diagnostic : result.diagnostics())
            diagnostics.append(columns.get(diagnostic.variableIndex())).append('\t')
                .append(diagnostic.type().name().toLowerCase(Locale.ROOT)).append('\t')
                .append(diagnostic.observedCount()).append('\t')
                .append(diagnostic.missingCount()).append('\t')
                .append(diagnostic.betweenChainVariance()).append('\t')
                .append(diagnostic.averageWithinChainVariance()).append('\t')
                .append(join(diagnostic.chainMeans())).append('\t')
                .append(join(diagnostic.chainVariances())).append('\n');
        files.put(Path.of(prefix + ".diagnostics.tsv"), diagnostics.toString());
        files.put(Path.of(prefix + ".metadata.tsv"),
            "key\tvalue\nimputations\t" + miceOptions.imputations()
            + "\niterations\t" + miceOptions.iterations()
            + "\npredictive_mean_donors\t" + miceOptions.predictiveMeanDonors()
            + "\nseed\t" + miceOptions.seed()
            + "\nconditional_models\tbootstrap PMM and baseline-category multinomial logistic"
            + "\nridge\t" + miceOptions.ridge()
            + "\ndonor_ties\treshuffled per missing cell within bootstrap donor pool"
            + "\nstreams\tindependent deterministic split streams\n");
        XwasFiles.publish(files);
        output.printf(Locale.ROOT, "rows=%d variables=%d imputations=%d seed=%d%n",
            rows, columns.size(), miceOptions.imputations(), miceOptions.seed());
    }

    private static void pool(String[] arguments, PrintStream output)
            throws IOException {
        Map<String, String> options = XwasFiles.options(arguments,
            "--input", "--parameter", "--imputation", "--estimate",
            "--variance", "--complete-df", "--out");
        DelimitedData input = DelimitedData.read(XwasFiles.path(options, "--input"));
        String parameterColumn = options.getOrDefault("--parameter", "parameter");
        String imputationColumn = options.getOrDefault("--imputation", "imputation");
        String estimateColumn = options.getOrDefault("--estimate", "estimate");
        String varianceColumn = options.getOrDefault("--variance", "variance");
        Map<String, Map<String, double[]>> grouped = new LinkedHashMap<>();
        for (String[] row : input.rows()) {
            String parameter = XwasFiles.id(row[input.column(parameterColumn)]);
            String imputation = XwasFiles.id(row[input.column(imputationColumn)]);
            double[] value = { XwasFiles.number(row[input.column(estimateColumn)]),
                XwasFiles.number(row[input.column(varianceColumn)]) };
            if (value[1] < 0.0) throw new IllegalArgumentException(
                "within-imputation variances must be nonnegative");
            if (grouped.computeIfAbsent(parameter, ignored -> new LinkedHashMap<>())
                    .put(imputation, value) != null)
                throw new IllegalArgumentException(
                    "duplicate parameter/imputation row: " + parameter + "/" + imputation);
        }
        double completeDf = options.containsKey("--complete-df")
            ? XwasFiles.number(options.get("--complete-df")) : Double.POSITIVE_INFINITY;
        StringBuilder text = new StringBuilder(
            "parameter\timputations\testimate\twithin_variance\tbetween_variance"
            + "\ttotal_variance\tse\tdf\trelative_increase\tfmi\tstatistic\tp\n");
        for (Map.Entry<String, Map<String, double[]>> entry : grouped.entrySet()) {
            double[] estimates = entry.getValue().values().stream()
                .mapToDouble(value -> value[0]).toArray();
            double[] variances = entry.getValue().values().stream()
                .mapToDouble(value -> value[1]).toArray();
            RubinPooling.Estimate pooled = RubinPooling.pool(estimates, variances, completeDf);
            text.append(entry.getKey()).append('\t').append(estimates.length).append('\t')
                .append(pooled.estimate()).append('\t').append(pooled.withinVariance()).append('\t')
                .append(pooled.betweenVariance()).append('\t').append(pooled.totalVariance()).append('\t')
                .append(pooled.standardError()).append('\t').append(pooled.degreesOfFreedom()).append('\t')
                .append(pooled.relativeIncreaseInVariance()).append('\t')
                .append(pooled.fractionMissingInformation()).append('\t')
                .append(pooled.statistic()).append('\t').append(pooled.pValue()).append('\n');
        }
        Path destination = XwasFiles.path(options, "--out");
        Map<Path, String> files = new LinkedHashMap<>();
        files.put(destination, text.toString());
        files.put(suffix(destination, ".metadata.tsv"),
            "key\tvalue\nmethod\tRubin rules with Barnard-Rubin finite-df correction"
            + "\ncomplete_df\t" + completeDf + "\nparameters\t" + grouped.size() + "\n");
        XwasFiles.publish(files);
        output.printf(Locale.ROOT, "parameters=%d complete_df=%s%n",
            grouped.size(), Double.toString(completeDf));
    }

    private static OmicsInput omics(Path path, DelimitedData phenotype,
            String idColumn) throws IOException {
        DelimitedMatrixSource source = DelimitedMatrixSource.open(path);
        Map<String, String[]> byId = new LinkedHashMap<>();
        int id = phenotype.column(idColumn);
        for (String[] row : phenotype.rows()) {
            String key = XwasFiles.id(row[id]);
            if (byId.put(key, row) != null)
                throw new IllegalArgumentException("duplicate phenotype ID: " + key);
        }
        List<Integer> order = new ArrayList<>();
        List<String> samples = new ArrayList<>();
        List<String[]> rows = new ArrayList<>();
        List<String> sourceSamples = source.metadata().sampleIds();
        for (int index = 0; index < sourceSamples.size(); index++) {
            String[] row = byId.get(sourceSamples.get(index));
            if (row != null) {
                order.add(index); samples.add(sourceSamples.get(index)); rows.add(row);
            }
        }
        if (samples.isEmpty()) throw new IllegalArgumentException(
            "omics and phenotype have no aligned sample IDs");
        List<String> featureIds = new ArrayList<>();
        List<double[]> values = new ArrayList<>();
        int[] selected = order.stream().mapToInt(Integer::intValue).toArray();
        try (NumericBlockReader reader = source.open(selected)) {
            for (NumericBlock block; (block = reader.read(1024)) != null;)
                for (var row : block.rows()) {
                    featureIds.add(row.id()); values.add(row.values());
                }
        }
        return new OmicsInput(featureIds, samples, values.toArray(double[][]::new), rows);
    }

    private static Design design(List<String[]> rows, DelimitedData phenotype,
            String groupColumn, List<String> covariates) {
        int group = phenotype.column(groupColumn);
        LinkedHashMap<String, Integer> levels = new LinkedHashMap<>();
        for (String[] row : rows) levels.computeIfAbsent(
            XwasFiles.id(row[group]), ignored -> levels.size());
        if (levels.size() != 2) throw new IllegalArgumentException(
            "--group must have exactly two aligned levels");
        String reference = levels.keySet().iterator().next();
        String tested = levels.keySet().stream().skip(1).findFirst().orElseThrow();
        int[] covariateIndices = covariates.stream().mapToInt(phenotype::column).toArray();
        double[][] matrix = new double[rows.size()][covariates.size() + 2];
        for (int row = 0; row < rows.size(); row++) {
            matrix[row][0] = 1.0;
            for (int column = 0; column < covariateIndices.length; column++)
                matrix[row][column + 1] = XwasFiles.number(
                    rows.get(row)[covariateIndices[column]]);
            matrix[row][matrix[row].length - 1] =
                rows.get(row)[group].equals(tested) ? 1.0 : 0.0;
        }
        double[] contrast = new double[matrix[0].length];
        contrast[contrast.length - 1] = 1.0;
        return new Design(matrix, contrast, reference, tested);
    }

    private static LinkedHashMap<String, VariableType> types(String specification) {
        LinkedHashMap<String, VariableType> result = new LinkedHashMap<>();
        for (String item : specification.split(",", -1)) {
            String[] fields = item.trim().split(":", -1);
            if (fields.length != 2 || fields[0].isBlank() || fields[1].isBlank())
                throw new IllegalArgumentException(
                    "--types uses column:continuous|binary|categorical entries");
            VariableType type;
            try { type = VariableType.valueOf(fields[1].trim().toUpperCase(Locale.ROOT)); }
            catch (IllegalArgumentException failure) {
                throw new IllegalArgumentException("unknown imputation type: " + fields[1], failure);
            }
            if (result.put(fields[0].trim(), type) != null)
                throw new IllegalArgumentException("duplicate imputation column: " + fields[0]);
        }
        if (result.size() < 2) throw new IllegalArgumentException(
            "MICE requires at least two columns in --types");
        return result;
    }

    private static boolean missing(String value) {
        return value.isEmpty() || Set.of("na", "nan", ".")
            .contains(value.toLowerCase(Locale.ROOT));
    }

    private static String label(Map<String, Integer> encoder, int code) {
        return encoder.entrySet().stream().filter(entry -> entry.getValue() == code)
            .map(Map.Entry::getKey).findFirst().orElseThrow(() ->
                new IllegalArgumentException("imputed categorical code is unknown: " + code));
    }

    private static double probability(String value, String name) {
        double result = XwasFiles.number(value);
        if (result < 0.0 || result > 1.0)
            throw new IllegalArgumentException(name + " must be in [0,1]");
        return result;
    }

    private static int integerOption(Map<String, String> options,
            String name, int fallback) {
        if (!options.containsKey(name)) return fallback;
        int value = Integer.parseInt(options.get(name));
        if (value < 1) throw new IllegalArgumentException(name + " must be positive");
        return value;
    }

    private static long longOption(Map<String, String> options,
            String name, long fallback) {
        if (!options.containsKey(name)) return fallback;
        long value = Long.parseLong(options.get(name));
        if (value < 1) throw new IllegalArgumentException(name + " must be positive");
        return value;
    }

    private static double doubleOption(Map<String, String> options,
            String name, double fallback) {
        return options.containsKey(name) ? XwasFiles.number(options.get(name)) : fallback;
    }

    private static long positiveLong(String value, String name) {
        long result = Long.parseLong(value);
        if (result < 1) throw new IllegalArgumentException(name + " must be positive");
        return result;
    }

    private static String join(double[] values) {
        StringBuilder result = new StringBuilder();
        for (int index = 0; index < values.length; index++) {
            if (index > 0) result.append(';');
            result.append(values[index]);
        }
        return result.toString();
    }

    private static Path suffix(Path path, String suffix) {
        String value = path.toString();
        if (value.endsWith(".tsv")) value = value.substring(0, value.length() - 4);
        return Path.of(value + suffix);
    }

    static String help(String command) {
        return switch (command) {
            case "differential" -> """
                Usage: jlinalg differential --method limma|voom|negative-binomial
                  --omics FEATURES_BY_SAMPLE.tsv --pheno SAMPLE_TABLE.tsv --id sample_id
                  --group condition [--covariates age,sex] --out results.tsv
                The first observed group is the reference and the second is tested. limma uses
                moderated Gaussian models; voom adds count precision weights; negative-binomial
                uses median-ratio size factors, dispersion shrinkage, and NB log-link IRLS.
                """;
            case "ewas-regions" -> """
                Usage: jlinalg ewas-regions --input probes.tsv --genome-build GRCh38 --out regions.tsv
                  [--probe probe_id --chromosome chromosome --position position --effect effect --p p]
                  [--max-gap 500 --min-probes 3 --correlation-length 200]
                Forms nonoverlapping coordinate runs and uses signed Stouffer aggregation with
                exp(-distance/correlation-length) covariance. Output includes probe coverage and BH.
                """;
            case "multiple-test" -> """
                Usage: jlinalg multiple-test --method ihw --input tests.tsv --id id --p p
                  --covariate mean --fold fold --null-independent true --out adjusted.tsv
                Usage: jlinalg multiple-test --method hierarchy --input tests.tsv --id id --parent parent
                  --p p [--status status --alpha 0.05] --out adjusted.tsv
                Failed status rows remain in the prespecified family with p=1. IHW weights are learned
                outside each fold; hierarchy uses gated leaf-weighted Bonferroni FWER control.
                """;
            case "multiple-impute" -> """
                Usage: jlinalg multiple-impute --input data.tsv [--id sample_id]
                  --types age:continuous,case:binary,site:categorical
                  --out completed [--imputations 5 --iterations 10 --donors 5 --seed 1]
                Writes OUT.impN.tsv, OUT.diagnostics.tsv, and OUT.metadata.tsv. Missing tokens are
                blank, NA, NaN, or '.'. Each chain uses an independent deterministic random stream.
                """;
            default -> """
                Usage: jlinalg mi-pool --input estimates.tsv --out pooled.tsv [--complete-df N]
                  [--parameter parameter --imputation imputation --estimate estimate --variance variance]
                Pools model estimates and within-imputation variances using Rubin rules, with the
                Barnard-Rubin finite-sample correction when --complete-df is supplied.
                """;
        };
    }

    private record OmicsInput(List<String> featureIds, List<String> sampleIds,
        double[][] values, List<String[]> phenotypeRows) { }
    private record Design(double[][] matrix, double[] contrast,
        String referenceGroup, String testedGroup) { }
}
