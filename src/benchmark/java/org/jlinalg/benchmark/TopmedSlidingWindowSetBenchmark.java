/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmLaplaceOptions;
import org.jlinalg.glmm.GlmmLaplaceResult;
import org.jlinalg.glmm.SparseGlmmLaplace;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.pedigree.SparsePedigreeMixedModel;
import org.jlinalg.pipeline.VariantBlock;
import org.jlinalg.pipeline.VariantBlockReader;
import org.jlinalg.pipeline.VariantFilterOptions;
import org.jlinalg.pipeline.VariantFilterResult;
import org.jlinalg.pipeline.VariantFilters;
import org.jlinalg.pipeline.VariantRecord;
import org.jlinalg.pipeline.VariantSource;
import org.jlinalg.pipeline.VariantSources;
import org.jlinalg.pipeline.VariantStatistics;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.settest.EffectAllele;
import org.jlinalg.settest.KnownCovarianceSetTestNullModel;
import org.jlinalg.settest.PreparedVariantSet;
import org.jlinalg.settest.SetTestMissingPolicy;
import org.jlinalg.settest.SetTestOptions;
import org.jlinalg.settest.SetTestResult;
import org.jlinalg.settest.SetTestScoreNullModel;
import org.jlinalg.settest.SetTestSuiteResult;
import org.jlinalg.settest.SetTests;
import org.jlinalg.settest.SkatOResult;
import org.jlinalg.settest.VariantSet;
import org.jlinalg.settest.VariantWeights;
import org.jlinalg.settest.WeightedVariant;

/** Real TOPMed pedigree-aware 2 kb / 500 bp sliding-window set benchmark. */
public final class TopmedSlidingWindowSetBenchmark {
    private static volatile double checksum;

    private TopmedSlidingWindowSetBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        quietLogging();
        Options options = Options.parse(arguments);
        Files.createDirectories(options.outputDirectory);
        PreparedData data = prepare(options);
        exportPrepared(data, options);
        if (options.prepareOnly) return;

        try (BackendContext context = BackendContext.select(options.backend)) {
            System.out.printf("backend requested=%s selected=%s device=%s "
                    + "accelerated=%s automatic=%s%n",
                context.provenance().requested(),
                context.provenance().selectedBackend(),
                context.provenance().deviceDescription(),
                context.provenance().accelerated(),
                context.provenance().automaticRouting());
        }

        SetTestOptions setOptions = new SetTestOptions(
            rareFilter(options.maximumMaf), SetTestMissingPolicy.MEAN_IMPUTE,
            new double[] {0, 0.25, 0.5, 0.75, 1},
            options.skatOSimulations, 20260903L);
        List<PreparedVariantSet> preparedWindows = data.windows.stream()
            .map(window -> SetTests.prepare(window.variantSet,
                data.samples.size(), setOptions)).toList();
        double[][] design = design(data.samples);
        double[] relationship = relationship(data.pedigree, data.samples);
        List<String> observationIds = data.samples.stream()
            .map(Sample::sabreid).toList();
        PedigreeRandomEffectTerm animal = PedigreeRandomEffectTerm.of(
            "pedigree", observationIds, data.pedigree);
        List<Timing> timings = new ArrayList<>();
        List<Result> results = new ArrayList<>();

        double[] bmi = data.samples.stream().mapToDouble(Sample::bmi).toArray();
        SparseLinearMixedModelResult bmiFit = SparsePedigreeMixedModel.fit(
            bmi, design, List.of(animal), List.of(), RemlOptions.defaults(),
            BackendPolicy.PREFERRED);
        double[] bmiVariances = bmiFit.varianceComponents();
        System.out.printf(Locale.ROOT,
            "BMI sparse null genetic_variance=%.9g residual_variance=%.9g "
                + "converged=%s%n",
            bmiVariances[0], bmiVariances[1], bmiFit.converged());
        try (KnownCovarianceSetTestNullModel bmiNull =
                KnownCovarianceSetTestNullModel.prepare(bmi, design,
                marginalCovariance(relationship, bmiVariances[0],
                    bmiVariances[1]), options.backend)) {
            runModel("BMI", bmiNull,
                preparedWindows, data.windows, setOptions, options,
                timings, results);
        }

        double[] obesity = data.samples.stream()
            .mapToDouble(Sample::obesity).toArray();
        GlmmLaplaceResult obesityFit = SparseGlmmLaplace.fitWithPrecision(
            obesity, MatrixOps.rowMajor(design, obesity.length),
            obesity.length, design[0].length, GlmFamilies.binomial(),
            List.of(animal.randomEffect()), List.of(animal.precision()),
            null, null, GlmmLaplaceOptions.defaults(),
            BackendPolicy.PREFERRED);
        double obesityVariance = obesityFit.varianceComponents()[0];
        System.out.printf(Locale.ROOT,
            "Obesity sparse null genetic_variance=%.9g converged=%s%n",
            obesityVariance, obesityFit.converged());
        WorkingBinary working = workingBinary(obesity,
            obesityFit.linearPredictor(), obesityFit.fittedMeans(),
            GlmFamilies.binomial());
        try (KnownCovarianceSetTestNullModel obesityNull =
                KnownCovarianceSetTestNullModel.prepare(
                working.response, design,
                marginalCovariance(relationship, obesityVariance,
                    working.residualDiagonal), options.backend)) {
            runModel("Obesity", obesityNull, preparedWindows, data.windows,
                setOptions, options, timings, results);
        }

        writeTimings(options.outputDirectory.resolve("jlinalg_timings.csv"),
            timings);
        writeResults(options.outputDirectory.resolve("jlinalg_results.csv"),
            results);
        if (!Double.isFinite(checksum))
            throw new IllegalStateException("non-finite benchmark checksum");
    }

    private static void runModel(
            String model, SetTestScoreNullModel nullModel,
            List<PreparedVariantSet> sets, List<Window> windows,
            SetTestOptions setOptions, Options options,
            List<Timing> timings, List<Result> results) {
        for (String method : List.of("burden", "skat", "skat-o", "suite")) {
            run(method, sets.get(0), nullModel, setOptions);
            for (int measurement = 1;
                    measurement <= options.measurements; measurement++) {
                System.gc();
                long started = System.nanoTime();
                Object last = null;
                for (PreparedVariantSet set : sets)
                    last = run(method, set, nullModel, setOptions);
                double seconds = (System.nanoTime() - started) / 1e9;
                consume(last);
                timings.add(new Timing("JLinAlg", model, method,
                    options.backend.name(), measurement, sets.size(), seconds,
                    sets.size() / seconds));
                System.out.printf(Locale.ROOT,
                    "JLinAlg model=%s method=%s measurement=%d windows=%d "
                        + "seconds=%.6f windows_per_second=%.3f%n",
                    model, method, measurement, sets.size(), seconds,
                    sets.size() / seconds);
            }
        }
        for (int index = 0; index < sets.size(); index++) {
            SetTestSuiteResult suite = SetTests.scoreSuite(
                sets.get(index), nullModel, setOptions);
            Window window = windows.get(index);
            results.add(result(model, window, "burden",
                suite.burden().pValue(), suite.burden().statistic()));
            results.add(result(model, window, "skat",
                suite.skat().pValue(), suite.skat().statistic()));
            results.add(result(model, window, "skat-o",
                suite.skatO().adjustedPValue(),
                suite.skatO().minimumComponentPValue()));
        }
    }

    private static Object run(
            String method, PreparedVariantSet set,
            SetTestScoreNullModel nullModel, SetTestOptions options) {
        return switch (method) {
            case "burden" -> SetTests.burden(set, nullModel);
            case "skat" -> SetTests.skat(set, nullModel);
            case "skat-o" -> SetTests.skatO(set, nullModel, options);
            case "suite" -> SetTests.scoreSuite(set, nullModel, options);
            default -> throw new IllegalArgumentException("unknown method: " + method);
        };
    }

    private static void consume(Object result) {
        if (result instanceof SetTestResult value) checksum += value.pValue();
        else if (result instanceof SkatOResult value)
            checksum += value.adjustedPValue();
        else if (result instanceof SetTestSuiteResult value)
            checksum += value.skatO().adjustedPValue();
    }

    private static PreparedData prepare(Options options) throws IOException {
        Map<String, PhenotypeRow> phenotype = readPhenotype(options.phenotype);
        VariantSource source = VariantSources.open(options.genotype);
        List<String> vcfIds = source.metadata().sampleIds();
        List<Sample> samples = new ArrayList<>();
        List<Integer> order = new ArrayList<>();
        for (int index = 0; index < vcfIds.size(); index++) {
            PhenotypeRow row = phenotype.get(vcfIds.get(index));
            if (row == null || !row.complete()) continue;
            order.add(index);
            samples.add(row.sample());
        }
        if (samples.isEmpty())
            throw new IllegalArgumentException("no complete matched framid samples");
        int[] requestedOrder = order.stream().mapToInt(Integer::intValue).toArray();
        List<Gene> genes = readGenes(options.annotation);
        List<WindowSeed> seeds = scanWindows(
            source, requestedOrder, genes, options);
        Map<String, VariantDatum> variants = new LinkedHashMap<>();
        List<Window> windows = new ArrayList<>();
        int variantNumber = 0;
        for (int index = 0; index < seeds.size(); index++) {
            WindowSeed seed = seeds.get(index);
            List<WeightedVariant> members = new ArrayList<>();
            List<String> keys = new ArrayList<>();
            for (VariantRecord record : seed.variants) {
                String canonical = canonical(record);
                VariantDatum datum = variants.get(canonical);
                if (datum == null) {
                    VariantStatistics statistics = VariantStatistics.of(record);
                    EffectAllele effect = statistics.alternateAlleleFrequency() <= 0.5
                        ? EffectAllele.ALTERNATE : EffectAllele.REFERENCE;
                    String key = "v" + (++variantNumber);
                    datum = new VariantDatum(key, record, statistics, effect,
                        VariantWeights.betaBurden(
                            statistics.minorAlleleFrequency(), 1, 25));
                    variants.put(canonical, datum);
                }
                keys.add(datum.key);
                members.add(new WeightedVariant(
                    datum.record, datum.effectAllele, datum.weight));
            }
            String id = String.format(Locale.ROOT, "window%02d", index + 1);
            String labels = genes.stream()
                .filter(gene -> sameChromosome(gene.chromosome, seed.chromosome)
                    && gene.start <= seed.end && gene.end >= seed.start)
                .map(gene -> gene.name + "|" + gene.id)
                .distinct().reduce((left, right) -> left + ";" + right)
                .orElse("");
            windows.add(new Window(id, seed.chromosome, seed.start, seed.end,
                labels, keys, new VariantSet(id, members)));
        }
        Pedigree pedigree = readPedigree(options.pedigree, samples);
        System.out.printf(Locale.ROOT,
            "prepared samples=%d windows=%d unique_rare_variants=%d "
                + "maf_max=%.4f pedigree_individuals=%d%n",
            samples.size(), windows.size(), variants.size(),
            options.maximumMaf, pedigree.size());
        return new PreparedData(List.copyOf(samples), List.copyOf(windows),
            List.copyOf(variants.values()), pedigree);
    }

    private static Map<String, PhenotypeRow> readPhenotype(Path path)
            throws IOException {
        Map<String, PhenotypeRow> result = new LinkedHashMap<>();
        try (BufferedReader reader = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int framid = required(header, "framid");
            int sabreid = required(header, "sabreid");
            int bmi = required(header, "BMI");
            int obesity = required(header, "Obesity");
            int sex = required(header, "Sex");
            int age = required(header, "Age");
            int wbc = required(header, "WBC_Pred");
            int ly = required(header, "LY_PER_Pred");
            int mo = required(header, "MO_PER_Pred");
            int eo = required(header, "EO_PER_Pred");
            int ba = required(header, "BA_PER_Pred");
            for (String line; (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                String id = row.get(framid).trim();
                if (id.isEmpty()) continue;
                Sample sample = new Sample(id, row.get(sabreid).trim(),
                    number(row.get(bmi)), number(row.get(obesity)),
                    number(row.get(sex)), number(row.get(age)),
                    number(row.get(wbc)), number(row.get(ly)),
                    number(row.get(mo)), number(row.get(eo)),
                    number(row.get(ba)));
                if (result.putIfAbsent(id, new PhenotypeRow(sample)) != null)
                    throw new IllegalArgumentException(
                        "duplicate phenotype framid: " + id);
            }
        }
        return result;
    }

    private static List<WindowSeed> scanWindows(
            VariantSource source, int[] order, List<Gene> genes,
            Options options)
            throws IOException {
        TreeMap<Long, List<VariantRecord>> active = new TreeMap<>();
        List<WindowSeed> result = new ArrayList<>();
        String chromosome = null;
        VariantFilterOptions filter = rareFilter(options.maximumMaf);
        boolean done = false;
        try (VariantBlockReader reader = source.open(order)) {
            for (VariantBlock block; !done && (block = reader.read(256)) != null;) {
                for (VariantRecord original : block.variants()) {
                    if (chromosome == null) chromosome = original.chromosome();
                    if (!sameChromosome(chromosome, original.chromosome())) {
                        finalizeWindows(active, chromosome, Long.MAX_VALUE,
                            genes, options, result);
                        done = true;
                        break;
                    }
                    finalizeWindows(active, chromosome, original.position(),
                        genes, options, result);
                    if (result.size() >= options.windowCount) {
                        done = true;
                        break;
                    }
                    VariantRecord record = canonicalRecord(original);
                    VariantFilterResult evaluated = VariantFilters.evaluate(
                        record, filter);
                    if (!evaluated.included()) continue;
                    long latest = ((record.position() - 1) / options.stride)
                        * options.stride + 1;
                    long earliest = Math.max(1,
                        record.position() - options.windowSize + 1);
                    for (long start = latest; start >= earliest;
                            start -= options.stride)
                        active.computeIfAbsent(start,
                            ignored -> new ArrayList<>()).add(record);
                }
            }
        }
        if (result.size() < options.windowCount)
            finalizeWindows(active, chromosome, Long.MAX_VALUE,
                genes, options, result);
        if (result.size() < options.windowCount)
            throw new IllegalArgumentException("VCF contains only "
                + result.size() + " testable rare-variant windows");
        return List.copyOf(result.subList(0, options.windowCount));
    }

    private static void finalizeWindows(
            TreeMap<Long, List<VariantRecord>> active, String chromosome,
            long nextPosition, List<Gene> genes, Options options,
            List<WindowSeed> output) {
        while (!active.isEmpty()
                && active.firstKey() + options.windowSize - 1 < nextPosition) {
            Map.Entry<Long, List<VariantRecord>> entry = active.pollFirstEntry();
            long end = entry.getKey() + options.windowSize - 1;
            boolean annotated = genes.stream().anyMatch(gene ->
                sameChromosome(gene.chromosome, chromosome)
                    && gene.start <= end && gene.end >= entry.getKey());
            if (entry.getValue().size() >= options.minimumVariants && annotated)
                output.add(new WindowSeed(chromosome, entry.getKey(),
                    end, List.copyOf(entry.getValue())));
            if (output.size() >= options.windowCount) return;
        }
    }

    private static Pedigree readPedigree(Path path, List<Sample> samples)
            throws IOException {
        Map<String, PedigreeIndividual> all = new LinkedHashMap<>();
        try (BufferedReader reader = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int id = required(header, "sabreid");
            int sire = required(header, "fid");
            int dam = required(header, "mid");
            for (String line; (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                String value = nullable(row.get(id));
                if (value == null) continue;
                PedigreeIndividual individual = new PedigreeIndividual(
                    value, nullable(row.get(sire)), nullable(row.get(dam)));
                if (all.putIfAbsent(value, individual) != null)
                    throw new IllegalArgumentException(
                        "duplicate pedigree sabreid: " + value);
            }
        }
        Map<String, PedigreeIndividual> selected = new LinkedHashMap<>();
        ArrayDeque<String> queue = new ArrayDeque<>();
        for (Sample sample : samples) queue.add(sample.sabreid);
        while (!queue.isEmpty()) {
            String id = queue.removeFirst();
            if (selected.containsKey(id)) continue;
            PedigreeIndividual individual = all.get(id);
            if (individual == null)
                individual = new PedigreeIndividual(id, null, null);
            selected.put(id, individual);
            String[] parents = {individual.sireId(), individual.damId()};
            for (String parent : parents) {
                if (parent == null) continue;
                if (!all.containsKey(parent))
                    all.put(parent, new PedigreeIndividual(parent, null, null));
                queue.add(parent);
            }
        }
        return Pedigree.of(new ArrayList<>(selected.values()));
    }

    private static List<Gene> readGenes(Path path) throws IOException {
        List<Gene> result = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int chromosome = required(header, "seqname");
            int start = required(header, "start");
            int end = required(header, "end");
            int id = required(header, "gene_id");
            int name = required(header, "gene_name");
            for (String line; (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                if (!sameChromosome("22", row.get(chromosome))) continue;
                result.add(new Gene(row.get(chromosome),
                    Long.parseLong(row.get(start)), Long.parseLong(row.get(end)),
                    row.get(id), row.get(name)));
            }
        }
        return result;
    }

    private static double[][] design(List<Sample> samples) {
        double[][] result = new double[samples.size()][8];
        for (int row = 0; row < samples.size(); row++) {
            Sample value = samples.get(row);
            result[row] = new double[] {1, value.sex, value.age, value.wbc,
                value.ly, value.mo, value.eo, value.ba};
        }
        return result;
    }

    private static double[] relationship(
            Pedigree pedigree, List<Sample> samples) {
        int rows = samples.size();
        double[] result = new double[rows * rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column <= row; column++) {
                double value = pedigree.relationship(
                    samples.get(row).sabreid, samples.get(column).sabreid);
                result[row * rows + column] = value;
                result[column * rows + row] = value;
            }
        }
        return result;
    }

    private static double[] marginalCovariance(
            double[] relationship, double geneticVariance,
            double residualVariance) {
        double[] residual = new double[(int) Math.sqrt(relationship.length)];
        Arrays.fill(residual, residualVariance);
        return marginalCovariance(relationship, geneticVariance, residual);
    }

    private static double[] marginalCovariance(
            double[] relationship, double geneticVariance,
            double[] residualDiagonal) {
        int rows = residualDiagonal.length;
        if (relationship.length != rows * rows)
            throw new IllegalArgumentException("covariance dimensions differ");
        double[] result = relationship.clone();
        for (int index = 0; index < result.length; index++)
            result[index] *= geneticVariance;
        for (int row = 0; row < rows; row++)
            result[row * rows + row] += residualDiagonal[row];
        return result;
    }

    private static WorkingBinary workingBinary(
            double[] response, double[] predictor, double[] means,
            GlmFamily family) {
        double[] working = new double[response.length];
        double[] residual = new double[response.length];
        for (int row = 0; row < response.length; row++) {
            double derivative = family.meanDerivative(predictor[row]);
            double variance = family.variance(means[row]);
            if (!Double.isFinite(derivative) || derivative == 0
                    || !(variance > 0) || !Double.isFinite(variance))
                throw new IllegalArgumentException(
                    "invalid binary working-model weights");
            working[row] = predictor[row]
                + (response[row] - means[row]) / derivative;
            residual[row] = variance / (derivative * derivative);
        }
        return new WorkingBinary(working, residual);
    }

    private static void exportPrepared(PreparedData data, Options options)
            throws IOException {
        Path directory = options.outputDirectory;
        try (BufferedWriter writer = Files.newBufferedWriter(
                directory.resolve("analysis.csv"), StandardCharsets.UTF_8)) {
            writer.write("framid,sabreid,BMI,Obesity,Sex,Age,WBC_Pred,"
                + "LY_PER_Pred,MO_PER_Pred,EO_PER_Pred,BA_PER_Pred\n");
            for (Sample value : data.samples)
                writer.write(String.format(Locale.ROOT,
                    "%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g%n",
                    csvField(value.framid), csvField(value.sabreid), value.bmi,
                    value.obesity, value.sex, value.age, value.wbc, value.ly,
                    value.mo, value.eo, value.ba));
        }
        try (BufferedWriter writer = Files.newBufferedWriter(
                directory.resolve("variants.csv"), StandardCharsets.UTF_8)) {
            writer.write("key,variant_id,chromosome,position,reference,alternate,"
                + "maf,effect_allele,weight\n");
            for (VariantDatum value : data.variants)
                writer.write(String.format(Locale.ROOT,
                    "%s,%s,%s,%d,%s,%s,%.17g,%s,%.17g%n", value.key,
                    csvField(value.record.id()), csvField(value.record.chromosome()),
                    value.record.position(), csvField(value.record.referenceAllele()),
                    csvField(value.record.alternateAllele()),
                    value.statistics.minorAlleleFrequency(), value.effectAllele,
                    value.weight));
        }
        try (BufferedWriter writer = Files.newBufferedWriter(
                directory.resolve("genotypes.csv"), StandardCharsets.UTF_8)) {
            writer.write("framid");
            for (VariantDatum value : data.variants)
                writer.write("," + value.key);
            writer.write("\n");
            List<double[]> dosage = data.variants.stream()
                .map(value -> value.record.dosages()).toList();
            for (int sample = 0; sample < data.samples.size(); sample++) {
                writer.write(csvField(data.samples.get(sample).framid));
                for (double[] values : dosage) {
                    writer.write(",");
                    if (Double.isFinite(values[sample]))
                        writer.write(String.format(Locale.ROOT,
                            "%.17g", values[sample]));
                }
                writer.write("\n");
            }
        }
        try (BufferedWriter writer = Files.newBufferedWriter(
                directory.resolve("windows.csv"), StandardCharsets.UTF_8)) {
            writer.write("window_id,chromosome,start,end,genes,variant_keys\n");
            for (Window value : data.windows)
                writer.write(String.format(Locale.ROOT, "%s,%s,%d,%d,%s,%s%n",
                    value.id, csvField(value.chromosome), value.start, value.end,
                    csvField(value.genes), csvField(String.join(";", value.keys))));
        }
        double[] matrix = relationship(data.pedigree, data.samples);
        try (DataOutputStream output = new DataOutputStream(
                Files.newOutputStream(directory.resolve("relationship.bin")))) {
            for (double value : matrix) output.writeDouble(value);
        }
        try (BufferedWriter writer = Files.newBufferedWriter(
                directory.resolve("manifest.csv"), StandardCharsets.UTF_8)) {
            writer.write("metric,value\n");
            writer.write("samples," + data.samples.size() + "\n");
            writer.write("windows," + data.windows.size() + "\n");
            writer.write("unique_rare_variants," + data.variants.size() + "\n");
            writer.write("window_size," + options.windowSize + "\n");
            writer.write("stride," + options.stride + "\n");
            writer.write("maximum_maf," + options.maximumMaf + "\n");
            writer.write("matching_id,framid\n");
            writer.write("pedigree_id,sabreid\n");
            writer.write("unlisted_pedigree_ids,singleton\n");
            writer.write("filter_and_matching_timed,false\n");
        }
    }

    private static void writeTimings(Path path, List<Timing> values)
            throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,method,backend,measurement,windows,"
                + "seconds,windows_per_second\n");
            for (Timing value : values)
                writer.write(String.format(Locale.ROOT,
                    "%s,%s,%s,%s,%d,%d,%.9f,%.9f%n", value.runtime,
                    value.model, value.method, value.backend, value.measurement,
                    value.windows, value.seconds, value.windowsPerSecond));
        }
    }

    private static void writeResults(Path path, List<Result> values)
            throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,window_id,chromosome,start,end,genes,"
                + "variants,method,p_value,statistic\n");
            for (Result value : values)
                writer.write(String.format(Locale.ROOT,
                    "%s,%s,%s,%s,%d,%d,%s,%d,%s,%.17g,%.17g%n",
                    value.runtime, value.model, value.windowId,
                    csvField(value.chromosome), value.start, value.end,
                    csvField(value.genes), value.variants, value.method,
                    value.pValue, value.statistic));
        }
    }

    private static Result result(
            String model, Window window, String method,
            double pValue, double statistic) {
        return new Result("JLinAlg", model, window.id, window.chromosome,
            window.start, window.end, window.genes, window.keys.size(), method,
            pValue, statistic);
    }

    private static VariantFilterOptions rareFilter(double maximumMaf) {
        return VariantFilterOptions.builder().maximumMaf(maximumMaf).build();
    }

    private static VariantRecord canonicalRecord(VariantRecord value) {
        return new VariantRecord(canonical(value), value.chromosome(),
            value.position(), value.referenceAllele(), value.alternateAllele(),
            value.dosages(), value.imputationQuality());
    }

    private static String canonical(VariantRecord value) {
        return value.chromosome() + ":" + value.position() + ":"
            + value.referenceAllele() + ":" + value.alternateAllele();
    }

    private static boolean sameChromosome(String left, String right) {
        return left.replaceFirst("(?i)^chr", "").equalsIgnoreCase(
            right.replaceFirst("(?i)^chr", ""));
    }

    private static Double number(String value) {
        if (value == null || value.isBlank() || value.equalsIgnoreCase("NA")
                || value.equals(".")) return null;
        try {
            double result = Double.parseDouble(value);
            return Double.isFinite(result) ? result : null;
        } catch (NumberFormatException ignored) {
            return null;
        }
    }

    private static String nullable(String value) {
        if (value == null) return null;
        String result = value.trim();
        return result.isEmpty() || result.equals("0")
            || result.equals(".") || result.equalsIgnoreCase("NA")
                ? null : result;
    }

    private static int required(List<String> header, String name) {
        int result = header.indexOf(name);
        if (result < 0)
            throw new IllegalArgumentException("missing column: " + name);
        return result;
    }

    private static List<String> csv(String line) {
        if (line == null) throw new IllegalArgumentException("empty CSV file");
        List<String> result = new ArrayList<>();
        StringBuilder field = new StringBuilder();
        boolean quoted = false;
        for (int index = 0; index < line.length(); index++) {
            char value = line.charAt(index);
            if (value == '"') {
                if (quoted && index + 1 < line.length()
                        && line.charAt(index + 1) == '"') {
                    field.append('"');
                    index++;
                } else quoted = !quoted;
            } else if (value == ',' && !quoted) {
                result.add(field.toString());
                field.setLength(0);
            } else field.append(value);
        }
        result.add(field.toString());
        return result;
    }

    private static String csvField(String value) {
        if (value.indexOf(',') < 0 && value.indexOf('"') < 0
                && value.indexOf('\n') < 0) return value;
        return '"' + value.replace("\"", "\"\"") + '"';
    }

    private static void quietLogging() {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        Logger.getLogger("").setLevel(Level.WARNING);
        Arrays.stream(Logger.getLogger("").getHandlers())
            .forEach(handler -> handler.setLevel(Level.WARNING));
    }

    private record Sample(
        String framid, String sabreid, Double bmi, Double obesity,
        Double sex, Double age, Double wbc, Double ly, Double mo,
        Double eo, Double ba) { }

    private record PhenotypeRow(Sample sample) {
        boolean complete() {
            return sample.sabreid != null && !sample.sabreid.isBlank()
                && sample.bmi != null && sample.obesity != null
                && sample.sex != null && sample.age != null
                && sample.wbc != null && sample.ly != null && sample.mo != null
                && sample.eo != null && sample.ba != null;
        }
    }

    private record WindowSeed(
        String chromosome, long start, long end, List<VariantRecord> variants) { }
    private record Gene(
        String chromosome, long start, long end, String id, String name) { }
    private record VariantDatum(
        String key, VariantRecord record, VariantStatistics statistics,
        EffectAllele effectAllele, double weight) { }
    private record Window(
        String id, String chromosome, long start, long end, String genes,
        List<String> keys, VariantSet variantSet) { }
    private record PreparedData(
        List<Sample> samples, List<Window> windows,
        List<VariantDatum> variants, Pedigree pedigree) { }
    private record WorkingBinary(
        double[] response, double[] residualDiagonal) { }
    private record Timing(
        String runtime, String model, String method, String backend,
        int measurement, int windows, double seconds,
        double windowsPerSecond) { }
    private record Result(
        String runtime, String model, String windowId, String chromosome,
        long start, long end, String genes, int variants, String method,
        double pValue, double statistic) { }

    private static final class Options {
        private Path genotype = Path.of("D:/Research/topmed/sqtl/sqtl-batch1234/"
            + "freeze.12c.chr22-fhs-framid.vcf.gz");
        private Path phenotype = Path.of("D:/Research/topmed/sqtl/sqtl-batch1234/"
            + "mastermat-batch1234-wcbc-ext.csv");
        private Path pedigree = Path.of("D:/Research/topmed/splicing-bmi/new/"
            + "sabre_ped_0407_v1_rj.csv");
        private Path annotation = Path.of(
            "D:/Research/topmed/gencode.v48.annotation.gene.csv");
        private Path outputDirectory = Path.of(
            "build/benchmarks/topmed-sliding-set");
        private BackendPolicy backend = BackendPolicy.ONEMKL;
        private int windowSize = 2000;
        private int stride = 500;
        private int windowCount = 10;
        private int minimumVariants = 2;
        private int measurements = 3;
        private int skatOSimulations = 10_000;
        private double maximumMaf = 0.01;
        private boolean prepareOnly;
        private final Set<String> pedigreeIds = new LinkedHashSet<>();

        private static Options parse(String[] arguments) {
            Options result = new Options();
            for (int index = 0; index < arguments.length; index++) {
                String key = arguments[index];
                if (key.equals("--prepare-only")) {
                    result.prepareOnly = true;
                    continue;
                }
                if (!key.startsWith("--") || index + 1 >= arguments.length)
                    throw new IllegalArgumentException("invalid option: " + key);
                String value = arguments[++index];
                switch (key.substring(2)) {
                    case "genotype" -> result.genotype = Path.of(value);
                    case "phenotype" -> result.phenotype = Path.of(value);
                    case "pedigree" -> result.pedigree = Path.of(value);
                    case "annotation" -> result.annotation = Path.of(value);
                    case "output-dir" -> result.outputDirectory = Path.of(value);
                    case "backend" -> result.backend = BackendPolicy.valueOf(
                        value.toUpperCase(Locale.ROOT));
                    case "window-size" -> result.windowSize = Integer.parseInt(value);
                    case "stride" -> result.stride = Integer.parseInt(value);
                    case "windows" -> result.windowCount = Integer.parseInt(value);
                    case "minimum-variants" ->
                        result.minimumVariants = Integer.parseInt(value);
                    case "measurements" -> result.measurements = Integer.parseInt(value);
                    case "skato-simulations" ->
                        result.skatOSimulations = Integer.parseInt(value);
                    case "maximum-maf" -> result.maximumMaf = Double.parseDouble(value);
                    default -> throw new IllegalArgumentException(
                        "unknown option: " + key);
                }
            }
            if (result.windowSize < 1 || result.stride < 1
                    || result.windowCount < 1 || result.minimumVariants < 1
                    || result.measurements < 1 || result.skatOSimulations < 1
                    || !(result.maximumMaf > 0 && result.maximumMaf <= 0.5))
                throw new IllegalArgumentException("benchmark options are invalid");
            return result;
        }
    }
}
