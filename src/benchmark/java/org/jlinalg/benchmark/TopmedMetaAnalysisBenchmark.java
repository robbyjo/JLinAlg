/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.meta.MetaAnalysis;
import org.jlinalg.meta.MetaAnalysisBatchResult;
import org.jlinalg.meta.MetaAnalysisOptions;
import org.jlinalg.meta.MetaAnalysisResult;
import org.jlinalg.meta.MetaStudy;

/** Four-cohort TOPMed BMI fixed- and random-effects benchmark against R. */
public final class TopmedMetaAnalysisBenchmark {
    private static final String[] COHORTS = {"CARDIA", "FHS", "JHS", "WHI"};
    private static final String[] FILES = {
        "splicing-bmi-cardia-geneadj.csv",
        "splicing-bmi-fhs-batch1234-geneadj.csv",
        "splicing-bmi-jhs-geneadj.csv",
        "splicing-bmi-whi-geneadj.csv"
    };
    private static volatile double checksum;

    private TopmedMetaAnalysisBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        Options options = Options.parse(arguments);
        Files.createDirectories(options.outputDirectory());
        Data data = read(options.inputDirectory(), options.limit());
        writePrepared(options.outputDirectory().resolve("analysis.csv"), data);
        System.out.printf(Locale.ROOT,
            "prepared analyses=%d studies=%d input=%s%n",
            data.analyses(), COHORTS.length, options.inputDirectory());
        if (options.prepareOnly()) return;

        List<Timing> timings = new ArrayList<>();
        List<NamedResults> outputs = new ArrayList<>();
        org.jlinalg.meta.PreparedMetaAnalysisBatch prepared =
            MetaAnalysis.prepareBatch(data.effects(), data.standardErrors(),
                data.analyses(), COHORTS.length);
        for (String model : options.models()) {
            MetaAnalysisOptions fitOptions = model.equals("fixed")
                ? MetaAnalysisOptions.fixedEffect()
                : MetaAnalysisOptions.randomEffects();
            run(data, prepared, model, fitOptions, options, timings, outputs);
        }
        writeTimings(options.outputDirectory().resolve(
            "jlinalg_" + options.engine() + "_t" + options.threads()
                + "_timings.csv"), timings);
        writeResults(options.outputDirectory().resolve(
            "jlinalg_" + options.engine() + "_t" + options.threads()
                + "_results.csv"), data.keys(), outputs);
        if (!Double.isFinite(checksum))
            throw new IllegalStateException("non-finite benchmark checksum");
    }

    private static void run(
            Data data, org.jlinalg.meta.PreparedMetaAnalysisBatch prepared,
            String model, MetaAnalysisOptions fitOptions, Options options,
            List<Timing> timings, List<NamedResults> outputs) {
        fit(data, prepared, fitOptions, options);
        Results first = null;
        for (int measurement = 1;
                measurement <= options.measurements(); measurement++) {
            System.gc();
            long started = System.nanoTime();
            Results result = fit(data, prepared, fitOptions, options);
            double seconds = (System.nanoTime() - started) / 1e9;
            if (first == null) first = result;
            checksum += Arrays.stream(result.pValues()).sum()
                + Arrays.stream(result.tauSquared()).sum();
            int usedThreads = options.engine().equals("batch")
                ? options.threads() : 1;
            timings.add(new Timing("JLinAlg", options.engine(), model,
                fitOptions.tauSquaredEstimator().name(), usedThreads,
                measurement, data.analyses(), seconds,
                data.analyses() / seconds));
            System.out.printf(Locale.ROOT,
                "JLinAlg engine=%s model=%s threads=%d measurement=%d "
                    + "analyses=%d seconds=%.6f analyses_per_second=%.3f%n",
                options.engine(), model, usedThreads, measurement,
                data.analyses(), seconds, data.analyses() / seconds);
        }
        outputs.add(new NamedResults(model, first));
    }

    private static Results fit(
            Data data, org.jlinalg.meta.PreparedMetaAnalysisBatch prepared,
            MetaAnalysisOptions options, Options benchmarkOptions) {
        if (benchmarkOptions.engine().equals("batch"))
            return results(prepared.fit(options, benchmarkOptions.threads()));
        int analyses = data.analyses();
        Results output = new Results(analyses);
        for (int analysis = 0; analysis < analyses; analysis++) {
            int offset = analysis * COHORTS.length;
            List<MetaStudy> studies = List.of(
                study(0, offset, data), study(1, offset, data),
                study(2, offset, data), study(3, offset, data));
            MetaAnalysisResult value = MetaAnalysis.fit(
                studies, options, BackendPolicy.CPU);
            output.pooledEffects()[analysis] = value.pooledEffectSize();
            output.standardErrors()[analysis] = value.standardError();
            output.pValues()[analysis] = value.pValue();
            output.cochranQ()[analysis] = value.cochranQ();
            output.cochranQPValues()[analysis] = value.cochranQPValue();
            output.tauSquared()[analysis] = value.tauSquared();
            output.iSquared()[analysis] = value.iSquared();
        }
        return output;
    }

    private static MetaStudy study(int study, int offset, Data data) {
        return new MetaStudy(COHORTS[study], data.effects()[offset + study],
            data.standardErrors()[offset + study]);
    }

    private static Results results(MetaAnalysisBatchResult value) {
        return new Results(value.pooledEffectSizes(), value.standardErrors(),
            value.pValues(), value.cochranQ(), value.cochranQPValues(),
            value.tauSquared(), value.iSquared());
    }

    private static Data read(Path inputDirectory, int limit) throws IOException {
        Map<String, double[]> common = new HashMap<>(350_000);
        for (int cohort = 0; cohort < FILES.length; cohort++) {
            Path path = inputDirectory.resolve(FILES[cohort]);
            readCohort(path, cohort, common);
            if (cohort > 0) {
                int column = cohort * 2;
                common.entrySet().removeIf(entry ->
                    !Double.isFinite(entry.getValue()[column]));
            }
        }
        List<Map.Entry<String, double[]>> rows = new ArrayList<>(common.entrySet());
        rows.sort(Comparator.comparing(Map.Entry::getKey));
        int analyses = limit == 0 ? rows.size() : Math.min(limit, rows.size());
        List<String> keys = new ArrayList<>(analyses);
        double[] effects = new double[analyses * COHORTS.length];
        double[] errors = new double[effects.length];
        for (int analysis = 0; analysis < analyses; analysis++) {
            Map.Entry<String, double[]> row = rows.get(analysis);
            keys.add(row.getKey());
            for (int cohort = 0; cohort < COHORTS.length; cohort++) {
                effects[analysis * COHORTS.length + cohort] =
                    row.getValue()[cohort * 2];
                errors[analysis * COHORTS.length + cohort] =
                    row.getValue()[cohort * 2 + 1];
            }
        }
        return new Data(keys, effects, errors);
    }

    private static void readCohort(
            Path path, int cohort, Map<String, double[]> common)
            throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            String header = reader.readLine();
            if (header == null || !header.startsWith(",Fx_BMI,SE_BMI,"))
                throw new IllegalArgumentException(
                    "unexpected cohort header in " + path);
            String line;
            while ((line = reader.readLine()) != null) {
                int first = line.indexOf(',');
                int second = line.indexOf(',', first + 1);
                int third = line.indexOf(',', second + 1);
                if (first <= 0 || second <= first + 1 || third <= second + 1)
                    continue;
                String key = line.substring(0, first);
                double effect = Double.parseDouble(
                    line.substring(first + 1, second));
                double error = Double.parseDouble(
                    line.substring(second + 1, third));
                if (!Double.isFinite(effect) || !(error > 0.0)
                        || !Double.isFinite(error)) continue;
                double[] values;
                if (cohort == 0) {
                    values = new double[COHORTS.length * 2];
                    Arrays.fill(values, Double.NaN);
                    common.put(key, values);
                } else {
                    values = common.get(key);
                    if (values == null) continue;
                }
                values[cohort * 2] = effect;
                values[cohort * 2 + 1] = error;
            }
        }
    }

    private static void writePrepared(Path path, Data data) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("feature");
            for (String cohort : COHORTS)
                writer.write(",effect_" + cohort + ",se_" + cohort);
            writer.newLine();
            for (int analysis = 0; analysis < data.analyses(); analysis++) {
                writer.write(data.keys().get(analysis));
                int offset = analysis * COHORTS.length;
                for (int study = 0; study < COHORTS.length; study++) {
                    writer.write("," + data.effects()[offset + study]);
                    writer.write("," + data.standardErrors()[offset + study]);
                }
                writer.newLine();
            }
        }
    }

    private static void writeTimings(Path path, List<Timing> values)
            throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,engine,model,tau_estimator,threads,measurement,"
                + "analyses,seconds,analyses_per_second\n");
            for (Timing value : values) writer.write(String.format(Locale.ROOT,
                "%s,%s,%s,%s,%d,%d,%d,%.9f,%.9f%n", value.runtime(),
                value.engine(), value.model(), value.tauEstimator(),
                value.threads(), value.measurement(), value.analyses(),
                value.seconds(), value.analysesPerSecond()));
        }
    }

    private static void writeResults(
            Path path, List<String> keys, List<NamedResults> values)
            throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,feature,beta,se,p_value,cochran_q,"
                + "cochran_q_p_value,tau_squared,i_squared\n");
            for (NamedResults named : values) {
                Results result = named.results();
                for (int index = 0; index < keys.size(); index++)
                    writer.write(String.format(Locale.ROOT,
                        "JLinAlg,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g%n",
                        named.model(), keys.get(index),
                        result.pooledEffects()[index],
                        result.standardErrors()[index],
                        result.pValues()[index], result.cochranQ()[index],
                        result.cochranQPValues()[index],
                        result.tauSquared()[index], result.iSquared()[index]));
            }
        }
    }

    private record Data(List<String> keys, double[] effects,
        double[] standardErrors) {
        int analyses() { return keys.size(); }
    }
    private record Results(double[] pooledEffects, double[] standardErrors,
            double[] pValues, double[] cochranQ, double[] cochranQPValues,
            double[] tauSquared, double[] iSquared) {
        Results(int analyses) {
            this(new double[analyses], new double[analyses],
                new double[analyses], new double[analyses],
                new double[analyses], new double[analyses],
                new double[analyses]);
        }
    }
    private record NamedResults(String model, Results results) { }
    private record Timing(String runtime, String engine, String model,
        String tauEstimator, int threads, int measurement, int analyses,
        double seconds, double analysesPerSecond) { }

    private record Options(Path inputDirectory, Path outputDirectory,
            int measurements, int threads, int limit, String engine,
            List<String> models, boolean prepareOnly) {
        static Options parse(String[] arguments) {
            Map<String, String> values = new LinkedHashMap<>();
            for (int index = 0; index < arguments.length; index += 2) {
                if (index + 1 >= arguments.length
                        || !arguments[index].startsWith("--"))
                    throw new IllegalArgumentException(
                        "arguments must be --name value pairs");
                values.put(arguments[index].substring(2), arguments[index + 1]);
            }
            Path input = Path.of(values.getOrDefault("input-dir",
                "D:/Research/topmed/splicing-bmi/new"));
            Path output = Path.of(values.getOrDefault("output-dir",
                "build/benchmarks/topmed-meta-analysis"));
            int measurements = positive(values, "measurements", 3);
            int threads = positive(values, "threads",
                Math.min(8, Runtime.getRuntime().availableProcessors()));
            int limit = Integer.parseInt(values.getOrDefault("limit", "0"));
            if (limit < 0)
                throw new IllegalArgumentException("--limit must not be negative");
            String engine = values.getOrDefault("engine", "batch");
            if (!engine.equals("batch") && !engine.equals("scalar"))
                throw new IllegalArgumentException(
                    "--engine must be batch or scalar");
            List<String> models = Arrays.asList(values.getOrDefault(
                "models", "fixed,reml").split(","));
            if (models.stream().anyMatch(value ->
                    !value.equals("fixed") && !value.equals("reml")))
                throw new IllegalArgumentException(
                    "--models values must be fixed or reml");
            boolean prepareOnly = Boolean.parseBoolean(
                values.getOrDefault("prepare-only", "false"));
            return new Options(input, output, measurements, threads, limit,
                engine, List.copyOf(models), prepareOnly);
        }

        private static int positive(
                Map<String, String> values, String name, int fallback) {
            int result = Integer.parseInt(values.getOrDefault(
                name, Integer.toString(fallback)));
            if (result < 1)
                throw new IllegalArgumentException(
                    "--" + name + " must be positive");
            return result;
        }
    }
}
