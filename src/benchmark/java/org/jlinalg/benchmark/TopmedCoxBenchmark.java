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
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.survival.CoxMixedModel;
import org.jlinalg.survival.CoxMixedOptions;
import org.jlinalg.survival.CoxMixedResult;
import org.jlinalg.survival.CoxOptions;
import org.jlinalg.survival.CoxRandomEffectTerm;
import org.jlinalg.survival.CoxRegression;
import org.jlinalg.survival.CoxResult;
import org.jlinalg.survival.CoxSurvivalData;

/** TOPMed fixed and Gaussian-frailty Cox benchmark against survival/coxme. */
public final class TopmedCoxBenchmark {
    private static volatile double checksum;

    private TopmedCoxBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        Options options = Options.parse(arguments);
        Data data = readAnalysis(options.preparedDirectory(), options.genes(),
            options.maximumRows()).variableGenes();
        try (BackendContext context = BackendContext.select(options.backend())) {
            System.out.printf("backend requested=%s selected=%s device=%s%n",
                context.provenance().requested(),
                context.provenance().selectedBackend(),
                context.provenance().deviceDescription());
        }
        System.out.printf(Locale.ROOT,
            "cohort rows=%d events=%d genes=%d%n",
            data.rows(), data.events(), data.genes());

        Map<String, Fitter> fitters = new LinkedHashMap<>();
        if (options.models().contains("cox"))
            fitters.put("cox", TopmedCoxBenchmark::fixed);
        CoxRandomEffectTerm batch = null;
        if (options.models().contains("coxme")
                || options.models().contains("pedigree")) {
            batch = CoxRandomEffectTerm.independent(
                RandomEffectTerm.randomIntercept("Levy_Set", data.batch()));
        }
        if (options.models().contains("coxme")) {
            CoxRandomEffectTerm retainedBatch = batch;
            fitters.put("coxme", (value, backend) -> mixed(
                value, List.of(retainedBatch),
                CoxMixedOptions.defaults(), backend));
        }
        if (options.models().contains("pedigree")) {
            CoxRandomEffectTerm genetic = pedigreeTerm(data, options);
            CoxRandomEffectTerm retainedBatch = batch;
            CoxMixedOptions mixedOptions = new CoxMixedOptions(
                CoxOptions.defaults(), new double[] {0.5, 0.5},
                30, 1e-4, 1e-8, 1e4);
            fitters.put("pedigree", (value, backend) -> mixed(value,
                List.of(genetic, retainedBatch), mixedOptions, backend));
        }

        List<Timing> timings = new ArrayList<>();
        List<Result> results = new ArrayList<>();
        for (String model : options.models()) {
            Fitter fitter = fitters.get(model);
            if (fitter == null)
                throw new IllegalArgumentException("unknown model: " + model);
            run(data, fitter, options.backend(), options.threads(), 1);
            for (int measurement = 1;
                    measurement <= options.measurements(); measurement++) {
                System.gc();
                long started = System.nanoTime();
                List<Fit> fitted = run(data, fitter, options.backend(),
                    options.threads(), data.genes());
                double seconds = (System.nanoTime() - started) / 1e9;
                timings.add(new Timing(model, options.backend().name(),
                    options.threads(), measurement, data.genes(),
                    data.rows(), seconds));
                for (Fit value : fitted) checksum += value.beta();
                if (measurement == 1)
                    for (int gene = 0; gene < fitted.size(); gene++) {
                        Fit value = fitted.get(gene);
                        results.add(new Result(model, options.threads(),
                            data.featureKeys().get(gene),
                            data.featureIds().get(gene), value.beta(),
                            value.standardError(), value.converged()));
                    }
                System.out.printf(Locale.ROOT,
                    "JLinAlg model=%s threads=%d measurement=%d rows=%d "
                        + "genes=%d seconds=%.6f genes_per_second=%.3f%n",
                    model, options.threads(), measurement, data.rows(),
                    data.genes(), seconds, data.genes() / seconds);
            }
        }
        writeTimings(options.outputPrefix() + "_timings.csv", timings);
        writeResults(options.outputPrefix() + "_results.csv", results);
        if (!Double.isFinite(checksum))
            throw new IllegalStateException("non-finite benchmark checksum");
    }

    private static Fit fixed(GeneData data, BackendPolicy backend) {
        CoxResult fit = CoxRegression.fit(data.survival(), data.design(),
            null, CoxOptions.defaults(), backend);
        return new Fit(fit.beta()[0], fit.standardErrors()[0], fit.converged());
    }

    private static Fit mixed(
            GeneData data, List<CoxRandomEffectTerm> randomEffects,
            CoxMixedOptions options, BackendPolicy backend) {
        CoxMixedResult fit = CoxMixedModel.fit(data.survival(), data.design(),
            randomEffects, null, options, backend);
        return new Fit(fit.beta()[0], fit.standardErrors()[0], fit.converged());
    }

    private static List<Fit> run(
            Data data, Fitter fitter, BackendPolicy backend,
            int threads, int genes) throws Exception {
        if (threads == 1) {
            List<Fit> result = new ArrayList<>(genes);
            for (int gene = 0; gene < genes; gene++)
                result.add(fitter.fit(data.gene(gene), backend));
            return result;
        }
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        try {
            List<Callable<Fit>> tasks = new ArrayList<>(genes);
            for (int gene = 0; gene < genes; gene++) {
                int selected = gene;
                tasks.add(() -> fitter.fit(data.gene(selected), backend));
            }
            List<Fit> result = new ArrayList<>(genes);
            for (Future<Fit> future : executor.invokeAll(tasks))
                result.add(future.get());
            return result;
        } finally {
            executor.shutdown();
        }
    }

    private static CoxRandomEffectTerm pedigreeTerm(
            Data data, Options options) throws IOException {
        List<PedigreeIndividual> all = readPedigreeIndividuals(
            options.preparedDirectory().resolve("pedigree.csv"));
        Map<String, PedigreeIndividual> byId = new LinkedHashMap<>();
        for (PedigreeIndividual value : all) byId.put(value.id(), value);
        Set<String> retained = new LinkedHashSet<>(data.animals());
        boolean changed;
        do {
            changed = false;
            for (String id : List.copyOf(retained)) {
                PedigreeIndividual value = byId.get(id);
                if (value == null)
                    throw new IllegalArgumentException(
                        "animal is absent from pedigree: " + id);
                if (value.sireId() != null)
                    changed |= retained.add(value.sireId());
                if (value.damId() != null)
                    changed |= retained.add(value.damId());
            }
        } while (changed);
        List<PedigreeIndividual> selected = all.stream()
            .filter(value -> retained.contains(value.id())).toList();
        Pedigree pedigree = Pedigree.of(selected);

        List<String> coefficientNames = new ArrayList<>();
        Map<String, Integer> coefficientById = new LinkedHashMap<>();
        for (String animal : data.animals())
            coefficientById.computeIfAbsent(animal, value -> {
                coefficientNames.add(value);
                return coefficientNames.size() - 1;
            });
        int coefficients = coefficientNames.size();
        double[][] incidence = new double[data.rows()][coefficients];
        for (int row = 0; row < data.rows(); row++)
            incidence[row][coefficientById.get(data.animals().get(row))] = 1.0;
        double[] covariance = new double[coefficients * coefficients];
        for (int row = 0; row < coefficients; row++)
            for (int column = 0; column <= row; column++) {
                double value = pedigree.relationship(
                    coefficientNames.get(row), coefficientNames.get(column));
                covariance[row * coefficients + column] = value;
                covariance[column * coefficients + row] = value;
            }
        System.out.printf(Locale.ROOT,
            "pedigree observed=%d ancestor-closure=%d%n",
            coefficients, pedigree.size());
        return CoxRandomEffectTerm.fromCovariance(
            "additive genetic", incidence, covariance, coefficientNames,
            1e-10, options.backend());
    }

    private static Data readAnalysis(
            Path directory, int requestedGenes, int maximumRows)
            throws IOException {
        List<String> featureKeys = new ArrayList<>();
        List<String> featureIds = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                directory.resolve("features.csv"), StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int key = required(header, "feature_key");
            int id = required(header, "feature_id");
            for (String line; featureKeys.size() < requestedGenes
                    && (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                featureKeys.add(row.get(key));
                featureIds.add(row.get(id));
            }
        }
        if (featureKeys.size() != requestedGenes)
            throw new IllegalArgumentException(
                "prepared feature count is too small");

        Map<String, SurvivalRow> survival = new LinkedHashMap<>();
        try (BufferedReader reader = Files.newBufferedReader(
                directory.resolve("survival.csv"), StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int sample = required(header, "SampleName");
            int time = required(header, "time");
            int event = required(header, "event");
            for (String line; (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                survival.put(row.get(sample), new SurvivalRow(
                    number(row.get(time)), number(row.get(event)) == 1.0));
            }
        }

        List<Double> times = new ArrayList<>();
        List<Boolean> events = new ArrayList<>();
        List<double[]> covariateRows = new ArrayList<>();
        List<double[]> predictorRows = new ArrayList<>();
        List<String> batches = new ArrayList<>();
        List<String> animals = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                directory.resolve("analysis.csv"), StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int sample = required(header, "SampleName");
            int batch = required(header, "Batch");
            int animal = required(header, "animal_id");
            int[] covariates = {
                required(header, "Sex"), required(header, "Age"),
                required(header, "WBC_Pred"), required(header, "LY_PER_Pred"),
                required(header, "MO_PER_Pred"), required(header, "EO_PER_Pred"),
                required(header, "BA_PER_Pred")};
            int[] features = featureKeys.stream().mapToInt(
                value -> required(header, value)).toArray();
            for (String line; (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                SurvivalRow outcome = survival.get(row.get(sample));
                if (outcome == null)
                    throw new IllegalArgumentException(
                        "survival row is missing for " + row.get(sample));
                double[] covariate = new double[covariates.length];
                for (int index = 0; index < covariates.length; index++)
                    covariate[index] = number(row.get(covariates[index]));
                double[] predictor = new double[features.length];
                for (int index = 0; index < features.length; index++)
                    predictor[index] = number(row.get(features[index]));
                times.add(outcome.time());
                events.add(outcome.event());
                covariateRows.add(covariate);
                predictorRows.add(predictor);
                batches.add(row.get(batch));
                animals.add(row.get(animal));
            }
        }
        int availableRows = times.size();
        int rows = Math.min(maximumRows, availableRows);
        double[] time = new double[rows];
        boolean[] event = new boolean[rows];
        double[][] covariates = new double[rows][];
        double[][] predictors = new double[rows][];
        List<String> selectedBatches = new ArrayList<>(rows);
        List<String> selectedAnimals = new ArrayList<>(rows);
        int eventCount = 0;
        for (int row = 0; row < rows; row++) {
            int source = rows == availableRows ? row
                : rows == 1 ? 0
                    : (int) Math.round(row * (availableRows - 1.0)
                        / (rows - 1.0));
            time[row] = times.get(source);
            event[row] = events.get(source);
            if (event[row]) eventCount++;
            covariates[row] = covariateRows.get(source);
            predictors[row] = predictorRows.get(source);
            selectedBatches.add(batches.get(source));
            selectedAnimals.add(animals.get(source));
        }
        return new Data(rows, requestedGenes, eventCount,
            CoxSurvivalData.rightCensored(time, event), covariates, predictors,
            List.copyOf(selectedBatches), List.copyOf(selectedAnimals),
            List.copyOf(featureKeys), List.copyOf(featureIds));
    }

    private static List<PedigreeIndividual> readPedigreeIndividuals(Path path)
            throws IOException {
        List<PedigreeIndividual> result = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int id = required(header, "id");
            int sire = required(header, "sire");
            int dam = required(header, "dam");
            for (String line; (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                result.add(new PedigreeIndividual(row.get(id),
                    nullable(row.get(sire)), nullable(row.get(dam))));
            }
        }
        return List.copyOf(result);
    }

    private static String nullable(String value) {
        return value == null || value.isBlank() ? null : value;
    }

    private static double number(String value) {
        double result = Double.parseDouble(value);
        if (!Double.isFinite(result))
            throw new IllegalArgumentException("non-finite prepared value");
        return result;
    }

    private static int required(List<String> header, String name) {
        int result = header.indexOf(name);
        if (result < 0)
            throw new IllegalArgumentException(
                "missing prepared column: " + name);
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

    private static void writeTimings(String filename, List<Timing> values)
            throws IOException {
        Path path = Path.of(filename);
        Files.createDirectories(path.toAbsolutePath().getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,backend,threads,measurement,genes,rows,seconds,genes_per_second\n");
            for (Timing value : values)
                writer.write(String.format(Locale.ROOT,
                    "JLinAlg,%s,%s,%d,%d,%d,%d,%.9f,%.9f%n",
                    value.model(), value.backend(), value.threads(),
                    value.measurement(), value.genes(), value.rows(),
                    value.seconds(), value.genes() / value.seconds()));
        }
    }

    private static void writeResults(String filename, List<Result> values)
            throws IOException {
        Path path = Path.of(filename);
        Files.createDirectories(path.toAbsolutePath().getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,threads,feature_key,feature_id,beta,standard_error,converged\n");
            for (Result value : values)
                writer.write(String.format(Locale.ROOT,
                    "JLinAlg,%s,%d,%s,%s,%.17g,%.17g,%s%n",
                    value.model(), value.threads(), value.featureKey(),
                    value.featureId(), value.beta(), value.standardError(),
                    value.converged()));
        }
    }

    private record SurvivalRow(double time, boolean event) { }
    private record Data(
            int rows, int genes, int events, CoxSurvivalData survival,
            double[][] covariates, double[][] predictors, List<String> batch,
            List<String> animals, List<String> featureKeys,
            List<String> featureIds) {
        GeneData gene(int gene) {
            double[][] design = new double[rows][covariates[0].length + 1];
            for (int row = 0; row < rows; row++) {
                design[row][0] = predictors[row][gene];
                System.arraycopy(covariates[row], 0, design[row], 1,
                    covariates[row].length);
            }
            return new GeneData(survival, design);
        }

        Data variableGenes() {
            List<Integer> retained = new ArrayList<>();
            for (int gene = 0; gene < genes; gene++) {
                double minimum = Double.POSITIVE_INFINITY;
                double maximum = Double.NEGATIVE_INFINITY;
                for (int row = 0; row < rows; row++) {
                    minimum = Math.min(minimum, predictors[row][gene]);
                    maximum = Math.max(maximum, predictors[row][gene]);
                }
                if (maximum > minimum) retained.add(gene);
            }
            if (retained.size() == genes) return this;
            double[][] selected = new double[rows][retained.size()];
            for (int row = 0; row < rows; row++)
                for (int gene = 0; gene < retained.size(); gene++)
                    selected[row][gene] = predictors[row][retained.get(gene)];
            List<String> keys = retained.stream()
                .map(featureKeys::get).toList();
            List<String> ids = retained.stream()
                .map(featureIds::get).toList();
            System.out.printf("excluded constant genes=%d%n",
                genes - retained.size());
            return new Data(rows, retained.size(), events, survival,
                covariates, selected, batch, animals, keys, ids);
        }
    }

    @FunctionalInterface
    private interface Fitter {
        Fit fit(GeneData data, BackendPolicy backend);
    }

    private record Fit(double beta, double standardError, boolean converged) { }
    private record GeneData(CoxSurvivalData survival, double[][] design) { }
    private record Timing(String model, String backend, int threads,
        int measurement, int genes, int rows, double seconds) { }
    private record Result(String model, int threads, String featureKey,
        String featureId, double beta, double standardError,
        boolean converged) { }

    private record Options(
            Path preparedDirectory, int genes, int maximumRows,
            int threads, int measurements, Set<String> models,
            BackendPolicy backend, String outputPrefix) {
        static Options parse(String[] arguments) {
            Map<String, String> values = new LinkedHashMap<>();
            for (int index = 0; index < arguments.length; index += 2) {
                if (index + 1 >= arguments.length
                        || !arguments[index].startsWith("--"))
                    throw new IllegalArgumentException(
                        "arguments must be --name value pairs");
                values.put(arguments[index].substring(2), arguments[index + 1]);
            }
            Path prepared = Path.of(values.getOrDefault(
                "prepared-dir", "build/benchmarks/topmed-obesity"));
            int genes = integer(values, "genes", 20);
            int maximumRows = integer(values, "max-rows", Integer.MAX_VALUE);
            int threads = integer(values, "threads", 1);
            int measurements = integer(values, "measurements", 3);
            Set<String> models = java.util.Collections.unmodifiableSet(
                new LinkedHashSet<>(Arrays.asList(values.getOrDefault(
                    "models", "cox,coxme,pedigree").split(","))));
            BackendPolicy backend = BackendPolicy.valueOf(values.getOrDefault(
                "backend", "preferred").toUpperCase(Locale.ROOT));
            String prefix = values.getOrDefault("output-prefix",
                prepared.resolve("jlinalg_cox_t" + threads).toString());
            return new Options(prepared, genes, maximumRows, threads,
                measurements, models, backend, prefix);
        }

        private static int integer(
                Map<String, String> values, String name, int defaultValue) {
            int value = Integer.parseInt(values.getOrDefault(
                name, Integer.toString(defaultValue)));
            if (value < 1)
                throw new IllegalArgumentException(
                    "--" + name + " must be positive");
            return value;
        }
    }
}
