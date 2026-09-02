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
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jlinalg.association.AssociationBatchResult;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.association.AssociationFailurePolicy;
import org.jlinalg.association.AssociationFitter;
import org.jlinalg.association.FastOlsAssociation;
import org.jlinalg.association.ParallelAssociationEngine;
import org.jlinalg.association.VariableMissingPolicy;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.pedigree.SparsePedigreeMixedModel;
import org.jlinalg.reml.RemlOptions;

/** Real-data 100-gene BMI association benchmark shared with the R harness. */
public final class Topmed100GeneBenchmark {
    private static volatile double checksum;

    private Topmed100GeneBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        Logger.getLogger("").setLevel(Level.WARNING);
        Arrays.stream(Logger.getLogger("").getHandlers())
            .forEach(handler -> handler.setLevel(Level.WARNING));
        Options options = Options.parse(arguments);
        Data data = readAnalysis(options.preparedDirectory, options.genes);
        AssociationEngineOptions engine = new AssociationEngineOptions(
            options.threads, 16, options.backend,
            AssociationFailurePolicy.FAIL_FAST, VariableMissingPolicy.ERROR);
        try (BackendContext context = BackendContext.select(options.backend)) {
            System.out.printf("backend requested=%s selected=%s device=%s "
                    + "accelerated=%s automatic=%s%n",
                context.provenance().requested(),
                context.provenance().selectedBackend(),
                context.provenance().deviceDescription(),
                context.provenance().accelerated(),
                context.provenance().automaticRouting());
        }
        RandomEffectTerm batch = RandomEffectTerm.randomIntercept(
            "Batch", data.batch());
        RemlOptions remlOptions = RemlOptions.defaults();
        Map<String, AssociationFitter> fitters = new LinkedHashMap<>();
        List<AutoCloseable> preparedModels = new ArrayList<>();
        FastOlsAssociation.PreparedPredictorScan preparedOls = null;
        if (options.models.contains("ols")) {
            preparedOls = FastOlsAssociation.preparePredictors(
                data.response(), data.baseDesign(), data.rows(), data.baseColumns(),
                null, null, OlsOptions.defaults(), engine);
            preparedModels.add(preparedOls);
        }
        if (options.models.contains("reml")) {
            SparseLinearMixedModel.Prepared preparedBatch =
                SparseLinearMixedModel.prepare(data.rows(), List.of(batch),
                    remlOptions, options.backend);
            preparedBatch.warmStart(preparedBatch.fit(
                data.response(), data.baseDesign(), data.baseColumns())
                .varianceComponents());
            preparedModels.add(preparedBatch);
            fitters.put("reml", (response, design, rows, columns, backend) ->
                preparedBatch.fit(response, design, columns)
                    .associationStatistics());
        }

        if (options.models.contains("pedigree")) {
            Pedigree pedigree = readPedigree(
                options.preparedDirectory.resolve("pedigree.csv"));
            PedigreeRandomEffectTerm animal = PedigreeRandomEffectTerm.of(
                "animal_id", data.animal(), pedigree);
            SparseLinearMixedModel.Prepared preparedPedigree =
                SparsePedigreeMixedModel.prepare(data.rows(), List.of(animal),
                    List.of(batch), remlOptions, options.backend);
            preparedPedigree.warmStart(preparedPedigree.fit(
                data.response(), data.baseDesign(), data.baseColumns())
                .varianceComponents());
            preparedModels.add(preparedPedigree);
            fitters.put("pedigree", (response, design, rows, columns, backend) ->
                preparedPedigree.fit(response, design, columns)
                    .associationStatistics());
            System.out.printf(Locale.ROOT,
                "loaded pedigree individuals=%d, observations=%d%n",
                pedigree.size(), data.rows());
        }

        List<Timing> timings = new ArrayList<>();
        List<Result> results = new ArrayList<>();
        for (String model : options.models) {
            if (!model.equals("ols") && !fitters.containsKey(model))
                throw new IllegalArgumentException("unknown model: " + model);
            run(model, data, firstPredictor(data), engine, fitters.get(model),
                preparedOls);
            for (int measurement = 1; measurement <= options.measurements;
                    measurement++) {
                System.gc();
                long started = System.nanoTime();
                AssociationBatchResult result = run(
                    model, data, data.predictors(), engine, fitters.get(model),
                    preparedOls);
                double seconds = (System.nanoTime() - started) / 1e9;
                consume(result);
                timings.add(new Timing("JLinAlg", model, options.backend.name(), options.threads,
                    measurement, data.genes(), seconds,
                    data.genes() / seconds));
                if (measurement == 1) addResults(
                    results, model, options.threads, data, result);
                System.out.printf(Locale.ROOT,
                    "JLinAlg model=%s threads=%d measurement=%d genes=%d "
                        + "seconds=%.6f genes_per_second=%.3f failures=%d%n",
                    model, options.threads, measurement, data.genes(), seconds,
                    data.genes() / seconds, result.failures().size());
            }
        }
        writeTimings(options.outputPrefix + "_timings.csv", timings);
        writeResults(options.outputPrefix + "_results.csv", results);
        for (AutoCloseable prepared : preparedModels)
            prepared.close();
        if (!Double.isFinite(checksum))
            throw new IllegalStateException("non-finite benchmark checksum");
    }

    private static AssociationBatchResult run(
            String model, Data data, double[] predictors,
            AssociationEngineOptions engine, AssociationFitter fitter,
            FastOlsAssociation.PreparedPredictorScan preparedOls) {
        int genes = predictors.length / data.rows();
        List<String> names = data.featureKeys().subList(0, genes);
        if (model.equals("ols")) {
            return preparedOls.scan(predictors, genes, names);
        }
        AssociationEngineOptions exactFitEngine = engine.withChunkSize(1);
        return ParallelAssociationEngine.scanPredictors(
            data.response(), data.baseDesign(), data.rows(), data.baseColumns(),
            predictors, genes, names, fitter, exactFitEngine);
    }

    private static double[] firstPredictor(Data data) {
        double[] result = new double[data.rows()];
        for (int row = 0; row < data.rows(); row++)
            result[row] = data.predictors()[row * data.genes()];
        return result;
    }

    private static void consume(AssociationBatchResult result) {
        double[] beta = result.beta();
        checksum += beta.length == 0 ? 0.0 : beta[beta.length - 1];
    }

    private static void addResults(
            List<Result> destination, String model, int threads,
            Data data, AssociationBatchResult source) {
        double[] beta = source.beta();
        double[] standardErrors = source.standardErrors();
        for (int index = 0; index < data.genes(); index++) {
            destination.add(new Result("JLinAlg", model, threads,
                data.featureKeys().get(index), data.featureIds().get(index),
                beta[index], standardErrors[index]));
        }
    }

    private static Data readAnalysis(Path directory, int requestedGenes)
            throws IOException {
        List<String> featureKeys = new ArrayList<>();
        List<String> featureIds = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                directory.resolve("features.csv"), StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int keyIndex = required(header, "feature_key");
            int idIndex = required(header, "feature_id");
            for (String line; featureKeys.size() < requestedGenes
                    && (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                featureKeys.add(row.get(keyIndex));
                featureIds.add(row.get(idIndex));
            }
        }
        if (featureKeys.size() != requestedGenes)
            throw new IllegalArgumentException("prepared feature count is too small");

        List<double[]> designs = new ArrayList<>();
        List<double[]> predictorRows = new ArrayList<>();
        List<Double> responses = new ArrayList<>();
        List<String> batches = new ArrayList<>();
        List<String> animals = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                directory.resolve("analysis.csv"), StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int response = required(header, "BMI");
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
                double[] design = new double[covariates.length + 1];
                design[0] = 1.0;
                for (int index = 0; index < covariates.length; index++)
                    design[index + 1] = number(row.get(covariates[index]));
                double[] predictors = new double[features.length];
                for (int index = 0; index < features.length; index++)
                    predictors[index] = number(row.get(features[index]));
                responses.add(number(row.get(response)));
                designs.add(design);
                predictorRows.add(predictors);
                batches.add(row.get(batch));
                animals.add(row.get(animal));
            }
        }
        int rows = responses.size();
        int baseColumns = designs.get(0).length;
        double[] response = new double[rows];
        double[] baseDesign = new double[rows * baseColumns];
        double[] predictors = new double[rows * requestedGenes];
        for (int row = 0; row < rows; row++) {
            response[row] = responses.get(row);
            System.arraycopy(designs.get(row), 0, baseDesign,
                row * baseColumns, baseColumns);
            System.arraycopy(predictorRows.get(row), 0, predictors,
                row * requestedGenes, requestedGenes);
        }
        return new Data(rows, requestedGenes, baseColumns, response,
            baseDesign, predictors, List.copyOf(batches), List.copyOf(animals),
            List.copyOf(featureKeys), List.copyOf(featureIds));
    }

    private static Pedigree readPedigree(Path path) throws IOException {
        List<PedigreeIndividual> individuals = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                path, StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int id = required(header, "id");
            int sire = required(header, "sire");
            int dam = required(header, "dam");
            for (String line; (line = reader.readLine()) != null;) {
                List<String> row = csv(line);
                individuals.add(new PedigreeIndividual(row.get(id),
                    nullable(row.get(sire)), nullable(row.get(dam))));
            }
        }
        return Pedigree.of(individuals);
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
        if (result < 0) throw new IllegalArgumentException(
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
            writer.write("runtime,model,backend,threads,measurement,genes,seconds,genes_per_second\n");
            for (Timing value : values) writer.write(String.format(Locale.ROOT,
                "%s,%s,%s,%d,%d,%d,%.9f,%.9f%n", value.runtime(), value.model(),
                value.backend(), value.threads(), value.measurement(), value.genes(),
                value.seconds(), value.genesPerSecond()));
        }
    }

    private static void writeResults(String filename, List<Result> values)
            throws IOException {
        Path path = Path.of(filename);
        Files.createDirectories(path.toAbsolutePath().getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,threads,feature_key,feature_id,beta,standard_error\n");
            for (Result value : values) writer.write(String.format(Locale.ROOT,
                "%s,%s,%d,%s,%s,%.17g,%.17g%n", value.runtime(), value.model(),
                value.threads(), value.featureKey(), value.featureId(),
                value.beta(), value.standardError()));
        }
    }

    private record Data(
        int rows, int genes, int baseColumns, double[] response,
        double[] baseDesign, double[] predictors, List<String> batch,
        List<String> animal, List<String> featureKeys,
        List<String> featureIds) { }
    private record Timing(String runtime, String model, String backend, int threads,
        int measurement, int genes, double seconds, double genesPerSecond) { }
    private record Result(String runtime, String model, int threads,
        String featureKey, String featureId, double beta, double standardError) { }

    private record Options(
            Path preparedDirectory, int genes, int threads, int measurements,
            Set<String> models, BackendPolicy backend, String outputPrefix) {
        static Options parse(String[] arguments) {
            Map<String, String> values = new LinkedHashMap<>();
            for (int index = 0; index < arguments.length; index += 2) {
                if (index + 1 >= arguments.length || !arguments[index].startsWith("--"))
                    throw new IllegalArgumentException("arguments must be --name value pairs");
                values.put(arguments[index].substring(2), arguments[index + 1]);
            }
            Path prepared = Path.of(values.getOrDefault(
                "prepared-dir", "build/benchmarks/topmed100"));
            int genes = integer(values, "genes", 100);
            int threads = integer(values, "threads", 1);
            int measurements = integer(values, "measurements", 1);
            Set<String> models = java.util.Collections.unmodifiableSet(
                new LinkedHashSet<>(Arrays.asList(values.getOrDefault(
                    "models", "ols,reml,pedigree").split(","))));
            BackendPolicy backend = BackendPolicy.valueOf(values.getOrDefault(
                "backend", "cpu").toUpperCase(Locale.ROOT));
            String prefix = values.getOrDefault("output-prefix",
                prepared.resolve("jlinalg_t" + threads).toString());
            return new Options(prepared, genes, threads, measurements,
                models, backend, prefix);
        }

        private static int integer(
                Map<String, String> values, String name, int defaultValue) {
            int value = Integer.parseInt(values.getOrDefault(
                name, Integer.toString(defaultValue)));
            if (value < 1) throw new IllegalArgumentException(
                "--" + name + " must be positive");
            return value;
        }
    }
}
