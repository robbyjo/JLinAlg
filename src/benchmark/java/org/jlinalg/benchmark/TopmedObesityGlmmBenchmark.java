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
import org.jlinalg.association.AssociationModels;
import org.jlinalg.association.ParallelAssociationEngine;
import org.jlinalg.association.VariableMissingPolicy;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glmm.GlmmLaplaceOptions;
import org.jlinalg.glmm.SparseGlmmLaplace;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

/** Exact logistic GLM and first-order Laplace GLMM TOPMed benchmark. */
public final class TopmedObesityGlmmBenchmark {
    private static volatile double checksum;

    private TopmedObesityGlmmBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        Options options = Options.parse(arguments);
        Data data = readAnalysis(options.preparedDirectory, options.genes);
        AssociationEngineOptions engine = new AssociationEngineOptions(
            options.threads, 1, options.backend,
            AssociationFailurePolicy.FAIL_FAST, VariableMissingPolicy.ERROR);
        try (BackendContext context = BackendContext.select(options.backend)) {
            System.out.printf("backend requested=%s selected=%s device=%s%n",
                context.provenance().requested(),
                context.provenance().selectedBackend(),
                context.provenance().deviceDescription());
        }

        Map<String, AssociationFitter> fitters = new LinkedHashMap<>();
        List<SparseGlmmLaplace.Prepared> preparedScans = new ArrayList<>();
        if (options.models.contains("glm")) {
            fitters.put("glm", AssociationModels.glm(
                GlmFamilies.binomial(), GlmOptions.defaults()));
        }
        RandomEffectTerm batch = null;
        if (options.models.contains("glmm") || options.models.contains("pedigree")) {
            batch = RandomEffectTerm.randomIntercept("Levy_Set", data.batch());
        }
        if (options.models.contains("glmm")) {
            SparseGlmmLaplace.Prepared prepared = SparseGlmmLaplace.prepare(
                data.rows(), GlmFamilies.binomial(), List.of(batch),
                GlmmLaplaceOptions.defaults(), options.backend);
            preparedScans.add(prepared);
            fitters.put("glmm", laplace(prepared));
        }
        if (options.models.contains("pedigree")) {
            Pedigree pedigree = readPedigree(
                options.preparedDirectory.resolve("pedigree.csv"));
            PedigreeRandomEffectTerm genetic = PedigreeRandomEffectTerm.of(
                "additive genetic", data.animal(), pedigree);
            SparseGlmmLaplace.Prepared prepared =
                SparseGlmmLaplace.prepareWithPrecision(data.rows(),
                    GlmFamilies.binomial(),
                    List.of(batch, genetic.randomEffect()),
                    List.of(SparsePrecisionMatrix.identity(
                        batch.coefficients()), genetic.precision()),
                    GlmmLaplaceOptions.defaults(), options.backend);
            preparedScans.add(prepared);
            fitters.put("pedigree", laplace(prepared));
            System.out.printf(Locale.ROOT,
                "loaded pedigree individuals=%d observations=%d%n",
                pedigree.size(), data.rows());
        }

        try {
            List<Timing> timings = new ArrayList<>();
            List<Result> results = new ArrayList<>();
            for (String model : options.models) {
                AssociationFitter fitter = fitters.get(model);
                if (fitter == null) throw new IllegalArgumentException(
                    "unknown model: " + model);
                scan(data, firstPredictor(data), fitter, engine);
                for (int measurement = 1; measurement <= options.measurements;
                        measurement++) {
                    System.gc();
                    long started = System.nanoTime();
                    AssociationBatchResult result = scan(
                        data, data.predictors(), fitter, engine);
                    double seconds = (System.nanoTime() - started) / 1e9;
                    consume(result);
                    timings.add(new Timing(model, options.backend.name(),
                        options.threads, measurement, data.genes(), seconds));
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
        } finally {
            for (int index = preparedScans.size() - 1; index >= 0; index--)
                preparedScans.get(index).close();
        }
        if (!Double.isFinite(checksum))
            throw new IllegalStateException("non-finite benchmark checksum");
    }

    private static AssociationFitter laplace(
            SparseGlmmLaplace.Prepared prepared) {
        return (response, design, rows, columns, ignoredBackend) ->
            prepared.fit(response, design, columns).associationStatistics();
    }

    private static AssociationBatchResult scan(
            Data data, double[] predictors, AssociationFitter fitter,
            AssociationEngineOptions engine) {
        int genes = predictors.length / data.rows();
        return ParallelAssociationEngine.scanPredictors(
            data.response(), data.baseDesign(), data.rows(), data.baseColumns(),
            predictors, genes, data.featureKeys().subList(0, genes),
            fitter, engine);
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
            destination.add(new Result(model, threads,
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
            throw new IllegalArgumentException("prepared feature count is too small");

        List<double[]> designs = new ArrayList<>();
        List<double[]> predictorRows = new ArrayList<>();
        List<Double> responses = new ArrayList<>();
        List<String> batches = new ArrayList<>();
        List<String> animals = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(
                directory.resolve("analysis.csv"), StandardCharsets.UTF_8)) {
            List<String> header = csv(reader.readLine());
            int response = required(header, "Obesity");
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
                "JLinAlg,%s,%s,%d,%d,%d,%.9f,%.9f%n", value.model(),
                value.backend(), value.threads(), value.measurement(), value.genes(),
                value.seconds(), value.genes() / value.seconds()));
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
                "JLinAlg,%s,%d,%s,%s,%.17g,%.17g%n", value.model(),
                value.threads(), value.featureKey(), value.featureId(),
                value.beta(), value.standardError()));
        }
    }

    private record Data(
        int rows, int genes, int baseColumns, double[] response,
        double[] baseDesign, double[] predictors, List<String> batch,
        List<String> animal, List<String> featureKeys, List<String> featureIds) { }
    private record Timing(String model, String backend, int threads,
        int measurement, int genes, double seconds) { }
    private record Result(String model, int threads, String featureKey,
        String featureId, double beta, double standardError) { }

    private record Options(
            Path preparedDirectory, int genes, int threads, int measurements,
            Set<String> models, BackendPolicy backend, String outputPrefix) {
        static Options parse(String[] arguments) {
            Map<String, String> values = new LinkedHashMap<>();
            for (int index = 0; index < arguments.length; index += 2) {
                if (index + 1 >= arguments.length || !arguments[index].startsWith("--"))
                    throw new IllegalArgumentException(
                        "arguments must be --name value pairs");
                values.put(arguments[index].substring(2), arguments[index + 1]);
            }
            Path prepared = Path.of(values.getOrDefault(
                "prepared-dir", "build/benchmarks/topmed100"));
            int genes = integer(values, "genes", 20);
            int threads = integer(values, "threads", 1);
            int measurements = integer(values, "measurements", 3);
            Set<String> models = java.util.Collections.unmodifiableSet(
                new LinkedHashSet<>(Arrays.asList(values.getOrDefault(
                    "models", "glm,glmm,pedigree").split(","))));
            BackendPolicy backend = BackendPolicy.valueOf(values.getOrDefault(
                "backend", "preferred").toUpperCase(Locale.ROOT));
            String prefix = values.getOrDefault("output-prefix",
                prepared.resolve("jlinalg_obesity_t" + threads).toString());
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
