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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.GammScanResult;
import org.jlinalg.gam.PreparedGammPredictorScan;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.reml.RemlOptions;

/** TOPMed 100-gene batch and additive-pedigree Gaussian GAMM benchmark. */
public final class TopmedGammBenchmark {
    private static volatile double checksum;

    private TopmedGammBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        Options options = Options.parse(arguments);
        Topmed100GeneBenchmark.Data data =
            Topmed100GeneBenchmark.readAnalysis(
                options.preparedDirectory(), options.genes());
        List<Timing> timings = new ArrayList<>();
        List<Result> results = new ArrayList<>();
        for (String model : options.models()) {
            try (PreparedGammPredictorScan scan =
                    prepare(model, data, options)) {
                GammScanResult warm = fit(scan, data, 0).fit();
                scan.warmStart(warm);
                System.out.printf(Locale.ROOT,
                    "prepared model=%s rows=%d genes=%d components=%s%n",
                    model, data.rows(), data.genes(),
                    warm.reml().componentNames());
                for (int measurement = 1;
                        measurement <= options.measurements();
                        measurement++) {
                    System.gc();
                    long started = System.nanoTime();
                    IndexedFit[] fitted =
                        run(scan, data, options.threads());
                    double seconds =
                        (System.nanoTime() - started) / 1e9;
                    timings.add(new Timing("JLinAlg", model,
                        options.backend().name(), options.threads(),
                        measurement, data.genes(), seconds,
                        data.genes() / seconds));
                    if (measurement == 1) {
                        for (IndexedFit value : fitted) {
                            results.add(value.result());
                        }
                    }
                    checksum += Arrays.stream(fitted)
                        .mapToDouble(value ->
                            value.result().fittedChecksum()).sum();
                    double meanEvaluations = Arrays.stream(fitted)
                        .mapToInt(value -> value.fit().reml()
                            .iterations())
                        .average().orElse(Double.NaN);
                    int maximumEvaluations = Arrays.stream(fitted)
                        .mapToInt(value -> value.fit().reml().iterations())
                        .max().orElse(0);
                    long nonconverged = Arrays.stream(fitted)
                        .filter(value -> !value.fit().reml().converged())
                        .count();
                    System.out.printf(Locale.ROOT,
                        "JLinAlg model=%s backend=%s threads=%d "
                            + "measurement=%d genes=%d seconds=%.6f "
                            + "genes_per_second=%.3f mean_evaluations=%.2f "
                            + "max_evaluations=%d nonconverged=%d%n",
                        model, options.backend(), options.threads(),
                        measurement, data.genes(), seconds,
                        data.genes() / seconds, meanEvaluations,
                        maximumEvaluations, nonconverged);
                }
            }
        }
        writeTimings(Path.of(
            options.outputPrefix() + "_timings.csv"), timings);
        writeResults(Path.of(
            options.outputPrefix() + "_results.csv"), results);
        if (!Double.isFinite(checksum)) {
            throw new IllegalStateException(
                "non-finite benchmark checksum");
        }
    }

    private static PreparedGammPredictorScan prepare(
            String model,
            Topmed100GeneBenchmark.Data data,
            Options options) throws IOException {
        if (model.equals("batch-gamm")) {
            return PreparedGammPredictorScan.batch(
                data.response(), data.baseDesign(), data.rows(),
                data.baseColumns(), "Batch", data.batch(), 10,
                RemlOptions.defaults(), options.backend());
        }
        if (model.equals("pedigree-gamm")) {
            PedigreeRandomEffectTerm pedigree =
                pedigreeTerm(data, options.preparedDirectory()
                    .resolve("pedigree.csv"));
            RandomEffectTerm batch =
                RandomEffectTerm.randomIntercept(
                    "Batch", data.batch());
            return new PreparedGammPredictorScan(
                data.response(), data.baseDesign(), data.rows(),
                data.baseColumns(),
                List.of(pedigree.randomEffect(), batch),
                List.of(pedigree.precision(),
                    SparsePrecisionMatrix.identity(
                        batch.coefficients())),
                10, RemlOptions.defaults(), options.backend());
        }
        throw new IllegalArgumentException(
            "unknown GAMM model: " + model);
    }

    private static IndexedFit[] run(
            PreparedGammPredictorScan scan,
            Topmed100GeneBenchmark.Data data,
            int threads) {
        IndexedFit[] result = new IndexedFit[data.genes()];
        if (threads == 1) {
            for (int gene = 0; gene < data.genes(); gene++) {
                result[gene] = fit(scan, data, gene);
            }
            return result;
        }
        ForkJoinPool pool =
            new ForkJoinPool(Math.min(threads, data.genes()));
        try {
            pool.submit(() -> IntStream.range(0, data.genes())
                .parallel().forEach(gene ->
                    result[gene] = fit(scan, data, gene))).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException(
                "GAMM benchmark interrupted", exception);
        } catch (ExecutionException exception) {
            Throwable cause = exception.getCause();
            if (cause instanceof RuntimeException runtime) {
                throw runtime;
            }
            throw new IllegalStateException(
                "GAMM benchmark failed", cause);
        } finally {
            pool.shutdownNow();
        }
        return result;
    }

    private static IndexedFit fit(
            PreparedGammPredictorScan scan,
            Topmed100GeneBenchmark.Data data,
            int gene) {
        double[] predictor = new double[data.rows()];
        for (int row = 0; row < data.rows(); row++) {
            predictor[row] =
                data.predictors()[row * data.genes() + gene];
        }
        GammScanResult fit =
            scan.fit("s(omics)", predictor);
        double[] variances = fit.reml().varianceComponents();
        Map<String, Double> byName = new LinkedHashMap<>();
        for (int index = 0;
                index < fit.reml().componentNames().size(); index++) {
            byName.put(fit.reml().componentNames().get(index),
                variances[index]);
        }
        Result result = new Result(
            "JLinAlg",
            fit.reml().componentNames().contains("additive genetic")
                ? "pedigree-gamm" : "batch-gamm",
            data.featureKeys().get(gene),
            data.featureIds().get(gene),
            fit.smoothTerm().effectiveDegreesOfFreedom(),
            fit.smoothTerm().smoothingParameter(),
            fit.reml().logLikelihood(),
            Arrays.stream(fit.fittedValues()).sum(),
            Arrays.stream(fit.residuals())
                .map(value -> value * value).sum(),
            byName.getOrDefault(
                "additive genetic", Double.NaN),
            byName.get("Batch"),
            byName.get("__scan_smooth__"),
            byName.get("residual"));
        return new IndexedFit(fit, result);
    }

    private static PedigreeRandomEffectTerm pedigreeTerm(
            Topmed100GeneBenchmark.Data data, Path path)
            throws IOException {
        List<PedigreeIndividual> all =
            readPedigreeIndividuals(path);
        Map<String, PedigreeIndividual> byId =
            new LinkedHashMap<>();
        for (PedigreeIndividual value : all) {
            byId.put(value.id(), value);
        }
        Set<String> retained =
            new LinkedHashSet<>(data.animal());
        boolean changed;
        do {
            changed = false;
            for (String id : List.copyOf(retained)) {
                PedigreeIndividual value = byId.get(id);
                if (value == null) {
                    throw new IllegalArgumentException(
                        "animal is absent from pedigree: " + id);
                }
                if (value.sireId() != null) {
                    changed |= retained.add(value.sireId());
                }
                if (value.damId() != null) {
                    changed |= retained.add(value.damId());
                }
            }
        } while (changed);
        List<PedigreeIndividual> selected = all.stream()
            .filter(value -> retained.contains(value.id())).toList();
        System.out.printf(Locale.ROOT,
            "pedigree observed=%d ancestor_closure=%d "
                + "inbreeding=assumed-zero%n",
            new LinkedHashSet<>(data.animal()).size(),
            selected.size());
        return PedigreeRandomEffectTerm.ofUninbred(
            "additive genetic", data.animal(), selected);
    }

    private static List<PedigreeIndividual> readPedigreeIndividuals(
            Path path) throws IOException {
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

    private static void writeTimings(
            Path path, List<Timing> values) throws IOException {
        Files.createDirectories(path.toAbsolutePath().getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,backend,threads,measurement,"
                + "genes,seconds,genes_per_second\n");
            for (Timing value : values) {
                writer.write(String.format(Locale.ROOT,
                    "%s,%s,%s,%d,%d,%d,%.9f,%.9f%n",
                    value.runtime(), value.model(), value.backend(),
                    value.threads(), value.measurement(), value.genes(),
                    value.seconds(), value.genesPerSecond()));
            }
        }
    }

    private static void writeResults(
            Path path, List<Result> values) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,feature_key,feature_id,edf,"
                + "smoothing_parameter,log_likelihood,fitted_checksum,"
                + "residual_sum_squares,pedigree_variance,batch_variance,"
                + "smooth_variance,residual_variance\n");
            for (Result value : values) {
                writer.write(String.format(Locale.ROOT,
                    "%s,%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,"
                        + "%.17g,%.17g,%.17g,%.17g%n",
                    value.runtime(), value.model(), value.featureKey(),
                    value.featureId(), value.edf(),
                    value.smoothingParameter(), value.logLikelihood(),
                    value.fittedChecksum(),
                    value.residualSumSquares(),
                    value.pedigreeVariance(), value.batchVariance(),
                    value.smoothVariance(), value.residualVariance()));
            }
        }
    }

    private static List<String> csv(String line) {
        if (line == null) {
            throw new IllegalArgumentException("empty CSV file");
        }
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
                } else {
                    quoted = !quoted;
                }
            } else if (value == ',' && !quoted) {
                result.add(field.toString());
                field.setLength(0);
            } else {
                field.append(value);
            }
        }
        result.add(field.toString());
        return result;
    }

    private static int required(
            List<String> header, String name) {
        int result = header.indexOf(name);
        if (result < 0) {
            throw new IllegalArgumentException(
                "missing prepared column: " + name);
        }
        return result;
    }

    private static String nullable(String value) {
        return value == null || value.isBlank() ? null : value;
    }

    private record IndexedFit(GammScanResult fit, Result result) { }
    private record Timing(
        String runtime, String model, String backend, int threads,
        int measurement, int genes, double seconds,
        double genesPerSecond) { }
    private record Result(
        String runtime, String model, String featureKey,
        String featureId, double edf, double smoothingParameter,
        double logLikelihood, double fittedChecksum,
        double residualSumSquares, double pedigreeVariance,
        double batchVariance, double smoothVariance,
        double residualVariance) { }

    private record Options(
            Path preparedDirectory,
            int genes,
            int threads,
            int measurements,
            Set<String> models,
            BackendPolicy backend,
            String outputPrefix) {
        static Options parse(String[] arguments) {
            Map<String, String> values = new LinkedHashMap<>();
            for (int index = 0;
                    index < arguments.length; index += 2) {
                if (index + 1 >= arguments.length
                        || !arguments[index].startsWith("--")) {
                    throw new IllegalArgumentException(
                        "arguments must be --name value pairs");
                }
                values.put(arguments[index].substring(2),
                    arguments[index + 1]);
            }
            Path prepared = Path.of(values.getOrDefault(
                "prepared-dir",
                "build/benchmarks/topmed100"));
            int genes = positive(values, "genes", 100);
            int threads = positive(values, "threads", 1);
            int measurements =
                positive(values, "measurements", 1);
            Set<String> models =
                java.util.Collections.unmodifiableSet(
                    new LinkedHashSet<>(Arrays.asList(
                        values.getOrDefault("models",
                            "batch-gamm,pedigree-gamm")
                            .split(","))));
            BackendPolicy backend = BackendPolicy.valueOf(
                values.getOrDefault("backend", "preferred")
                    .toUpperCase(Locale.ROOT));
            String prefix = values.getOrDefault(
                "output-prefix",
                prepared.resolve(
                    "jlinalg_gamm_t" + threads).toString());
            return new Options(prepared, genes, threads,
                measurements, models, backend, prefix);
        }

        private static int positive(
                Map<String, String> values,
                String name, int fallback) {
            int result = Integer.parseInt(values.getOrDefault(
                name, Integer.toString(fallback)));
            if (result < 1) {
                throw new IllegalArgumentException(
                    "--" + name + " must be positive");
            }
            return result;
        }
    }
}
