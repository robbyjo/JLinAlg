/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.GamResult;
import org.jlinalg.gam.PreparedGamPredictorScan;
import org.jlinalg.reml.RemlOptions;

/** Real-data TOPMed BMI P-spline benchmark shared with the mgcv harness. */
public final class TopmedGamBenchmark {
    private static volatile double checksum;

    private TopmedGamBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        Options options = Options.parse(arguments);
        Topmed100GeneBenchmark.Data data =
            Topmed100GeneBenchmark.readAnalysis(
                options.preparedDirectory(), options.genes());
        List<Timing> timings = new ArrayList<>();
        Result[] firstResults = null;
        try (PreparedGamPredictorScan prepared = new PreparedGamPredictorScan(
                data.response(), data.baseDesign(), data.rows(),
                data.baseColumns(), RemlOptions.defaults(), options.backend())) {
            fit(prepared, data, 0);
            for (int measurement = 1;
                    measurement <= options.measurements(); measurement++) {
                System.gc();
                long started = System.nanoTime();
                Result[] results = run(prepared, data, options.threads());
                double seconds = (System.nanoTime() - started) / 1e9;
                timings.add(new Timing("JLinAlg", options.backend().name(),
                    options.threads(), measurement, data.genes(), seconds,
                    data.genes() / seconds));
                if (measurement == 1) firstResults = results;
                checksum += Arrays.stream(results)
                    .mapToDouble(Result::fittedChecksum).sum();
                System.out.printf(Locale.ROOT,
                    "JLinAlg model=gam backend=%s threads=%d measurement=%d "
                        + "genes=%d seconds=%.6f genes_per_second=%.3f%n",
                    options.backend(), options.threads(), measurement,
                    data.genes(), seconds, data.genes() / seconds);
            }
        }
        writeTimings(Path.of(options.outputPrefix() + "_timings.csv"), timings);
        writeResults(Path.of(options.outputPrefix() + "_results.csv"), firstResults);
        if (!Double.isFinite(checksum))
            throw new IllegalStateException("non-finite benchmark checksum");
    }

    private static Result[] run(
            PreparedGamPredictorScan prepared,
            Topmed100GeneBenchmark.Data data,
            int threads) {
        Result[] results = new Result[data.genes()];
        if (threads == 1) {
            for (int gene = 0; gene < data.genes(); gene++)
                results[gene] = fit(prepared, data, gene);
            return results;
        }
        ForkJoinPool pool = new ForkJoinPool(Math.min(threads, data.genes()));
        try {
            pool.submit(() -> IntStream.range(0, data.genes()).parallel()
                .forEach(gene -> results[gene] =
                    fit(prepared, data, gene))).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("GAM benchmark interrupted", exception);
        } catch (ExecutionException exception) {
            Throwable cause = exception.getCause();
            if (cause instanceof RuntimeException runtime) throw runtime;
            throw new IllegalStateException("GAM benchmark failed", cause);
        } finally {
            pool.shutdownNow();
        }
        return results;
    }

    private static Result fit(
            PreparedGamPredictorScan prepared,
            Topmed100GeneBenchmark.Data data,
            int gene) {
        double[] predictor = new double[data.rows()];
        for (int row = 0; row < data.rows(); row++)
            predictor[row] = data.predictors()[row * data.genes() + gene];
        GamResult result = prepared.fit("s(omics)", predictor, 10);
        return new Result(data.featureKeys().get(gene),
            data.featureIds().get(gene),
            result.smoothTerms().get(0).effectiveDegreesOfFreedom(),
            result.smoothTerms().get(0).smoothingParameter(),
            result.mixedModel().reml().logLikelihood(),
            Arrays.stream(result.fittedValues()).sum(),
            Arrays.stream(result.residuals())
                .map(value -> value * value).sum());
    }

    private static void writeTimings(Path path, List<Timing> timings)
            throws IOException {
        Files.createDirectories(path.toAbsolutePath().getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,backend,threads,measurement,genes,"
                + "seconds,genes_per_second\n");
            for (Timing value : timings) writer.write(String.format(Locale.ROOT,
                "%s,gam,%s,%d,%d,%d,%.9f,%.9f%n", value.runtime(),
                value.backend(), value.threads(), value.measurement(),
                value.genes(), value.seconds(), value.genesPerSecond()));
        }
    }

    private static void writeResults(Path path, Result[] results)
            throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,feature_key,feature_id,edf,"
                + "smoothing_parameter,log_likelihood,fitted_checksum,"
                + "residual_sum_squares\n");
            for (Result value : results) writer.write(String.format(Locale.ROOT,
                "JLinAlg,gam,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g%n",
                value.featureKey(), value.featureId(), value.edf(),
                value.smoothingParameter(), value.logLikelihood(),
                value.fittedChecksum(), value.residualSumSquares()));
        }
    }

    private record Timing(String runtime, String backend, int threads,
        int measurement, int genes, double seconds, double genesPerSecond) { }
    private record Result(String featureKey, String featureId, double edf,
        double smoothingParameter, double logLikelihood, double fittedChecksum,
        double residualSumSquares) { }

    private record Options(Path preparedDirectory, int genes, int threads,
            int measurements, BackendPolicy backend, String outputPrefix) {
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
                "prepared-dir", "build/benchmarks/topmed100"));
            int genes = positive(values, "genes", 100);
            int threads = positive(values, "threads", 1);
            int measurements = positive(values, "measurements", 1);
            BackendPolicy backend = BackendPolicy.valueOf(values.getOrDefault(
                "backend", "preferred").toUpperCase(Locale.ROOT));
            String prefix = values.getOrDefault("output-prefix",
                prepared.resolve("jlinalg_gam_t" + threads).toString());
            return new Options(prepared, genes, threads, measurements,
                backend, prefix);
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
