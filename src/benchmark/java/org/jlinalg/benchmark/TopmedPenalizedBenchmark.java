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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import org.jlinalg.penalized.ElasticNetOptions;
import org.jlinalg.penalized.PenalizedRegression;
import org.jlinalg.penalized.PenalizedRegressionPath;
import org.jlinalg.penalized.PenalizedRegressionResult;

/** Real-data 100-gene BMI penalized-regression benchmark shared with R. */
public final class TopmedPenalizedBenchmark {
    private static volatile double checksum;
    private static final List<String> MODEL_ORDER =
        List.of("ridge", "lasso", "elastic-net");

    private TopmedPenalizedBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        Options options = Options.parse(arguments);
        Topmed100GeneBenchmark.Data source =
            Topmed100GeneBenchmark.readAnalysis(
                options.preparedDirectory(), options.genes());
        Data data = prepare(source);
        List<Timing> timings = new ArrayList<>();
        List<Coefficient> coefficients = new ArrayList<>();
        ElasticNetOptions commonOptions = ElasticNetOptions.builder()
            .penaltyFactors(data.penaltyFactors())
            .parallelism(options.threads())
            .relativeTolerance(options.tolerance()).build();
        long prepareStarted = System.nanoTime();
        PenalizedRegression.Prepared prepared = PenalizedRegression.prepare(
            data.response(), data.design(), data.rows(), data.columns(),
            commonOptions);
        double prepareSeconds = (System.nanoTime() - prepareStarted) / 1e9;

        System.out.printf(Locale.ROOT,
            "samples=%d genes=%d covariates=%d lambdas=%d threads=%d%n",
            data.rows(), data.genes(), data.columns() - data.genes(),
            options.lambdaCount(), options.threads());
        System.out.printf(Locale.ROOT, "preprocessing_seconds=%.6f%n",
            prepareSeconds);
        for (String model : MODEL_ORDER) {
            if (!options.models().contains(model)) continue;
            double alpha = alpha(model);
            for (int warmup = 0; warmup < options.warmups(); warmup++) {
                consume(prepared.automaticPath(options.lambdaCount(),
                    options.minimumRatio(), alpha));
            }
            for (int measurement = 1;
                    measurement <= options.measurements(); measurement++) {
                System.gc();
                long started = System.nanoTime();
                PenalizedRegressionPath path = null;
                for (int repetition = 0;
                        repetition < options.repetitions(); repetition++) {
                    path = prepared.automaticPath(options.lambdaCount(),
                        options.minimumRatio(), alpha);
                }
                double seconds = (System.nanoTime() - started) / 1e9
                    / options.repetitions();
                consume(path);
                PenalizedRegressionPath measuredPath = path;
                double meanIterations = java.util.stream.IntStream.range(
                    0, measuredPath.size()).map(
                        index -> measuredPath.fit(index).iterations())
                    .average().orElse(Double.NaN);
                int maximumIterations = java.util.stream.IntStream.range(
                    0, measuredPath.size()).map(
                        index -> measuredPath.fit(index).iterations())
                    .max().orElse(0);
                timings.add(new Timing("JLinAlg", model, options.threads(),
                    measurement, data.rows(), data.genes(), path.size(), seconds));
                if (measurement == 1) {
                    addCoefficients(coefficients, model, data, path);
                }
                System.out.printf(Locale.ROOT,
                    "JLinAlg model=%s threads=%d measurement=%d "
                        + "seconds=%.6f lambda_fits_per_second=%.1f "
                        + "mean_iterations=%.1f max_iterations=%d%n",
                    model, options.threads(), measurement, seconds,
                    path.size() / seconds, meanIterations, maximumIterations);
            }
        }
        runSuite(options, prepared, timings);
        writeTimings(options.outputPrefix() + "_timings.csv", timings);
        writeCoefficients(options.outputPrefix() + "_coefficients.csv", coefficients);
        System.out.printf(Locale.ROOT, "checksum=%.17g%n", checksum);
    }

    private static void runSuite(Options options,
            PenalizedRegression.Prepared prepared, List<Timing> timings) {
        List<String> models = MODEL_ORDER.stream()
            .filter(options.models()::contains).toList();
        for (int warmup = 0; warmup < options.warmups(); warmup++) {
            models.forEach(model -> consume(prepared.automaticPath(
                options.lambdaCount(), options.minimumRatio(), alpha(model))));
        }
        for (int measurement = 1;
                measurement <= options.measurements(); measurement++) {
            long started = System.nanoTime();
            List<PenalizedRegressionPath> paths = List.of();
            ForkJoinPool pool = options.threads() == 1 || models.size() == 1
                ? null : new ForkJoinPool(
                    Math.min(options.threads(), models.size()));
            try {
                for (int repetition = 0;
                        repetition < options.repetitions(); repetition++) {
                    if (pool == null) {
                        paths = models.stream().map(model ->
                            prepared.automaticPath(options.lambdaCount(),
                                options.minimumRatio(), alpha(model))).toList();
                    } else {
                        paths = pool.submit(() -> models.parallelStream()
                            .map(model -> prepared.automaticPath(
                                options.lambdaCount(), options.minimumRatio(),
                                alpha(model))).toList()).join();
                    }
                }
            } finally {
                if (pool != null) pool.shutdown();
            }
            double seconds = (System.nanoTime() - started) / 1e9
                / options.repetitions();
            paths.forEach(TopmedPenalizedBenchmark::consume);
            timings.add(new Timing("JLinAlg", "all", options.threads(),
                measurement, prepared.rows(), options.genes(),
                options.lambdaCount() * paths.size(), seconds));
            System.out.printf(Locale.ROOT,
                "JLinAlg model=all threads=%d measurement=%d seconds=%.6f "
                    + "lambda_fits_per_second=%.1f%n", options.threads(),
                measurement, seconds,
                options.lambdaCount() * paths.size() / seconds);
        }
    }

    private static Data prepare(Topmed100GeneBenchmark.Data source) {
        int covariates = source.baseColumns() - 1;
        int columns = source.genes() + covariates;
        double[] design = new double[source.rows() * columns];
        for (int row = 0; row < source.rows(); row++) {
            System.arraycopy(source.predictors(), row * source.genes(), design,
                row * columns, source.genes());
            System.arraycopy(source.baseDesign(), row * source.baseColumns() + 1,
                design, row * columns + source.genes(), covariates);
        }
        double[] penaltyFactors = new double[columns];
        Arrays.fill(penaltyFactors, 0, source.genes(),
            (double) columns / source.genes());
        List<String> names = new ArrayList<>(source.featureKeys());
        names.addAll(List.of("Sex", "Age", "WBC_Pred", "LY_PER_Pred",
            "MO_PER_Pred", "EO_PER_Pred", "BA_PER_Pred"));
        return new Data(source.rows(), source.genes(), columns,
            source.response(), design, penaltyFactors, List.copyOf(names));
    }

    private static double alpha(String model) {
        return switch (model) {
            case "ridge" -> 0.0;
            case "lasso" -> 1.0;
            case "elastic-net" -> 0.5;
            default -> throw new IllegalArgumentException("unknown model: " + model);
        };
    }

    private static void consume(PenalizedRegressionPath path) {
        PenalizedRegressionResult fit = path.fit(path.size() - 1);
        checksum += fit.intercept() + fit.objective()
            + fit.coefficients()[0] + path.lambdas()[0];
    }

    private static void addCoefficients(List<Coefficient> output, String model,
            Data data, PenalizedRegressionPath path) {
        PenalizedRegressionResult fit = path.fit(path.size() - 1);
        output.add(new Coefficient("JLinAlg", model, optionsAlpha(model),
            path.lambdas()[path.size() - 1], "(Intercept)", fit.intercept()));
        double[] beta = fit.coefficients();
        for (int column = 0; column < beta.length; column++) {
            output.add(new Coefficient("JLinAlg", model, optionsAlpha(model),
                path.lambdas()[path.size() - 1], data.names().get(column),
                beta[column]));
        }
    }

    private static double optionsAlpha(String model) { return alpha(model); }

    private static void writeTimings(String filename, List<Timing> values)
            throws IOException {
        Path path = Path.of(filename);
        Files.createDirectories(path.toAbsolutePath().getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,threads,measurement,samples,genes,"
                + "lambdas,seconds,lambda_fits_per_second\n");
            for (Timing value : values) writer.write(String.format(Locale.ROOT,
                "%s,%s,%d,%d,%d,%d,%d,%.9f,%.9f%n", value.runtime(),
                value.model(), value.threads(), value.measurement(), value.samples(),
                value.genes(), value.lambdas(), value.seconds(),
                value.lambdas() / value.seconds()));
        }
    }

    private static void writeCoefficients(String filename,
            List<Coefficient> values) throws IOException {
        Path path = Path.of(filename);
        Files.createDirectories(path.toAbsolutePath().getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(
                path, StandardCharsets.UTF_8)) {
            writer.write("runtime,model,alpha,lambda,term,beta\n");
            for (Coefficient value : values) writer.write(String.format(Locale.ROOT,
                "%s,%s,%.17g,%.17g,%s,%.17g%n", value.runtime(), value.model(),
                value.alpha(), value.lambda(), value.term(), value.beta()));
        }
    }

    private record Data(int rows, int genes, int columns, double[] response,
        double[] design, double[] penaltyFactors, List<String> names) { }
    private record Timing(String runtime, String model, int threads,
        int measurement, int samples, int genes, int lambdas, double seconds) { }
    private record Coefficient(String runtime, String model, double alpha,
        double lambda, String term, double beta) { }

    private record Options(Path preparedDirectory, int genes, int threads,
            int measurements, int warmups, int repetitions, int lambdaCount,
            double minimumRatio,
            double tolerance, Set<String> models, String outputPrefix) {
        static Options parse(String[] arguments) {
            Map<String, String> values = new LinkedHashMap<>();
            for (int index = 0; index < arguments.length; index += 2) {
                if (index + 1 >= arguments.length
                        || !arguments[index].startsWith("--")) {
                    throw new IllegalArgumentException(
                        "arguments must be --name value pairs");
                }
                values.put(arguments[index].substring(2), arguments[index + 1]);
            }
            Path prepared = Path.of(values.getOrDefault(
                "prepared-dir", "build/benchmarks/topmed100"));
            int genes = integer(values, "genes", 100);
            int threads = integer(values, "threads", 1);
            int measurements = integer(values, "measurements", 3);
            int warmups = integer(values, "warmups", 5);
            int repetitions = integer(values, "repetitions", 10);
            int lambdas = integer(values, "lambdas", 100);
            double ratio = Double.parseDouble(values.getOrDefault(
                "minimum-ratio", "1e-4"));
            double tolerance = Double.parseDouble(values.getOrDefault(
                "tolerance", "1e-8"));
            Set<String> models = java.util.Collections.unmodifiableSet(
                new LinkedHashSet<>(Arrays.asList(values.getOrDefault("models",
                    "ridge,lasso,elastic-net").split(","))));
            String prefix = values.getOrDefault("output-prefix",
                prepared.resolve("jlinalg_penalized_t" + threads).toString());
            return new Options(prepared, genes, threads, measurements, warmups,
                repetitions, lambdas, ratio, tolerance, models, prefix);
        }

        private static int integer(Map<String, String> values,
                String name, int defaultValue) {
            int value = Integer.parseInt(values.getOrDefault(
                name, Integer.toString(defaultValue)));
            if (value < 1) throw new IllegalArgumentException(
                "--" + name + " must be positive");
            return value;
        }
    }
}
