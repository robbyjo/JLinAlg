/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.timeseries.Arima;
import org.jlinalg.timeseries.ArimaOptions;
import org.jlinalg.timeseries.ArimaOrder;
import org.jlinalg.timeseries.ArimaResult;
import org.jlinalg.timeseries.ArimaSelectionResult;
import org.jlinalg.timeseries.AutomaticArima;
import org.jlinalg.timeseries.ExactArma;
import org.jlinalg.timeseries.ExactArmaResult;
import org.jlinalg.timeseries.SeasonalArimaOrder;

/** Real-data macrobenchmarks for conditional ARIMA and exact ARMA fitting. */
public final class TimeSeriesBenchmark {
    private static final String RESOURCE_ROOT = "/timeseries/";
    private static final Logger BACKEND_LOGGER =
        Logger.getLogger("jdistlib.accelerator.ComputeBackends");
    private static volatile double checksum;

    private TimeSeriesBenchmark() { }

    public static void main(String[] arguments) {
        BACKEND_LOGGER.setLevel(Level.WARNING);
        int warmups = integerProperty("jlinalg.benchmark.warmups", 1);
        int measurements = integerProperty("jlinalg.benchmark.measurements", 3);
        Set<String> selected = selectedCases();
        Map<String, Series> data = loadData();
        List<BenchmarkCase> cases = cases(data);
        if (!selected.isEmpty()) {
            cases = cases.stream()
                .filter(value -> selected.contains(value.name()))
                .toList();
        }
        if (cases.isEmpty()) {
            throw new IllegalArgumentException(
                "no benchmark cases matched jlinalg.benchmark.cases; available cases: "
                    + cases(data).stream().map(BenchmarkCase::name)
                        .collect(Collectors.joining(",")));
        }

        System.out.println("benchmark,dataset,model,observations,period,"
            + "warmups,measurements,median_seconds,function_evaluations,converged,result_metric");
        for (BenchmarkCase benchmark : cases) {
            measure(benchmark, warmups, measurements);
        }
    }

    private static List<BenchmarkCase> cases(Map<String, Series> data) {
        Series nile = data.get("Nile");
        Series internet = data.get("WWWusage");
        Series airline = data.get("AirPassengers");
        Series gas = data.get("UKgas");
        Series temperature = data.get("nottem");
        Series sunspots = data.get("sunspots");
        ArimaOptions airlineOptions = seasonal(0, 1, 1, 12);
        ArimaOptions gasOptions = seasonal(1, 1, 1, 4);
        ArimaOptions temperatureOptions = seasonal(1, 0, 0, 12);

        return List.of(
            conditional("conditional_ar2_nile", nile, "AR(2)",
                nile.values(), ArimaOrder.ar(2), ArimaOptions.defaults()),
            conditional("conditional_arima310_wwwusage", internet,
                "ARIMA(3,1,0)", internet.values(),
                ArimaOrder.arima(3, 1, 0), ArimaOptions.defaults()),
            conditional("conditional_sarima_airline", airline,
                "log SARIMA(0,1,1)(0,1,1)[12]", log(airline.values()),
                ArimaOrder.arima(0, 1, 1), airlineOptions),
            conditional("conditional_sarima_ukgas", gas,
                "log SARIMA(1,1,1)(1,1,1)[4]", log(gas.values()),
                ArimaOrder.arima(1, 1, 1), gasOptions),
            conditional("conditional_seasonal_ar_nottem", temperature,
                "AR(2) SAR(1)[12]", temperature.values(),
                ArimaOrder.ar(2), temperatureOptions),
            conditional("conditional_arma21_sunspots", sunspots,
                "ARMA(2,1)", sunspots.values(), ArimaOrder.arma(2, 1),
                ArimaOptions.defaults()),
            new BenchmarkCase("automatic_nile", nile, "auto p<=3,d<=1,q<=3",
                () -> summarize(AutomaticArima.select(
                    nile.values(), 3, 1, 3))),
            new BenchmarkCase("exact_ar2_nile", nile, "exact AR(2)",
                () -> summarize(ExactArma.fit(nile.values(), ArimaOrder.ar(2),
                    true, BackendPolicy.CPU)))
        );
    }

    private static BenchmarkCase conditional(
            String name, Series series, String model, double[] values,
            ArimaOrder order, ArimaOptions options) {
        return new BenchmarkCase(name, series, model,
            () -> summarize(Arima.fit(values, order, options)));
    }

    private static ArimaOptions seasonal(int ar, int differences, int ma,
            int period) {
        return ArimaOptions.builder()
            .seasonalOrder(SeasonalArimaOrder.of(ar, differences, ma, period))
            .build();
    }

    private static FitSummary summarize(ArimaResult result) {
        return new FitSummary(result.aicc() + result.innovationVariance(),
            result.functionEvaluations(), result.converged());
    }

    private static FitSummary summarize(ArimaSelectionResult result) {
        ArimaResult best = result.bestModel();
        return new FitSummary(best.aicc() + result.candidates().size(), -1,
            result.candidates().stream().allMatch(candidate -> candidate.converged()));
    }

    private static FitSummary summarize(ExactArmaResult result) {
        return new FitSummary(result.aic() + result.innovationVariance(),
            result.functionEvaluations(), result.converged());
    }

    private static void measure(
            BenchmarkCase benchmark, int warmups, int measurements) {
        for (int iteration = 0; iteration < warmups; iteration++) {
            consume(benchmark.operation().run());
        }
        double[] seconds = new double[measurements];
        FitSummary result = null;
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime();
            result = benchmark.operation().run();
            seconds[iteration] = (System.nanoTime() - started) / 1e9;
            consume(result);
        }
        Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        Series series = benchmark.series();
        System.out.printf(Locale.ROOT,
            "%s,%s,%s,%d,%d,%d,%d,%.6f,%d,%s,%.12g%n",
            benchmark.name(), series.name(), csv(benchmark.model()),
            series.values().length, series.period(), warmups, measurements,
            median, result.functionEvaluations(), result.converged(),
            result.checksum());
    }

    private static void consume(FitSummary result) {
        checksum += result.checksum();
    }

    private static String csv(String value) {
        return "\"" + value.replace("\"", "\"\"") + "\"";
    }

    private static Map<String, Series> loadData() {
        Map<String, Series> result = new LinkedHashMap<>();
        add(result, read("AirPassengers", 12, 144,
            "bdb98adbd418a6de6842a742e0602f363c3b841c26677e7e61a1e055e9509bd8"));
        add(result, read("Nile", 1, 100,
            "d0452bea38c61e796a4eeb950bf91d20fb5c7f13d5822eadf5990fe54f9c8d07"));
        add(result, read("nottem", 12, 240,
            "634e2e374c5d59c1332e6e6de207d691cb4ebdcd007f7f0e2993ac2a19c90473"));
        add(result, read("sunspots", 12, 2_820,
            "8db04db1d4406f8a005c00fc3a4647de1845fd3e7e7bbde4ad176a19037c9e60"));
        add(result, read("UKgas", 4, 108,
            "69ec013fbeeaf956a42603aaf8f68d943a3aab5f90b6ed89c5de2ed02eba59ff"));
        add(result, read("WWWusage", 1, 100,
            "245774828aa5bc64292f29db15434004b7472022b73b81a3d5d1c02f55c42e08"));
        return result;
    }

    private static void add(Map<String, Series> destination, Series series) {
        destination.put(series.name(), series);
    }

    private static Series read(
            String name, int period, int expectedSize, String expectedSha256) {
        String resource = RESOURCE_ROOT + name + ".csv";
        InputStream input = TimeSeriesBenchmark.class.getResourceAsStream(resource);
        if (input == null) {
            throw new IllegalStateException("missing benchmark resource " + resource);
        }
        List<Double> values = new ArrayList<>();
        double previousTime = Double.NEGATIVE_INFINITY;
        try (input) {
            byte[] bytes = input.readAllBytes();
            String actualSha256 = HexFormat.of().formatHex(
                MessageDigest.getInstance("SHA-256").digest(bytes));
            if (!expectedSha256.equals(actualSha256)) {
                throw new IllegalStateException(resource + " SHA-256 is "
                    + actualSha256 + "; expected " + expectedSha256);
            }
            BufferedReader reader = new BufferedReader(new InputStreamReader(
                new ByteArrayInputStream(bytes), StandardCharsets.UTF_8));
            String header = reader.readLine();
            if (!"rownames,time,value".equals(header)) {
                throw new IllegalStateException("unexpected header in " + resource);
            }
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split(",", -1);
                if (fields.length != 3) {
                    throw new IllegalStateException("malformed row in " + resource);
                }
                int row = Integer.parseInt(fields[0]);
                double time = Double.parseDouble(fields[1]);
                double value = Double.parseDouble(fields[2]);
                if (row != values.size() + 1 || !(time > previousTime)
                        || !Double.isFinite(value)) {
                    throw new IllegalStateException("invalid row in " + resource);
                }
                previousTime = time;
                values.add(value);
            }
        } catch (IOException | NoSuchAlgorithmException
                | NumberFormatException exception) {
            throw new IllegalStateException("cannot read " + resource, exception);
        }
        if (values.size() != expectedSize) {
            throw new IllegalStateException(resource + " has " + values.size()
                + " rows; expected " + expectedSize);
        }
        double[] array = values.stream().mapToDouble(Double::doubleValue).toArray();
        return new Series(name, period, array);
    }

    private static double[] log(double[] values) {
        double[] result = values.clone();
        for (int index = 0; index < result.length; index++) {
            if (!(result[index] > 0.0)) {
                throw new IllegalArgumentException("log series must be positive");
            }
            result[index] = Math.log(result[index]);
        }
        return result;
    }

    private static Set<String> selectedCases() {
        String property = System.getProperty("jlinalg.benchmark.cases", "").trim();
        if (property.isEmpty()) return Set.of();
        return Arrays.stream(property.split(","))
            .map(String::trim).filter(value -> !value.isEmpty())
            .collect(Collectors.toUnmodifiableSet());
    }

    private static int integerProperty(String name, int defaultValue) {
        int value = Integer.getInteger(name, defaultValue);
        if (value < 1) throw new IllegalArgumentException(name + " must be positive");
        return value;
    }

    @FunctionalInterface
    private interface Operation { FitSummary run(); }

    private record BenchmarkCase(
        String name, Series series, String model, Operation operation) { }

    private record FitSummary(
        double checksum, int functionEvaluations, boolean converged) { }

    private record Series(String name, int period, double[] values) {
        private Series {
            values = values.clone();
        }
    }
}
