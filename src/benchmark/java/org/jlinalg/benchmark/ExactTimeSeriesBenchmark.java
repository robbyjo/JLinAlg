/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.timeseries.*;

/** Same-estimand full ML fits; reports accuracy and consumes every fitted result.
 * Run from repository root with test resources available on disk. */
public final class ExactTimeSeriesBenchmark {
    private static volatile double checksum;
    private ExactTimeSeriesBenchmark() { }
    public static void main(String[] args) throws Exception {
        int warmups = Integer.getInteger("jlinalg.benchmark.warmups", 3);
        int repetitions = Integer.getInteger("jlinalg.benchmark.measurements", 7);
        System.out.println("case,java_median_ms,loglik_error,max_coefficient_error,max_forecast_error,max_se_error,converged,checksum");
        for (String name : new String[] {"stationary", "stationary-missing", "integrated-ar", "integrated-ma",
                "integrated-missing", "seasonal", "seasonal-missing", "seasonal-only", "twice-integrated", "drift",
                "seasonal-drift", "seasonal-ar-ma", "monthly-missing"}) {
            Path root = Path.of("src/test/resources/timeseries");
            double[] y = Files.readAllLines(root.resolve(name + "-series.csv")).stream().skip(1)
                .mapToDouble(Double::parseDouble).toArray();
            Map<String, Double> ref = new LinkedHashMap<>();
            Files.readAllLines(root.resolve(name + "-reference.csv")).stream().skip(1).forEach(line -> {
                String[] fields = line.split(","); ref.put(fields[0], Double.parseDouble(fields[1]));
            });
            ArimaOrder order = new ArimaOrder(ref.get("p").intValue(), ref.get("d").intValue(), ref.get("q").intValue());
            SeasonalArimaOrder seasonal = new SeasonalArimaOrder(ref.get("P").intValue(), ref.get("D").intValue(),
                ref.get("Q").intValue(), ref.get("period").intValue());
            ArimaOptions options = ArimaOptions.builder().seasonalOrder(seasonal).includeDrift(ref.get("drift") == 1)
                .optimizationTolerance(1e-9).build();
            double[] times = new double[repetitions];
            double ll = 0, variance = 0;
            double[] coefficients = null;
            ArimaForecast forecast = null;
            boolean converged = true;
            for (int iteration = -warmups; iteration < repetitions; iteration++) {
                long start = System.nanoTime();
                if (order.differences() + seasonal.differences() == 0) {
                    ExactArmaResult f = ExactArma.fit(y, order, false, BackendPolicy.CPU);
                    ll = f.logLikelihood(); variance = f.innovationVariance();
                    coefficients = join(f.autoregressive(), f.movingAverage());
                    converged &= f.converged();
                } else {
                    ArimaResult f = DiffuseArima.fit(y, order, options).fit();
                    ll = f.logLikelihood(); variance = f.innovationVariance();
                    coefficients = join(f.autoregressive(), f.movingAverage(), f.seasonalAutoregressive(), f.seasonalMovingAverage(),
                        f.drift() ? new double[] {f.location() / (seasonal.differences() > 0 ? seasonal.period() : 1)} : new double[0]);
                    forecast = f.forecast(8); converged &= f.converged();
                }
                if (iteration >= 0) times[iteration] = (System.nanoTime() - start) / 1e6;
                checksum += ll + variance + Arrays.stream(coefficients).sum();
            }
            Arrays.sort(times);
            double coefError = error(ref, "coef", coefficients);
            double forecastError = forecast == null ? Double.NaN : error(ref, "forecast", forecast.means());
            double seError = forecast == null ? Double.NaN : error(ref, "se", forecast.standardErrors());
            System.out.printf(Locale.ROOT, "%s,%.6f,%.9g,%.9g,%.9g,%.9g,%s,%.9g%n", name,
                times[repetitions / 2], ll - ref.get("loglik"), coefError, forecastError, seError, converged, checksum);
            if (!converged || Math.abs(ll - ref.get("loglik")) > 1e-3 || coefError > 1e-3
                    || (forecast != null && (forecastError > .02 || seError > .002)))
                throw new AssertionError("accuracy/convergence gate failed for " + name);
        }
    }
    private static double[] join(double[]... blocks) {
        return Arrays.stream(blocks).flatMapToDouble(Arrays::stream).toArray();
    }
    private static double error(Map<String, Double> ref, String prefix, double[] values) {
        double error = 0;
        for (int i = 0; i < values.length; i++) error = Math.max(error,
            Math.abs(values[i] - ref.get(values.length == 1 ? prefix : prefix + (i + 1))));
        return error;
    }
}
