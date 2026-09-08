/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

/** Fixed-parameter Gaussian likelihood including state initialization, excluding
 * optimization, Hessians, forecasts and I/O. Paired with exact_arima_likelihood_benchmark.R. */
public final class ExactArimaLikelihoodBenchmark {
    private static volatile double checksum;
    private ExactArimaLikelihoodBenchmark() { }
    public static void main(String[] args) throws Exception {
        int batch = Integer.getInteger("jlinalg.benchmark.batch", 1000);
        int repetitions = Integer.getInteger("jlinalg.benchmark.measurements", 11);
        System.out.println("case,java_likelihood_ms,loglik_error,kappa10000_error,kappa10000000_error,checksum");
        for (String name : new String[] {"stationary", "stationary-missing", "integrated-ar", "integrated-ma",
                "integrated-missing", "seasonal", "seasonal-missing", "seasonal-only", "twice-integrated", "drift",
                "seasonal-drift", "seasonal-ar-ma", "monthly-missing"}) {
            Path root = Path.of("src/test/resources/timeseries");
            double[] y = Files.readAllLines(root.resolve(name + "-series.csv")).stream().skip(1).mapToDouble(Double::parseDouble).toArray();
            Map<String, Double> ref = new LinkedHashMap<>();
            Files.readAllLines(root.resolve(name + "-reference.csv")).stream().skip(1).forEach(line -> {
                String[] fields = line.split(","); ref.put(fields[0], Double.parseDouble(fields[1]));
            });
            ArimaOrder order = new ArimaOrder(ref.get("p").intValue(), ref.get("d").intValue(), ref.get("q").intValue());
            SeasonalArimaOrder seasonal = new SeasonalArimaOrder(ref.get("P").intValue(), ref.get("D").intValue(),
                ref.get("Q").intValue(), ref.get("period").intValue());
            double[] coefficients = ref.entrySet().stream().filter(e -> e.getKey().startsWith("coef")).mapToDouble(Map.Entry::getValue).toArray();
            int position = 0;
            double[] ar = Arrays.copyOfRange(coefficients, position, position += order.autoregressive());
            double[] ma = Arrays.copyOfRange(coefficients, position, position += order.movingAverage());
            double[] sar = Arrays.copyOfRange(coefficients, position, position += seasonal.autoregressive());
            double[] sma = Arrays.copyOfRange(coefficients, position, position += seasonal.movingAverage());
            double slope = ref.get("drift") == 1 ? coefficients[position] : 0;
            ar = ArimaMath.multiplyAr(ar, sar, seasonal.period()); ma = ArimaMath.multiplyMa(ma, sma, seasonal.period());
            double[] delta = ArimaMath.differencingPolynomial(order, seasonal);
            double[] times = new double[repetitions];
            double ll = 0;
            for (int repetition = -3; repetition < repetitions; repetition++) {
                long start = System.nanoTime();
                for (int iteration = 0; iteration < batch; iteration++) {
                    ll = -new ArimaStateSpace(ar, ma, delta).filter(y, 0, slope, false).negativeLogLikelihood();
                    checksum += ll;
                }
                if (repetition >= 0) times[repetition] = (System.nanoTime() - start) / (1e6 * batch);
            }
            if (Math.abs(ll - ref.get("loglik")) > 1e-4) throw new AssertionError(name + " likelihood accuracy");
            Arrays.sort(times);
            System.out.printf(Locale.ROOT, "%s,%.6f,%.9g,%.9g,%.9g,%.9g%n", name, times[repetitions / 2],
                ll - ref.get("loglik"), ll - ref.getOrDefault("kappa10000", ll),
                ll - ref.getOrDefault("kappa10000000", ll), checksum);
        }
    }
}
