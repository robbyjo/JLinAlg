/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.susie.Susie;
import org.jlinalg.susie.SusieOptions;
import org.jlinalg.susie.SusieResult;

/** End-to-end benchmark using susieR's N3finemapping vignette data. */
public final class SusieBenchmark {
    private static final byte[] MAGIC = "JLSUSIE1".getBytes(StandardCharsets.US_ASCII);
    private static final String RESOURCE = "/susie/N3finemapping.bin.gz";
    private static volatile double consumedChecksum;

    private SusieBenchmark() { }

    public static void main(String[] arguments) {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        int warmups = integerProperty("jlinalg.benchmark.susie.warmups", 3);
        int measurements = integerProperty(
            "jlinalg.benchmark.susie.measurements", 7);
        BackendPolicy policy = BackendPolicy.valueOf(System.getProperty(
            "jlinalg.benchmark.susie.backend", "CPU").toUpperCase(Locale.ROOT));
        Data data = readData();
        List<String> names = new ArrayList<>(data.variables());
        for (int index = 0; index < data.variables(); index++) {
            names.add("variable" + (index + 1));
        }
        SusieOptions options = new SusieOptions(
            10, 200, 1e-6, 0.2, true, 0.95, 0.5);

        for (int iteration = 0; iteration < warmups; iteration++) {
            consume(Susie.fit(data.response(), data.design(), names,
                options, policy));
        }
        double[] seconds = new double[measurements];
        SusieResult result = null;
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime();
            result = Susie.fit(data.response(), data.design(), names,
                options, policy);
            seconds[iteration] = (System.nanoTime() - started) / 1e9;
            consume(result);
        }
        Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        System.out.println("runtime,dataset,observations,variables,effects,warmups,"
            + "measurements,median_seconds,iterations,converged,residual_variance,"
            + "credible_sets,backend,checksum");
        System.out.printf(Locale.ROOT,
            "Java,N3finemapping,%d,%d,%d,%d,%d,%.9f,%d,%s,%.12g,%d,%s,%.12g%n",
            data.observations(), data.variables(), result.effects(), warmups,
            measurements, median, result.iterations(), result.converged(),
            result.residualVariance(), result.credibleSets().size(),
            csv(result.backend().selectedBackend()), checksum(result));
    }

    private static void consume(SusieResult result) {
        consumedChecksum += checksum(result);
    }

    private static double checksum(SusieResult result) {
        double value = result.residualVariance() + result.intercept();
        double[] pip = result.pip();
        double[] mean = result.posteriorMean();
        for (int index = 0; index < pip.length; index++) {
            value += pip[index] * (index + 1.0) + mean[index];
        }
        return value;
    }

    private static Data readData() {
        InputStream resource = SusieBenchmark.class.getResourceAsStream(RESOURCE);
        if (resource == null) {
            throw new IllegalStateException("missing benchmark resource " + RESOURCE);
        }
        try (DataInputStream input = new DataInputStream(new BufferedInputStream(
                new GZIPInputStream(resource)))) {
            byte[] magic = input.readNBytes(MAGIC.length);
            if (!Arrays.equals(MAGIC, magic)) {
                throw new IllegalStateException("invalid SuSiE benchmark header");
            }
            int observations = input.readInt();
            int variables = input.readInt();
            double[][] design = new double[observations][variables];
            for (double[] row : design) {
                for (int column = 0; column < variables; column++) {
                    row[column] = input.readDouble();
                }
            }
            double[] response = new double[observations];
            for (int row = 0; row < observations; row++) {
                response[row] = input.readDouble();
            }
            double[] trueCoefficient = new double[variables];
            for (int column = 0; column < variables; column++) {
                trueCoefficient[column] = input.readDouble();
            }
            if (input.read() != -1) {
                throw new IllegalStateException("trailing SuSiE benchmark data");
            }
            return new Data(design, response, trueCoefficient);
        } catch (IOException exception) {
            throw new IllegalStateException("cannot read " + RESOURCE, exception);
        }
    }

    private static int integerProperty(String name, int fallback) {
        int value = Integer.getInteger(name, fallback);
        if (value < 1) throw new IllegalArgumentException(name + " must be positive");
        return value;
    }

    private static String csv(String value) {
        return '"' + value.replace("\"", "\"\"") + '"';
    }

    private record Data(
            double[][] design, double[] response, double[] trueCoefficient) {
        int observations() { return design.length; }
        int variables() { return design[0].length; }
    }
}
