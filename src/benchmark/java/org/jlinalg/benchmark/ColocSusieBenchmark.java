/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import org.jlinalg.coloc.ColocOptions;
import org.jlinalg.coloc.ColocSusie;
import org.jlinalg.coloc.ColocSusieInput;
import org.jlinalg.coloc.ColocSusieResult;

/** Deterministic throughput benchmark for the coloc.susie combination step. */
public final class ColocSusieBenchmark {
    private ColocSusieBenchmark() { }

    public static void main(String[] arguments) {
        int variants = Integer.getInteger(
            "jlinalg.benchmark.coloc.variants", 10_000);
        int signals = Integer.getInteger(
            "jlinalg.benchmark.coloc.signals", 10);
        int warmups = Integer.getInteger(
            "jlinalg.benchmark.coloc.warmups", 3);
        int measurements = Integer.getInteger(
            "jlinalg.benchmark.coloc.measurements", 7);
        if (variants < 1 || signals < 1 || warmups < 0 || measurements < 1) {
            throw new IllegalArgumentException("invalid benchmark dimensions");
        }

        List<String> names = new ArrayList<>(variants);
        for (int variant = 0; variant < variants; variant++) {
            names.add("rs" + (variant + 1));
        }
        Random random = new Random(34_903L);
        double[][] first = logBayesFactors(random, signals, variants);
        double[][] second = logBayesFactors(random, signals, variants);
        ColocSusieInput input1 = new ColocSusieInput(names, first);
        ColocSusieInput input2 = new ColocSusieInput(names, second);
        ColocOptions options = new ColocOptions(
            1e-4, 1e-4, 5e-6, 0.5, false, null, null);

        double checksum = 0.0;
        for (int iteration = 0; iteration < warmups; iteration++) {
            checksum += run(input1, input2, options);
        }
        double[] seconds = new double[measurements];
        for (int iteration = 0; iteration < measurements; iteration++) {
            long start = System.nanoTime();
            checksum += run(input1, input2, options);
            seconds[iteration] = (System.nanoTime() - start) / 1e9;
        }
        Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        System.out.println("benchmark,variants,signals_per_trait,signal_pairs,"
            + "warmups,measurements,median_seconds,variant_pair_updates_per_second,checksum");
        System.out.printf(java.util.Locale.ROOT,
            "coloc_susie,%d,%d,%d,%d,%d,%.9f,%.3f,%.9f%n",
            variants, signals, signals * signals, warmups, measurements,
            median, variants * (double) signals * signals / median, checksum);
    }

    private static double run(
            ColocSusieInput first, ColocSusieInput second,
            ColocOptions options) {
        ColocSusieResult result = ColocSusie.analyze(first, second, options);
        double[] posterior = result.sharedVariantPosterior();
        return result.signalPairs().get(0).posteriorH4()
            + posterior[posterior.length - 1];
    }

    private static double[][] logBayesFactors(
            Random random, int signals, int variants) {
        double[][] result = new double[signals][variants];
        for (int signal = 0; signal < signals; signal++) {
            int lead = (signal * 997 + 101) % variants;
            for (int variant = 0; variant < variants; variant++) {
                result[signal][variant] = random.nextGaussian() * 0.5;
            }
            result[signal][lead] += 12.0;
        }
        return result;
    }
}
