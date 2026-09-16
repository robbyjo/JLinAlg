/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.benchmark;

import java.util.Arrays;
import java.util.Random;
import java.util.function.DoubleSupplier;
import org.jlinalg.settest.*;

/** Reproducible fused-versus-composed ACAT-O score-summary benchmark. */
public final class AcatBenchmark {
    private static volatile double consumed;

    private AcatBenchmark() { }

    public static void main(String[] args) {
        int variants = Integer.getInteger("jlinalg.benchmark.acat.variants", 80);
        int batch = Integer.getInteger("jlinalg.benchmark.acat.batch", 8);
        int warmups = Integer.getInteger("jlinalg.benchmark.acat.warmups", 3);
        int measurements = Integer.getInteger("jlinalg.benchmark.acat.measurements", 7);
        if (variants < 2 || batch < 1 || warmups < 0 || measurements < 1)
            throw new IllegalArgumentException("invalid ACAT benchmark controls");
        Input input = input(variants);

        double[] fused = measure(warmups, measurements, batch,
            () -> AcatTests.acatO("benchmark", input.state(), input.mafs(),
                input.macs()).pValue());
        double[] composed = measure(warmups, measurements, batch,
            () -> composed(input));
        double fusedP = AcatTests.acatO("benchmark", input.state(),
            input.mafs(), input.macs()).pValue();
        double composedP = composed(input);
        System.out.println("variants,batch,warmups,measurements,fused_median_ms,"
            + "composed_median_ms,speedup,absolute_p_difference,checksum");
        System.out.printf(java.util.Locale.ROOT,
            "%d,%d,%d,%d,%.6f,%.6f,%.3f,%.9g,%.17g%n",
            variants, batch, warmups, measurements, fused[measurements / 2],
            composed[measurements / 2],
            composed[measurements / 2] / fused[measurements / 2],
            Math.abs(fusedP - composedP), consumed);
    }

    private static double composed(Input input) {
        double[] w125 = weights(input.mafs(), 25);
        double[] w11 = weights(input.mafs(), 1);
        double[] p = {
            SummarySetTests.skat("benchmark", input.state(), w125).pValue(),
            SummarySetTests.skat("benchmark", input.state(), w11).pValue(),
            SummarySetTests.burden("benchmark", input.state(), w125).pValue(),
            SummarySetTests.burden("benchmark", input.state(), w11).pValue(),
            SummarySetTests.acatV("benchmark", input.state(), input.mafs(),
                input.macs(), 1, 25, AcatTests.DEFAULT_MAC_THRESHOLD).pValue(),
            SummarySetTests.acatV("benchmark", input.state(), input.mafs(),
                input.macs(), 1, 1, AcatTests.DEFAULT_MAC_THRESHOLD).pValue()
        };
        return Acat.combine(p).pValue();
    }

    private static double[] weights(double[] mafs, double shape2) {
        double[] result = new double[mafs.length];
        for (int index = 0; index < result.length; index++)
            result[index] = VariantWeights.betaBurden(
                mafs[index], 1, shape2);
        return result;
    }

    private static double[] measure(
            int warmups, int measurements, int batch,
            DoubleSupplier operation) {
        for (int iteration = 0; iteration < warmups; iteration++)
            for (int item = 0; item < batch; item++) consumed += operation.getAsDouble();
        double[] elapsed = new double[measurements];
        for (int iteration = 0; iteration < measurements; iteration++) {
            long start = System.nanoTime();
            for (int item = 0; item < batch; item++) consumed += operation.getAsDouble();
            elapsed[iteration] = (System.nanoTime() - start) / 1e6 / batch;
        }
        Arrays.sort(elapsed);
        return elapsed;
    }

    private static Input input(int variants) {
        Random random = new Random(20260916L);
        double[] scores = new double[variants];
        double[] information = new double[variants * variants];
        double[] mafs = new double[variants], macs = new double[variants];
        for (int row = 0; row < variants; row++) {
            scores[row] = random.nextGaussian();
            mafs[row] = 0.00025 + 0.04975 * (row + 1.0) / variants;
            macs[row] = Math.max(1, Math.rint(20_000 * mafs[row]));
            for (int column = 0; column < variants; column++)
                information[row * variants + column] =
                    Math.pow(0.15, Math.abs(row - column));
        }
        return new Input(new SetTestScoreState(scores, information, variants),
            mafs, macs);
    }

    private record Input(
        SetTestScoreState state, double[] mafs, double[] macs) { }
}
