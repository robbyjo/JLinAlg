/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.benchmark;

import java.util.Arrays;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gee.Gee;
import org.jlinalg.gee.GeeCorrelation;
import org.jlinalg.gee.GeeOptions;
import org.jlinalg.gee.GeeResult;
import org.jlinalg.glm.GlmFamilies;

/** Reproducible cluster-scaling smoke benchmark for the GEE kernel. */
public final class GeeBenchmark {
    private GeeBenchmark() { }

    /** Runs independence, exchangeable, AR(1), and unstructured fits. */
    public static void main(String[] arguments) {
        int clusters = arguments.length == 0 ? 10_000
            : Integer.parseInt(arguments[0]);
        int size = arguments.length < 2 ? 5 : Integer.parseInt(arguments[1]);
        int rows = clusters * size;
        double[] response = new double[rows];
        double[][] design = new double[rows][3];
        int[] id = new int[rows];
        int[] repeated = new int[rows];
        long state = 0x9e3779b97f4a7c15L;
        for (int cluster = 0; cluster < clusters; cluster++) {
            double clusterEffect = ((cluster % 17) - 8) * 0.025;
            for (int visit = 0; visit < size; visit++) {
                int row = cluster * size + visit;
                state ^= state << 13;
                state ^= state >>> 7;
                state ^= state << 17;
                double noise = ((state >>> 11) * 0x1.0p-53 - 0.5) * 0.4;
                design[row][0] = 1.0;
                design[row][1] = visit;
                design[row][2] = cluster % 2;
                response[row] = 1.0 + 0.4 * visit + 0.25 * design[row][2]
                    + clusterEffect + noise;
                id[row] = cluster;
                repeated[row] = visit;
            }
        }
        for (GeeCorrelation correlation : new GeeCorrelation[] {
                GeeCorrelation.INDEPENDENCE, GeeCorrelation.EXCHANGEABLE,
                GeeCorrelation.AR1, GeeCorrelation.UNSTRUCTURED}) {
            long start = System.nanoTime();
            GeeResult result = Gee.fit(response, design, id, repeated,
                GlmFamilies.gaussian(), null, null,
                GeeOptions.builder().correlation(correlation).build(),
                BackendPolicy.CPU);
            double seconds = (System.nanoTime() - start) / 1e9;
            System.out.printf(java.util.Locale.ROOT,
                "%s rows=%d clusters=%d seconds=%.6f beta=%s alpha=%s%n",
                correlation, rows, clusters, seconds,
                Arrays.toString(result.coefficients()),
                Arrays.toString(result.associationParameters()));
        }
    }
}
