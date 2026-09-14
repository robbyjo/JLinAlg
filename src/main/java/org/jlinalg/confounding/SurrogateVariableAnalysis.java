/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Iteratively reweighted SVA following Bioconductor {@code sva::irwsva.build}. */
public final class SurrogateVariableAnalysis {
    private SurrogateVariableAnalysis() { }

    /** Standard SVA settings; Bioconductor uses five iterations. */
    public record Options(int factors, int iterations) {
        public Options {
            if (factors < 1) throw new IllegalArgumentException("factors must be positive");
            if (iterations < 1) throw new IllegalArgumentException("iterations must be positive");
        }
        public static Options defaults(int factors) { return new Options(factors, 5); }
    }

    public static SvaResult fit(double[][] data, double[][] fullDesign,
            double[][] nullDesign, Options options) {
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            return fit(data, fullDesign, nullDesign, options, context.backend());
        }
    }

    public static SvaResult fit(double[][] data, double[][] fullDesign,
            double[][] nullDesign, Options options, ComputeBackend backend) {
        int n = ConfounderMath.samples(data);
        ConfounderMath.validateDesign(fullDesign, n, "full design");
        if (nullDesign == null) nullDesign = firstColumn(fullDesign);
        ConfounderMath.validateDesign(nullDesign, n, "null design");
        if (options.factors() + fullDesign[0].length >= n)
            throw new IllegalArgumentException("SVA factors leave no residual degrees of freedom");

        double[][] residual = ConfounderMath.residualize(data, fullDesign);
        double[][] sv = ConfounderMath.sampleSvd(residual, options.factors(), backend).factors();
        double[] pprobB = new double[data.length];
        double[] pprobGam = new double[data.length];
        double[] weights = new double[data.length];
        double[][] weighted = null;
        for (int iteration = 0; iteration < options.iterations(); iteration++) {
            double[][] fullWithSv = ConfounderMath.append(fullDesign, sv);
            double[][] nullWithSv = ConfounderMath.append(nullDesign, sv);
            double[] biologicalP = ConfounderMath.fPValues(data, fullWithSv, nullWithSv);
            double[] biologicalLfdr = ConfounderMath.edgeLocalFdr(biologicalP);
            for (int i = 0; i < data.length; i++) pprobB[i] = 1.0 - biologicalLfdr[i];

            double[] heterogeneityP = ConfounderMath.fPValues(data, nullWithSv, nullDesign);
            double[] heterogeneityLfdr = ConfounderMath.edgeLocalFdr(heterogeneityP);
            weighted = new double[data.length][n];
            for (int i = 0; i < data.length; i++) {
                pprobGam[i] = 1.0 - heterogeneityLfdr[i];
                weights[i] = pprobGam[i] * (1.0 - pprobB[i]);
                double mean = Arrays.stream(data[i]).average().orElseThrow();
                for (int j = 0; j < n; j++) weighted[i][j] = (data[i][j] - mean) * weights[i];
            }
            sv = ConfounderMath.sampleSvd(weighted, options.factors(), backend).factors();
        }
        ConfounderMath.Svd finalSvd = ConfounderMath.sampleSvd(weighted, options.factors(), backend);
        return new SvaResult(finalSvd.factors(), pprobGam, pprobB, weights,
            finalSvd.varianceExplained(), ConfounderMath.removeFactors(data, finalSvd.factors()),
            options.factors(), options.iterations(), true, Double.NaN, List.of());
    }

    /** Buja-Eyuboglu estimate used by {@code sva::num.sv(method="be")}. */
    public static int estimateFactorCount(double[][] data, double[][] design,
            int permutations, long seed, ComputeBackend backend) {
        int n = ConfounderMath.samples(data), features = data.length;
        if (permutations < 1) throw new IllegalArgumentException("permutations must be positive");
        int rank = ConfounderMath.rank(design, 1e-12);
        if (rank != design[0].length) throw new IllegalArgumentException("SVA design is rank deficient");
        double[][] residual = ConfounderMath.residualize(data, design);
        ConfounderMath.Svd observed = ConfounderMath.sampleSvd(residual,
            Math.min(features, n), backend);
        int residualComponents = Math.min(features, n) - rank;
        if (residualComponents <= 0) return 0;
        double[] observedFraction = fractions(observed.singularValues(), residualComponents);
        int[] exceed = new int[residualComponents];
        Random random = new Random(seed);
        for (int permutation = 0; permutation < permutations; permutation++) {
            double[][] shuffled = ConfounderMath.copy(residual);
            for (double[] row : shuffled) for (int j = row.length - 1; j > 0; j--) {
                int replacement = random.nextInt(j + 1);
                double swap = row[j]; row[j] = row[replacement]; row[replacement] = swap;
            }
            shuffled = ConfounderMath.residualize(shuffled, design);
            ConfounderMath.Svd nullSvd = ConfounderMath.sampleSvd(shuffled,
                Math.min(features, n), backend);
            double[] nullFraction = fractions(nullSvd.singularValues(), residualComponents);
            for (int k = 0; k < residualComponents; k++)
                if (nullFraction[k] >= observedFraction[k]) exceed[k]++;
        }
        double previous = 0.0; int selected = 0;
        for (int k = 0; k < residualComponents; k++) {
            double p = Math.max(previous, (double) exceed[k] / permutations);
            if (p <= 0.10) selected++;
            previous = p;
        }
        return selected;
    }

    private static double[][] firstColumn(double[][] design) {
        double[][] result = new double[design.length][1];
        for (int i = 0; i < design.length; i++) result[i][0] = design[i][0];
        return result;
    }

    private static double[] fractions(double[] singular, int count) {
        double total = 0.0;
        for (int i = 0; i < count; i++) total += singular[i] * singular[i];
        double[] result = new double[count];
        for (int i = 0; i < count; i++) result[i] = singular[i] * singular[i] / total;
        return result;
    }
}
