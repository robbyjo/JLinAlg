/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Deterministic port of Roby Joehanes' AutoSVA adaptation. */
public final class AutoSva {
    private AutoSva() { }

    /** AutoSVA settings; the source adaptation uses B=20 and MRSE=0.001. */
    public record Options(int factors, int iterations, double mrseTolerance,
                          int maximumFactors, int increment) {
        public Options {
            if (factors < 0) throw new IllegalArgumentException("factors must be nonnegative");
            if (iterations < 1) throw new IllegalArgumentException("iterations must be positive");
            if (!(mrseTolerance >= 0.0) || !Double.isFinite(mrseTolerance))
                throw new IllegalArgumentException("MRSE tolerance must be finite and nonnegative");
            if (maximumFactors < 0 || increment < 1)
                throw new IllegalArgumentException("AutoSVA search bounds are invalid");
        }
        public static Options fixed(int factors) { return new Options(factors, 20, 0.001, 0, 5); }
        public static Options automatic(int maximumFactors) { return new Options(0, 20, 0.001, maximumFactors, 5); }
    }

    public static SvaResult fit(double[][] data, double[][] preserveDesign,
            double[][] nullDesign, Options options) {
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            return fit(data, preserveDesign, nullDesign, preserveDesign, options, context.backend());
        }
    }

    public static SvaResult fit(double[][] data, double[][] preserveDesign,
            double[][] nullDesign, double[][] selectionDesign, Options options,
            ComputeBackend backend) {
        if (options.factors() > 0)
            return fitFixed(data, preserveDesign, nullDesign, options.factors(), options, backend);
        int n = ConfounderMath.samples(data);
        int maximum = options.maximumFactors() > 0 ? options.maximumFactors() : n / 2;
        maximum = Math.min(maximum, Math.min(data.length, n - preserveDesign[0].length - 1));
        if (maximum < 1) throw new IllegalArgumentException("no valid AutoSVA factor count");
        List<SvaResult.FactorSelection> path = new ArrayList<>();
        SvaResult best = null;
        double bestRatio = preservedSignalRatio(data, selectionDesign);
        path.add(new SvaResult.FactorSelection(0, bestRatio, true));
        boolean plateau = false;
        int lastImprovement = 0;
        outer: for (int start = 1; start <= maximum; start += options.increment()) {
            int end = Math.min(maximum, start + options.increment() - 1);
            boolean improved = false;
            for (int k = start; k <= end; k++) {
                try {
                    SvaResult candidate = fitFixed(data, preserveDesign, nullDesign, k, options, backend);
                    double ratio = preservedSignalRatio(candidate.adjusted(), selectionDesign);
                    path.add(new SvaResult.FactorSelection(k, ratio, true));
                    if (ratio > bestRatio) {
                        bestRatio = ratio;
                        best = candidate;
                        lastImprovement = path.size() - 1;
                        improved = true;
                    }
                } catch (RuntimeException failure) {
                    path.add(new SvaResult.FactorSelection(k, 0.0, false));
                }
            }
            double last = path.get(path.size() - 1).fRatio();
            if (!improved) {
                if (plateau || last < 0.75 * bestRatio || path.size() - 1 - lastImprovement > 5)
                    break outer;
                plateau = true;
            } else plateau = false;
        }
        if (best == null) throw new IllegalArgumentException("AutoSVA found no improvement over zero factors");
        return new SvaResult(best.surrogateVariables(), best.probabilityHeterogeneity(),
            best.probabilityBiological(), best.featureWeights(), best.varianceExplained(),
            best.adjusted(), best.factorCount(), best.iterations(), best.converged(),
            best.finalWeightMrse(), path);
    }

    private static SvaResult fitFixed(double[][] data, double[][] preserve,
            double[][] nullDesign, int factors, Options options, ComputeBackend backend) {
        int n = ConfounderMath.samples(data);
        ConfounderMath.validateDesign(preserve, n, "preserve design");
        if (nullDesign == null) nullDesign = firstColumn(preserve);
        ConfounderMath.validateDesign(nullDesign, n, "null design");
        if (factors + preserve[0].length >= n)
            throw new IllegalArgumentException("AutoSVA factors leave no residual degrees of freedom");
        double[][] residual = ConfounderMath.residualize(data, preserve);
        double[][] sv = ConfounderMath.sampleSvd(residual, factors, backend).factors();
        double[] rssNull = ConfounderMath.residualSums(data, nullDesign);
        double[] means = new double[data.length];
        for (int i = 0; i < data.length; i++) means[i] = Arrays.stream(data[i]).average().orElseThrow();
        double[] old = new double[data.length];
        Arrays.fill(old, 1.0);
        double[] pprobGam = null, pprobB = null, weights = null;
        double[][] weighted = null;
        double mrse = Double.POSITIVE_INFINITY;
        int completed = 0;
        boolean converged = false;
        for (int iteration = 1; iteration <= options.iterations(); iteration++) {
            double[][] fullSv = ConfounderMath.append(preserve, sv);
            double[][] nullSv = ConfounderMath.append(nullDesign, sv);
            double[] biologicalP = ConfounderMath.fPValues(data, fullSv, nullSv);
            double[] biologicalLfdr = ConfounderMath.edgeLocalFdr(biologicalP);
            pprobB = new double[data.length];
            for (int i = 0; i < data.length; i++) pprobB[i] = 1.0 - biologicalLfdr[i];

            double[] rss1 = ConfounderMath.residualSums(data, nullSv);
            int df1 = nullSv[0].length, df0 = nullDesign[0].length;
            double[] heterogeneityP = new double[data.length];
            for (int i = 0; i < data.length; i++) {
                double f = ((rssNull[i] - rss1[i]) / (df1 - df0)) / (rss1[i] / (n - df1));
                if (!Double.isFinite(f) || f <= 0.0) f = 1e-12;
                heterogeneityP[i] = jdistlib.F.cumulative(f, df1 - df0, n - df1, false, false);
            }
            double[] heterogeneityLfdr = ConfounderMath.edgeLocalFdr(heterogeneityP);
            pprobGam = new double[data.length];
            weights = new double[data.length];
            weighted = new double[data.length][n];
            double squared = 0.0;
            for (int i = 0; i < data.length; i++) {
                pprobGam[i] = 1.0 - heterogeneityLfdr[i];
                weights[i] = pprobGam[i] * biologicalLfdr[i];
                squared += (old[i] - weights[i]) * (old[i] - weights[i]);
                for (int j = 0; j < n; j++) weighted[i][j] = (data[i][j] - means[i]) * weights[i];
            }
            ConfounderMath.Svd decomposition = ConfounderMath.sampleSvd(weighted, factors, backend);
            sv = decomposition.factors();
            mrse = Math.sqrt(squared / data.length);
            completed = iteration;
            if (mrse <= options.mrseTolerance()) { converged = true; break; }
            old = weights.clone();
        }
        ConfounderMath.Svd finalSvd = ConfounderMath.sampleSvd(weighted, factors, backend);
        return new SvaResult(finalSvd.factors(), pprobGam, pprobB, weights,
            finalSvd.varianceExplained(), ConfounderMath.removeFactors(data, finalSvd.factors()),
            factors, completed, converged, mrse, List.of());
    }

    private static double preservedSignalRatio(double[][] data, double[][] design) {
        int n = ConfounderMath.samples(data), rank = ConfounderMath.rank(design, 1e-12);
        if (rank < 2 || n <= rank) throw new IllegalArgumentException(
            "AutoSVA selection design needs an intercept and preserved terms");
        double[][] centered = ConfounderMath.prepare(data, true, false);
        double[][] q = ConfounderMath.orthonormalColumns(design, 1e-12);
        double model = 0.0, residual = 0.0;
        for (double[] row : centered) {
            double fitted = 0.0;
            for (int k = 0; k < q[0].length; k++) {
                double coefficient = 0.0;
                for (int j = 0; j < n; j++) coefficient += row[j] * q[j][k];
                if (k > 0) fitted += coefficient * coefficient;
            }
            double total = 0.0;
            for (double value : row) total += value * value;
            model += fitted;
            residual += Math.max(0.0, total - fitted);
        }
        return (model / (rank - 1.0)) / (residual / (n - rank));
    }

    private static double[][] firstColumn(double[][] design) {
        double[][] result = new double[design.length][1];
        for (int i = 0; i < design.length; i++) result[i][0] = design[i][0];
        return result;
    }
}
