/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Empirical-Bayes ComBat for normalized continuous feature-by-sample data.
 * The equations and edge-case behavior follow Bioconductor sva 3.60.0.
 */
public final class ComBat {
    private ComBat() { }

    /** ComBat settings. */
    public record Options(boolean parametric, boolean meanOnly, String referenceBatch,
                          double convergenceTolerance, int maximumIterations) {
        public Options {
            if (!(convergenceTolerance > 0.0) || !Double.isFinite(convergenceTolerance))
                throw new IllegalArgumentException("convergence tolerance must be positive and finite");
            if (maximumIterations < 1)
                throw new IllegalArgumentException("maximum iterations must be positive");
        }
        public static Options defaults() { return new Options(true, false, null, 1e-4, 10000); }
    }

    /** Adjusted matrix and fitted location/scale effects (batch by feature). */
    public record Result(double[][] adjusted, double[][] gamma, double[][] delta,
                         List<String> batches, boolean meanOnly,
                         int[] unadjustedFeatures) {
        public Result { batches = List.copyOf(batches); }
    }

    public static Result adjust(double[][] data, String[] batch, double[][] covariates,
            Options options) {
        int samples = validate(data, batch);
        if (covariates != null) ConfounderMath.validateDesign(covariates, samples, "covariates");
        LinkedHashMap<String, List<Integer>> groups = groups(batch);
        int batchCount = groups.size(), features = data.length;
        if (batchCount < 2) throw new IllegalArgumentException("ComBat needs at least two batches");
        List<String> levels = List.copyOf(groups.keySet());
        int reference = options.referenceBatch() == null ? -1 : levels.indexOf(options.referenceBatch());
        if (options.referenceBatch() != null && reference < 0)
            throw new IllegalArgumentException("reference batch is absent: " + options.referenceBatch());

        boolean meanOnly = options.meanOnly() || groups.values().stream().anyMatch(g -> g.size() == 1);
        boolean[] unadjusted = new boolean[features];
        for (int f = 0; f < features; f++) for (List<Integer> group : groups.values()) {
            if (group.size() > 1 && variance(data[f], group) == 0.0) unadjusted[f] = true;
        }
        List<Integer> keep = new ArrayList<>();
        List<Integer> excluded = new ArrayList<>();
        for (int f = 0; f < features; f++) (unadjusted[f] ? excluded : keep).add(f);
        if (keep.isEmpty()) return new Result(ConfounderMath.copy(data),
            new double[batchCount][features], ones(batchCount, features), levels,
            meanOnly, excluded.stream().mapToInt(Integer::intValue).toArray());

        double[][] batchDesign = new double[samples][batchCount];
        for (int j = 0; j < samples; j++) batchDesign[j][levels.indexOf(batch[j])] = 1.0;
        if (reference >= 0) for (int j = 0; j < samples; j++) batchDesign[j][reference] = 1.0;
        double[][] design = combineAndDropIntercept(batchDesign, covariates, reference);
        if (ConfounderMath.rank(design, 1e-12) < design[0].length)
            throw confounded(design, batchCount);

        int retained = keep.size();
        if (!options.parametric() && retained < 2)
            throw new IllegalArgumentException(
                "nonparametric ComBat requires at least two adjustable features");
        double[][] values = new double[retained][];
        for (int i = 0; i < retained; i++) values[i] = data[keep.get(i)].clone();
        double[][] beta = new double[retained][design[0].length];
        for (int f = 0; f < retained; f++) beta[f] = leastSquares(design, values[f]);
        int[] batchSizes = groups.values().stream().mapToInt(List::size).toArray();

        double[] grandMean = new double[retained];
        for (int f = 0; f < retained; f++) {
            if (reference >= 0) grandMean[f] = beta[f][reference];
            else for (int b = 0; b < batchCount; b++)
                grandMean[f] += ((double) batchSizes[b] / samples) * beta[f][b];
        }
        double[][] standMean = new double[retained][samples];
        for (int f = 0; f < retained; f++) for (int j = 0; j < samples; j++) {
            double value = grandMean[f];
            for (int c = batchCount; c < design[0].length; c++) value += design[j][c] * beta[f][c];
            standMean[f][j] = value;
        }
        double[] pooledVariance = new double[retained];
        for (int f = 0; f < retained; f++) {
            List<Integer> selected = reference >= 0 ? groups.get(levels.get(reference)) : null;
            if (selected == null) {
                int count = 0;
                for (int j = 0; j < samples; j++) if (Double.isFinite(values[f][j])) {
                    double fitted = dot(design[j], beta[f]);
                    pooledVariance[f] += square(values[f][j] - fitted); count++;
                }
                pooledVariance[f] /= count;
            } else {
                int count = 0;
                for (int j : selected) if (Double.isFinite(values[f][j])) {
                    double fitted = dot(design[j], beta[f]);
                    pooledVariance[f] += square(values[f][j] - fitted); count++;
                }
                pooledVariance[f] /= count;
            }
            if (!(pooledVariance[f] > 0.0) || !Double.isFinite(pooledVariance[f]))
                throw new IllegalArgumentException("ComBat pooled variance is not positive for feature " + keep.get(f));
        }
        double[][] standardized = new double[retained][samples];
        for (int f = 0; f < retained; f++) for (int j = 0; j < samples; j++)
            standardized[f][j] = Double.isNaN(values[f][j]) ? Double.NaN
                : (values[f][j] - standMean[f][j]) / Math.sqrt(pooledVariance[f]);

        double[][] gammaHatFeatureBatch = new double[retained][batchCount];
        for (int f = 0; f < retained; f++) gammaHatFeatureBatch[f] = leastSquares(batchDesign, standardized[f]);
        double[][] gammaHat = transpose(gammaHatFeatureBatch);
        double[][] deltaHat = new double[batchCount][retained];
        int bIndex = 0;
        for (List<Integer> group : groups.values()) {
            for (int f = 0; f < retained; f++) deltaHat[bIndex][f] = meanOnly ? 1.0 : variance(standardized[f], group);
            bIndex++;
        }
        double[] gammaBar = new double[batchCount], t2 = new double[batchCount];
        double[] aPrior = new double[batchCount], bPrior = new double[batchCount];
        for (int b = 0; b < batchCount; b++) {
            gammaBar[b] = mean(gammaHat[b]);
            t2[b] = sampleVariance(gammaHat[b]);
            double dm = mean(deltaHat[b]), dv = sampleVariance(deltaHat[b]);
            aPrior[b] = (2.0 * dv + dm * dm) / dv;
            bPrior[b] = (dm * dv + dm * dm * dm) / dv;
        }

        double[][] gammaStar = new double[batchCount][retained];
        double[][] deltaStar = new double[batchCount][retained];
        bIndex = 0;
        for (List<Integer> group : groups.values()) {
            if (options.parametric()) {
                if (meanOnly) {
                    for (int f = 0; f < retained; f++) {
                        gammaStar[bIndex][f] = postMean(gammaHat[bIndex][f], gammaBar[bIndex],
                            1, 1.0, t2[bIndex]);
                        deltaStar[bIndex][f] = 1.0;
                    }
                } else {
                    iterateParametric(standardized, group, gammaHat[bIndex], deltaHat[bIndex],
                        gammaBar[bIndex], t2[bIndex], aPrior[bIndex], bPrior[bIndex],
                        options, gammaStar[bIndex], deltaStar[bIndex]);
                }
            } else {
                nonparametric(standardized, group, gammaHat[bIndex], deltaHat[bIndex],
                    meanOnly, gammaStar[bIndex], deltaStar[bIndex]);
            }
            bIndex++;
        }
        if (reference >= 0) {
            Arrays.fill(gammaStar[reference], 0.0);
            Arrays.fill(deltaStar[reference], 1.0);
        }

        double[][] adjustedKept = ConfounderMath.copy(standardized);
        bIndex = 0;
        for (List<Integer> group : groups.values()) {
            for (int f = 0; f < retained; f++) for (int j : group) if (Double.isFinite(adjustedKept[f][j]))
                adjustedKept[f][j] = (adjustedKept[f][j] - gammaStar[bIndex][f])
                    / Math.sqrt(deltaStar[bIndex][f]);
            bIndex++;
        }
        for (int f = 0; f < retained; f++) for (int j = 0; j < samples; j++)
            if (Double.isFinite(adjustedKept[f][j])) adjustedKept[f][j] = adjustedKept[f][j]
                * Math.sqrt(pooledVariance[f]) + standMean[f][j];
        if (reference >= 0) for (int j : groups.get(levels.get(reference)))
            for (int f = 0; f < retained; f++) adjustedKept[f][j] = values[f][j];

        double[][] adjusted = ConfounderMath.copy(data);
        double[][] gamma = new double[batchCount][features], delta = ones(batchCount, features);
        for (int f = 0; f < retained; f++) {
            int original = keep.get(f); adjusted[original] = adjustedKept[f];
            for (int b = 0; b < batchCount; b++) {
                gamma[b][original] = gammaStar[b][f]; delta[b][original] = deltaStar[b][f];
            }
        }
        return new Result(adjusted, gamma, delta, levels, meanOnly,
            excluded.stream().mapToInt(Integer::intValue).toArray());
    }

    private static void iterateParametric(double[][] s, List<Integer> group,
            double[] gammaHat, double[] deltaHat, double gammaBar, double t2,
            double a, double b, Options options, double[] gammaOut, double[] deltaOut) {
        int features = s.length;
        double[] oldGamma = gammaHat.clone(), oldDelta = deltaHat.clone();
        for (int iteration = 0; iteration < options.maximumIterations(); iteration++) {
            double change = 0.0;
            for (int f = 0; f < features; f++) {
                int count = 0;
                for (int j : group) if (Double.isFinite(s[f][j])) count++;
                gammaOut[f] = postMean(gammaHat[f], gammaBar, count, oldDelta[f], t2);
                double sum = 0.0;
                for (int j : group) if (Double.isFinite(s[f][j])) sum += square(s[f][j] - gammaOut[f]);
                deltaOut[f] = (0.5 * sum + b) / (count / 2.0 + a - 1.0);
                change = Math.max(change, relative(gammaOut[f], oldGamma[f]));
                change = Math.max(change, relative(deltaOut[f], oldDelta[f]));
            }
            if (change <= options.convergenceTolerance()) return;
            System.arraycopy(gammaOut, 0, oldGamma, 0, features);
            System.arraycopy(deltaOut, 0, oldDelta, 0, features);
        }
        throw new IllegalArgumentException("ComBat parametric adjustment did not converge");
    }

    private static void nonparametric(double[][] s, List<Integer> group,
            double[] gammaHat, double[] deltaHat, boolean meanOnly,
            double[] gammaOut, double[] deltaOut) {
        int features = s.length;
        for (int target = 0; target < features; target++) {
            double maxLog = Double.NEGATIVE_INFINITY;
            double[] logLikelihood = new double[features];
            for (int candidate = 0; candidate < features; candidate++) {
                if (candidate == target) { logLikelihood[candidate] = Double.NEGATIVE_INFINITY; continue; }
                double value = 0.0;
                for (int j : group) if (Double.isFinite(s[target][j])) {
                    double d = meanOnly ? 1.0 : deltaHat[candidate];
                    value += -0.5 * Math.log(2.0 * Math.PI * d)
                        - square(s[target][j] - gammaHat[candidate]) / (2.0 * d);
                }
                logLikelihood[candidate] = value;
                maxLog = Math.max(maxLog, value);
            }
            double denominator = 0.0;
            for (int candidate = 0; candidate < features; candidate++) if (candidate != target)
                denominator += Math.exp(logLikelihood[candidate] - maxLog);
            for (int candidate = 0; candidate < features; candidate++) if (candidate != target) {
                double weight = Math.exp(logLikelihood[candidate] - maxLog) / denominator;
                gammaOut[target] += gammaHat[candidate] * weight;
                deltaOut[target] += (meanOnly ? 1.0 : deltaHat[candidate]) * weight;
            }
        }
    }

    private static double[][] combineAndDropIntercept(double[][] batch, double[][] covariates, int reference) {
        int rows = batch.length, bc = batch[0].length, cc = covariates == null ? 0 : covariates[0].length;
        boolean[] keep = new boolean[bc + cc]; int columns = 0;
        for (int c = 0; c < keep.length; c++) {
            boolean allOne = true;
            for (int r = 0; r < rows; r++) {
                double value = c < bc ? batch[r][c] : covariates[r][c - bc];
                allOne &= value == 1.0;
            }
            keep[c] = !allOne || c == reference;
            if (keep[c]) columns++;
        }
        double[][] result = new double[rows][columns];
        for (int r = 0; r < rows; r++) { int out = 0; for (int c = 0; c < keep.length; c++) if (keep[c])
            result[r][out++] = c < bc ? batch[r][c] : covariates[r][c - bc]; }
        return result;
    }

    private static IllegalArgumentException confounded(double[][] design, int batches) {
        if (design[0].length == batches + 1) return new IllegalArgumentException("the covariate is confounded with batch");
        if (design[0].length > batches + 1) return new IllegalArgumentException("at least one covariate is confounded with batch");
        return new IllegalArgumentException("batch design is rank deficient");
    }

    private static double[] leastSquares(double[][] design, double[] response) {
        int p = design[0].length;
        double[][] gram = new double[p][p]; double[] rhs = new double[p]; int count = 0;
        for (int i = 0; i < design.length; i++) if (Double.isFinite(response[i])) {
            count++;
            for (int a = 0; a < p; a++) { rhs[a] += design[i][a] * response[i];
                for (int b = 0; b < p; b++) gram[a][b] += design[i][a] * design[i][b]; }
        }
        if (count < p) throw new IllegalArgumentException("too few observed values for ComBat design");
        return solve(gram, rhs);
    }

    private static double[] solve(double[][] matrix, double[] rhs) {
        int n = rhs.length; double[][] a = new double[n][n + 1];
        for (int i = 0; i < n; i++) { System.arraycopy(matrix[i], 0, a[i], 0, n); a[i][n] = rhs[i]; }
        for (int c = 0; c < n; c++) {
            int pivot = c; for (int r = c + 1; r < n; r++) if (Math.abs(a[r][c]) > Math.abs(a[pivot][c])) pivot = r;
            if (!(Math.abs(a[pivot][c]) > 1e-12)) throw new IllegalArgumentException("ComBat design is rank deficient");
            double[] swap = a[c]; a[c] = a[pivot]; a[pivot] = swap;
            double d = a[c][c]; for (int j = c; j <= n; j++) a[c][j] /= d;
            for (int r = 0; r < n; r++) if (r != c) { double m = a[r][c]; for (int j = c; j <= n; j++) a[r][j] -= m * a[c][j]; }
        }
        double[] result = new double[n]; for (int i = 0; i < n; i++) result[i] = a[i][n]; return result;
    }

    private static int validate(double[][] data, String[] batch) {
        if (data == null || data.length == 0 || data[0] == null || data[0].length == 0)
            throw new IllegalArgumentException("data must contain features and samples");
        int n = data[0].length;
        if (batch == null || batch.length != n) throw new IllegalArgumentException("batch length must equal samples");
        for (String value : batch) if (value == null || value.isBlank()) throw new IllegalArgumentException("batch values must be nonblank");
        for (double[] row : data) { if (row == null || row.length != n) throw new IllegalArgumentException("data must be rectangular");
            for (double value : row) if (Double.isInfinite(value)) throw new IllegalArgumentException("data cannot contain infinities"); }
        return n;
    }

    private static LinkedHashMap<String, List<Integer>> groups(String[] batch) {
        LinkedHashMap<String, List<Integer>> result = new LinkedHashMap<>();
        for (int i = 0; i < batch.length; i++) result.computeIfAbsent(batch[i], ignored -> new ArrayList<>()).add(i);
        return result;
    }
    private static double variance(double[] row, List<Integer> indices) { double mean = 0; int n = 0; for (int i : indices) if (Double.isFinite(row[i])) { mean += row[i]; n++; } mean /= n; double ss = 0; for (int i : indices) if (Double.isFinite(row[i])) ss += square(row[i] - mean); return ss / (n - 1.0); }
    private static double mean(double[] x) { return Arrays.stream(x).average().orElseThrow(); }
    private static double sampleVariance(double[] x) { double m = mean(x), ss = 0; for (double v : x) ss += square(v - m); return ss / (x.length - 1.0); }
    private static double postMean(double estimate, double mean, int n, double variance, double t2) { return (t2 * n * estimate + variance * mean) / (t2 * n + variance); }
    private static double relative(double value, double old) { return Math.abs(value - old) / Math.max(Math.abs(old), 1e-12); }
    private static double dot(double[] a, double[] b) { double result = 0; for (int i = 0; i < a.length; i++) result += a[i] * b[i]; return result; }
    private static double square(double x) { return x * x; }
    private static double[][] transpose(double[][] x) { double[][] result = new double[x[0].length][x.length]; for (int i = 0; i < x.length; i++) for (int j = 0; j < x[0].length; j++) result[j][i] = x[i][j]; return result; }
    private static double[][] ones(int rows, int columns) { double[][] result = new double[rows][columns]; for (double[] row : result) Arrays.fill(row, 1.0); return result; }
}
