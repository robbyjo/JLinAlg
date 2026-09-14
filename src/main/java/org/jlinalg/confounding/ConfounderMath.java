/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.Arrays;
import jdistlib.Normal;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import jdistlib.accelerator.SingularValueDecomposition;

/** Package-private dense operations shared by confounder estimators. */
final class ConfounderMath {
    private ConfounderMath() { }

    record Svd(double[][] factors, double[][] loadings,
               double[] singularValues, double[] varianceExplained) { }

    static int samples(double[][] matrix) {
        if (matrix == null || matrix.length == 0 || matrix[0] == null
                || matrix[0].length == 0) {
            throw new IllegalArgumentException("data must contain features and samples");
        }
        int samples = matrix[0].length;
        for (double[] row : matrix) {
            if (row == null || row.length != samples)
                throw new IllegalArgumentException("data must be rectangular");
            for (double value : row) if (!Double.isFinite(value))
                throw new IllegalArgumentException("data must contain only finite values");
        }
        return samples;
    }

    static double[][] prepare(double[][] data, boolean center, boolean scale) {
        int n = samples(data);
        double[][] result = copy(data);
        for (double[] row : result) {
            double mean = center ? Arrays.stream(row).average().orElseThrow() : 0.0;
            double sumSquares = 0.0;
            for (int j = 0; j < n; j++) {
                row[j] -= mean;
                sumSquares += row[j] * row[j];
            }
            if (scale) {
                double standardDeviation = Math.sqrt(sumSquares / Math.max(1, n - 1));
                if (!(standardDeviation > 0.0))
                    throw new IllegalArgumentException("cannot scale a constant feature");
                for (int j = 0; j < n; j++) row[j] /= standardDeviation;
            }
        }
        return result;
    }

    static Svd sampleSvd(double[][] data, int requested, ComputeBackend backend) {
        int features = data.length, samples = samples(data);
        int maximum = Math.min(features, samples);
        if (requested > maximum)
            throw new IllegalArgumentException("factors exceed matrix rank bound");
        double[] flat = flatten(data);
        if (features > samples * 2L)
            return sampleGramSvd(flat, features, samples, requested, backend);
        SingularValueDecomposition svd = backend.dgesvd(flat, features, samples);
        int components = svd.components();
        double[] singular = Arrays.copyOf(svd.singularValues(), requested);
        double[] rightT = svd.rightSingularVectorsTransposed();
        double[] left = svd.leftSingularVectors();
        double[][] factors = new double[samples][requested];
        double[][] loadings = new double[features][requested];
        double total = 0.0;
        for (double value : svd.singularValues()) total += value * value;
        double[] fraction = new double[requested];
        for (int k = 0; k < requested; k++) {
            fraction[k] = total == 0.0 ? 0.0 : singular[k] * singular[k] / total;
            for (int j = 0; j < samples; j++)
                factors[j][k] = rightT[k * samples + j];
            for (int i = 0; i < features; i++)
                loadings[i][k] = left[i * components + k] * singular[k];
            orient(factors, loadings, k);
        }
        return new Svd(factors, loadings, singular, fraction);
    }

    /** Tall-matrix path used by the AutoSVA source: diagonalize X'X, then Xv. */
    private static Svd sampleGramSvd(double[] data, int features, int samples,
            int requested, ComputeBackend backend) {
        double[] gram = new double[samples * samples];
        backend.dgemm(MatrixTranspose.TRANSPOSE, MatrixTranspose.NONE,
            samples, samples, features, 1.0, data, data, 0.0, gram);
        SingularValueDecomposition decomposition = backend.dgesvd(gram, samples, samples);
        double[] vectors = decomposition.leftSingularVectors();
        int components = decomposition.components();
        double[][] factors = new double[samples][requested];
        double[][] loadings = new double[features][requested];
        double[] singular = new double[requested];
        double total = 0.0;
        for (double eigenvalue : decomposition.singularValues()) total += Math.max(0.0, eigenvalue);
        double[] fraction = new double[requested];
        for (int k = 0; k < requested; k++) {
            double eigenvalue = Math.max(0.0, decomposition.singularValues()[k]);
            singular[k] = Math.sqrt(eigenvalue);
            fraction[k] = total == 0.0 ? 0.0 : eigenvalue / total;
            for (int j = 0; j < samples; j++) factors[j][k] = vectors[j * components + k];
            for (int i = 0; i < features; i++) {
                double value = 0.0;
                for (int j = 0; j < samples; j++) value += data[i * samples + j] * factors[j][k];
                loadings[i][k] = value;
            }
            orient(factors, loadings, k);
        }
        return new Svd(factors, loadings, singular, fraction);
    }

    private static void orient(double[][] factors, double[][] loadings, int factor) {
        int pivot = 0;
        for (int i = 1; i < loadings.length; i++)
            if (Math.abs(loadings[i][factor]) > Math.abs(loadings[pivot][factor])) pivot = i;
        if (loadings[pivot][factor] >= 0.0) return;
        for (double[] row : factors) row[factor] = -row[factor];
        for (double[] row : loadings) row[factor] = -row[factor];
    }

    static double[][] residualize(double[][] data, double[][] design) {
        int n = samples(data);
        validateDesign(design, n, "design");
        double[][] q = orthonormalColumns(design, 1e-12);
        double[][] residual = copy(data);
        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < q[0].length; k++) {
                double coefficient = 0.0;
                for (int j = 0; j < n; j++) coefficient += data[i][j] * q[j][k];
                for (int j = 0; j < n; j++) residual[i][j] -= coefficient * q[j][k];
            }
        }
        return residual;
    }

    static double[][] removeFactors(double[][] data, double[][] factors) {
        if (factors == null || factors.length == 0 || factors[0].length == 0)
            return copy(data);
        return residualize(data, factors);
    }

    static double[] residualSums(double[][] data, double[][] design) {
        double[][] residual = residualize(data, design);
        double[] result = new double[data.length];
        for (int i = 0; i < residual.length; i++)
            for (double value : residual[i]) result[i] += value * value;
        return result;
    }

    static double[] fPValues(double[][] data, double[][] full, double[][] reduced) {
        int n = samples(data);
        int fullRank = rank(full, 1e-12), reducedRank = rank(reduced, 1e-12);
        if (fullRank != full[0].length || reducedRank != reduced[0].length)
            throw new IllegalArgumentException("SVA design matrices must have full column rank");
        if (fullRank <= reducedRank || n <= fullRank)
            throw new IllegalArgumentException("SVA models require nested positive degrees of freedom");
        double[] rss1 = residualSums(data, full), rss0 = residualSums(data, reduced);
        double[] result = new double[data.length];
        for (int i = 0; i < result.length; i++) {
            double f = ((rss0[i] - rss1[i]) / (fullRank - reducedRank))
                / (rss1[i] / (n - fullRank));
            if (!Double.isFinite(f) || f < 0.0) f = 0.0;
            result[i] = jdistlib.F.cumulative(f, fullRank - reducedRank,
                n - fullRank, false, false);
            result[i] = Math.max(0.0, Math.min(1.0, result[i]));
        }
        return result;
    }

    static double[] edgeLocalFdr(double[] probabilities) {
        int n = probabilities.length;
        if (n < 2) throw new IllegalArgumentException("local FDR needs at least two features");
        double pi0 = 0.0;
        for (double value : probabilities) if (value >= 0.8) pi0++;
        pi0 = Math.min(1.0, pi0 / n / 0.2);
        double[] p = probabilities.clone(), x = new double[n];
        for (int i = 0; i < n; i++) {
            p[i] = Math.max(1e-8, Math.min(1.0 - 1e-8, p[i]));
            x[i] = Normal.quantile(p[i], 0.0, 1.0, true, false);
        }
        double bandwidth = 1.5 * bandwidthNrd0(x);
        double[] lfdr = new double[n];
        double normalizer = n * bandwidth * Math.sqrt(2.0 * Math.PI);
        for (int i = 0; i < n; i++) {
            double density = 0.0;
            for (double value : x) {
                double z = (x[i] - value) / bandwidth;
                if (Math.abs(z) < 38.0) density += Math.exp(-0.5 * z * z);
            }
            density /= normalizer;
            double nullDensity = Math.exp(-0.5 * x[i] * x[i]) / Math.sqrt(2.0 * Math.PI);
            lfdr[i] = Math.min(1.0, pi0 * nullDensity / Math.max(density, Double.MIN_NORMAL));
        }
        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) order[i] = i;
        Arrays.sort(order, (a, b) -> Double.compare(p[a], p[b]));
        double maximum = 0.0;
        for (int index : order) {
            maximum = Math.max(maximum, lfdr[index]);
            lfdr[index] = maximum;
        }
        return lfdr;
    }

    private static double bandwidthNrd0(double[] values) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        double mean = Arrays.stream(values).average().orElseThrow(), ss = 0.0;
        for (double value : values) ss += (value - mean) * (value - mean);
        double sd = Math.sqrt(ss / Math.max(1, values.length - 1));
        double iqr = quantile(sorted, 0.75) - quantile(sorted, 0.25);
        double scale = Math.min(sd, iqr / 1.34);
        if (!(scale > 0.0)) scale = Math.max(sd, Math.abs(sorted[sorted.length - 1] - sorted[0]));
        if (!(scale > 0.0)) scale = 1.0;
        return 0.9 * scale * Math.pow(values.length, -0.2);
    }

    private static double quantile(double[] sorted, double probability) {
        double position = (sorted.length - 1) * probability;
        int lower = (int) Math.floor(position), upper = (int) Math.ceil(position);
        return sorted[lower] + (position - lower) * (sorted[upper] - sorted[lower]);
    }

    static double[][] append(double[][] left, double[][] right) {
        if (left.length != right.length) throw new IllegalArgumentException("design rows differ");
        double[][] result = new double[left.length][left[0].length + right[0].length];
        for (int i = 0; i < left.length; i++) {
            System.arraycopy(left[i], 0, result[i], 0, left[i].length);
            System.arraycopy(right[i], 0, result[i], left[i].length, right[i].length);
        }
        return result;
    }

    static int rank(double[][] matrix, double tolerance) {
        return orthonormalColumns(matrix, tolerance)[0].length;
    }

    static double[][] orthonormalColumns(double[][] matrix, double tolerance) {
        validateDesign(matrix, matrix.length, "matrix");
        int rows = matrix.length, columns = matrix[0].length, rank = 0;
        double[][] q = new double[rows][columns];
        for (int column = 0; column < columns; column++) {
            double[] value = new double[rows];
            double originalNorm = 0.0;
            for (int row = 0; row < rows; row++) {
                value[row] = matrix[row][column];
                originalNorm += value[row] * value[row];
            }
            for (int pass = 0; pass < 2; pass++) for (int k = 0; k < rank; k++) {
                double projection = 0.0;
                for (int row = 0; row < rows; row++) projection += q[row][k] * value[row];
                for (int row = 0; row < rows; row++) value[row] -= projection * q[row][k];
            }
            double norm = 0.0;
            for (double entry : value) norm += entry * entry;
            norm = Math.sqrt(norm);
            if (norm <= tolerance * Math.max(1.0, Math.sqrt(originalNorm))) continue;
            for (int row = 0; row < rows; row++) q[row][rank] = value[row] / norm;
            rank++;
        }
        double[][] compact = new double[rows][rank];
        for (int row = 0; row < rows; row++) System.arraycopy(q[row], 0, compact[row], 0, rank);
        return compact;
    }

    static void validateDesign(double[][] design, int rows, String name) {
        if (design == null || design.length != rows || design[0] == null
                || design[0].length == 0) throw new IllegalArgumentException(name + " dimensions are invalid");
        int columns = design[0].length;
        for (double[] row : design) {
            if (row == null || row.length != columns) throw new IllegalArgumentException(name + " must be rectangular");
            for (double value : row) if (!Double.isFinite(value))
                throw new IllegalArgumentException(name + " must be finite");
        }
    }

    static double[] flatten(double[][] matrix) {
        int columns = matrix[0].length;
        double[] result = new double[matrix.length * columns];
        for (int row = 0; row < matrix.length; row++)
            System.arraycopy(matrix[row], 0, result, row * columns, columns);
        return result;
    }

    static double[][] copy(double[][] matrix) {
        double[][] result = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) result[i] = matrix[i].clone();
        return result;
    }
}
