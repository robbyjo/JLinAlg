/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.Normal;
import jdistlib.T;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.SymmetricEigenDecomposition;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;

/** Correlated-effect GLS with cluster sandwiches and explicit finite-sample corrections. */
public final class MetaClusterRobust {
    private MetaClusterRobust() { }

    public static MetaClusterRobustResult fit(List<MetaStudy> studies, String[] clusters,
            double[][] moderators, List<String> names, MetaAnalysisOptions options, BackendPolicy policy) {
        return fit(studies, clusters, moderators, names, options, MetaRobustCorrection.CR0, policy);
    }
    public static MetaClusterRobustResult fit(List<MetaStudy> studies, String[] clusters,
            double[][] moderators, List<String> names, MetaAnalysisOptions options,
            MetaRobustCorrection correction, BackendPolicy policy) {
        MetaMath.Data data = MetaMath.data(studies); int n = studies.size();
        double[][] sampling = new double[n][n];
        for (int i = 0; i < n; i++) sampling[i][i] = data.variances()[i];
        return fit(studies, clusters, moderators, names, sampling, options, correction, policy);
    }

    /** Known correlated sampling covariance plus an estimated independent tau-squared. */
    public static MetaClusterRobustResult fit(List<MetaStudy> studies, String[] clusters,
            double[][] moderators, List<String> names, double[][] samplingCovariance,
            MetaAnalysisOptions options, MetaRobustCorrection correction, BackendPolicy policy) {
        double[] y = MetaMath.data(studies).effects();
        double[] x = MetaGls.design(moderators, names, y.length, true), v = MetaGls.matrix(samplingCovariance, y.length);
        validateOptions(options, correction);
        try (BackendContext context = BackendContext.select(policy)) {
            double tau = MetaGls.tau(y, x, names.size() + 1, v, options, context.backend());
            return calculate(y, x, names, clusters, MetaGls.addDiagonal(v, y.length, tau), tau, options, correction, context);
        }
    }

    /** Apply a cluster sandwich to the jointly estimated hierarchical covariance.
     * The supplied design must match the hierarchical fit, including its intercept. */
    public static MetaClusterRobustResult fit(List<MetaStudy> studies, String[] clusters,
            double[][] moderators, List<String> names, MetaHierarchicalResult hierarchical,
            MetaAnalysisOptions options, MetaRobustCorrection correction, BackendPolicy policy) {
        if (hierarchical == null || !hierarchical.converged()) throw new IllegalArgumentException("a converged hierarchical fit is required");
        if (!hierarchical.model().coefficientNames().equals(MetaMultilevelRegression.names(names, true)))
            throw new IllegalArgumentException("hierarchical fixed design must include the same intercept and moderators");
        double[] y = MetaMath.data(studies).effects(), x = MetaGls.design(moderators, names, studies.size(), true);
        validateOptions(options, correction);
        try (BackendContext context = BackendContext.select(policy)) {
            MetaClusterRobustResult result = calculate(y, x, names, clusters,
                MetaGls.matrix(hierarchical.totalCovariance(), y.length), Double.NaN, options, correction, context);
            double[] beta = result.beta(), expected = hierarchical.beta();
            for (int j = 0; j < beta.length; j++) if (Math.abs(beta[j] - expected[j]) > 1e-8 * (1 + Math.abs(expected[j])))
                throw new IllegalArgumentException("studies or moderators do not match the hierarchical fit");
            return result;
        }
    }

    private static void validateOptions(MetaAnalysisOptions options, MetaRobustCorrection correction) {
        if (options == null || correction == null) throw new IllegalArgumentException("options and correction required");
        if (options.inferenceMethod() == MetaInferenceMethod.HARTUNG_KNAPP
                || options.inferenceMethod() == MetaInferenceMethod.MODIFIED_HARTUNG_KNAPP)
            throw new IllegalArgumentException("Hartung-Knapp is not a cluster sandwich correction; select NORMAL or STUDENT_T and CR0/CR1/CR2");
    }

    private static MetaClusterRobustResult calculate(double[] y, double[] x, List<String> names,
            String[] clusters, double[] v, double tau, MetaAnalysisOptions options,
            MetaRobustCorrection correction, BackendContext context) {
        int n = y.length, p = names.size() + 1;
        if (clusters == null || clusters.length != n) throw new IllegalArgumentException("one cluster label per effect required");
        Map<String, List<Integer>> groups = new LinkedHashMap<>();
        for (int i = 0; i < n; i++) {
            if (clusters[i] == null || clusters[i].isBlank()) throw new IllegalArgumentException("cluster labels must not be blank");
            groups.computeIfAbsent(clusters[i], ignored -> new ArrayList<>()).add(i);
        }
        int g = groups.size();
        if (g < 2) throw new IllegalArgumentException("cluster covariance is not estimable from one independent cluster");
        if (correction == MetaRobustCorrection.CR1P && g <= p)
            throw new IllegalArgumentException("CR1P requires more clusters than coefficients");
        if (options.inferenceMethod() != MetaInferenceMethod.NORMAL && correction != MetaRobustCorrection.CR2 && g <= p)
            throw new IllegalArgumentException("residual cluster t inference requires more clusters than coefficients");
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            if (!clusters[i].equals(clusters[j]) && Math.abs(v[i * n + j]) > 1e-12 * Math.sqrt(v[i * n + i] * v[j * n + j]))
                throw new IllegalArgumentException("working covariance must be block diagonal across independent clusters");
        ComputeBackend backend = context.backend();
        MetaGls.Fit fit = MetaGls.fit(y, x, p, v, backend);
        // Whiten independent clusters for score and Satterthwaite contractions.
        double[] xw = new double[n * p], ew = new double[n];
        List<int[]> indices = new ArrayList<>();
        for (List<Integer> group : groups.values()) {
            int[] ix = group.stream().mapToInt(Integer::intValue).toArray(); indices.add(ix); int m = ix.length;
            double[] block = new double[m * m];
            for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) block[a * m + b] = v[ix[a] * n + ix[b]];
            double[] root = inverseRoot(block, m, backend);
            for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) {
                ew[ix[a]] += root[a * m + b] * fit.residual()[ix[b]];
                for (int j = 0; j < p; j++) xw[ix[a] * p + j] += root[a * m + b] * x[ix[b] * p + j];
            }
        }
        double[] cov = new double[p * p];
        double[][][] dfVectors = correction == MetaRobustCorrection.CR2 ? new double[p][g][n] : null;
        int clusterIndex = 0;
        for (int[] ix : indices) {
            int m = ix.length; double[] a = MetaGls.identity(m);
            if (correction == MetaRobustCorrection.CR2) {
                double[] vi = new double[m * m], residualCov = new double[m * m];
                for (int r = 0; r < m; r++) for (int c = 0; c < m; c++) {
                    vi[r * m + c] = v[ix[r] * n + ix[c]];
                    residualCov[r * m + c] = vi[r * m + c];
                    for (int j = 0; j < p; j++) for (int k = 0; k < p; k++)
                        residualCov[r * m + c] -= x[ix[r] * p + j] * fit.bread()[j * p + k] * x[ix[c] * p + k];
                }
                // Symmetric Bell-McCaffrey adjustment A = V^.5 (V^.5 D V^.5)^-.5 V^.5.
                // In whitened score coordinates its transpose is V (V^.5 D V^.5)^-.5.
                double[] root = power(vi, m, 0.5, backend);
                double[] middle = multiply(multiply(root, residualCov, m), root, m);
                a = multiply(vi, inverseRoot(middle, m, backend), m);
            }
            // u = K_i' X_i B; each coefficient's adjusted score is u' e_i.
            double[] u = new double[m * p];
            for (int r = 0; r < m; r++) for (int s = 0; s < m; s++)
                for (int j = 0; j < p; j++) for (int k = 0; k < p; k++)
                    u[r * p + j] += a[r * m + s] * xw[ix[s] * p + k] * fit.bread()[k * p + j];
            double[] score = new double[p];
            for (int r = 0; r < m; r++) for (int j = 0; j < p; j++) score[j] += u[r * p + j] * ew[ix[r]];
            for (int j = 0; j < p; j++) for (int k = 0; k < p; k++) cov[j * p + k] += score[j] * score[k];
            if (correction == MetaRobustCorrection.CR2) for (int j = 0; j < p; j++) {
                double[] vector = dfVectors[j][clusterIndex], xtu = new double[p], bxtu = new double[p];
                for (int r = 0; r < m; r++) {
                    vector[ix[r]] = u[r * p + j];
                    for (int k = 0; k < p; k++) xtu[k] += xw[ix[r] * p + k] * u[r * p + j];
                }
                for (int k = 0; k < p; k++) for (int l = 0; l < p; l++) bxtu[k] += fit.bread()[k * p + l] * xtu[l];
                for (int i = 0; i < n; i++) for (int k = 0; k < p; k++) vector[i] -= xw[i * p + k] * bxtu[k];
            }
            clusterIndex++;
        }
        double multiplier = switch (correction) {
            case CR0, CR2 -> 1;
            case CR1 -> (double) g / (g - 1);
            case CR1P -> (double) g / (g - p);
            case CR1S -> (double) g * (n - 1) / ((g - 1.0) * (n - p));
        };
        for (int j = 0; j < cov.length; j++) cov[j] *= multiplier;
        double[] df = new double[p], se = new double[p], lo = new double[p], hi = new double[p]; Arrays.fill(df, g - p);
        for (int j = 0; j < p; j++) {
            se[j] = Math.sqrt(cov[j * p + j]);
            if (correction == MetaRobustCorrection.CR2) {
                double trace = 0, squareTrace = 0;
                for (int a = 0; a < g; a++) for (int b = 0; b < g; b++) {
                    double dot = 0; for (int i = 0; i < n; i++) dot += dfVectors[j][a][i] * dfVectors[j][b][i];
                    if (a == b) trace += dot; squareTrace += dot * dot;
                }
                if (!(squareTrace > 0)) throw new IllegalArgumentException("CR2 coefficient degrees of freedom are unidentifiable");
                df[j] = trace * trace / squareTrace;
            }
            double probability = (1 + options.confidenceLevel()) / 2;
            double critical = options.inferenceMethod() == MetaInferenceMethod.NORMAL
                ? Normal.quantile(probability, 0, 1, true, false) : T.quantile(probability, df[j], true, false);
            lo[j] = fit.beta()[j] - critical * se[j]; hi[j] = fit.beta()[j] + critical * se[j];
        }
        AssociationStatistics stats = options.inferenceMethod() == MetaInferenceMethod.NORMAL
            ? AssociationStatistics.normal(fit.beta(), se) : AssociationStatistics.studentT(fit.beta(), se, df,
                correction == MetaRobustCorrection.CR2 ? DegreesOfFreedomMethod.SATTERTHWAITE : DegreesOfFreedomMethod.CLUSTER);
        return new MetaClusterRobustResult(MetaMultilevelRegression.names(names, true), stats, cov, tau, g,
            context.provenance(), correction, lo, hi);
    }

    private static double[] inverseRoot(double[] matrix, int n, ComputeBackend backend) {
        return power(matrix, n, -0.5, backend);
    }
    private static double[] power(double[] matrix, int n, double exponent, ComputeBackend backend) {
        SymmetricEigenDecomposition eigen = backend.dsyev(matrix, n);
        double[] values = eigen.eigenvalues(), vectors = eigen.eigenvectors(), out = new double[n * n];
        double max = Arrays.stream(values).max().orElseThrow();
        for (int k = 0; k < n; k++) {
            if (!(values[k] > 1e-12 * max)) throw new IllegalArgumentException("covariance or CR2 residual block is singular; cluster leverage is too high");
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
                out[i * n + j] += vectors[i * n + k] * vectors[j * n + k] * Math.pow(values[k], exponent);
        }
        return out;
    }
    private static double[] multiply(double[] a, double[] b, int n) {
        double[] out = new double[n * n];
        for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) for (int j = 0; j < n; j++)
            out[i * n + j] += a[i * n + k] * b[k * n + j];
        return out;
    }
}
