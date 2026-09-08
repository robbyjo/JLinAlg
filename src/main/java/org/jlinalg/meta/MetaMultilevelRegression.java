/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.function.ToDoubleFunction;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Known-covariance GLS and joint ML/REML estimation of grouped random coefficients. */
public final class MetaMultilevelRegression {
    private MetaMultilevelRegression() { }
    public enum Estimation { ML, REML }

    /** Fit with a supplied, known total covariance (no covariance estimation). */
    public static MetaMultilevelRegressionResult fit(List<MetaStudy> studies, double[][] moderators,
            List<String> names, double[][] covariance, boolean includeIntercept, BackendPolicy policy) {
        double[] y = MetaMath.data(studies).effects();
        double[] x = MetaGls.design(moderators, names, y.length, includeIntercept);
        double[] v = MetaGls.matrix(covariance, y.length); int p = x.length / y.length;
        try (BackendContext context = BackendContext.select(policy)) {
            MetaGls.Fit f = MetaGls.fit(y, x, p, v, context.backend());
            return new MetaMultilevelRegressionResult(names(names, includeIntercept), f.beta(), f.bread(),
                f.q(), y.length - p, f.ml(), true);
        }
    }

    /** Estimate all random covariance matrices jointly; known sampling covariance is not rescaled.
     * The optimizer budget is maximumIterations times the number of Cholesky coordinates per start.
     * Inference is plug-in z/t or residual Hartung-Knapp; it does not include covariance-parameter uncertainty.
     * tauSquaredEstimator must be REML here; use estimation to select ML versus REML. */
    public static MetaHierarchicalResult fit(List<MetaStudy> studies, double[][] moderators,
            List<String> moderatorNames, double[][] samplingCovariance, boolean includeIntercept,
            List<MetaRandomEffect> randomEffects, Estimation estimation,
            MetaAnalysisOptions options, BackendPolicy policy) {
        double[] y = MetaMath.data(studies).effects(); int n = y.length;
        double[] x = MetaGls.design(moderators, moderatorNames, n, includeIntercept), v = MetaGls.matrix(samplingCovariance, n);
        if (randomEffects == null || randomEffects.isEmpty() || estimation == null || options == null || policy == null)
            throw new IllegalArgumentException("random effects, estimation and controls required");
        if (options.method() != MetaAnalysisMethod.RANDOM_EFFECT || options.tauSquaredEstimator() != TauSquaredEstimator.REML)
            throw new IllegalArgumentException("hierarchical covariance uses joint ML/REML; scalar DL/PM and fixed-effect options do not apply");
        List<MetaRandomEffect> terms = List.copyOf(randomEffects);
        HashSet<String> unique = new HashSet<>(); int d = 0, p = x.length / n;
        for (MetaRandomEffect term : terms) {
            if (term.rows() != n || !unique.add(term.name())) throw new IllegalArgumentException("random terms need unique names and matching rows");
            d += term.parameters();
        }
        try (BackendContext context = BackendContext.select(policy)) {
            ComputeBackend backend = context.backend();
            MetaGls.Fit initial = MetaGls.fit(y, x, p, v, backend);
            checkIdentifiable(terms, y, x, p, estimation, backend);
            double variance = 0;
            for (double e : initial.residual()) variance += e * e / (n - p);
            double scale = Math.sqrt(Math.max(1e-8, variance) / terms.size());
            double[] scales = new double[d], start = new double[d]; int k = 0;
            for (MetaRandomEffect term : terms) for (int a = 0; a < term.columns(); a++) {
                double norm = 0; for (int i = 0; i < n; i++) norm += term.z(i, a) * term.z(i, a) / n;
                for (int b = 0; b <= a; b++) if (term.structure() == MetaRandomEffect.Structure.UNSTRUCTURED || a == b) {
                    scales[k] = scale / Math.sqrt(norm); start[k++] = a == b ? 0.5 : 0;
                }
            }
            ToDoubleFunction<double[]> objective = point -> {
                double[] total = total(v, terms, physical(point, scales));
                for (double value : total) if (!Double.isFinite(value)) return Double.POSITIVE_INFINITY;
                MetaGls.Fit fit = MetaGls.fit(y, x, p, total, backend);
                return -(estimation == Estimation.REML ? fit.reml() : fit.ml());
            };
            MetaCovarianceOptimizer.Result best = MetaCovarianceOptimizer.minimize(objective, start,
                options.maximumIterations() * d, options.tolerance());
            for (int j = 0; j < d; j++) start[j] *= 3;
            MetaCovarianceOptimizer.Result second = MetaCovarianceOptimizer.minimize(objective, start,
                options.maximumIterations() * d, options.tolerance());
            if (second.objective() < best.objective() || (!best.converged() && second.converged()
                    && second.objective() <= best.objective() + options.tolerance())) best = second;
            double[] theta = physical(best.point(), scales), total = total(v, terms, theta);
            MetaGls.Fit f = MetaGls.fit(y, x, p, total, backend);
            List<double[]> covariances = new ArrayList<>(); k = 0;
            for (MetaRandomEffect term : terms) { covariances.add(term.covariance(theta, k)); k += term.parameters(); }
            double[][] totalMatrix = new double[n][n];
            for (int i = 0; i < n; i++) System.arraycopy(total, i * n, totalMatrix[i], 0, n);
            double[] cov = f.bread().clone(), se = new double[p], lo = new double[p], hi = new double[p];
            double inflation = MetaAnalysis.inferenceScale(options, f.q(), n - p), critical = MetaAnalysis.critical(options, n - p);
            for (int j = 0; j < cov.length; j++) cov[j] *= inflation;
            for (int j = 0; j < p; j++) {
                se[j] = Math.sqrt(cov[j * p + j]); lo[j] = f.beta()[j] - critical * se[j]; hi[j] = f.beta()[j] + critical * se[j];
            }
            return new MetaHierarchicalResult(new MetaMultilevelRegressionResult(names(moderatorNames, includeIntercept),
                f.beta(), cov, f.q(), n - p, estimation == Estimation.REML ? f.reml() : f.ml(), best.converged()),
                terms.stream().map(MetaRandomEffect::name).toList(), covariances, totalMatrix,
                MetaAnalysis.statistics(options, f.beta(), se, n - p), lo, hi, estimation, best.iterations());
        }
    }

    private static double[] physical(double[] point, double[] scales) {
        double[] result = point.clone(); for (int i = 0; i < result.length; i++) result[i] *= scales[i]; return result;
    }

    private static double[] total(double[] sampling, List<MetaRandomEffect> terms, double[] theta) {
        double[] v = sampling.clone(); int offset = 0, n = terms.get(0).rows();
        for (MetaRandomEffect term : terms) {
            int q = term.columns(); double[] g = term.covariance(theta, offset); offset += term.parameters();
            for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) if (term.sameGroup(i, j)) {
                double value = 0;
                for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) value += term.z(i, a) * g[a * q + b] * term.z(j, b);
                v[i * n + j] += value; if (i != j) v[j * n + i] += value;
            }
        }
        return v;
    }

    /** Reject covariance bases that are confounded with each other or the REML fixed space. */
    private static void checkIdentifiable(List<MetaRandomEffect> terms, double[] y, double[] x,
            int p, Estimation estimation, ComputeBackend backend) {
        int n = y.length; List<double[]> basis = new ArrayList<>();
        double[] projection = MetaGls.identity(n);
        if (estimation == Estimation.REML) {
            double[] bread = MetaGls.fit(y, x, p, MetaGls.identity(n), backend).bread();
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
                for (int a = 0; a < p; a++) for (int b = 0; b < p; b++)
                    projection[i * n + j] -= x[i * p + a] * bread[a * p + b] * x[j * p + b];
        }
        for (MetaRandomEffect term : terms) for (int a = 0; a < term.columns(); a++)
            for (int b = 0; b <= a; b++) if (a == b || term.structure() == MetaRandomEffect.Structure.UNSTRUCTURED) {
                double[] component = new double[n * n];
                for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (term.sameGroup(i, j))
                    component[i * n + j] = term.z(i, a) * term.z(j, b) + (a == b ? 0 : term.z(i, b) * term.z(j, a));
                double[] projected = multiply(multiply(projection, component, n), projection, n);
                double norm = 0; for (double value : projected) norm += value * value;
                if (!(norm > 1e-20)) throw new IllegalArgumentException("random covariance is unidentifiable after fixed effects");
                for (int j = 0; j < projected.length; j++) projected[j] /= Math.sqrt(norm);
                basis.add(projected);
            }
        int d = basis.size(); double[] gram = new double[d * d];
        for (int a = 0; a < d; a++) for (int b = 0; b < d; b++)
            for (int j = 0; j < n * n; j++) gram[a * d + b] += basis.get(a)[j] * basis.get(b)[j];
        for (double value : backend.dsyev(gram, d).eigenvalues()) if (value < 1e-9)
            throw new IllegalArgumentException("random covariance components are linearly dependent");
    }
    private static double[] multiply(double[] a, double[] b, int n) {
        double[] result = new double[n * n];
        for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) for (int j = 0; j < n; j++) result[i * n + j] += a[i * n + k] * b[k * n + j];
        return result;
    }
    static List<String> names(List<String> names, boolean intercept) {
        List<String> result = new ArrayList<>(); if (intercept) result.add("(Intercept)"); result.addAll(names); return result;
    }
}
