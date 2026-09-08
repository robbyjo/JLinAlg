/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Projected BFGS on scaled fixed coefficients and a nonnegative variance coordinate.
 * Numerical derivatives always differentiate the complete adaptive integral. */
final class QuadratureOptimizer {
    private QuadratureOptimizer() { }
    @FunctionalInterface interface Objective {
        double value(double[] parameters);
        default double[] exactGradient(double[] parameters) { return null; }
    }
    record Fit(double[] parameters, double value) { }

    static Fit minimize(Objective f, double[] start, int variance, boolean fixedVariance,
            GlmmQuadratureOptions options) {
        double[] x = start.clone(); double fx = f.value(x);
        if (!Double.isFinite(fx)) return new Fit(x, fx);
        int n = x.length;
        double[][] inverse = identity(n);
        double[] g = gradient(f, x, variance, fixedVariance);
        for (int iteration = 0; iteration < options.maximumIterations(); iteration++) {
            if (projectedNorm(g, x, variance, fixedVariance) <= options.gradientTolerance() * .25) break;
            double[] direction = new double[n];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) direction[i] -= inverse[i][j] * g[j];
            if (fixedVariance || (x[variance] == 0 && direction[variance] < 0)) direction[variance] = 0;
            if (!(dot(g, direction) < 0)) {
                inverse = identity(n);
                for (int i = 0; i < n; i++) direction[i] = -g[i];
                if (fixedVariance || (x[variance] == 0 && direction[variance] < 0)) direction[variance] = 0;
            }
            double max = 0;
            for (double d : direction) max = Math.max(max, Math.abs(d));
            if (max > 5) for (int i = 0; i < n; i++) direction[i] *= 5 / max;
            double[] candidate = x.clone(); double next = Double.POSITIVE_INFINITY;
            boolean accepted = false;
            for (double alpha = 1; alpha >= 1e-12; alpha *= .5) {
                for (int i = 0; i < n; i++) candidate[i] = x[i] + alpha * direction[i];
                candidate[variance] = fixedVariance ? 0 : Math.max(0, candidate[variance]);
                double slope = 0;
                for (int i = 0; i < n; i++) slope += g[i] * (candidate[i] - x[i]);
                next = f.value(candidate);
                if (Double.isFinite(next) && slope < 0 && next <= fx + 1e-4 * slope) {
                    accepted = true; break;
                }
            }
            if (!accepted) break; // No success based merely on an unchanged objective.
            double[] ng = gradient(f, candidate, variance, fixedVariance);
            double[] s = new double[n], y = new double[n], hy = new double[n];
            for (int i = 0; i < n; i++) { s[i] = candidate[i] - x[i]; y[i] = ng[i] - g[i]; }
            double sy = dot(s, y);
            if (sy > 1e-12 * Math.sqrt(dot(s, s) * dot(y, y))) {
                for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) hy[i] += inverse[i][j] * y[j];
                double factor = (sy + dot(y, hy)) / (sy * sy);
                for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
                    inverse[i][j] += factor * s[i] * s[j] - (hy[i] * s[j] + s[i] * hy[j]) / sy;
            } else inverse = identity(n);
            x = candidate.clone(); g = ng; fx = next;
        }
        return new Fit(x, fx);
    }

    static double[] gradient(Objective f, double[] x, int variance, boolean fixedVariance) {
        double[] exact = f.exactGradient(x);
        if (exact != null) {
            if (fixedVariance) exact[variance] = 0;
            return exact;
        }
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            if (fixedVariance && i == variance) continue;
            double h = 1e-5 * Math.max(1, Math.abs(x[i]));
            double[] a = x.clone(), b = x.clone();
            if (i == variance && x[i] < h) {
                a[i] += h; b[i] += 2 * h;
                result[i] = (-3 * f.value(x) + 4 * f.value(a) - f.value(b)) / (2 * h);
            } else {
                a[i] += h; b[i] -= h;
                result[i] = (f.value(a) - f.value(b)) / (2 * h);
            }
        }
        return result;
    }
    static double projectedNorm(double[] gradient, double[] x, int variance, boolean fixedVariance) {
        double result = 0;
        for (int i = 0; i < x.length; i++) {
            if (i == variance && (fixedVariance || (x[i] == 0 && gradient[i] > 0))) continue;
            result = Math.max(result, Math.abs(gradient[i]));
        }
        return result;
    }

    static double[][] covariance(Objective f, double[] x, int dimension) {
        double[][] hessian = new double[dimension][dimension];
        double center = f.value(x);
        double[] steps = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            steps[i] = 2e-4 * Math.max(1, Math.abs(x[i]));
            if (i == x.length - 1) steps[i] = Math.min(steps[i], x[i] / 4);
            if (!(steps[i] > 1e-8)) return null;
        }
        for (int i = 0; i < dimension; i++) {
            double[] a = x.clone(), b = x.clone(); a[i] += steps[i]; b[i] -= steps[i];
            hessian[i][i] = (f.value(a) - 2 * center + f.value(b)) / (steps[i] * steps[i]);
            for (int j = 0; j < i; j++) {
                double[] pp = x.clone(), pm = x.clone(), mp = x.clone(), mm = x.clone();
                pp[i] += steps[i]; pp[j] += steps[j]; pm[i] += steps[i]; pm[j] -= steps[j];
                mp[i] -= steps[i]; mp[j] += steps[j]; mm[i] -= steps[i]; mm[j] -= steps[j];
                hessian[i][j] = (f.value(pp) - f.value(pm) - f.value(mp) + f.value(mm)) / (4 * steps[i] * steps[j]);
                hessian[j][i] = hessian[i][j];
            }
        }
        return inversePositiveDefinite(hessian);
    }

    static double[][] inversePositiveDefinite(double[][] a) {
        int n = a.length; double[][] l = new double[n][n]; double scale = 0;
        for (int i = 0; i < n; i++) scale = Math.max(scale, Math.abs(a[i][i]));
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
            double value = a[i][j];
            for (int k = 0; k < j; k++) value -= l[i][k] * l[j][k];
            if (i == j) {
                if (!Double.isFinite(value) || value <= Math.max(1e-10 * scale, 1e-10)) return null;
                l[i][j] = Math.sqrt(value);
            } else l[i][j] = value / l[j][j];
        }
        double[][] inverse = new double[n][n];
        for (int column = 0; column < n; column++) {
            double[] v = new double[n];
            for (int i = 0; i < n; i++) {
                double sum = i == column ? 1 : 0;
                for (int j = 0; j < i; j++) sum -= l[i][j] * v[j];
                v[i] = sum / l[i][i];
            }
            for (int i = n - 1; i >= 0; i--) {
                double sum = v[i];
                for (int j = i + 1; j < n; j++) sum -= l[j][i] * inverse[j][column];
                inverse[i][column] = sum / l[i][i];
            }
        }
        return inverse;
    }
    private static double[][] identity(int n) {
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) result[i][i] = 1;
        return result;
    }
    private static double dot(double[] a, double[] b) {
        double result = 0;
        for (int i = 0; i < a.length; i++) result += a[i] * b[i];
        return result;
    }
}
