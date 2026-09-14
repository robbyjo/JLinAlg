/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.special.Gamma;

/**
 * Variational Bayesian factor analysis ported from official PMBio PEER.
 * Input data are feature by sample; PEER internally models sample by feature.
 */
public final class Peer {
    private Peer() { }

    /** Official PEER dense-VBFA defaults, with an explicit reproducibility seed. */
    public record Options(int factors, int maximumIterations, double tolerance,
                          double varianceTolerance, boolean addMean, long seed,
                          double alphaA, double alphaB, double epsilonA,
                          double epsilonB, double covariatePrecision) {
        public Options {
            if (factors < 1 || maximumIterations < 1)
                throw new IllegalArgumentException("PEER factors and iterations must be positive");
            if (tolerance < 0.0 || varianceTolerance < 0.0
                    || !Double.isFinite(tolerance) || !Double.isFinite(varianceTolerance))
                throw new IllegalArgumentException("PEER tolerances must be finite and nonnegative");
            if (!(alphaA > 0.0 && alphaB > 0.0 && epsilonA > 0.0
                    && epsilonB > 0.0 && covariatePrecision > 0.0))
                throw new IllegalArgumentException("PEER prior parameters must be positive");
        }
        public static Options defaults(int factors) {
            return new Options(factors, 1000, 1e-3, 1e-8, false, 1L,
                0.001, 0.1, 0.1, 10.0, 100.0);
        }
    }

    public static PeerResult fit(double[][] featuresBySamples,
            double[][] covariates, Options options) {
        return fit(featuresBySamples, null, covariates, options);
    }

    /** Fits with optional feature-by-sample measurement variances. */
    public static PeerResult fit(double[][] featuresBySamples,
            double[][] measurementVariance, double[][] covariates, Options options) {
        int n = ConfounderMath.samples(featuresBySamples), d = featuresBySamples.length;
        if (measurementVariance != null) {
            if (ConfounderMath.samples(measurementVariance) != n
                    || measurementVariance.length != d)
                throw new IllegalArgumentException("measurement variance dimensions must equal data");
            for (double[] row : measurementVariance) for (double value : row)
                if (value < 0.0) throw new IllegalArgumentException("measurement variances must be nonnegative");
        }
        if (covariates != null) ConfounderMath.validateDesign(covariates, n, "covariates");
        int suppliedCovariates = covariates == null ? 0 : covariates[0].length;
        int nc = suppliedCovariates + (options.addMean() ? 1 : 0);
        int k = nc + options.factors();
        if (k >= n || k >= d)
            throw new IllegalArgumentException("PEER total factors must be less than samples and features");

        double[][] y = transpose(featuresBySamples);
        double[][] yVariance = measurementVariance == null
            ? filled(n, d, 0.01) : transpose(measurementVariance);
        double[][] observed = new double[n][nc];
        for (int i = 0; i < n; i++) {
            if (suppliedCovariates > 0)
                System.arraycopy(covariates[i], 0, observed[i], 0, suppliedCovariates);
            if (options.addMean()) observed[i][nc - 1] = 1.0;
        }

        Random random = new Random(options.seed());
        double[][] x = gaussian(n, k, random), w = gaussian(d, k, random);
        if (nc > 0) {
            double[][] coefficients = leastSquares(observed, y);
            for (int i = 0; i < n; i++) System.arraycopy(observed[i], 0, x[i], 0, nc);
            for (int feature = 0; feature < d; feature++)
                for (int c = 0; c < nc; c++) w[feature][c] = coefficients[c][feature];
        }
        double[] priorPrecision = new double[k]; Arrays.fill(priorPrecision, 1.0);
        for (int c = 0; c < nc; c++) priorPrecision[c] = options.covariatePrecision();
        double[][] priorMean = new double[n][k];
        for (int i = 0; i < n; i++) if (nc > 0)
            System.arraycopy(observed[i], 0, priorMean[i], 0, nc);

        double[] alpha = new double[k], lnAlpha = new double[k];
        double[] epsilon = new double[d], lnEpsilon = new double[d];
        Arrays.fill(epsilon, options.epsilonA() / options.epsilonB());
        Arrays.fill(lnEpsilon, Gamma.digamma(options.epsilonA()) - Math.log(options.epsilonB()));
        double[][] xCov = identity(k), x2 = crossPlus(x, xCov, n);
        double[][] w2 = cross(w), wXPrecision = cross(w);
        updateAlpha(alpha, lnAlpha, w2, d, options);

        List<Double> bounds = new ArrayList<>(), residualVariances = new ArrayList<>();
        double lastBound = Double.NEGATIVE_INFINITY, lastResidual = Double.POSITIVE_INFINITY;
        boolean converged = false;
        int completed = 0;
        for (int iteration = 0; iteration < options.maximumIterations(); iteration++) {
            WUpdate wu = updateW(y, x, x2, alpha, epsilon);
            w = wu.mean(); w2 = wu.secondMomentSum(); wXPrecision = wu.xPrecision();
            updateAlpha(alpha, lnAlpha, w2, d, options);
            XUpdate xu = updateX(y, w, wXPrecision, epsilon, priorMean, priorPrecision, nc);
            x = xu.mean(); xCov = xu.covariance(); x2 = xu.secondMomentSum();
            updateEpsilon(epsilon, lnEpsilon, y, yVariance, x, x2, w,
                wu.alphaAtUpdate(), wu.xSecondAtUpdate(), wu.epsilonAtUpdate(), options);

            double[][] residual = subtract(y, multiply(x, transpose(w)));
            double residualVariance = meanSquares(residual);
            double bound = evidenceBound(y, yVariance, x, xCov, x2, w, w2,
                alpha, lnAlpha, epsilon, lnEpsilon, priorMean, priorPrecision,
                wu, options);
            bounds.add(bound); residualVariances.add(residualVariance); completed = iteration + 1;
            double deltaBound = bound - lastBound;
            double deltaResidual = lastResidual - residualVariance;
            if (iteration > 0 && (Math.abs(deltaBound) < options.tolerance()
                    || Math.abs(deltaResidual) < options.varianceTolerance())) {
                converged = true; break;
            }
            lastBound = bound; lastResidual = residualVariance;
        }
        double[][] residual = subtract(y, multiply(x, transpose(w)));
        double[][] hiddenX = sliceColumns(x, nc, k), hiddenW = sliceColumns(w, nc, k);
        return new PeerResult(hiddenX, hiddenW, transpose(residual),
            Arrays.copyOfRange(alpha, nc, k), epsilon.clone(), completed, converged,
            bounds, residualVariances);
    }

    private static WUpdate updateW(double[][] y, double[][] x, double[][] x2,
            double[] alpha, double[] epsilon) {
        int n = y.length, d = y[0].length, k = alpha.length;
        double[][] w = new double[d][k], w2 = new double[k][k], xPrecision = new double[k][k];
        double logDetCovariance = 0.0;
        double[][] alphaMatrix = diagonal(alpha);
        for (int feature = 0; feature < d; feature++) {
            double[][] precision = add(alphaMatrix, scale(x2, epsilon[feature]));
            double[][] covariance = inversePositiveDefinite(precision);
            logDetCovariance += logDetPositiveDefinite(covariance);
            double[] xty = new double[k];
            for (int i = 0; i < n; i++) for (int factor = 0; factor < k; factor++)
                xty[factor] += x[i][factor] * y[i][feature];
            double[] mean = multiply(covariance, xty);
            for (int factor = 0; factor < k; factor++) w[feature][factor] = epsilon[feature] * mean[factor];
            double[][] second = add(covariance, outer(w[feature], w[feature]));
            addInPlace(w2, second); addScaledInPlace(xPrecision, second, epsilon[feature]);
        }
        return new WUpdate(w, w2, xPrecision, logDetCovariance,
            alpha.clone(), copy(x2), epsilon.clone());
    }

    private static XUpdate updateX(double[][] y, double[][] w,
            double[][] wXPrecision, double[] epsilon, double[][] priorMean,
            double[] priorPrecision, int nc) {
        int n = y.length, d = y[0].length, k = priorPrecision.length;
        double[][] covariance = inversePositiveDefinite(add(wXPrecision, diagonal(priorPrecision)));
        double[][] x = new double[n][k];
        for (int i = 0; i < n; i++) for (int feature = 0; feature < d; feature++)
            for (int factor = 0; factor < k; factor++)
                x[i][factor] += epsilon[feature] * y[i][feature] * w[feature][factor];
        for (int i = 0; i < n; i++) for (int factor = 0; factor < k; factor++)
            x[i][factor] += priorMean[i][factor] * priorPrecision[factor];
        x = multiply(x, covariance);
        for (int i = 0; i < n; i++) for (int c = 0; c < nc; c++) x[i][c] = priorMean[i][c];
        return new XUpdate(x, covariance, crossPlus(x, covariance, n));
    }

    private static void updateAlpha(double[] alpha, double[] lnAlpha,
            double[][] w2, int features, Options options) {
        double a = options.alphaA() + 0.5 * features;
        for (int factor = 0; factor < alpha.length; factor++) {
            double b = options.alphaB() + 0.5 * w2[factor][factor];
            alpha[factor] = a / b; lnAlpha[factor] = Gamma.digamma(a) - Math.log(b);
        }
    }

    private static void updateEpsilon(double[] epsilon, double[] lnEpsilon,
            double[][] y, double[][] yVariance, double[][] x, double[][] x2,
            double[][] w, double[] alphaLast, double[][] x2Last,
            double[] epsilonLast, Options options) {
        int n = y.length, d = y[0].length;
        double a = options.epsilonA() + 0.5 * n;
        for (int feature = 0; feature < d; feature++) {
            double b1 = 0.0, b2 = 0.0;
            for (int i = 0; i < n; i++) {
                b1 += yVariance[i][feature] + y[i][feature] * y[i][feature];
                b2 += y[i][feature] * dot(x[i], w[feature]);
            }
            double[][] oldCovariance = inversePositiveDefinite(add(diagonal(alphaLast),
                scale(x2Last, epsilonLast[feature])));
            double[][] wSecond = add(oldCovariance, outer(w[feature], w[feature]));
            double b3 = elementProductSum(x2, wSecond);
            double b = options.epsilonB() + 0.5 * b1 - b2 + 0.5 * b3;
            if (!(b > 0.0) || !Double.isFinite(b))
                throw new IllegalArgumentException("PEER noise posterior became non-positive");
            epsilon[feature] = a / b; lnEpsilon[feature] = Gamma.digamma(a) - Math.log(b);
        }
    }

    private static double evidenceBound(double[][] y, double[][] yVariance,
            double[][] x, double[][] xCov, double[][] x2, double[][] w,
            double[][] w2, double[] alpha, double[] lnAlpha, double[] epsilon,
            double[] lnEpsilon, double[][] priorMean, double[] priorPrecision,
            WUpdate wu, Options options) {
        int n = y.length, d = y[0].length, k = alpha.length;
        double logProb = 0.5 * n * (sum(lnEpsilon) - d * Math.log(2.0 * Math.PI));
        for (int feature = 0; feature < d; feature++) {
            double b1 = 0.0, b2 = 0.0;
            for (int i = 0; i < n; i++) {
                b1 += yVariance[i][feature] + y[i][feature] * y[i][feature];
                b2 += y[i][feature] * dot(x[i], w[feature]);
            }
            double[][] covariance = inversePositiveDefinite(add(diagonal(wu.alphaAtUpdate()),
                scale(wu.xSecondAtUpdate(), wu.epsilonAtUpdate()[feature])));
            double b3 = elementProductSum(x2, add(covariance, outer(w[feature], w[feature])));
            logProb += -0.5 * epsilon[feature] * b1 + epsilon[feature] * b2
                - 0.5 * epsilon[feature] * b3;
        }
        double wBound = -0.5 * d * (k * Math.log(2.0 * Math.PI) - sum(lnAlpha))
            - 0.5 * diagonalProduct(alpha, w2)
            + 0.5 * d * k * Math.log(2.0 * Math.PI) + 0.5 * wu.logDetCovariance();
        double xBound = -0.5 * n * k * Math.log(2.0 * Math.PI)
            + 0.5 * n * Arrays.stream(priorPrecision).map(Math::log).sum();
        for (int factor = 0; factor < k; factor++) {
            double quadratic = n * xCov[factor][factor];
            for (int i = 0; i < n; i++) quadratic += x[i][factor] * x[i][factor]
                - 2.0 * x[i][factor] * priorMean[i][factor]
                + priorMean[i][factor] * priorMean[i][factor];
            xBound -= 0.5 * priorPrecision[factor] * quadratic;
        }
        xBound += 0.5 * n * (k * Math.log(2.0 * Math.PI) + k
            + logDetPositiveDefinite(xCov));
        return logProb + wBound + xBound
            + gammaBound(alpha, lnAlpha, options.alphaA(), options.alphaB(), d, w2)
            + gammaBound(epsilon, lnEpsilon, options.epsilonA(), options.epsilonB(), n, null);
    }

    private static double gammaBound(double[] mean, double[] lnMean,
            double pa, double pb, int plate, double[][] second) {
        double a = pa + 0.5 * plate, result = 0.0;
        for (int i = 0; i < mean.length; i++) {
            double b = a / mean[i];
            double entropy = a - Math.log(b) + Gamma.logGamma(a) + (1.0 - a) * Gamma.digamma(a);
            result += pa * Math.log(pb) + (pa - 1.0) * lnMean[i]
                - pb * mean[i] - Gamma.logGamma(pa) + entropy;
        }
        return result;
    }

    private record WUpdate(double[][] mean, double[][] secondMomentSum,
            double[][] xPrecision, double logDetCovariance,
            double[] alphaAtUpdate, double[][] xSecondAtUpdate,
            double[] epsilonAtUpdate) { }
    private record XUpdate(double[][] mean, double[][] covariance,
                           double[][] secondMomentSum) { }

    private static double[][] leastSquares(double[][] x, double[][] y) {
        return multiply(inversePositiveDefinite(cross(x)), multiply(transpose(x), y));
    }
    private static double[][] gaussian(int rows, int columns, Random random) { double[][] r = new double[rows][columns]; for (double[] row : r) for (int j = 0; j < columns; j++) row[j] = random.nextGaussian(); return r; }
    private static double[][] filled(int rows, int columns, double value) { double[][] r = new double[rows][columns]; for (double[] row : r) Arrays.fill(row, value); return r; }
    private static double[][] identity(int n) { double[][] r = new double[n][n]; for (int i = 0; i < n; i++) r[i][i] = 1.0; return r; }
    private static double[][] diagonal(double[] values) { double[][] r = new double[values.length][values.length]; for (int i = 0; i < values.length; i++) r[i][i] = values[i]; return r; }
    private static double[][] cross(double[][] x) { return multiply(transpose(x), x); }
    private static double[][] crossPlus(double[][] x, double[][] covariance, int multiplier) { double[][] r = cross(x); addScaledInPlace(r, covariance, multiplier); return r; }
    private static double[][] transpose(double[][] x) { double[][] r = new double[x[0].length][x.length]; for (int i = 0; i < x.length; i++) for (int j = 0; j < x[0].length; j++) r[j][i] = x[i][j]; return r; }
    private static double[][] multiply(double[][] a, double[][] b) { if (a[0].length != b.length) throw new IllegalArgumentException("matrix dimensions differ"); double[][] r = new double[a.length][b[0].length]; for (int i = 0; i < a.length; i++) for (int k = 0; k < b.length; k++) { double v = a[i][k]; for (int j = 0; j < b[0].length; j++) r[i][j] += v * b[k][j]; } return r; }
    private static double[] multiply(double[][] a, double[] x) { double[] r = new double[a.length]; for (int i = 0; i < a.length; i++) r[i] = dot(a[i], x); return r; }
    private static double[][] add(double[][] a, double[][] b) { double[][] r = copy(a); addInPlace(r, b); return r; }
    private static double[][] scale(double[][] a, double value) { double[][] r = copy(a); for (double[] row : r) for (int j = 0; j < row.length; j++) row[j] *= value; return r; }
    private static void addInPlace(double[][] a, double[][] b) { addScaledInPlace(a, b, 1.0); }
    private static void addScaledInPlace(double[][] a, double[][] b, double scale) { for (int i = 0; i < a.length; i++) for (int j = 0; j < a[i].length; j++) a[i][j] += scale * b[i][j]; }
    private static double[][] subtract(double[][] a, double[][] b) { double[][] r = copy(a); for (int i = 0; i < a.length; i++) for (int j = 0; j < a[i].length; j++) r[i][j] -= b[i][j]; return r; }
    private static double[][] outer(double[] a, double[] b) { double[][] r = new double[a.length][b.length]; for (int i = 0; i < a.length; i++) for (int j = 0; j < b.length; j++) r[i][j] = a[i] * b[j]; return r; }
    private static double[][] sliceColumns(double[][] x, int from, int to) { double[][] r = new double[x.length][to - from]; for (int i = 0; i < x.length; i++) System.arraycopy(x[i], from, r[i], 0, to - from); return r; }
    private static double[][] copy(double[][] x) { double[][] r = new double[x.length][]; for (int i = 0; i < x.length; i++) r[i] = x[i].clone(); return r; }
    private static double dot(double[] a, double[] b) { double r = 0.0; for (int i = 0; i < a.length; i++) r += a[i] * b[i]; return r; }
    private static double sum(double[] x) { double r = 0.0; for (double v : x) r += v; return r; }
    private static double elementProductSum(double[][] a, double[][] b) { double r = 0.0; for (int i = 0; i < a.length; i++) for (int j = 0; j < a[i].length; j++) r += a[i][j] * b[i][j]; return r; }
    private static double diagonalProduct(double[] diagonal, double[][] matrix) { double r = 0.0; for (int i = 0; i < diagonal.length; i++) r += diagonal[i] * matrix[i][i]; return r; }
    private static double meanSquares(double[][] x) { double r = 0.0; for (double[] row : x) for (double v : row) r += v * v; return r / (x.length * x[0].length); }

    private static double[][] inversePositiveDefinite(double[][] matrix) {
        int n = matrix.length; double[][] l = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
            double value = matrix[i][j]; for (int k = 0; k < j; k++) value -= l[i][k] * l[j][k];
            if (i == j) { if (!(value > 0.0) || !Double.isFinite(value)) throw new IllegalArgumentException("PEER precision is not positive definite"); l[i][j] = Math.sqrt(value); }
            else l[i][j] = value / l[j][j];
        }
        double[][] inverse = new double[n][n];
        for (int c = 0; c < n; c++) { double[] y = new double[n], x = new double[n];
            for (int i = 0; i < n; i++) { double value = i == c ? 1.0 : 0.0; for (int k = 0; k < i; k++) value -= l[i][k] * y[k]; y[i] = value / l[i][i]; }
            for (int i = n - 1; i >= 0; i--) { double value = y[i]; for (int k = i + 1; k < n; k++) value -= l[k][i] * x[k]; x[i] = value / l[i][i]; }
            for (int i = 0; i < n; i++) inverse[i][c] = x[i]; }
        return inverse;
    }
    private static double logDetPositiveDefinite(double[][] matrix) { int n = matrix.length; double[][] l = new double[n][n]; double result = 0.0; for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) { double value = matrix[i][j]; for (int k = 0; k < j; k++) value -= l[i][k] * l[j][k]; if (i == j) { if (!(value > 0.0)) throw new IllegalArgumentException("matrix is not positive definite"); l[i][j] = Math.sqrt(value); result += 2.0 * Math.log(l[i][j]); } else l[i][j] = value / l[j][j]; } return result; }
}
