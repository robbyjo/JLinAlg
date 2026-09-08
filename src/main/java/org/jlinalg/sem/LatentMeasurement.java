/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

/** Joint confirmatory measurement fitting, plus a legacy descriptive PCA extractor. */
public final class LatentMeasurement {
    private LatentMeasurement() { }

    /** Joint latent RAM ML; specify factor scales and loadings in the model. */
    public static SemFitResult fit(double[][] data, SemModel model) {
        if(model.latentVariables().isEmpty())throw new IllegalArgumentException("measurement model needs a latent variable");
        return Sem.fit(data,model);
    }

    /** Descriptive PCA extraction, not a fitted latent measurement/structural model. */
    public static Result fit(double[][] data, int factorCount, int maximumIterations,
                             double tolerance) {
        if (data == null || data.length < 3 || data[0] == null || factorCount < 1
                || factorCount >= data[0].length || maximumIterations < 1 || !(tolerance > 0.0))
            throw new IllegalArgumentException("latent measurement inputs are invalid");
        int observations = data.length, indicators = data[0].length;
        double[] means = new double[indicators];
        for (double[] row : data) { if (row == null || row.length != indicators) throw new IllegalArgumentException("indicator data must be rectangular"); for (int column = 0; column < indicators; column++) { if (!Double.isFinite(row[column])) throw new IllegalArgumentException("latent measurement requires complete finite data"); means[column] += row[column]; } }
        for (int column = 0; column < indicators; column++) means[column] /= observations;
        double[] covariance = new double[indicators * indicators];
        for (double[] row : data) for (int left = 0; left < indicators; left++) for (int right = 0; right < indicators; right++) covariance[left * indicators + right] += (row[left] - means[left]) * (row[right] - means[right]) / observations;
        double[] loading = new double[indicators * factorCount]; double[] residual = covariance.clone();
        for (int factor = 0; factor < factorCount; factor++) {
            double[] vector = new double[indicators];
            int seed=0;for(int column=1;column<indicators;column++)if(residual[column*indicators+column]>residual[seed*indicators+seed])seed=column;
            vector[seed]=1;
            for (int iteration = 0; iteration < maximumIterations; iteration++) { double[] next = multiply(residual, vector, indicators); double norm = norm(next); if (!(norm > 0.0)) break; for (int column = 0; column < indicators; column++) next[column] /= norm; double change = 0.0; for (int column = 0; column < indicators; column++) change = Math.max(change, Math.abs(next[column] - vector[column])); vector = next; if (change <= tolerance) break; }
            double[] product = multiply(residual, vector, indicators); double eigenvalue = dot(vector, product); double scale = Math.sqrt(Math.max(0.0, eigenvalue)); for (int column = 0; column < indicators; column++) loading[column * factorCount + factor] = vector[column] * scale; for (int left = 0; left < indicators; left++) for (int right = 0; right < indicators; right++) residual[left * indicators + right] -= loading[left * factorCount + factor] * loading[right * factorCount + factor];
        }
        double[] residualVariances = new double[indicators]; for (int column = 0; column < indicators; column++) residualVariances[column] = Math.max(1e-10, residual[column * indicators + column]);
        return new Result(means, loading, residualVariances, observations, indicators, factorCount);
    }

    public static Result fit(double[][] data, int factorCount) { return fit(data, factorCount, 500, 1e-8); }
    private static double[] multiply(double[] matrix, double[] vector, int dimension) { double[] result = new double[dimension]; for (int row = 0; row < dimension; row++) for (int column = 0; column < dimension; column++) result[row] += matrix[row * dimension + column] * vector[column]; return result; }
    private static double dot(double[] left, double[] right) { double result = 0.0; for (int index = 0; index < left.length; index++) result += left[index] * right[index]; return result; }
    private static double norm(double[] values) { return Math.sqrt(dot(values, values)); }

    public record Result(double[] means, double[] loadings, double[] residualVariances,
                         int observations, int indicatorCount, int factorCount) {
        public Result { means = means.clone(); loadings = loadings.clone(); residualVariances = residualVariances.clone(); }
        public double[] means() { return means.clone(); }
        public double[] loadings() { return loadings.clone(); }
        public double[] residualVariances() { return residualVariances.clone(); }
    }
}
