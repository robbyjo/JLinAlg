/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import org.jlinalg.compute.BackendPolicy;

/** Gaussian FIML saturated moment estimation followed by covariance-structure fitting. */
public final class SemFiml {
    private SemFiml() { }

    public static SemFimlResult fit(double[][] data, SemModel model) {
        return fit(data, model, 200, 1e-8, BackendPolicy.PREFERRED);
    }

    public static SemFimlResult fit(double[][] data, SemModel model, int maximumIterations,
                                    double tolerance, BackendPolicy backendPolicy) {
        if (data == null || data.length < 3 || model == null || maximumIterations < 1
                || !(tolerance > 0.0) || backendPolicy == null)
            throw new IllegalArgumentException("FIML inputs are invalid");
        int variables = model.variables().size();
        for (double[] row : data) if (row == null || row.length != variables) throw new IllegalArgumentException("FIML data must be rectangular");
        double[] means = initialMeans(data, variables), covariance = initialCovariance(data, means, variables);
        Set<String> patterns = new HashSet<>(); boolean converged = false; int iterations = 0;
        for (int iteration = 1; iteration <= maximumIterations; iteration++) {
            iterations = iteration; double[] first = new double[variables], second = new double[variables * variables]; patterns.clear();
            for (double[] row : data) {
                int[] observed = observed(row); StringBuilder pattern = new StringBuilder(); for (int index : observed) pattern.append(index).append(','); patterns.add(pattern.toString());
                double[] conditionalMean = means.clone(), conditionalCovariance = covariance.clone();
                int[] missing = missing(row);
                if (observed.length > 0 && missing.length > 0) {
                    double[] oo = submatrix(covariance, observed, observed, variables), om = submatrix(covariance, missing, observed, variables);
                    double[] inv = inverse(oo, observed.length); double[] difference = new double[observed.length];
                    for (int index = 0; index < observed.length; index++) difference[index] = row[observed[index]] - means[observed[index]];
                    double[] conditional = multiply(om, inverseTimes(inv, difference), missing.length, observed.length);
                    for (int index = 0; index < missing.length; index++) conditionalMean[missing[index]] = means[missing[index]] + conditional[index];
                    double[] mm = submatrix(covariance, missing, missing, variables);
                    double[] mo = transpose(om, missing.length, observed.length);
                    double[] conditionalVariance = subtract(mm, multiply(multiply(om, inv, missing.length, observed.length), mo, missing.length, missing.length), missing.length);
                    conditionalCovariance = new double[variables * variables]; for (int left = 0; left < missing.length; left++) for (int right = 0; right < missing.length; right++) conditionalCovariance[missing[left] * variables + missing[right]] = conditionalVariance[left * missing.length + right];
                }
                for (int variable = 0; variable < variables; variable++) { if (Double.isFinite(row[variable])) conditionalMean[variable] = row[variable]; first[variable] += conditionalMean[variable]; }
                for (int left = 0; left < variables; left++) for (int right = 0; right < variables; right++) second[left * variables + right] += conditionalMean[left] * conditionalMean[right] + conditionalCovariance[left * variables + right];
            }
            double[] nextMeans = first.clone(), nextCovariance = second.clone(); for (int variable = 0; variable < variables; variable++) nextMeans[variable] /= data.length;
            for (int left = 0; left < variables; left++) for (int right = 0; right < variables; right++) nextCovariance[left * variables + right] = nextCovariance[left * variables + right] / data.length - nextMeans[left] * nextMeans[right];
            for (int variable = 0; variable < variables; variable++) nextCovariance[variable * variables + variable] = Math.max(1e-8, nextCovariance[variable * variables + variable]);
            double change = 0.0; for (int index = 0; index < means.length; index++) change = Math.max(change, Math.abs(nextMeans[index] - means[index])); for (int index = 0; index < covariance.length; index++) change = Math.max(change, Math.abs(nextCovariance[index] - covariance[index]));
            means = nextMeans; covariance = nextCovariance; if (change <= tolerance * (1.0 + maxAbs(covariance))) { converged = true; break; }
        }
        SemFitResult fit = Sem.fitCovariance(covariance, data.length, model, SemOptions.defaults(), backendPolicy);
        return new SemFimlResult(fit, means, covariance, patterns.size(), iterations, converged);
    }

    private static double[] initialMeans(double[][] data, int variables) { double[] result = new double[variables]; int[] counts = new int[variables]; for (double[] row : data) for (int variable = 0; variable < variables; variable++) if (Double.isFinite(row[variable])) { result[variable] += row[variable]; counts[variable]++; } for (int variable = 0; variable < variables; variable++) result[variable] /= Math.max(1, counts[variable]); return result; }
    private static double[] initialCovariance(double[][] data, double[] means, int variables) { double[] result = new double[variables * variables]; int[] counts = new int[result.length]; for (double[] row : data) for (int left = 0; left < variables; left++) if (Double.isFinite(row[left])) for (int right = 0; right < variables; right++) if (Double.isFinite(row[right])) { result[left * variables + right] += (row[left] - means[left]) * (row[right] - means[right]); counts[left * variables + right]++; } for (int index = 0; index < result.length; index++) result[index] /= Math.max(1, counts[index]); for (int variable = 0; variable < variables; variable++) result[variable * variables + variable] = Math.max(1e-4, result[variable * variables + variable]); return result; }
    private static int[] observed(double[] row) { return java.util.stream.IntStream.range(0, row.length).filter(i -> Double.isFinite(row[i])).toArray(); }
    private static int[] missing(double[] row) { return java.util.stream.IntStream.range(0, row.length).filter(i -> !Double.isFinite(row[i])).toArray(); }
    private static double[] submatrix(double[] matrix, int[] rows, int[] columns, int dimension) { double[] result = new double[rows.length * columns.length]; for (int row = 0; row < rows.length; row++) for (int column = 0; column < columns.length; column++) result[row * columns.length + column] = matrix[rows[row] * dimension + columns[column]]; return result; }
    private static double[] inverseTimes(double[] matrix, double[] vector) { double[] result = new double[vector.length]; for (int row = 0; row < vector.length; row++) for (int column = 0; column < vector.length; column++) result[row] += matrix[row * vector.length + column] * vector[column]; return result; }
    private static double[] multiply(double[] left, double[] right, int rows, int columns) { int shared = right.length / columns; double[] result = new double[rows * columns]; for (int row = 0; row < rows; row++) for (int column = 0; column < columns; column++) for (int k = 0; k < shared; k++) result[row * columns + column] += left[row * shared + k] * right[k * columns + column]; return result; }
    private static double[] transpose(double[] matrix, int rows, int columns) { double[] result = new double[matrix.length]; for (int row = 0; row < rows; row++) for (int column = 0; column < columns; column++) result[column * rows + row] = matrix[row * columns + column]; return result; }
    private static double[] subtract(double[] left, double[] right, int dimension) { double[] result = left.clone(); for (int index = 0; index < result.length; index++) result[index] -= right[index]; return result; }
    private static double[] inverse(double[] matrix, int dimension) { double[] result = new double[dimension * dimension]; for (int index = 0; index < dimension; index++) result[index * dimension + index] = 1.0; double[] a = matrix.clone(); for (int column = 0; column < dimension; column++) { int pivot = column; for (int row = column + 1; row < dimension; row++) if (Math.abs(a[row * dimension + column]) > Math.abs(a[pivot * dimension + column])) pivot = row; if (Math.abs(a[pivot * dimension + column]) < 1e-12) throw new IllegalArgumentException("FIML observed covariance is singular"); swapRows(a, result, dimension, column, pivot); double scale = a[column * dimension + column]; for (int j = 0; j < dimension; j++) { a[column * dimension + j] /= scale; result[column * dimension + j] /= scale; } for (int row = 0; row < dimension; row++) if (row != column) { double factor = a[row * dimension + column]; for (int j = 0; j < dimension; j++) { a[row * dimension + j] -= factor * a[column * dimension + j]; result[row * dimension + j] -= factor * result[column * dimension + j]; } } } return result; }
    private static void swapRows(double[] left, double[] right, int dimension, int first, int second) { if (first == second) return; for (int column = 0; column < dimension; column++) { double value = left[first * dimension + column]; left[first * dimension + column] = left[second * dimension + column]; left[second * dimension + column] = value; value = right[first * dimension + column]; right[first * dimension + column] = right[second * dimension + column]; right[second * dimension + column] = value; } }
    private static double maxAbs(double[] values) { double result = 0.0; for (double value : values) result = Math.max(result, Math.abs(value)); return result; }
}
