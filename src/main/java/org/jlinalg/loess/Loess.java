/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.loess;

import java.util.Arrays;

/** Fast one-dimensional direct-surface local polynomial regression. */
public final class Loess {
    private static final double SINGULAR_TOLERANCE = 1e-12;
    private Loess() { }

    public static LoessResult fit(double[] x, double[] response) {
        return fit(x, response, null, LoessOptions.defaults());
    }

    public static LoessResult fit(double[] x, double[] response,
            LoessOptions options) {
        return fit(x, response, null, options);
    }

    public static LoessResult fit(double[] x, double[] response,
            double[] weights, LoessOptions options) {
        return prepare(x, options).fit(response, weights);
    }

    /** Sorts and validates the predictor once for repeated response fits. */
    public static Prepared prepare(double[] x) {
        return prepare(x, LoessOptions.defaults());
    }

    /** Sorts and validates the predictor once for repeated response fits. */
    public static Prepared prepare(double[] x, LoessOptions options) {
        return new Prepared(x, options);
    }

    public static final class Prepared {
        private final double[] x;
        private final int[] order;
        private final double[] sortedX;
        private final LoessOptions options;
        private final int neighborhoodSize;

        private Prepared(double[] x, LoessOptions options) {
            if (x == null || x.length < 2 || options == null)
                throw new IllegalArgumentException(
                    "LOESS requires at least two predictors and options");
            this.x = x.clone();
            Integer[] boxed = new Integer[x.length];
            for (int index = 0; index < x.length; index++) {
                if (!Double.isFinite(x[index]))
                    throw new IllegalArgumentException(
                        "LOESS predictors must be finite");
                boxed[index] = index;
            }
            Arrays.sort(boxed, (left, right) -> {
                int compared = Double.compare(x[left], x[right]);
                return compared != 0 ? compared : Integer.compare(left, right);
            });
            this.order = Arrays.stream(boxed).mapToInt(Integer::intValue).toArray();
            this.sortedX = new double[x.length];
            for (int position = 0; position < x.length; position++)
                sortedX[position] = x[order[position]];
            if (!(sortedX[x.length - 1] > sortedX[0]))
                throw new IllegalArgumentException(
                    "LOESS requires predictor variation");
            this.options = options;
            int requested = (int) Math.floor(options.span() * x.length + 1e-7);
            this.neighborhoodSize = Math.min(x.length,
                Math.max(options.degree() + 1, requested));
        }

        public LoessOptions options() { return options; }
        public int observations() { return x.length; }
        public int neighborhoodSize() { return neighborhoodSize; }

        public LoessResult fit(double[] response, double[] weights) {
            validateResponse(response);
            double responseScale = responseScale(response);
            double[] normalizedResponse = response.clone();
            for (int i = 0; i < response.length; i++) normalizedResponse[i] /= responseScale;
            double[] prior = weights(weights);
            double[] robust = new double[x.length];
            Arrays.fill(robust, 1.0);
            int robustUpdates = options.family() == LoessFamily.SYMMETRIC
                ? Math.max(0, options.robustnessIterations() - 1) : 0;
            Fit fit = null;
            for (int iteration = 0; iteration <= robustUpdates; iteration++) {
                fit = fitted(normalizedResponse, prior, robust, true);
                if (iteration < robustUpdates)
                    robust = robustness(fit.residuals());
            }
            for (int i = 0; i < response.length; i++) {
                fit.fitted()[i] *= responseScale;
                fit.residuals()[i] = response[i] - fit.fitted()[i];
            }
            return new LoessResult(this, response, prior, robust,
                fit.fitted(), fit.residuals(), fit.leverage(),
                Arrays.stream(fit.leverage()).sum(), robustUpdates + 1);
        }

        double[] predict(double[] response, double[] prior,
                double[] robust, double[] queries) {
            if (queries == null)
                throw new IllegalArgumentException(
                    "LOESS prediction points are required");
            double[] result = new double[queries.length];
            double responseScale = responseScale(response);
            double[] normalizedResponse = response.clone();
            for (int i = 0; i < response.length; i++) normalizedResponse[i] /= responseScale;
            for (int index = 0; index < queries.length; index++) {
                if (!Double.isFinite(queries[index]))
                    throw new IllegalArgumentException(
                        "LOESS prediction points must be finite");
                result[index] = local(normalizedResponse, prior, robust,
                    queries[index], -1).value() * responseScale;
            }
            return result;
        }

        private Fit fitted(double[] response, double[] prior,
                double[] robust, boolean leverage) {
            double[] fitted = new double[x.length];
            double[] hat = new double[x.length];
            for (int row = 0; row < x.length; row++) {
                Local local = local(response, prior, robust, x[row],
                    leverage ? row : -1);
                fitted[row] = local.value();
                hat[row] = local.leverage();
            }
            double[] residuals = new double[x.length];
            for (int row = 0; row < x.length; row++)
                residuals[row] = response[row] - fitted[row];
            return new Fit(fitted, residuals, hat);
        }

        private Local local(double[] response, double[] prior,
                double[] robust, double query, int targetRow) {
            Window window = window(query);
            double radius = Math.max(
                Math.abs(sortedX[window.first()] - query),
                Math.abs(sortedX[window.last() - 1] - query));
            double coordinateScale = 1.0;
            if (!Double.isFinite(radius)) {
                coordinateScale = Math.max(Math.abs(query),
                    Math.max(Math.abs(sortedX[window.first()]), Math.abs(sortedX[window.last() - 1])));
                radius = Math.max(Math.abs(sortedX[window.first()] / coordinateScale - query / coordinateScale),
                    Math.abs(sortedX[window.last() - 1] / coordinateScale - query / coordinateScale));
            }
            if (!(radius > 0.0)) radius = Math.ulp(Math.abs(query) + 1.0);
            int maximumDegree = options.degree();
            double[] normal = new double[5];
            double[] rhs = new double[3];
            double targetWeight = 0.0;
            for (int position = window.first(); position < window.last(); position++) {
                int row = order[position];
                double scaled = coordinateScale == 1.0 ? (x[row] - query) / radius
                    : (x[row] / coordinateScale - query / coordinateScale) / radius;
                double distance = Math.abs(scaled);
                double kernel = distance >= 1.0 ? 0.0
                    : cube(1.0 - cube(distance));
                double weight = prior[row] * robust[row] * kernel;
                if (!(weight > 0.0)) continue;
                double power = 1.0;
                for (int exponent = 0; exponent <= 2 * maximumDegree; exponent++) {
                    normal[exponent] += weight * power;
                    if (exponent <= maximumDegree)
                        rhs[exponent] += weight * power * response[row];
                    power *= scaled;
                }
                if (row == targetRow) targetWeight = weight;
            }
            // Normalize each local system; common weight units must not alter
            // rank detection or overflow the determinant used for leverage.
            double mass = normal[0];
            if (mass > 0) {
                for (int i = 0; i < normal.length; i++) normal[i] /= mass;
                for (int i = 0; i < rhs.length; i++) rhs[i] /= mass;
                targetWeight /= mass;
            }
            for (int degree = maximumDegree; degree >= 0; degree--) {
                Solve solve = solve(normal, rhs, degree);
                if (solve != null)
                    return new Local(solve.coefficients()[0],
                        targetWeight * solve.inverse00());
            }
            throw new IllegalArgumentException(
                "LOESS neighborhood has no positive-weight observations");
        }

        private Window window(double query) {
            int lower = 0, upper = sortedX.length - neighborhoodSize;
            while (lower < upper) {
                int middle = lower + (upper - lower) / 2;
                double left = query - sortedX[middle];
                double right = sortedX[middle + neighborhoodSize] - query;
                if (!Double.isFinite(left) || !Double.isFinite(right)) {
                    double scale = Math.max(Math.abs(query), Math.max(Math.abs(sortedX[middle]),
                        Math.abs(sortedX[middle + neighborhoodSize])));
                    left = query / scale - sortedX[middle] / scale;
                    right = sortedX[middle + neighborhoodSize] / scale - query / scale;
                }
                if (left > right) lower = middle + 1;
                else upper = middle;
            }
            return new Window(lower, lower + neighborhoodSize);
        }

        private void validateResponse(double[] response) {
            if (response == null || response.length != x.length)
                throw new IllegalArgumentException(
                    "LOESS response must match predictors");
            for (double value : response) if (!Double.isFinite(value))
                throw new IllegalArgumentException(
                    "LOESS response must be finite");
        }

        private double[] weights(double[] supplied) {
            if (supplied == null) {
                double[] result = new double[x.length];
                Arrays.fill(result, 1.0);
                return result;
            }
            if (supplied.length != x.length)
                throw new IllegalArgumentException(
                    "LOESS weights must match predictors");
            double[] result = supplied.clone();
            double positive = 0.0;
            for (double value : result) {
                if (!(value >= 0.0) || !Double.isFinite(value))
                    throw new IllegalArgumentException(
                        "LOESS weights must be finite and nonnegative");
                positive = Math.max(positive, value);
            }
            if (!(positive > 0.0))
                throw new IllegalArgumentException(
                    "LOESS requires a positive weight");
            for (int i = 0; i < result.length; i++) result[i] /= positive;
            return result;
        }

        private static double responseScale(double[] response) {
            double scale = 0.0;
            for (double value : response) scale = Math.max(scale, Math.abs(value));
            return scale > 0.0 ? scale : 1.0;
        }
    }

    private static Solve solve(double[] moments, double[] rhs, int degree) {
        int dimension = degree + 1;
        double[][] augmented = new double[dimension][dimension + 1];
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++)
                augmented[row][column] = moments[row + column];
            augmented[row][dimension] = rhs[row];
        }
        double scale = Math.max(1.0, moments[0]);
        for (int pivot = 0; pivot < dimension; pivot++) {
            int best = pivot;
            for (int row = pivot + 1; row < dimension; row++)
                if (Math.abs(augmented[row][pivot])
                        > Math.abs(augmented[best][pivot])) best = row;
            if (Math.abs(augmented[best][pivot]) <= SINGULAR_TOLERANCE * scale)
                return null;
            double[] swap = augmented[pivot];
            augmented[pivot] = augmented[best];
            augmented[best] = swap;
            double divisor = augmented[pivot][pivot];
            for (int column = pivot; column <= dimension; column++)
                augmented[pivot][column] /= divisor;
            for (int row = 0; row < dimension; row++) {
                if (row == pivot) continue;
                double multiplier = augmented[row][pivot];
                for (int column = pivot; column <= dimension; column++)
                    augmented[row][column] -= multiplier
                        * augmented[pivot][column];
            }
        }
        double[] coefficients = new double[dimension];
        for (int row = 0; row < dimension; row++)
            coefficients[row] = augmented[row][dimension];
        double inverse00 = inverse00(moments, degree);
        return new Solve(coefficients, inverse00);
    }

    private static double inverse00(double[] moments, int degree) {
        if (degree == 0) return 1.0 / moments[0];
        if (degree == 1) {
            double determinant = moments[0] * moments[2]
                - moments[1] * moments[1];
            return moments[2] / determinant;
        }
        double a = moments[0], b = moments[1], c = moments[2];
        double d = moments[3], e = moments[4];
        double determinant = a * (c * e - d * d)
            - b * (b * e - c * d) + c * (b * d - c * c);
        return (c * e - d * d) / determinant;
    }

    private static double[] robustness(double[] residuals) {
        double[] absolute = Arrays.stream(residuals).map(Math::abs).sorted()
            .toArray();
        int middle = absolute.length / 2;
        double median = absolute.length % 2 == 0
            ? 0.5 * (absolute[middle - 1] + absolute[middle])
            : absolute[middle];
        double cutoff = 6.0 * median;
        double[] result = new double[residuals.length];
        if (!(cutoff > 0.0)) {
            Arrays.fill(result, 1.0);
            return result;
        }
        for (int row = 0; row < residuals.length; row++) {
            double ratio = Math.abs(residuals[row]) / cutoff;
            result[row] = ratio >= 1.0 ? 0.0
                : square(1.0 - square(ratio));
        }
        return result;
    }

    private static double square(double value) { return value * value; }
    private static double cube(double value) { return value * value * value; }
    private record Window(int first, int last) { }
    private record Local(double value, double leverage) { }
    private record Solve(double[] coefficients, double inverse00) { }
    private record Fit(double[] fitted, double[] residuals, double[] leverage) { }
}
