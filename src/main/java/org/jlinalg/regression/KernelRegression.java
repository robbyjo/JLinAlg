/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** One-dimensional Nadaraya-Watson Gaussian-kernel regression. */
public final class KernelRegression {
    private KernelRegression() { }

    public static Result fit(double[] predictor, double[] response, double bandwidth) {
        validate(predictor, response, bandwidth);
        double[] fitted = predict(predictor, response, predictor, bandwidth);
        return new Result(predictor.clone(), response.clone(), fitted, bandwidth);
    }

    public static Result fit(double[] predictor, double[] response) {
        validate(predictor, response, 1.0);
        double mean = 0.0; for (double value : predictor) mean += value;
        mean /= predictor.length; double scale = 0.0;
        for (double value : predictor) scale += (value - mean) * (value - mean);
        double bandwidth = 1.06 * Math.sqrt(scale / Math.max(1, predictor.length - 1))
            * Math.pow(predictor.length, -0.2);
        return fit(predictor, response, Math.max(bandwidth, 1e-8));
    }

    public static double[] predict(double[] predictor, double[] response, double[] query, double bandwidth) {
        validate(predictor, response, bandwidth);
        if (query == null) throw new IllegalArgumentException("query is required");
        double[] result = new double[query.length];
        for (int point = 0; point < query.length; point++) {
            double numerator = 0.0, denominator = 0.0;
            for (int row = 0; row < predictor.length; row++) {
                double distance = (query[point] - predictor[row]) / bandwidth;
                double weight = Math.exp(-0.5 * distance * distance);
                numerator += weight * response[row]; denominator += weight;
            }
            result[point] = numerator / denominator;
        }
        return result;
    }

    private static void validate(double[] predictor, double[] response, double bandwidth) {
        if (predictor == null || response == null || predictor.length != response.length || predictor.length < 2
                || !(bandwidth > 0) || !Double.isFinite(bandwidth)) throw new IllegalArgumentException("kernel inputs are invalid");
        for (double value : predictor) if (!Double.isFinite(value)) throw new IllegalArgumentException("predictors must be finite");
        for (double value : response) if (!Double.isFinite(value)) throw new IllegalArgumentException("response must be finite");
    }

    /** Stored training data, fitted values, and bandwidth. */
    public record Result(double[] predictor, double[] response, double[] fittedValues, double bandwidth) {
        public Result { predictor = predictor.clone(); response = response.clone(); fittedValues = fittedValues.clone(); }
        public double[] predictor() { return predictor.clone(); }
        public double[] response() { return response.clone(); }
        public double[] fittedValues() { return fittedValues.clone(); }
        public double[] predict(double[] query) { return KernelRegression.predict(predictor, response, query, bandwidth); }
    }
}
