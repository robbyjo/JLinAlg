/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Deterministic span-selection smoother inspired by R's stats::supsmu. */
public final class SuperSmoother {
    private SuperSmoother() { }

    public static Result fit(double[] predictor, double[] response) {
        return fit(predictor, response, new double[] {0.05, 0.20, 0.50});
    }

    public static Result fit(double[] predictor, double[] response, double[] spans) {
        if (predictor == null || response == null || predictor.length != response.length || predictor.length < 4
                || spans == null || spans.length == 0) throw new IllegalArgumentException("supersmoother inputs are invalid");
        for (double value : predictor) if (!Double.isFinite(value)) throw new IllegalArgumentException("predictors must be finite");
        for (double value : response) if (!Double.isFinite(value)) throw new IllegalArgumentException("response must be finite");
        double[] fitted = new double[predictor.length], selected = new double[predictor.length];
        for (int row = 0; row < predictor.length; row++) {
            double bestError = Double.POSITIVE_INFINITY, bestFit = response[row], bestSpan = spans[0];
            for (double span : spans) {
                if (!(span > 0) || !Double.isFinite(span)) throw new IllegalArgumentException("spans must be positive");
                double prediction = localLinear(predictor, response, predictor[row], span, row);
                double error = response[row] - prediction;
                if (error * error < bestError) { bestError = error * error; bestFit = prediction; bestSpan = span; }
            }
            fitted[row] = bestFit; selected[row] = bestSpan;
        }
        return new Result(predictor.clone(), response.clone(), fitted, selected);
    }

    private static double localLinear(double[] x, double[] y, double query, double span, int excluded) {
        int neighbors = Math.max(3, (int) Math.ceil(span * x.length));
        double[] distances = new double[x.length]; for (int i = 0; i < x.length; i++) distances[i] = Math.abs(x[i] - query);
        double radius = kth(distances, Math.min(x.length - 1, neighbors - 1));
        if (!(radius > 0)) radius = 1.0;
        double sw = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
        for (int i = 0; i < x.length; i++) if (i != excluded) {
            double u = Math.abs(x[i] - query) / radius; if (u >= 1.0) continue;
            double weight = Math.pow(1.0 - u * u * u, 3);
            sw += weight; sx += weight * x[i]; sy += weight * y[i]; sxx += weight * x[i] * x[i]; sxy += weight * x[i] * y[i];
        }
        double determinant = sw * sxx - sx * sx;
        if (Math.abs(determinant) < 1e-14) return sw == 0.0 ? y[excluded] : sy / sw;
        return (sy * sxx - sxy * sx + query * (sw * sxy - sx * sy)) / determinant;
    }
    private static double kth(double[] values, int index) { double[] copy = values.clone(); java.util.Arrays.sort(copy); return copy[index]; }

    /** Smoother output and per-observation selected spans. */
    public record Result(double[] predictor, double[] response, double[] fittedValues, double[] selectedSpans) {
        public Result { predictor = predictor.clone(); response = response.clone(); fittedValues = fittedValues.clone(); selectedSpans = selectedSpans.clone(); }
        public double[] predictor() { return predictor.clone(); }
        public double[] response() { return response.clone(); }
        public double[] fittedValues() { return fittedValues.clone(); }
        public double[] selectedSpans() { return selectedSpans.clone(); }
        public double[] predict(double[] query) {
            if (query == null) throw new IllegalArgumentException("query is required");
            double[] result = new double[query.length];
            for (int i = 0; i < query.length; i++) result[i] = localLinear(predictor, response, query[i], selectedSpans[Math.min(i, selectedSpans.length - 1)], -1);
            return result;
        }
    }
}
