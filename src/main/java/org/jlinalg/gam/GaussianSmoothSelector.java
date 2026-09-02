/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Direct fixed, GCV, UBRE, and AIC fits for arbitrary quadratic smooths. */
public final class GaussianSmoothSelector {
    private GaussianSmoothSelector() { }

    /** Selects every smoothing parameter with default GCV controls. */
    public static GaussianSmoothSelectionResult fit(
            double[] response,
            double[][] parametricDesign,
            List<QuadraticSmoothTerm> smoothTerms) {
        return fit(response, parametricDesign, smoothTerms, null,
            SmoothingSelectionOptions.gcv(), BackendPolicy.PREFERRED);
    }

    /** Selects all marginal penalties by cyclic log-scale coordinate search. */
    public static GaussianSmoothSelectionResult fit(
            double[] response,
            double[][] parametricDesign,
            List<QuadraticSmoothTerm> smoothTerms,
            List<double[]> initialSmoothingParameters,
            SmoothingSelectionOptions options,
            BackendPolicy backendPolicy) {
        validate(response, parametricDesign, smoothTerms, options, backendPolicy);
        Layout layout = layout(smoothTerms);
        double[] logSmoothing = initial(initialSmoothingParameters, layout);
        int evaluations = 0;
        Candidate best = evaluate(response, parametricDesign, smoothTerms,
            expand(logSmoothing, layout), options, backendPolicy);
        evaluations++;
        double step = options.initialLogStep();
        for (int sweep = 0; sweep < options.maximumSweeps(); sweep++) {
            boolean improved = false;
            for (int parameter = 0; parameter < logSmoothing.length; parameter++) {
                double original = logSmoothing[parameter];
                Candidate coordinateBest = best;
                double selected = original;
                for (double direction : new double[] {-1.0, 1.0}) {
                    double trial = Math.max(options.minimumLogSmoothing(),
                        Math.min(options.maximumLogSmoothing(),
                            original + direction * step));
                    if (trial == original) continue;
                    logSmoothing[parameter] = trial;
                    Candidate candidate = evaluate(response, parametricDesign,
                        smoothTerms, expand(logSmoothing, layout), options,
                        backendPolicy);
                    evaluations++;
                    if (candidate.score() < coordinateBest.score()) {
                        coordinateBest = candidate;
                        selected = trial;
                    }
                }
                logSmoothing[parameter] = selected;
                if (coordinateBest.score() < best.score()
                        - options.tolerance() * (1.0 + Math.abs(best.score()))) {
                    best = coordinateBest;
                    improved = true;
                }
            }
            if (!improved) step *= 0.5;
            if (step <= options.tolerance()) break;
        }
        List<double[]> smoothing = expand(logSmoothing, layout);
        // Re-evaluate so the returned fit always corresponds to the final coordinate vector.
        best = evaluate(response, parametricDesign, smoothTerms,
            smoothing, options, backendPolicy);
        evaluations++;
        return best.result(smoothing, evaluations);
    }

    /** Fits once with explicit smoothing parameters and no outer optimization. */
    public static GaussianSmoothSelectionResult fitFixed(
            double[] response,
            double[][] parametricDesign,
            List<QuadraticSmoothTerm> smoothTerms,
            List<double[]> smoothingParameters,
            BackendPolicy backendPolicy) {
        SmoothingSelectionOptions options = SmoothingSelectionOptions.gcv();
        validate(response, parametricDesign, smoothTerms, options, backendPolicy);
        Candidate candidate = evaluate(response, parametricDesign, smoothTerms,
            smoothingParameters, options, backendPolicy);
        return candidate.result(smoothingParameters, 1);
    }

    private static Candidate evaluate(
            double[] response,
            double[][] parametricDesign,
            List<QuadraticSmoothTerm> smoothTerms,
            List<double[]> smoothingParameters,
            SmoothingSelectionOptions options,
            BackendPolicy backendPolicy) {
        PenalizedPredictor predictor = QuadraticPenalizedPredictor.compile(
            parametricDesign, smoothTerms, smoothingParameters, backendPolicy);
        double[] design = predictor.designView();
        double[] penalty = predictor.penaltyDiagonalView();
        int rows = response.length;
        int columns = predictor.columns();
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] cross = MatrixOps.transposeMultiply(
                backend, design, rows, columns, design, columns);
            for (int column = 0; column < columns; column++) {
                cross[column * columns + column] += penalty[column];
            }
            double[] responseMatrix = response.clone();
            double[] right = MatrixOps.transposeMultiply(
                backend, design, rows, columns, responseMatrix, 1);
            CholeskyFactor factor = backend.dpotrf(cross, columns);
            double[] coefficients = factor.solve(right);
            double[] fitted = MatrixOps.multiply(
                backend, design, rows, columns, coefficients);
            double[] residuals = MatrixOps.subtract(response, fitted);
            double rss = dot(residuals, residuals);
            double[] inverse = factor.solve(MatrixOps.identity(columns), columns);
            double penaltyTrace = 0.0;
            for (int column = 0; column < columns; column++) {
                penaltyTrace += inverse[column * columns + column]
                    * penalty[column];
            }
            double edf = Math.max(0.0, Math.min(columns,
                columns - penaltyTrace));
            double denominator = Math.max(1e-12, rows - edf);
            double scale = rss / denominator;
            double score = switch (options.criterion()) {
                case GCV -> rows * rss / (denominator * denominator);
                case UBRE -> rss / rows - options.knownScale()
                    + 2.0 * options.knownScale() * edf / rows;
                case AIC -> rows * Math.log(Math.max(rss / rows, 1e-300))
                    + 2.0 * edf;
            };
            double[] covariance = inverse.clone();
            for (int index = 0; index < covariance.length; index++) {
                covariance[index] *= scale;
            }
            return new Candidate(predictor, coefficients, fitted, residuals,
                scale, covariance, edf, score);
        }
    }

    private static void validate(
            double[] response,
            double[][] parametricDesign,
            List<QuadraticSmoothTerm> smoothTerms,
            SmoothingSelectionOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || response.length < 2 || parametricDesign == null
                || parametricDesign.length != response.length
                || smoothTerms == null || smoothTerms.isEmpty()
                || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "response, parametric design, smooths, controls, and backend are required");
        }
        MatrixOps.requireFinite(response, "response");
    }

    private static Layout layout(List<QuadraticSmoothTerm> smoothTerms) {
        int[] starts = new int[smoothTerms.size() + 1];
        for (int term = 0; term < smoothTerms.size(); term++) {
            starts[term + 1] = starts[term] + smoothTerms.get(term).penaltyCount();
        }
        return new Layout(starts);
    }

    private static double[] initial(List<double[]> supplied, Layout layout) {
        double[] result = new double[layout.starts()[layout.starts().length - 1]];
        if (supplied == null) return result;
        if (supplied.size() != layout.starts().length - 1) {
            throw new IllegalArgumentException("one smoothing vector is required per term");
        }
        for (int term = 0; term < supplied.size(); term++) {
            double[] values = supplied.get(term);
            int count = layout.starts()[term + 1] - layout.starts()[term];
            if (values == null || values.length != count) {
                throw new IllegalArgumentException(
                    "one positive smoothing parameter is required per penalty");
            }
            for (int penalty = 0; penalty < count; penalty++) {
                double value = values[penalty];
                if (!(value > 0.0) || !Double.isFinite(value)) {
                    throw new IllegalArgumentException(
                        "smoothing parameters must be finite and positive");
                }
                result[layout.starts()[term] + penalty] = Math.log(value);
            }
        }
        return result;
    }

    private static List<double[]> expand(double[] logarithms, Layout layout) {
        List<double[]> result = new ArrayList<>(layout.starts().length - 1);
        for (int term = 0; term < layout.starts().length - 1; term++) {
            double[] values = new double[
                layout.starts()[term + 1] - layout.starts()[term]];
            for (int penalty = 0; penalty < values.length; penalty++) {
                values[penalty] = Math.exp(
                    logarithms[layout.starts()[term] + penalty]);
            }
            result.add(values);
        }
        return List.copyOf(result);
    }

    private static double dot(double[] first, double[] second) {
        double result = 0.0;
        for (int index = 0; index < first.length; index++) {
            result += first[index] * second[index];
        }
        return result;
    }

    private record Layout(int[] starts) {
        private Layout { starts = Arrays.copyOf(starts, starts.length); }
        @Override public int[] starts() { return starts.clone(); }
    }

    private record Candidate(
            PenalizedPredictor predictor,
            double[] coefficients,
            double[] fitted,
            double[] residuals,
            double scale,
            double[] covariance,
            double edf,
            double score) {
        private GaussianSmoothSelectionResult result(
                List<double[]> smoothing, int evaluations) {
            return new GaussianSmoothSelectionResult(predictor, smoothing,
                coefficients, fitted, residuals, scale, covariance,
                edf, score, evaluations);
        }
    }
}
