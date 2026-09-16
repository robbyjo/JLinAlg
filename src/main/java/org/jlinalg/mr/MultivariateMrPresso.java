/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/**
 * Covariance-aware multivariate MR-PRESSO using leave-one-instrument-out
 * Mahalanobis residuals and a reproducible parametric bootstrap.
 */
public final class MultivariateMrPresso {
    private MultivariateMrPresso() { }

    public static MultivariateMrPressoResult analyze(
            List<MultivariateInstrument> instruments,
            List<String> exposureNames, List<String> outcomeNames,
            double[][] outcomeCorrelation, double significanceThreshold,
            int simulations, long seed, BackendPolicy backendPolicy) {
        if (!(significanceThreshold > 0 && significanceThreshold < 1)
                || simulations < 20)
            throw new IllegalArgumentException(
                "outlier threshold must be in (0,1) and simulations at least 20");
        int exposures = exposureNames.size(), outcomes = outcomeNames.size();
        if (instruments.size() <= exposures + 1)
            throw new IllegalArgumentException(
                "multivariate MR-PRESSO requires at least two residual instruments");
        MultivariateMrResult raw = MultivariateMendelianRandomization.fit(
            instruments, exposureNames, outcomeNames, outcomeCorrelation,
            false, backendPolicy);
        double[] correlation = MatrixOps.rowMajor(outcomeCorrelation, outcomes);
        LooDiagnostics observedDiagnostics = looDiagnostics(instruments,
            exposureNames, outcomeNames, outcomeCorrelation, correlation,
            backendPolicy);
        double[] observed = observedDiagnostics.scores();
        double observedGlobal = sum(observed);
        int globalExceedances = 0;
        int[] instrumentExceedances = new int[instruments.size()];
        Random random = new Random(seed);
        double[][] correlationCholesky = cholesky(outcomeCorrelation);
        for (int simulation = 0; simulation < simulations; simulation++) {
            List<MultivariateInstrument> generated = simulate(instruments,
                observedDiagnostics.coefficients(), exposures, outcomes,
                correlationCholesky, random);
            double[] score = looDiagnostics(generated, exposureNames,
                outcomeNames, outcomeCorrelation, correlation, backendPolicy)
                .scores();
            if (sum(score) >= observedGlobal) globalExceedances++;
            for (int instrument = 0; instrument < score.length; instrument++)
                if (score[instrument] >= observed[instrument])
                    instrumentExceedances[instrument]++;
        }
        List<MultivariateMrOutlier> tests = new ArrayList<>();
        List<String> outlierIds = new ArrayList<>();
        List<MultivariateInstrument> retained = new ArrayList<>();
        for (int instrument = 0; instrument < instruments.size(); instrument++) {
            // Include the observed statistic in the Monte Carlo reference set.
            // This prevents impossible zero p-values and gives the usual valid
            // finite-simulation calibration.
            double p = (instrumentExceedances[instrument] + 1.0)
                / (simulations + 1.0);
            double adjusted = Math.min(1, p * instruments.size());
            boolean outlier = adjusted <= significanceThreshold;
            String id = instruments.get(instrument).variantId();
            tests.add(new MultivariateMrOutlier(id, observed[instrument], p,
                adjusted, outlier));
            if (outlier) outlierIds.add(id); else retained.add(instruments.get(instrument));
        }
        MultivariateMrResult corrected = null;
        if (!outlierIds.isEmpty() && retained.size() > exposures)
            corrected = MultivariateMendelianRandomization.fit(retained,
                exposureNames, outcomeNames, outcomeCorrelation, false,
                backendPolicy);
        double[] distortion = new double[raw.beta().length];
        java.util.Arrays.fill(distortion, Double.NaN);
        if (corrected != null) {
            double[] before = raw.beta(), after = corrected.beta();
            for (int index = 0; index < distortion.length; index++)
                distortion[index] = after[index] == 0
                    ? Double.NaN : 100 * (before[index] - after[index])
                        / Math.abs(after[index]);
        }
        return new MultivariateMrPressoResult(raw, corrected, observedGlobal,
            (globalExceedances + 1.0) / (simulations + 1.0), tests, outlierIds,
            distortion, simulations, seed);
    }

    private static LooDiagnostics looDiagnostics(
            List<MultivariateInstrument> instruments,
            List<String> exposureNames, List<String> outcomeNames,
            double[][] outcomeCorrelation, double[] correlation,
            BackendPolicy backendPolicy) {
        int exposures = exposureNames.size(), outcomes = outcomeNames.size();
        double[] result = new double[instruments.size()];
        double[][] coefficients = new double[instruments.size()][];
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            for (int omitted = 0; omitted < instruments.size(); omitted++) {
                List<MultivariateInstrument> training = new ArrayList<>(
                    instruments.size() - 1);
                for (int index = 0; index < instruments.size(); index++)
                    if (index != omitted) training.add(instruments.get(index));
                double[] beta = MultivariateMendelianRandomization.fit(training,
                    exposureNames, outcomeNames, outcomeCorrelation, false,
                    backendPolicy).beta();
                coefficients[omitted] = beta;
                MultivariateInstrument heldOut = instruments.get(omitted);
                double[] residual = heldOut.outcomeEffects();
                double[] x = heldOut.exposureEffects();
                for (int outcome = 0; outcome < outcomes; outcome++)
                    for (int exposure = 0; exposure < exposures; exposure++)
                        residual[outcome] -= x[exposure]
                            * beta[outcome * exposures + exposure];
                double[] inverse = MultivariateMendelianRandomization
                    .inverseOutcomeCovariance(heldOut, correlation, outcomes,
                        backend);
                double[] solved = MatrixOps.multiply(backend, inverse,
                    outcomes, outcomes, residual);
                result[omitted] = Math.max(0, backend.ddot(outcomes,
                    residual, 0, 1, solved, 0, 1));
            }
        }
        return new LooDiagnostics(result, coefficients);
    }

    private static List<MultivariateInstrument> simulate(
            List<MultivariateInstrument> observed, double[][] looCoefficients,
            int exposures, int outcomes, double[][] correlationCholesky,
            Random random) {
        List<MultivariateInstrument> result = new ArrayList<>(observed.size());
        for (int instrumentIndex = 0; instrumentIndex < observed.size();
                instrumentIndex++) {
            MultivariateInstrument instrument = observed.get(instrumentIndex);
            double[] observedX = instrument.exposureEffects();
            double[] generatedX = observedX.clone();
            double[] sx = instrument.exposureStandardErrors();
            for (int exposure = 0; exposure < exposures; exposure++)
                generatedX[exposure] += sx[exposure] * random.nextGaussian();
            double[] independent = new double[outcomes];
            for (int outcome = 0; outcome < outcomes; outcome++)
                independent[outcome] = random.nextGaussian();
            double[] y = new double[outcomes];
            double[] sy = instrument.outcomeStandardErrors();
            for (int outcome = 0; outcome < outcomes; outcome++) {
                for (int exposure = 0; exposure < exposures; exposure++)
                    y[outcome] += observedX[exposure]
                        * looCoefficients[instrumentIndex][
                            outcome * exposures + exposure];
                double error = 0;
                for (int column = 0; column <= outcome; column++)
                    error += correlationCholesky[outcome][column]
                        * independent[column];
                y[outcome] += sy[outcome] * error;
            }
            result.add(new MultivariateInstrument(instrument.variantId(),
                generatedX, sx, y, sy));
        }
        return result;
    }

    private static double[][] cholesky(double[][] correlation) {
        int size = correlation.length;
        double[][] lower = new double[size][size];
        for (int row = 0; row < size; row++) {
            for (int column = 0; column <= row; column++) {
                double value = correlation[row][column];
                for (int k = 0; k < column; k++)
                    value -= lower[row][k] * lower[column][k];
                if (row == column) {
                    if (!(value > 0) || !Double.isFinite(value))
                        throw new IllegalArgumentException(
                            "outcome correlation must be positive definite");
                    lower[row][column] = Math.sqrt(value);
                } else {
                    lower[row][column] = value / lower[column][column];
                }
            }
        }
        return lower;
    }

    private static double sum(double[] values) {
        double result = 0;
        for (double value : values) result += value;
        return result;
    }

    private record LooDiagnostics(double[] scores, double[][] coefficients) { }
}
