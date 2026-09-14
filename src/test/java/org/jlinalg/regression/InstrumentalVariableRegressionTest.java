/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.inference.StatisticDistribution;
import org.jlinalg.ols.Ols;
import org.junit.jupiter.api.Test;

class InstrumentalVariableRegressionTest {
    @Test
    void fitsEndogenousRegressionAndDoesNotUseNaiveSecondStageResiduals() {
        Data data = data(240);
        InstrumentalVariableResult fit = InstrumentalVariableRegression.fit(
            data.response(), data.exogenous(), data.endogenous(),
            data.instruments(), new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.HOMOSKEDASTIC, 0.95),
            BackendPolicy.CPU);

        assertEquals(2.0, fit.coefficients()[2], 0.04);
        assertEquals(3, fit.parameters());
        assertEquals(2, fit.exogenousParameters());
        assertEquals(1, fit.endogenousParameters());
        assertEquals(2, fit.excludedInstruments());
        assertEquals(1, fit.overidentificationDegreesOfFreedom());
        assertEquals(237, fit.residualDegreesOfFreedom());
        assertEquals(DegreesOfFreedomMethod.RESIDUAL,
            fit.associationStatistics().degreesOfFreedomMethod());
        assertTrue(fit.minimumFirstStageFStatistic() > 100.0);
        assertTrue(fit.strengthDiagnostics().get(0).partialRSquared() > 0.7);

        double[][] naiveSecondStage = new double[data.response().length][3];
        double[][] instrumented = fit.instrumentedEndogenous();
        for (int row = 0; row < naiveSecondStage.length; row++) {
            naiveSecondStage[row][0] = data.exogenous()[row][0];
            naiveSecondStage[row][1] = data.exogenous()[row][1];
            naiveSecondStage[row][2] = instrumented[row][0];
        }
        double naiveVariance = Ols.fit(data.response(), naiveSecondStage)
            .residualVariance();
        assertNotEquals(naiveVariance, fit.residualVariance(), 1e-3);
        assertEquals(sumSquares(fit.residuals())
            / fit.residualDegreesOfFreedom(), fit.residualVariance(), 1e-14);
    }

    @Test
    void supportsHeteroskedasticAndClusterRobustInference() {
        Data data = data(240);
        InstrumentalVariableResult hc1 = InstrumentalVariableRegression.fit(
            data.response(), data.exogenous(), data.endogenous(),
            data.instruments(), new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.HC1, 0.90),
            BackendPolicy.CPU);
        assertEquals(StatisticDistribution.STANDARD_NORMAL,
            hc1.associationStatistics().statisticDistribution());
        assertEquals(DegreesOfFreedomMethod.ASYMPTOTIC,
            hc1.associationStatistics().degreesOfFreedomMethod());
        assertEquals(StatisticDistribution.CHI_SQUARE,
            hc1.strengthDiagnostics().get(0).referenceDistribution());
        assertTrue(Double.isInfinite(
            hc1.strengthDiagnostics().get(0).denominatorDegreesOfFreedom()));

        int[] clusters = new int[240];
        for (int row = 0; row < clusters.length; row++) clusters[row] = row / 6;
        InstrumentalVariableResult clustered = InstrumentalVariableRegression.fit(
            data.response(), data.exogenous(), data.endogenous(),
            data.instruments(), clusters, new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.CLUSTER_CR1, 0.95),
            BackendPolicy.CPU);
        assertEquals(40, clustered.clusters());
        assertEquals(37.0, clustered.inferenceDegreesOfFreedom());
        assertEquals(DegreesOfFreedomMethod.CLUSTER,
            clustered.associationStatistics().degreesOfFreedomMethod());
        assertEquals(StatisticDistribution.F,
            clustered.strengthDiagnostics().get(0).referenceDistribution());
        assertEquals(36.0, clustered.strengthDiagnostics().get(0)
            .denominatorDegreesOfFreedom());
        assertEquals(1, clustered.testContrast(
            new double[][] {{0.0, 0.0, 1.0}})
            .numeratorDegreesOfFreedom());
    }

    @Test
    void projectionsAreInvariantToExtremeInstrumentUnits() {
        Data data = data(180);
        double[][] scaled = new double[180][2];
        for (int row = 0; row < scaled.length; row++) {
            scaled[row][0] = data.instruments()[row][0] * 1e-10;
            scaled[row][1] = data.instruments()[row][1] * 1e10;
        }
        InstrumentalVariableResult baseline = InstrumentalVariableRegression.fit(
            data.response(), data.exogenous(), data.endogenous(),
            data.instruments(), new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.HC1, 0.95),
            BackendPolicy.CPU);
        InstrumentalVariableResult changedUnits = InstrumentalVariableRegression.fit(
            data.response(), data.exogenous(), data.endogenous(), scaled,
            new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.HC1, 0.95),
            BackendPolicy.CPU);

        assertArrayEquals(baseline.coefficients(),
            changedUnits.coefficients(), 2e-11);
        assertArrayEquals(baseline.covariance(),
            changedUnits.covariance(), 2e-11);
    }

    @Test
    void perfectFirstStageHasInfiniteRobustAndClusterDiagnostics() {
        int observations = 12;
        double[] response = new double[observations];
        double[][] exogenous = new double[observations][1];
        double[][] endogenous = new double[observations][1];
        double[][] instruments = new double[observations][1];
        int[] clusters = new int[observations];
        for (int row = 0; row < observations; row++) {
            double instrument = row - 5.5;
            exogenous[row][0] = 1.0;
            instruments[row][0] = instrument;
            endogenous[row][0] = instrument;
            response[row] = 0.7 + 1.8 * instrument
                + 0.2 * Math.sin(0.9 * (row + 1));
            clusters[row] = row / 2;
        }

        InstrumentalVariableResult robust = InstrumentalVariableRegression.fit(
            response, exogenous, endogenous, instruments,
            new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.HC1, 0.95),
            BackendPolicy.CPU);
        InstrumentStrengthDiagnostic robustDiagnostic =
            robust.strengthDiagnostics().get(0);
        assertEquals(Double.POSITIVE_INFINITY,
            robustDiagnostic.effectiveFStatistic());
        assertEquals(0.0, robustDiagnostic.pValue(), 0.0);
        assertEquals(StatisticDistribution.CHI_SQUARE,
            robustDiagnostic.referenceDistribution());

        InstrumentalVariableResult clustered = InstrumentalVariableRegression.fit(
            response, exogenous, endogenous, instruments, clusters,
            new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.CLUSTER_CR1, 0.95),
            BackendPolicy.CPU);
        InstrumentStrengthDiagnostic clusterDiagnostic =
            clustered.strengthDiagnostics().get(0);
        assertEquals(Double.POSITIVE_INFINITY,
            clusterDiagnostic.effectiveFStatistic());
        assertEquals(0.0, clusterDiagnostic.pValue(), 0.0);
        assertEquals(StatisticDistribution.F,
            clusterDiagnostic.referenceDistribution());
    }

    @Test
    void rejectsOrderAndRankIdentificationFailures() {
        double[] y = {1, 2, 3, 4, 5, 6};
        double[][] w = {{1}, {1}, {1}, {1}, {1}, {1}};
        double[][] twoEndogenous = {
            {1, 2}, {2, 1}, {3, 4}, {4, 3}, {5, 6}, {6, 5}
        };
        double[][] oneInstrument = {{-2}, {-1}, {0}, {1}, {2}, {3}};
        IllegalArgumentException order = assertThrows(
            IllegalArgumentException.class,
            () -> InstrumentalVariableRegression.fit(y, w, twoEndogenous,
                oneInstrument));
        assertTrue(order.getMessage().contains("underidentified"));

        double[][] oneEndogenous = {{1}, {2}, {1}, {3}, {2}, {4}};
        double[][] redundant = {
            {-2, -4}, {-1, -2}, {0, 0}, {1, 2}, {2, 4}, {3, 6}
        };
        IllegalArgumentException instrumentRank = assertThrows(
            IllegalArgumentException.class,
            () -> InstrumentalVariableRegression.fit(y, w, oneEndogenous,
                redundant));
        assertTrue(instrumentRank.getMessage().contains("instrument design"));

        double[][] orthogonalEndogenous = {
            {1}, {-2}, {1}, {0}, {0}, {0}
        };
        double[][] z = {{-1}, {0}, {1}, {2}, {3}, {4}};
        IllegalArgumentException relevance = assertThrows(
            IllegalArgumentException.class,
            () -> InstrumentalVariableRegression.fit(y, w,
                orthogonalEndogenous, z));
        assertTrue(relevance.getMessage().contains("not identified"));
    }

    @Test
    void validatesClusterModeAndReturnsDefensiveCopies() {
        Data data = data(60);
        InstrumentalVariableOptions clusterOptions =
            new InstrumentalVariableOptions(
                InstrumentalVariableCovariance.CLUSTER_CR0, 0.95);
        assertThrows(IllegalArgumentException.class,
            () -> InstrumentalVariableRegression.fit(data.response(),
                data.exogenous(), data.endogenous(), data.instruments(),
                clusterOptions, BackendPolicy.CPU));
        assertThrows(IllegalArgumentException.class,
            () -> InstrumentalVariableRegression.fit(data.response(),
                data.exogenous(), data.endogenous(), data.instruments(),
                new int[60], clusterOptions, BackendPolicy.CPU));

        InstrumentalVariableResult fit = InstrumentalVariableRegression.fit(
            data.response(), data.exogenous(), data.endogenous(),
            data.instruments());
        double coefficient = fit.coefficients()[0];
        double fitted = fit.instrumentedEndogenous()[0][0];
        double[] coefficients = fit.coefficients();
        double[][] instrumented = fit.instrumentedEndogenous();
        coefficients[0] = -100;
        instrumented[0][0] = -100;
        assertEquals(coefficient, fit.coefficients()[0]);
        assertEquals(fitted, fit.instrumentedEndogenous()[0][0]);
    }

    private static Data data(int observations) {
        double[] response = new double[observations];
        double[][] exogenous = new double[observations][2];
        double[][] endogenous = new double[observations][1];
        double[][] instruments = new double[observations][2];
        for (int row = 0; row < observations; row++) {
            double w = ((row % 17) - 8) / 5.0;
            double z1 = Math.sin((row + 1) * 0.37) + 0.15 * w;
            double z2 = Math.cos((row + 2) * 0.23) - 0.10 * w;
            double u = 0.45 * Math.sin((row + 3) * 1.17)
                + ((row / 6) % 5 - 2) * 0.08;
            double x = 0.9 * z1 - 0.55 * z2 + 0.25 * w
                + 0.8 * u + 0.12 * Math.cos((row + 1) * 0.71);
            response[row] = 0.7 - 0.3 * w + 2.0 * x + u;
            exogenous[row][0] = 1.0;
            exogenous[row][1] = w;
            endogenous[row][0] = x;
            instruments[row][0] = z1;
            instruments[row][1] = z2;
        }
        return new Data(response, exogenous, endogenous, instruments);
    }

    private static double sumSquares(double[] values) {
        double result = 0.0;
        for (double value : values) result += value * value;
        return result;
    }

    private record Data(
            double[] response,
            double[][] exogenous,
            double[][] endogenous,
            double[][] instruments) { }
}
