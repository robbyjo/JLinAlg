/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class ExactDiffuseArimaTest {
    private static final String[] CASES = {"stationary", "stationary-missing", "integrated-ar", "integrated-ma",
        "integrated-missing", "seasonal", "seasonal-missing", "seasonal-only", "twice-integrated", "drift",
        "seasonal-drift", "seasonal-ar-ma", "monthly-missing"};

    @Test void fixedParametersMatchRKalmanLikelihoodAndForecasts() throws IOException {
        for (String name : CASES) {
            Fixture f = fixture(name);
            int position = 0;
            double[] ar = Arrays.copyOfRange(f.coefficients, position, position += f.order.autoregressive());
            double[] ma = Arrays.copyOfRange(f.coefficients, position, position += f.order.movingAverage());
            double[] sar = Arrays.copyOfRange(f.coefficients, position, position += f.seasonal.autoregressive());
            double[] sma = Arrays.copyOfRange(f.coefficients, position, position += f.seasonal.movingAverage());
            double slope = f.drift ? f.coefficients[position] : 0;
            ArimaStateSpace.Evaluation filtered = new ArimaStateSpace(
                ArimaMath.multiplyAr(ar, sar, f.seasonal.period()), ArimaMath.multiplyMa(ma, sma, f.seasonal.period()),
                ArimaMath.differencingPolynomial(f.order, f.seasonal)).filter(f.series, 0, slope, true);
            assertEquals(f.ref.get("loglik"), -filtered.negativeLogLikelihood(), 1e-4, name);
            assertEquals(f.ref.get("variance"), filtered.variance(), 2e-6, name);
            assertEquals(f.ref.get("nobs").intValue(), filtered.observations(), name);
            assertEquals(0, filtered.remainingDiffuse(), name);
            ArimaForecast forecast = filtered.state().forecast(8, .95, filtered.variance());
            for (int h = 0; h < 8; h++) {
                assertEquals(f.ref.get("forecast" + (h + 1)), forecast.means()[h], 2e-5, name);
                assertEquals(f.ref.get("se" + (h + 1)), forecast.standardErrors()[h], 2e-5, name);
            }
        }
    }

    @Test void optimizedFitsMatchRIncludingSeasonalOnlyAndMissingDrift() throws IOException {
        for (String name : CASES) {
            Fixture f = fixture(name);
            double[] coefficients;
            double likelihood, variance;
            if (f.order.differences() + f.seasonal.differences() == 0) {
                ExactArmaResult fit = ExactArma.fit(f.series, f.order, false, BackendPolicy.CPU);
                assertTrue(fit.converged(), name);
                assertTrue(fit.coefficientInferenceAvailable(), name);
                likelihood = fit.logLikelihood(); variance = fit.innovationVariance();
                coefficients = join(fit.autoregressive(), fit.movingAverage());
            } else {
                DiffuseArima.Result result = DiffuseArima.fit(f.series, f.order, ArimaOptions.builder()
                    .seasonalOrder(f.seasonal).includeDrift(f.drift).optimizationTolerance(1e-9).build());
                ArimaResult fit = result.fit();
                assertTrue(fit.converged(), name + ": " + fit.convergenceMessage());
                assertTrue(result.diffuseLikelihood());
                assertEquals(f.order.differences() + f.seasonal.differences() * f.seasonal.period(), result.diffuseStateCount());
                assertEquals(f.ref.get("nobs").intValue(), result.likelihoodObservations());
                likelihood = fit.logLikelihood(); variance = fit.innovationVariance();
                coefficients = join(fit.autoregressive(), fit.movingAverage(), fit.seasonalAutoregressive(),
                    fit.seasonalMovingAverage(), f.drift ? new double[] {
                        fit.location() / (f.seasonal.differences() > 0 ? f.seasonal.period() : 1)} : new double[0]);
                ArimaForecast forecast = fit.forecast(8);
                for (int h = 0; h < 8; h++) {
                    assertEquals(f.ref.get("forecast" + (h + 1)), forecast.means()[h], 5e-5, name);
                    assertEquals(f.ref.get("se" + (h + 1)), forecast.standardErrors()[h], 5e-5, name);
                }
                assertArrayEquals(forecast.means(), fit.forecast(8).means());
            }
            assertEquals(f.ref.get("loglik"), likelihood, 1e-4, name);
            assertEquals(f.ref.get("variance"), variance, 3e-6, name);
            assertArrayEquals(f.coefficients, coefficients, 1e-5, name);
        }
    }

    @Test void randomWalkMarginalizesGapsAndPropagatesTerminalUncertainty() {
        double[] y = {Double.NaN, 10, Double.NaN, Double.NaN, 14, 13, Double.NaN, Double.NaN};
        ArimaResult fit = DiffuseArima.fit(y, new ArimaOrder(0, 1, 0), null).fit();
        double variance = (16.0 / 3 + 1) / 2;
        // Two observed increments with variances 3*sigma2 and sigma2.
        double loglik = -.5 * (2 * (Math.log(2 * Math.PI * variance) + 1) + Math.log(3));
        assertEquals(variance, fit.innovationVariance(), 1e-12);
        assertEquals(loglik, fit.logLikelihood(), 1e-12);
        assertArrayEquals(new double[] {13, 13, 13}, fit.forecast(3).means(), 1e-12);
        assertEquals(Math.sqrt(3 * variance), fit.forecast(3).standardErrors()[0], 1e-12);
        assertEquals(Math.sqrt(5 * variance), fit.forecast(3).standardErrors()[2], 1e-12);
        assertThrows(IllegalStateException.class, () -> fit.ljungBox(1));
    }

    @Test void completeDiffuseLikelihoodEqualsExactLikelihoodOfDifferences() {
        double[] y = new Random(91).doubles(150, -2, 3).toArray();
        for (ArimaOrder order : List.of(new ArimaOrder(1, 1, 1), new ArimaOrder(1, 2, 1))) {
            for (SeasonalArimaOrder seasonal : List.of(SeasonalArimaOrder.none(), SeasonalArimaOrder.of(0, 1, 0, 4))) {
                double[] ar = {.71}, ma = {.23};
                ArimaStateSpace.Evaluation diffuse = new ArimaStateSpace(ar, ma,
                    ArimaMath.differencingPolynomial(order, seasonal)).filter(y, 0, 0, false);
                ArimaStateSpace.Evaluation stationary = new ArimaStateSpace(ar, ma, new double[] {1})
                    .filter(ArimaMath.difference(y, order, seasonal), 0, 0, false);
                assertEquals(stationary.negativeLogLikelihood(), diffuse.negativeLogLikelihood(), 1e-8);
                assertEquals(0, diffuse.diffuseLogDeterminant(), 1e-8);
            }
        }
    }

    @Test void stationaryMissingLikelihoodMatchesIndependentDenseMarginalization() {
        double[] ar = {.62, -.17}, ma = {.24};
        double[] y = {1, Double.NaN, -.3, .7, Double.NaN, 2, -.4};
        int[] indices = {0, 2, 3, 5, 6};
        double[] covariance = ArimaMath.correlationMatrix(y.length, ar, ma);
        double scale = ArimaMath.marginalVariancePerInnovation(ar, ma);
        int n = indices.length;
        double[][] l = new double[n][n];
        double[] solved = new double[n];
        double quadratic = 0, logDet = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double v = scale * covariance[indices[i] * y.length + indices[j]];
                for (int k = 0; k < j; k++) v -= l[i][k] * l[j][k];
                l[i][j] = i == j ? Math.sqrt(v) : v / l[j][j];
            }
            double v = y[indices[i]];
            for (int j = 0; j < i; j++) v -= l[i][j] * solved[j];
            solved[i] = v / l[i][i]; quadratic += solved[i] * solved[i];
            logDet += 2 * Math.log(l[i][i]);
        }
        ArimaStateSpace.Evaluation filtered = new ArimaStateSpace(ar, ma, new double[] {1}).filter(y, 0, 0, false);
        assertEquals(quadratic, filtered.quadratic(), 1e-11);
        assertEquals(logDet, filtered.logDeterminant(), 1e-11);
    }

    @Test void nearUnitStationaryCovarianceDoesNotTruncateImpulseResponse() {
        double phi = .999999;
        double[] y = {2, Double.NaN, 1, -.5};
        ArimaStateSpace.Evaluation f = new ArimaStateSpace(new double[] {phi}, new double[0], new double[] {1})
            .filter(y, 0, 0, false);
        double q = 4 * (1 - phi * phi) + Math.pow(1 - 2 * phi * phi, 2) / (1 + phi * phi) + Math.pow(-.5 - phi, 2);
        assertEquals(q, f.quadratic(), 1e-9);
        assertEquals(-Math.log(1 - phi * phi) + Math.log(1 + phi * phi), f.logDeterminant(), 1e-9);
    }

    @Test void longMissingSeriesUsesBoundedStateMemory() {
        Random random = new Random(7634);
        double[] y = new double[50_000];
        double previous = 0;
        for (int i = 0; i < y.length; i++) {
            previous = .6 * previous + random.nextGaussian();
            y[i] = i % 7 == 0 ? Double.NaN : previous;
        }
        SparseMissingSeries.Result fit = SparseMissingSeries.fit(y, ArimaOrder.ar(1), false, BackendPolicy.CPU);
        assertTrue(fit.fit().converged());
        assertEquals(.6, fit.fit().autoregressive()[0], .02);
        assertTrue(fit.likelihoodPath().contains("state-space"));
        assertEquals(y.length - fit.missingCount(), fit.fit().observedValues());
    }

    @Test void unidentifiedDiffuseDirectionsAndInvalidInputsAreRejected() {
        double[] y = new Random(12).doubles(60).toArray();
        for (int i = 0; i < y.length; i += 2) y[i] = Double.NaN;
        assertThrows(IllegalArgumentException.class, () -> DiffuseArima.fit(y, new ArimaOrder(0, 0, 0),
            ArimaOptions.builder().seasonalOrder(SeasonalArimaOrder.of(0, 1, 0, 2)).build()));
        assertThrows(IllegalArgumentException.class, () -> DiffuseArima.fit(new double[] {1, 2, Double.POSITIVE_INFINITY, 4},
            new ArimaOrder(0, 1, 0), null));
        assertThrows(IllegalArgumentException.class, () -> DiffuseArima.fit(new double[] {1, 1, 1, 1}, new ArimaOrder(0, 1, 0), null));
    }

    @Test void optimizerBudgetExhaustionIsNotConvergence() {
        BoundedOptimizer.Result f = BoundedOptimizer.minimize(new double[] {-1, 1}, new double[] {-3, -3},
            new double[] {3, 3}, p -> 100 * Math.pow(p[1] - p[0] * p[0], 2) + Math.pow(1 - p[0], 2), 20, 1e-10);
        assertFalse(f.converged());
        BoundedOptimizer.Result invalid = BoundedOptimizer.minimize(new double[0], new double[0], new double[0],
            p -> Double.POSITIVE_INFINITY, 20, 1e-8);
        assertFalse(invalid.converged());
    }

    @Test void driftCanBeEstimatedWhenEveryExplicitDifferenceIsMissing() {
        double[] y = new double[90];
        Arrays.fill(y, Double.NaN);
        for (int i = 0; i < y.length; i += 3) y[i] = 100 * i + Math.sin(i);
        ArimaResult fit = DiffuseArima.fit(y, new ArimaOrder(0, 1, 0),
            ArimaOptions.builder().includeDrift(true).build()).fit();
        assertTrue(fit.converged());
        // Profile ML drift for a random walk with equally spaced observed gaps.
        assertEquals((y[87] - y[0]) / 87, fit.location(), 1e-5);
        assertTrue(Arrays.stream(fit.differencedSeries()).allMatch(Double::isNaN));
        assertTrue(Double.isFinite(fit.forecast(1).means()[0]));
    }

    private static Fixture fixture(String name) throws IOException {
        Map<String, Double> ref = new LinkedHashMap<>();
        List<String> rows = lines(name + "-reference.csv");
        for (String line : rows.subList(1, rows.size())) {
            String[] fields = line.split(","); ref.put(fields[0], Double.parseDouble(fields[1]));
        }
        double[] series = lines(name + "-series.csv").stream().skip(1).mapToDouble(Double::parseDouble).toArray();
        ArimaOrder order = new ArimaOrder(ref.get("p").intValue(), ref.get("d").intValue(), ref.get("q").intValue());
        SeasonalArimaOrder seasonal = new SeasonalArimaOrder(ref.get("P").intValue(), ref.get("D").intValue(),
            ref.get("Q").intValue(), ref.get("period").intValue());
        double[] coefficients = ref.entrySet().stream().filter(e -> e.getKey().startsWith("coef"))
            .mapToDouble(Map.Entry::getValue).toArray();
        return new Fixture(series, ref, order, seasonal, ref.get("drift") == 1, coefficients);
    }
    private static List<String> lines(String name) throws IOException {
        try (var input = ExactDiffuseArimaTest.class.getResourceAsStream("/timeseries/" + name)) {
            if (input == null) throw new IOException("missing fixture " + name);
            return new String(input.readAllBytes(), StandardCharsets.UTF_8).lines().toList();
        }
    }
    private static double[] join(double[]... arrays) { return Arrays.stream(arrays).flatMapToDouble(Arrays::stream).toArray(); }
    private record Fixture(double[] series, Map<String, Double> ref, ArimaOrder order,
                           SeasonalArimaOrder seasonal, boolean drift, double[] coefficients) { }
}
