/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class TimeSeriesReviewRegressionTest {
    @Test void leadingMissingPrefixDoesNotChangeDiffuseLikelihoodOrForecast() {
        double[] base = new double[40];
        for (int i = 0; i < base.length; i++) base[i] = Math.sin(.71 * i) + .03 * i * i;
        for (SeasonalArimaOrder seasonal : List.of(SeasonalArimaOrder.none(), SeasonalArimaOrder.of(0, 1, 0, 4))) {
            ArimaOptions options = ArimaOptions.builder().seasonalOrder(seasonal).build();
            ArimaResult original = DiffuseArima.fit(base, new ArimaOrder(0, 2, 0), options).fit();
            for (int leading : new int[] {10_000, 30_000, 100_000}) {
                double[] padded = new double[leading + base.length];
                Arrays.fill(padded, Double.NaN); System.arraycopy(base, 0, padded, leading, base.length);
                DiffuseArima.Result result = DiffuseArima.fit(padded, new ArimaOrder(0, 2, 0), options);
                assertTrue(result.fit().converged());
                assertEquals(original.logLikelihood(), result.fit().logLikelihood(), 1e-11);
                assertEquals(original.innovationVariance(), result.fit().innovationVariance(), 1e-12);
                assertEquals(original.effectiveObservations(), result.likelihoodObservations());
                assertArrayEquals(original.forecast(6).means(), result.fit().forecast(6).means(), 1e-10);
                assertArrayEquals(original.forecast(6).standardErrors(), result.fit().forecast(6).standardErrors(), 1e-10);
            }
        }
    }

    @Test void missingOrdinaryDriftIsUnrestrictedGls() throws IOException {
        double[] y = new double[101]; Arrays.fill(y, Double.NaN);
        y[0] = 0; y[1] = 1; y[50] = 5000; y[100] = 10000;
        ArimaResult fit = DiffuseArima.fit(y, new ArimaOrder(0, 1, 0),
            ArimaOptions.builder().includeDrift(true).build()).fit();
        assertTrue(fit.converged());
        assertEquals(100, fit.location(), 1e-11);
        assertEquals(reference().get("drift_loglik"), fit.logLikelihood(), 1e-7);
        assertEquals(10100, fit.forecast(1).means()[0], 1e-9);
        assertEquals(3, fit.effectiveObservations());
        assertEquals(-2 * fit.logLikelihood() + 4, fit.aic(), 1e-12); // drift + variance
    }

    @Test void missingSeasonalDriftIsUnrestrictedGls() throws IOException {
        double[] y = new double[104]; Arrays.fill(y, Double.NaN);
        for (int r = 0; r < 4; r++) { y[r] = 0; y[4 + r] = 1; y[100 + r] = 10000; }
        ArimaResult fit = DiffuseArima.fit(y, new ArimaOrder(0, 0, 0),
            ArimaOptions.builder().includeDrift(true).seasonalOrder(SeasonalArimaOrder.of(0, 1, 0, 4)).build()).fit();
        assertTrue(fit.converged());
        assertEquals(400, fit.location(), 1e-10);
        assertEquals(reference().get("seasonal_drift_loglik"), fit.logLikelihood(), 1e-7);
        assertArrayEquals(new double[] {10400, 10400, 10400, 10400}, fit.forecast(4).means(), 1e-9);
    }

    @Test void scalarSearchKeepsBetterBasinAndHonorsBudget() throws IOException {
        double[] y = {-1.5699824440879728,1.7418872190064827,1.0916093731684835,
            .2963556882787251,.8255498200370827,-.9177455894496253};
        ExactArmaResult fit = ExactArma.fit(y, ArimaOrder.ma(1), false, BackendPolicy.CPU);
        assertTrue(fit.converged());
        assertEquals(reference().get("ma1"), fit.movingAverage()[0], 1e-6);
        assertEquals(reference().get("ma_loglik"), fit.logLikelihood(), 1e-10);
        BoundedOptimizer.Result capped = BoundedOptimizer.minimize(new double[] {0}, new double[] {-3.8}, new double[] {3.8},
            x -> Math.pow(x[0] - .071, 2), 20, 1e-12);
        assertTrue(capped.evaluations() <= 20);
        assertFalse(capped.converged());
        assertTrue(capped.objective() <= .071 * .071);
    }

    @Test void exactWhiteNoiseRetainsSmallPositiveVariance() throws IOException {
        double[] base = {1, -2, 3, -1, 2};
        ExactArmaResult original = ExactArma.fit(base, ArimaOrder.arma(0, 0), false, BackendPolicy.CPU);
        for (double scale : new double[] {1e-8, 1e-50, 1e50}) {
            double[] y = Arrays.stream(base).map(value -> value * scale).toArray();
            ExactArmaResult fit = ExactArma.fit(y, ArimaOrder.arma(0, 0), false, BackendPolicy.CPU);
            assertTrue(fit.converged());
            assertEquals(3.8, fit.innovationVariance() / (scale * scale), 1e-12);
            assertEquals(original.logLikelihood() - base.length * Math.log(scale), fit.logLikelihood(), 1e-10);
        }
        assertEquals(reference().get("small_loglik"), original.logLikelihood() - 5 * Math.log(1e-8), 1e-12);
        assertThrows(IllegalArgumentException.class, () -> ExactArma.fit(new double[5], ArimaOrder.arma(0, 0), false, BackendPolicy.CPU));
    }

    @Test void shortExactArFitDoesNotRequireAConditionalFit() throws IOException {
        ExactArmaResult fit = ExactArma.fit(new double[] {1, -2, 3, -1}, ArimaOrder.ar(2), false, BackendPolicy.CPU);
        assertTrue(fit.converged());
        assertEquals(reference().get("ar_loglik"), fit.logLikelihood(), 1e-9);
        assertEquals(reference().get("ar1"), fit.autoregressive()[0], 1e-5);
        assertEquals(reference().get("ar2"), fit.autoregressive()[1], 1e-5);
    }

    @Test void unidentifiedArMaCancellationHasNoManufacturedCovariance() {
        List<double[]> panels = new ArrayList<>();
        for (int bits = 0; bits < 8; bits++) {
            double[] y = new double[3];
            for (int j = 0; j < 3; j++) y[j] = (bits & (1 << j)) == 0 ? -1 : 1;
            panels.add(y);
        }
        ExactArmaResult fit = ExactArma.fitPanel(panels, ArimaOrder.arma(1, 1), false, BackendPolicy.CPU);
        assertEquals(-12 * (Math.log(2 * Math.PI) + 1), fit.logLikelihood(), 1e-9);
        assertFalse(fit.coefficientInferenceAvailable());
        assertTrue(Arrays.stream(fit.standardErrors()).allMatch(Double::isNaN));
        assertTrue(Arrays.stream(fit.coefficientCovariance()).allMatch(Double::isNaN));
    }

    @Test void exactAndMissingWrapperRejectBothInfinities() {
        for (double infinity : new double[] {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY}) {
            double[] y = {1, infinity, -2, 3, -1, 2};
            assertThrows(IllegalArgumentException.class, () -> ExactArma.fit(y, ArimaOrder.arma(0, 0), false, BackendPolicy.CPU));
            assertThrows(IllegalArgumentException.class, () -> SparseMissingSeries.fit(y, ArimaOrder.arma(0, 0), false, BackendPolicy.CPU));
        }
        ExactArmaResult fit = ExactArma.fit(new double[] {1, Double.NaN, -2, 3, -1, 2}, ArimaOrder.arma(0, 0), false, BackendPolicy.CPU);
        assertTrue(fit.converged()); assertEquals(5, fit.observedValues());
    }

    private static Map<String, Double> reference() throws IOException {
        try (var input = TimeSeriesReviewRegressionTest.class.getResourceAsStream("/timeseries/review-reference.csv")) {
            if (input == null) throw new IOException("missing review reference");
            return new String(input.readAllBytes(), StandardCharsets.UTF_8).lines().skip(1)
                .map(line -> line.split(",")).collect(Collectors.toMap(row -> row[0], row -> Double.parseDouble(row[1])));
        }
    }
}
