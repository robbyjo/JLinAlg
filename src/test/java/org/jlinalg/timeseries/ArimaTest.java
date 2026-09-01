/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;
import org.junit.jupiter.api.Test;

class ArimaTest {
    @Test
    void arOneRecoversSimulatedCoefficientAndForecasts() {
        double[] series = simulateArma(500, 200, 2.0,
            new double[] {0.7}, new double[0], 112358L);

        ArimaResult result = Arima.fit(series, ArimaOrder.ar(1));

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.7, result.autoregressive()[0], 0.07);
        assertEquals(2.0, result.location(), 0.2);
        assertTrue(result.innovationVariance() > 0.7);
        assertTrue(result.innovationVariance() < 1.3);
        ArimaForecast forecast = result.forecast(4);
        assertEquals(4, forecast.means().length);
        assertTrue(forecast.standardErrors()[3] > forecast.standardErrors()[0]);
        assertEquals(7, result.residualAutocorrelation(6).length);
        LjungBoxResult diagnostic = result.ljungBox(10);
        assertEquals(9, diagnostic.degreesOfFreedom());
        assertTrue(diagnostic.pValue() >= 0.0 && diagnostic.pValue() <= 1.0);
    }

    @Test
    void movingAverageOneUsesRCompatiblePositiveSign() {
        double[] series = simulateArma(600, 200, -1.0,
            new double[0], new double[] {0.55}, 271828L);

        ArimaResult result = Arima.fit(series, ArimaOrder.ma(1));

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.55, result.movingAverage()[0], 0.08);
        assertEquals(-1.0, result.location(), 0.15);
    }

    @Test
    void higherOrderArAndArmaUseTheRequestedOrders() {
        double[] arTwoSeries = simulateArma(800, 300, 0.5,
            new double[] {0.55, -0.2}, new double[0], 161803L);
        ArimaResult arTwo = Arima.fit(arTwoSeries, ArimaOrder.ar(2));
        assertEquals(2, arTwo.autoregressive().length);
        assertEquals(0.55, arTwo.autoregressive()[0], 0.1);
        assertEquals(-0.2, arTwo.autoregressive()[1], 0.1);

        double[] armaSeries = simulateArma(1_000, 300, 0.0,
            new double[] {0.55}, new double[] {0.35}, 141421L);
        ArimaResult arma = Arima.fit(armaSeries,
            ArimaOrder.arma(1, 1), ArimaOptions.builder()
                .includeMean(false)
                .build());
        assertEquals(0.55, arma.autoregressive()[0], 0.12);
        assertEquals(0.35, arma.movingAverage()[0], 0.12);
    }

    @Test
    void arimaRandomWalkWithDriftForecastsOnOriginalScale() {
        double[] series = new double[80];
        series[0] = 10.0;
        for (int index = 1; index < series.length; index++) {
            series[index] = series[index - 1] + 2.0
                + (index % 2 == 0 ? 0.1 : -0.1);
        }
        ArimaOptions options = ArimaOptions.builder()
            .includeDrift(true)
            .build();

        ArimaResult result = Arima.fit(
            series, ArimaOrder.arima(0, 1, 0), options);
        ArimaForecast forecast = result.forecast(3);

        assertEquals(2.0, result.location(), 0.01);
        assertArrayEquals(new double[] {
            series[series.length - 1] + 2.0,
            series[series.length - 1] + 4.0,
            series[series.length - 1] + 6.0
        }, forecast.means(), 0.03);
        assertTrue(forecast.standardErrors()[2] > forecast.standardErrors()[0]);
    }

    @Test
    void seasonalAutoregressionRecoversLagFourCoefficient() {
        double[] innovations = gaussian(700, 314159L);
        double[] all = new double[innovations.length];
        for (int time = 0; time < all.length; time++) {
            all[time] = innovations[time]
                + (time >= 4 ? 0.6 * all[time - 4] : 0.0);
        }
        double[] series = java.util.Arrays.copyOfRange(all, 200, all.length);
        ArimaOptions options = ArimaOptions.builder()
            .seasonalOrder(SeasonalArimaOrder.of(1, 0, 0, 4))
            .includeMean(false)
            .build();

        ArimaResult result = Arima.fit(
            series, ArimaOrder.arma(0, 0), options);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.6, result.seasonalAutoregressive()[0], 0.08);
    }

    @Test
    void arOneCorrelationIsExactToeplitzPowerSequence() {
        double[] correlation = Arima.correlationMatrix(
            4, new double[] {0.5}, new double[0]);

        assertArrayEquals(new double[] {
            1.0, 0.5, 0.25, 0.125,
            0.5, 1.0, 0.5, 0.25,
            0.25, 0.5, 1.0, 0.5,
            0.125, 0.25, 0.5, 1.0
        }, correlation, 1e-12);
    }

    @Test
    void rejectsDriftWithoutExactlyOneDifference() {
        ArimaOptions options = ArimaOptions.builder()
            .includeDrift(true)
            .build();
        assertThrows(IllegalArgumentException.class,
            () -> Arima.fit(new double[] {1, 2, 3, 4, 5},
                ArimaOrder.arma(1, 0), options));
    }

    private static double[] simulateArma(
            int length,
            int burn,
            double mean,
            double[] ar,
            double[] ma,
            long seed) {
        int total = length + burn;
        double[] innovations = gaussian(total, seed);
        double[] values = new double[total];
        for (int time = 0; time < total; time++) {
            double value = mean + innovations[time];
            for (int lag = 1; lag <= ar.length && lag <= time; lag++) {
                value += ar[lag - 1] * (values[time - lag] - mean);
            }
            for (int lag = 1; lag <= ma.length && lag <= time; lag++) {
                value += ma[lag - 1] * innovations[time - lag];
            }
            values[time] = value;
        }
        return java.util.Arrays.copyOfRange(values, burn, total);
    }

    private static double[] gaussian(int length, long seed) {
        Random random = new Random(seed);
        double[] result = new double[length];
        for (int index = 0; index < length; index++) {
            result[index] = random.nextGaussian();
        }
        return result;
    }
}
