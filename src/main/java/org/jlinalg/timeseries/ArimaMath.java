/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import java.util.Arrays;
import jdistlib.Normal;

/** Internal polynomial, innovation, covariance, and forecast calculations. */
final class ArimaMath {
    private static final int MINIMUM_PSI_TERMS = 2_000;
    private static final int MAXIMUM_PSI_TERMS = 100_000;
    private static final double PSI_TOLERANCE = 1e-12;

    private ArimaMath() { }

    static double[] difference(
            double[] series, ArimaOrder order, SeasonalArimaOrder seasonal) {
        double[] result = series.clone();
        for (int count = 0; count < seasonal.differences(); count++) {
            result = differenceOnce(result, seasonal.period());
        }
        for (int count = 0; count < order.differences(); count++) {
            result = differenceOnce(result, 1);
        }
        return result;
    }

    static double[] differencingPolynomial(
            ArimaOrder order, SeasonalArimaOrder seasonal) {
        double[] result = {1.0};
        for (int count = 0; count < seasonal.differences(); count++) {
            double[] factor = new double[seasonal.period() + 1];
            factor[0] = 1.0;
            factor[seasonal.period()] = -1.0;
            result = convolve(result, factor);
        }
        for (int count = 0; count < order.differences(); count++) {
            result = convolve(result, new double[] {1.0, -1.0});
        }
        return result;
    }

    static Coefficients decode(
            double[] parameters,
            ArimaOrder order,
            SeasonalArimaOrder seasonal,
            boolean locationIncluded,
            boolean drift) {
        int offset = 0;
        double[] ar = stableCoefficients(parameters, offset,
            order.autoregressive(), false);
        offset += order.autoregressive();
        double[] ma = stableCoefficients(parameters, offset,
            order.movingAverage(), true);
        offset += order.movingAverage();
        double[] seasonalAr = stableCoefficients(parameters, offset,
            seasonal.autoregressive(), false);
        offset += seasonal.autoregressive();
        double[] seasonalMa = stableCoefficients(parameters, offset,
            seasonal.movingAverage(), true);
        offset += seasonal.movingAverage();
        double location = locationIncluded ? parameters[offset] : 0.0;
        double[] effectiveAr = multiplyAr(ar, seasonalAr, seasonal.period());
        double[] effectiveMa = multiplyMa(ma, seasonalMa, seasonal.period());
        return new Coefficients(ar, ma, seasonalAr, seasonalMa,
            effectiveAr, effectiveMa, location, drift);
    }

    static Innovations innovations(double[] series, Coefficients coefficients) {
        double[] errors = new double[series.length];
        double[] ar = coefficients.effectiveAr();
        double[] ma = coefficients.effectiveMa();
        double location = coefficients.location();
        for (int time = 0; time < series.length; time++) {
            double prediction = location;
            for (int lag = 1; lag <= ar.length && lag <= time; lag++) {
                prediction += ar[lag - 1]
                    * (series[time - lag] - location);
            }
            for (int lag = 1; lag <= ma.length && lag <= time; lag++) {
                prediction += ma[lag - 1] * errors[time - lag];
            }
            errors[time] = series[time] - prediction;
        }
        int warmup = ar.length;
        int used = series.length - warmup;
        double rss = 0.0;
        for (int time = warmup; time < series.length; time++) {
            rss += errors[time] * errors[time];
        }
        return new Innovations(errors, rss, used, warmup);
    }

    static double[] correlationMatrix(
            int observations, double[] ar, double[] ma) {
        if (observations < 1) {
            throw new IllegalArgumentException("observations must be positive");
        }
        double[] psi = psiWeights(ar, ma,
            Math.max(MINIMUM_PSI_TERMS, observations + 500));
        double[] autocovariance = new double[observations];
        for (int lag = 0; lag < observations; lag++) {
            double sum = 0.0;
            for (int index = 0; index + lag < psi.length; index++) {
                sum += psi[index] * psi[index + lag];
            }
            autocovariance[lag] = sum;
        }
        if (!(autocovariance[0] > 0.0)
                || !Double.isFinite(autocovariance[0])) {
            throw new IllegalArgumentException(
                "ARMA parameters do not define a finite stationary covariance");
        }
        double[] result = new double[observations * observations];
        for (int row = 0; row < observations; row++) {
            for (int column = 0; column <= row; column++) {
                double value = autocovariance[row - column] / autocovariance[0];
                result[row * observations + column] = value;
                result[column * observations + row] = value;
            }
        }
        return result;
    }

    static double marginalVariancePerInnovation(double[] ar, double[] ma) {
        double[] psi = psiWeights(ar, ma, MINIMUM_PSI_TERMS);
        double result = 0.0;
        for (double value : psi) result += value * value;
        return result;
    }

    static ArimaForecast forecast(
            ArimaResult result,
            double[] original,
            int horizon,
            double confidenceLevel) {
        if (horizon < 1) {
            throw new IllegalArgumentException("forecast horizon must be positive");
        }
        if (!(confidenceLevel > 0.0) || !(confidenceLevel < 1.0)
                || !Double.isFinite(confidenceLevel)) {
            throw new IllegalArgumentException(
                "confidenceLevel must lie strictly between zero and one");
        }
        double[] differenced = result.differencedSeries();
        double[] innovations = result.innovations();
        double[] ar = multiplyAr(result.autoregressive(),
            result.seasonalAutoregressive(), result.seasonalOrder().period());
        double[] ma = multiplyMa(result.movingAverage(),
            result.seasonalMovingAverage(), result.seasonalOrder().period());
        double[] futureDifferenced = new double[differenced.length + horizon];
        System.arraycopy(differenced, 0, futureDifferenced, 0, differenced.length);
        for (int step = 0; step < horizon; step++) {
            int time = differenced.length + step;
            double value = result.location();
            for (int lag = 1; lag <= ar.length; lag++) {
                int index = time - lag;
                if (index >= 0) {
                    value += ar[lag - 1]
                        * (futureDifferenced[index] - result.location());
                }
            }
            for (int lag = 1; lag <= ma.length; lag++) {
                int index = time - lag;
                if (index >= 0 && index < innovations.length) {
                    value += ma[lag - 1] * innovations[index];
                }
            }
            futureDifferenced[time] = value;
        }

        double[] delta = differencingPolynomial(result.order(),
            result.seasonalOrder());
        double[] extended = Arrays.copyOf(original, original.length + horizon);
        for (int step = 0; step < horizon; step++) {
            int time = original.length + step;
            double value = futureDifferenced[differenced.length + step];
            for (int lag = 1; lag < delta.length; lag++) {
                value -= delta[lag] * extended[time - lag];
            }
            extended[time] = value;
        }
        double[] stationaryPsi = psiWeights(ar, ma, horizon);
        double[] integratedPsi = new double[horizon];
        for (int index = 0; index < horizon; index++) {
            double value = stationaryPsi[index];
            for (int lag = 1; lag < delta.length && lag <= index; lag++) {
                value -= delta[lag] * integratedPsi[index - lag];
            }
            integratedPsi[index] = value;
        }
        double critical = Normal.quantile(
            0.5 + confidenceLevel / 2.0, 0.0, 1.0, true, false);
        double[] means = new double[horizon];
        double[] standardErrors = new double[horizon];
        double[] lower = new double[horizon];
        double[] upper = new double[horizon];
        double cumulative = 0.0;
        for (int step = 0; step < horizon; step++) {
            means[step] = extended[original.length + step];
            cumulative += integratedPsi[step] * integratedPsi[step];
            standardErrors[step] = Math.sqrt(
                result.innovationVariance() * cumulative);
            lower[step] = means[step] - critical * standardErrors[step];
            upper[step] = means[step] + critical * standardErrors[step];
        }
        return new ArimaForecast(
            means, standardErrors, lower, upper, confidenceLevel);
    }

    static double[] multiplyAr(
            double[] ar, double[] seasonalAr, int period) {
        double[] ordinary = new double[ar.length + 1];
        ordinary[0] = 1.0;
        for (int index = 0; index < ar.length; index++) {
            ordinary[index + 1] = -ar[index];
        }
        double[] seasonal = new double[seasonalAr.length * period + 1];
        seasonal[0] = 1.0;
        for (int index = 0; index < seasonalAr.length; index++) {
            seasonal[(index + 1) * period] = -seasonalAr[index];
        }
        double[] polynomial = convolve(ordinary, seasonal);
        double[] result = new double[polynomial.length - 1];
        for (int index = 1; index < polynomial.length; index++) {
            result[index - 1] = -polynomial[index];
        }
        return trim(result);
    }

    static double[] multiplyMa(
            double[] ma, double[] seasonalMa, int period) {
        double[] ordinary = new double[ma.length + 1];
        ordinary[0] = 1.0;
        System.arraycopy(ma, 0, ordinary, 1, ma.length);
        double[] seasonal = new double[seasonalMa.length * period + 1];
        seasonal[0] = 1.0;
        for (int index = 0; index < seasonalMa.length; index++) {
            seasonal[(index + 1) * period] = seasonalMa[index];
        }
        double[] polynomial = convolve(ordinary, seasonal);
        return trim(Arrays.copyOfRange(polynomial, 1, polynomial.length));
    }

    private static double[] stableCoefficients(
            double[] parameters, int offset, int count, boolean movingAverage) {
        double[] coefficients = new double[count];
        for (int order = 0; order < count; order++) {
            double reflection = Math.tanh(parameters[offset + order]);
            double[] updated = coefficients.clone();
            updated[order] = reflection;
            for (int index = 0; index < order; index++) {
                updated[index] = coefficients[index]
                    - reflection * coefficients[order - index - 1];
            }
            coefficients = updated;
        }
        if (movingAverage) {
            for (int index = 0; index < coefficients.length; index++) {
                coefficients[index] = -coefficients[index];
            }
        }
        return coefficients;
    }

    private static double[] psiWeights(double[] ar, double[] ma, int requested) {
        int minimumLength = Math.max(requested, MINIMUM_PSI_TERMS);
        int length = MAXIMUM_PSI_TERMS;
        double[] psi = new double[length];
        psi[0] = 1.0;
        int quiet = 0;
        int used = length;
        for (int index = 1; index < length; index++) {
            double value = index <= ma.length ? ma[index - 1] : 0.0;
            for (int lag = 1; lag <= ar.length && lag <= index; lag++) {
                value += ar[lag - 1] * psi[index - lag];
            }
            psi[index] = value;
            if (index >= minimumLength && Math.abs(value) < PSI_TOLERANCE) {
                quiet++;
                if (quiet >= 64) {
                    used = index + 1;
                    break;
                }
            } else {
                quiet = 0;
            }
        }
        return used == length ? psi : Arrays.copyOf(psi, used);
    }

    private static double[] differenceOnce(double[] values, int lag) {
        if (values.length <= lag) {
            throw new IllegalArgumentException(
                "differencing order leaves no observations");
        }
        double[] result = new double[values.length - lag];
        for (int index = lag; index < values.length; index++) {
            result[index - lag] = values[index] - values[index - lag];
        }
        return result;
    }

    private static double[] convolve(double[] left, double[] right) {
        double[] result = new double[left.length + right.length - 1];
        for (int first = 0; first < left.length; first++) {
            for (int second = 0; second < right.length; second++) {
                result[first + second] += left[first] * right[second];
            }
        }
        return result;
    }

    private static double[] trim(double[] values) {
        int length = values.length;
        while (length > 0 && values[length - 1] == 0.0) {
            length--;
        }
        return Arrays.copyOf(values, length);
    }

    record Coefficients(
            double[] ar,
            double[] ma,
            double[] seasonalAr,
            double[] seasonalMa,
            double[] effectiveAr,
            double[] effectiveMa,
            double location,
            boolean drift) { }

    record Innovations(double[] values, double rss, int used, int warmup) { }
}
