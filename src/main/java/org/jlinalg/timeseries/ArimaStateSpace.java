/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import jdistlib.Normal;

/** Innovations form with an exact (finite + kappa * diffuse) covariance split.
 * All covariance workspaces are state-sized, never observation-sized. */
final class ArimaStateSpace {
    private final int size;
    private final double[] transition, noise, observation, initial;
    private final int[][] nonzero;
    private final int[] observedState;
    private final int stationary, diffuse;

    ArimaStateSpace(double[] ar, double[] ma, double[] difference) {
        stationary = Math.max(ar.length, ma.length + 1);
        diffuse = difference.length - 1;
        size = stationary + diffuse;
        transition = new double[size * size];
        observation = new double[size];
        observation[0] = 1;
        for (int i = 0; i < ar.length; i++) transition[i * size] = ar[i];
        for (int i = 1; i < stationary; i++) transition[(i - 1) * size + i] = 1;
        for (int i = 0; i < diffuse; i++) observation[stationary + i] = -difference[i + 1];
        int observedCount = 0;
        for (double weight : observation) if (weight != 0) observedCount++;
        observedState = new int[observedCount];
        observedCount = 0;
        for (int i = 0; i < size; i++) if (observation[i] != 0) observedState[observedCount++] = i;
        if (diffuse > 0) {
            System.arraycopy(observation, 0, transition, stationary * size, size);
            for (int i = stationary + 1; i < size; i++) transition[i * size + i - 1] = 1;
        }
        double[] loading = new double[size];
        loading[0] = 1;
        System.arraycopy(ma, 0, loading, 1, ma.length);
        noise = new double[size * size];
        for (int i = 0; i < size; i++) for (int j = 0; j < size; j++)
            noise[i * size + j] = loading[i] * loading[j];
        nonzero = new int[size][];
        for (int i = 0; i < size; i++) {
            int count = 0;
            for (int j = 0; j < size; j++) if (transition[i * size + j] != 0) count++;
            nonzero[i] = new int[count];
            count = 0;
            for (int j = 0; j < size; j++) if (transition[i * size + j] != 0)
                nonzero[i][count++] = j;
        }
        initial = stationaryCovariance();
    }

    /** Solves P = T P T' + R R' by doubling, including near-unit stationary roots. */
    private double[] stationaryCovariance() {
        int r = stationary;
        double[] a = new double[r * r], p = new double[r * r];
        for (int i = 0; i < r; i++) for (int j = 0; j < r; j++) {
            a[i * r + j] = transition[i * size + j];
            p[i * r + j] = noise[i * size + j];
        }
        for (int iteration = 0; iteration < 64; iteration++) {
            double[] ap = multiply(a, p, r), increment = new double[r * r];
            double change = 0, scale = 0;
            for (int i = 0; i < r; i++) for (int j = 0; j < r; j++) {
                for (int k = 0; k < r; k++) increment[i * r + j] += ap[i * r + k] * a[j * r + k];
                change = Math.max(change, Math.abs(increment[i * r + j]));
                p[i * r + j] += increment[i * r + j];
                scale = Math.max(scale, Math.abs(p[i * r + j]));
            }
            if (!Double.isFinite(scale)) break;
            if (change <= 2e-15 * scale) {
                double[] result = new double[size * size];
                for (int i = 0; i < r; i++) System.arraycopy(p, i * r, result, i * size, r);
                return result;
            }
            a = multiply(a, a, r);
        }
        throw new IllegalArgumentException("stationary covariance failed to converge");
    }

    private static double[] multiply(double[] a, double[] b, int n) {
        double[] c = new double[n * n];
        for (int i = 0; i < n; i++) for (int k = 0; k < n; k++)
            if (a[i * n + k] != 0) for (int j = 0; j < n; j++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
        return c;
    }

    Evaluation filter(double[] values, double location, double slope, boolean retain) {
        return filter(values, location, slope, retain, false);
    }

    /** Unrestricted GLS profiling of a linear trend through the same innovations. */
    Evaluation profileSlope(double[] values, boolean retain) {
        return filter(values, 0, 0, retain, true);
    }

    private Evaluation filter(double[] values, double location, double slope, boolean retain, boolean profile) {
        double[] a = new double[size], p = initial.clone(), pi = new double[size * size];
        for (int i = stationary; i < size; i++) pi[i * size + i] = 1;
        double[] u = new double[size], ui = new double[size], nextA = new double[size];
        double[] temp = new double[size * size], next = new double[size * size];
        double[] residuals = retain ? new double[values.length] : new double[0];
        double[] trendState = profile ? new double[size] : null;
        double[] trendResiduals = profile && retain ? new double[values.length] : null;
        double trendInformation = 0, estimatedSlope = 0;
        if (retain) Arrays.fill(residuals, Double.NaN);
        double quadratic = 0, logDet = 0, diffuseLogDet = 0;
        int remaining = diffuse, count = 0;
        // With no observations, the diffuse initial levels are arbitrary. Move
        // the time origin to the first observation instead of accumulating huge
        // cancelling finite/diffuse covariances across an unobserved prefix.
        // Keep original t for regression trends and the returned forecast grid.
        int first = 0;
        while (first < values.length && Double.isNaN(values[first])) first++;
        for (int t = first; t < values.length; t++) {
            if (Double.isInfinite(values[t])) throw new IllegalArgumentException("infinite observation; use NaN for missing");
            if (t > first) {
                predictMean(a, nextA);
                System.arraycopy(nextA, 0, a, 0, size);
                if (profile) {
                    predictMean(trendState, nextA);
                    System.arraycopy(nextA, 0, trendState, 0, size);
                }
                propagate(p, temp, next, true);
                System.arraycopy(next, 0, p, 0, p.length);
                if (remaining > 0) {
                    propagate(pi, temp, next, false);
                    System.arraycopy(next, 0, pi, 0, pi.length);
                }
            }
            if (!Double.isFinite(values[t])) continue;
            double v = values[t] - location - slope * (t + 1) - dot(observation, a);
            double trendV = profile ? t + 1 - dot(observation, trendState) : 0;
            project(p, u);
            double f = dot(observation, u);
            double fi = 0;
            if (remaining > 0) {
                project(pi, ui);
                fi = dot(observation, ui);
            }
            if (remaining > 0 && fi > 1e-9) {
                // Coefficients of kappa and 1 in P - P Z' (Z P Z')^-1 Z P.
                for (int i = 0; i < size; i++) {
                    a[i] += ui[i] * v / fi;
                    if (profile) trendState[i] += ui[i] * trendV / fi;
                    for (int j = 0; j <= i; j++) {
                        double finite = p[i * size + j] - (ui[i] * u[j] + u[i] * ui[j]) / fi
                            + f * (ui[i] / fi) * (ui[j] / fi);
                        p[i * size + j] = p[j * size + i] = finite;
                        double infinite = pi[i * size + j] - ui[i] * ui[j] / fi;
                        pi[i * size + j] = pi[j * size + i] = infinite;
                    }
                }
                diffuseLogDet += Math.log(fi);
                if (--remaining == 0) Arrays.fill(pi, 0);
            } else {
                if (!(f > 0) || !Double.isFinite(f))
                    throw new IllegalArgumentException("nonpositive innovation variance");
                double error = v / Math.sqrt(f);
                double trendError = trendV / Math.sqrt(f);
                if (profile) {
                    // Recursive least squares avoids Qyy-Qxy^2/Qxx cancellation
                    // when the series has a large drift and small random errors.
                    double nextInformation = trendInformation + trendError * trendError;
                    if (nextInformation > 0) {
                        double difference = error - trendError * estimatedSlope;
                        quadratic += trendInformation / nextInformation * difference * difference;
                        estimatedSlope += trendError / nextInformation * difference;
                        trendInformation = nextInformation;
                    } else quadratic += error * error;
                } else quadratic += error * error;
                logDet += Math.log(f);
                count++;
                if (retain) {
                    residuals[t] = error;
                    if (profile) trendResiduals[t] = trendError;
                }
                for (int i = 0; i < size; i++) {
                    a[i] += u[i] * v / f;
                    if (profile) trendState[i] += u[i] * trendV / f;
                    for (int j = 0; j <= i; j++) {
                        double updated = p[i * size + j] - u[i] * u[j] / f;
                        p[i * size + j] = p[j * size + i] = updated;
                    }
                }
            }
        }
        if (profile) {
            if (!(trendInformation > 0)) throw new IllegalArgumentException("trend is unidentified in the observed pattern");
            slope = estimatedSlope;
            for (int i = 0; i < size; i++) a[i] -= slope * trendState[i];
            if (retain) for (int t = 0; t < values.length; t++) residuals[t] -= slope * trendResiduals[t];
        }
        return new Evaluation(quadratic, logDet, diffuseLogDet, count, remaining,
            residuals, retain ? new ForecastState(this, a, p, location, slope, values.length) : null, slope);
    }

    private void project(double[] p, double[] u) {
        Arrays.fill(u, 0);
        for (int i = 0; i < size; i++) for (int j : observedState)
            u[i] += p[i * size + j] * observation[j];
    }

    private void predictMean(double[] a, double[] result) {
        Arrays.fill(result, 0);
        for (int i = 0; i < size; i++) for (int k : nonzero[i])
            result[i] += transition[i * size + k] * a[k];
    }

    private void propagate(double[] p, double[] temp, double[] result, boolean addNoise) {
        Arrays.fill(temp, 0);
        for (int i = 0; i < size; i++) for (int k : nonzero[i])
            for (int j = 0; j < size; j++) temp[i * size + j] += transition[i * size + k] * p[k * size + j];
        for (int i = 0; i < size; i++) for (int j = 0; j <= i; j++) {
            double value = addNoise ? noise[i * size + j] : 0;
            for (int k : nonzero[j]) value += temp[i * size + k] * transition[j * size + k];
            result[i * size + j] = result[j * size + i] = value;
        }
    }

    private static double dot(double[] a, double[] b) {
        double value = 0;
        for (int i = 0; i < a.length; i++) value += a[i] * b[i];
        return value;
    }

    record Evaluation(double quadratic, double logDeterminant, double diffuseLogDeterminant,
                      int observations, int remainingDiffuse, double[] innovations, ForecastState state,
                      double fittedSlope) {
        double variance() { return quadratic / observations; }
        // R convention: omit all diffuse innovation terms (including their determinants).
        double negativeLogLikelihood() {
            if (remainingDiffuse != 0 || observations == 0 || !(variance() > 0)
                    || !Double.isFinite(variance())) return Double.POSITIVE_INFINITY;
            return 0.5 * (observations * (Math.log(2 * Math.PI) + 1 + Math.log(variance())) + logDeterminant);
        }
    }

    /** A private immutable snapshot; each forecast uses fresh workspaces. */
    static final class ForecastState {
        private final ArimaStateSpace model;
        private final double[] a, p;
        private final double location, slope;
        private final int length;
        ForecastState(ArimaStateSpace model, double[] a, double[] p, double location, double slope, int length) {
            this.model = model; this.a = a.clone(); this.p = p.clone();
            this.location = location; this.slope = slope; this.length = length;
        }
        ArimaForecast forecast(int horizon, double level, double variance) {
            if (horizon < 1 || !(level > 0 && level < 1))
                throw new IllegalArgumentException("positive horizon and confidence level in (0,1) required");
            int n = model.size;
            double[] mean = new double[horizon], se = new double[horizon];
            double[] lower = new double[horizon], upper = new double[horizon];
            double[] state = a.clone(), covariance = p.clone(), nextA = new double[n];
            double[] temp = new double[n * n], next = new double[n * n], u = new double[n];
            double critical = Normal.quantile(0.5 + level / 2, 0, 1, true, false);
            for (int h = 0; h < horizon; h++) {
                model.predictMean(state, nextA);
                System.arraycopy(nextA, 0, state, 0, n);
                model.propagate(covariance, temp, next, true);
                System.arraycopy(next, 0, covariance, 0, next.length);
                model.project(covariance, u);
                mean[h] = dot(model.observation, state) + location + slope * (length + h + 1);
                se[h] = Math.sqrt(Math.max(0, variance * dot(model.observation, u)));
                lower[h] = mean[h] - critical * se[h]; upper[h] = mean[h] + critical * se[h];
            }
            return new ArimaForecast(mean, se, lower, upper, level);
        }
    }
}
