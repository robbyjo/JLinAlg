/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.enrichment;

/** Deterministic two-colour successive biased sampling, without replacement.
 * The dynamic program evaluates Wallenius' distribution, not Fisher's noncentral
 * hypergeometric distribution. No normal approximation or central-tail fallback is used.
 */
public final class Wallenius {
    private Wallenius() { }

    /** Inclusive upper tail. Work is bounded to prevent accidental unbounded calculations. */
    public static double upperTail(int red, int white, int draws, int observed, double odds) {
        return upperTail(red, white, draws, observed, odds, 100_000_000L);
    }

    /** Inclusive upper tail with a caller-specified maximum number of state updates. */
    public static double upperTail(int red, int white, int draws, int observed,
            double odds, long maximumUpdates) {
        if (red < 0 || white < 0 || draws < 0 || draws > (long) red + white
                || Double.isNaN(odds) || odds < 0 || maximumUpdates < 1)
            throw new IllegalArgumentException("invalid Wallenius parameters");
        int low = Math.max(0, draws - white), high = Math.min(draws, red);
        if (observed <= low) return 1;
        if (observed > high) return 0;
        if (odds == 0) return 0;
        if (Double.isInfinite(odds)) return 1;
        if (odds == 1) return jdistlib.HyperGeometric.cumulative(
            observed - 1.0, red, white, draws, false, false);
        // States at/above observed are absorbing: accumulate their probability
        // directly, avoiding cancellation in 1 - CDF for small upper tails.
        int width = Math.min(observed, red + 1);
        if ((long) draws * width > maximumUpdates)
            throw new IllegalArgumentException("Wallenius work limit exceeded: draws*threshold="
                + (long) draws * width + "; increase --wallenius-max-updates explicitly");
        double[] state = new double[width];
        state[0] = 1;
        double tail = 0, compensation = 0;
        for (int t = 0; t < draws; t++) {
            for (int x = Math.min(t, width - 1); x >= Math.max(0, t - white); x--) {
                double mass = state[x];
                if (mass == 0) continue;
                int r = red - x, w = white - (t - x);
                double probability = r == 0 ? 0 : w == 0 ? 1
                    : odds <= 1 ? (r * odds) / (r * odds + w)
                    : r / (r + w / odds);
                double success = mass * probability;
                state[x] = mass * (1 - probability);
                if (x + 1 == observed) {
                    double y = success - compensation;
                    double sum = tail + y;
                    compensation = (sum - tail) - y;
                    tail = sum;
                } else if (x + 1 < width) state[x + 1] += success;
            }
        }
        if (tail == 0 || !Double.isFinite(tail))
            throw new ArithmeticException("Wallenius tail underflow or numerical failure; no fallback used");
        return Math.min(1, tail);
    }
}
