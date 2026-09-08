/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.concurrent.ConcurrentHashMap;

/** Physicists' Hermite rules: integral exp(-x*x) f(x) dx. */
final class AdaptiveHermiteRule {
    private static final ConcurrentHashMap<Integer, AdaptiveHermiteRule> CACHE = new ConcurrentHashMap<>();
    final double[] nodes;
    final double[] logWeights;
    static AdaptiveHermiteRule of(int n) { return CACHE.computeIfAbsent(n, AdaptiveHermiteRule::new); }

    private AdaptiveHermiteRule(int n) {
        nodes = new double[n]; logWeights = new double[n];
        // Golub-Welsch via implicit QL on the symmetric tridiagonal Jacobi
        // matrix. Only the first eigenvector row is needed for the weights.
        double[] e = new double[n], first = new double[n]; first[0] = 1;
        for (int i = 0; i < n - 1; i++) e[i] = Math.sqrt((i + 1.0) / 2);
        for (int l = 0; l < n; l++) {
            int iterations = 0;
            while (true) {
                int m = l;
                while (m < n - 1 && Math.abs(e[m]) > 2e-16 * (Math.abs(nodes[m]) + Math.abs(nodes[m + 1]))) m++;
                if (m == l) break;
                if (++iterations > 100) throw new ArithmeticException("Hermite eigensolver failed to converge");
                double g = (nodes[l + 1] - nodes[l]) / (2 * e[l]);
                double r = Math.hypot(g, 1);
                g = nodes[m] - nodes[l] + e[l] / (g + Math.copySign(r, g));
                double s = 1, c = 1, p = 0;
                int i;
                for (i = m - 1; i >= l; i--) {
                    double f = s * e[i], b = c * e[i];
                    r = Math.hypot(f, g); e[i + 1] = r;
                    if (r == 0) { nodes[i + 1] -= p; e[m] = 0; break; }
                    s = f / r; c = g / r; g = nodes[i + 1] - p;
                    r = (nodes[i] - g) * s + 2 * c * b;
                    p = s * r; nodes[i + 1] = g + p; g = c * r - b;
                    f = first[i + 1]; first[i + 1] = s * first[i] + c * f;
                    first[i] = c * first[i] - s * f;
                }
                if (r == 0 && i >= l) continue;
                nodes[l] -= p; e[l] = g; e[m] = 0;
            }
        }
        for (int i = 0; i < n; i++) logWeights[i] = .5 * Math.log(Math.PI) + 2 * Math.log(Math.abs(first[i]));
    }
}
