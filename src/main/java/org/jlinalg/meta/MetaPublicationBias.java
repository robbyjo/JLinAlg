/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;
import java.util.List;
import jdistlib.Normal;

/** Egger's multiplicative-dispersion t test and an approximate Spearman rank diagnostic.
 * The rank diagnostic is not Begg's Kendall test (metafor ranktest). */
public final class MetaPublicationBias {
    private MetaPublicationBias() { }

    public static MetaPublicationBiasResult diagnose(List<MetaStudy> studies) {
        MetaMath.Data data = MetaMath.data(studies);
        int n = studies.size();
        if (n < 3) throw new IllegalArgumentException("Egger inference requires at least three studies");
        double[] effects = data.effects();
        double[] errors = new double[n];
        for (int i = 0; i < n; i++) errors[i] = Math.sqrt(data.variances()[i]);
        // WLS effect~SE is algebraically the standardized-effect Egger regression,
        // with columns reversed. Use the shared QR fit instead of difference-of-sums normal equations.
        MetaBiasCorrections.Result pet = MetaBiasCorrections.pet(studies);
        double intercept = pet.slope(), interceptSe = pet.standardErrors()[1];
        double statistic = pet.associationStatistics().statistics()[1], p = pet.pValues()[1];
        double rank = spearman(effects, errors);
        double rankStatistic = rank * Math.sqrt(Math.max(0.0, (n - 2.0) / Math.max(1e-15, 1.0 - rank * rank)));
        double rankP = 2.0 * Normal.cumulative(-Math.abs(rankStatistic), 0.0, 1.0, true, false);
        return new MetaPublicationBiasResult(intercept, interceptSe, statistic, p,
            rank, rankStatistic, rankP);
    }

    private static double spearman(double[] x, double[] y) {
        double[] rx = ranks(x), ry = ranks(y); double sx = 0, sy = 0, sxy = 0;
        for (int i = 0; i < x.length; i++) { sx += rx[i]; sy += ry[i]; sxy += rx[i] * ry[i]; }
        double mx = sx / x.length, my = sy / x.length, xx = 0, yy = 0, xy = 0;
        for (int i = 0; i < x.length; i++) { double a = rx[i] - mx, b = ry[i] - my; xx += a*a; yy += b*b; xy += a*b; }
        return xy / Math.sqrt(xx * yy);
    }

    private static double[] ranks(double[] values) {
        Integer[] order = new Integer[values.length];
        for (int i = 0; i < order.length; i++) order[i] = i;
        Arrays.sort(order, (a, b) -> Double.compare(values[a], values[b]));
        double[] result = new double[values.length];
        for (int i = 0; i < order.length;) { int j = i + 1; while (j < order.length && values[order[j]] == values[order[i]]) j++; double rank = 0.5 * (i + j - 1) + 1.0; for (int k = i; k < j; k++) result[order[k]] = rank; i = j; }
        return result;
    }
}
