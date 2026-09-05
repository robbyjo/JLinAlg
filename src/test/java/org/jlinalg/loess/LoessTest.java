/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.loess;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class LoessTest {
    private static final double[] X = {-2, -1.7, -1.2, -0.9, -0.4, -0.1,
        0.15, 0.3, 0.55, 0.8, 1.05, 1.3, 1.55, 1.9, 2.2, 2.6, 3, 3.5,
        4, 4.6};
    private static final double[] Y = {-1.1892974268256817,
        -1.2766648104524686, -1.1020390859672262, -0.8783269096274834,
        -0.4694183423086505, -0.08483341664682816, 0.1619381324735992,
        0.36052020666133955, 0.5651872289306592, 0.8473560908995228,
        1.054923225594017, 1.138558185417193, 1.242283764189357,
        1.2013000876874145, 1.1784964038195902, 0.8955013718214642,
        0.6111200080598672, 0.15421677231038016, -0.1468024953079282,
        -0.3136910036334644};

    @Test
    void weightedDirectFitAndPredictionMatchRStatsLoess() {
        double[] weights = {1, 2, 1, 1, 1.5, 1, 1, 2, 1, 1, 1, 1, 1.5,
            1, 1, 1, 2, 1, 1, 1};
        LoessResult fit = Loess.fit(X, Y, weights,
            new LoessOptions(0.6, 2, LoessFamily.GAUSSIAN, 4));
        double[] expected = {-1.2681971672704513, -1.233835715479825,
            -1.0487740567025654, -0.8632119757772718, -0.424938764055021,
            -0.09486557649006683, 0.1806014493964511, 0.34306401488420324,
            0.6034905990743569, 0.8308088232892525, 1.0196853178173417,
            1.1503286897419351, 1.2255564344925807, 1.2297105554740047,
            1.1409466081647144, 0.909944832319317, 0.6285897952404388,
            0.3103072505647049, -0.019017646385423796, -0.42373906733944733};
        assertArrayEquals(expected, fit.fittedValues(), 2e-12);
        assertArrayEquals(new double[] {-1.2577014098041226,
            -0.9312719810278446, 0.23354286230536264, 0.9839383020868171,
            1.0322877218176159, -0.557256827032348},
            fit.predict(new double[] {-2.2, -1, 0.2, 1, 2.4, 4.8}), 2e-12);
        assertEquals(6.413160526789568,
            fit.effectiveDegreesOfFreedom(), 2e-12);
    }

    @Test
    void symmetricRobustFitMatchesRAndRejectsInvalidInputs() {
        double[] contaminated = Y.clone();
        contaminated[9] += 3.0;
        LoessResult fit = Loess.fit(X, contaminated, null,
            new LoessOptions(0.6, 2, LoessFamily.SYMMETRIC, 4));
        double[] expected = {-1.2454577724060452, -1.220832143350492,
            -1.0477717106478122, -0.8669455630400865, -0.42296156149513453,
            -0.09342083119577294, 0.17087477357034966, 0.32941844085252214,
            0.5935218303339806, 0.8291191744866804, 1.0172103042969307,
            1.1484806729977661, 1.2251410953871473, 1.2291195144269553,
            1.1473418048850523, 0.9322130400789683, 0.7270784314730654,
            0.4388220081416173, 0.10614251023139602, -0.33941633183539544};
        assertArrayEquals(expected, fit.fittedValues(), 2e-10);
        assertEquals(0.0, fit.robustnessWeights()[9], 0.0);
        assertTrue(fit.robustnessWeights()[0] > 0.0);
        assertArrayEquals(new double[] {0.8846113930058513,
            0.8747985633639853, 0.8854478632115335, 0.995451906274354,
            0.9198309656064613, 0.9945094395653256, 0.9984877782668423,
            0.9576424972019014, 0.9693464317234725, 0.0,
            0.9422598921498084, 0.9961961340232988, 0.9875885284900472,
            0.9731996862286254, 0.9480978368955658, 0.9656011674101068,
            0.7855963453304391, 0.0, 0.06311268136331559,
            0.7556674916389116}, fit.robustnessWeights(), 2e-10);
        assertThrows(IllegalArgumentException.class,
            () -> Loess.fit(new double[] {1, 1}, new double[] {1, 2}));
        assertThrows(IllegalArgumentException.class,
            () -> Loess.fit(X, new double[X.length - 1]));
        assertThrows(IllegalArgumentException.class,
            () -> Loess.prepare(X, LoessOptions.defaults().withSpan(1.1)));
    }
}
