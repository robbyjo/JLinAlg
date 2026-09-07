/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

final class BetaRegressionLinkTest {
    @Test
    void allMeanLinksRoundTripAndExposeAnalyticDerivatives() {
        for (BetaMeanLink link : BetaMeanLink.values()) {
            for (double mean : new double[] {0.1, 0.35, 0.8}) {
                double predictor = link.link(mean);
                assertEquals(mean, link.inverse(predictor), 2e-14,
                    link.name());
                double step = 1e-6;
                double numerical = (link.inverse(predictor + step)
                    - link.inverse(predictor - step)) / (2.0 * step);
                assertEquals(numerical,
                    link.derivative(predictor, mean), 2e-8, link.name());
            }
        }
    }

    @Test
    void allPrecisionLinksRoundTripAndExposeAnalyticDerivatives() {
        for (BetaPrecisionLink link : BetaPrecisionLink.values()) {
            for (double precision : new double[] {0.5, 4.0, 100.0}) {
                double predictor = link.link(precision);
                assertEquals(precision, link.inverse(predictor), 2e-13,
                    link.name());
                double step = 1e-6;
                double numerical = (link.inverse(predictor + step)
                    - link.inverse(predictor - step)) / (2.0 * step);
                double tolerance = 2e-7 * Math.max(1.0, Math.abs(numerical));
                assertEquals(numerical,
                    link.derivative(predictor, precision), tolerance,
                    link.name());
            }
        }
    }

    @Test
    void identityPrecisionRejectsNonpositivePredictors() {
        assertTrue(Double.isNaN(BetaPrecisionLink.IDENTITY.inverse(0.0)));
        assertTrue(Double.isNaN(BetaPrecisionLink.IDENTITY.inverse(-1.0)));
    }
}
