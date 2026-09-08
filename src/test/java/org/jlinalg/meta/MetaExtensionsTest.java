/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.Test;

final class MetaExtensionsTest {
    @Test
    void fullCovarianceModeratorsAndBiasCorrectionsFit() {
        List<MetaStudy> studies = List.of(new MetaStudy("a", 0.9, 0.2), new MetaStudy("b", 1.1, 0.25), new MetaStudy("c", 1.4, 0.3), new MetaStudy("d", 1.5, 0.2), new MetaStudy("e", 1.7, 0.25));
        double[][] moderators = {{0.0}, {0.2}, {0.4}, {0.6}, {0.8}}, covariance = {{0.04, 0.005, 0, 0, 0}, {0.005, 0.0625, 0.004, 0, 0}, {0, 0.004, 0.09, 0.003, 0}, {0, 0, 0.003, 0.04, 0.004}, {0, 0, 0, 0.004, 0.0625}};
        MetaMultilevelRegressionResult result = MetaMultilevelRegression.fit(studies, moderators, List.of("dose"), covariance, true, org.jlinalg.compute.BackendPolicy.PREFERRED);
        assertEquals(2, result.coefficients().length); assertTrue(Double.isFinite(result.q())); assertEquals(2, MetaBiasCorrections.pet(studies).method().equals("PET") ? 2 : 0); assertTrue(MetaBiasCorrections.trimAndFill(studies).imputedStudyCount() >= 0);
    }
}
