/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.multipletesting;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.Arrays;
import org.junit.jupiter.api.Test;

final class MultipleTestingTest {
    @Test void weightedBhMatchesKnownUnweightedExample() {
        double[] adjusted = WeightedBenjaminiHochberg.adjust(
            new double[] {0.01, 0.04, 0.03, 0.20});
        assertArrayEquals(new double[] {0.04, 0.05333333333333334,
            0.05333333333333334, 0.20}, adjusted, 1e-15);
    }

    @Test void ihwWeightsForHeldOutFoldDoNotUseItsPValues() {
        double[] p = {0.001, 0.9, 0.01, 0.8, 0.2, 0.7, 0.3, 0.6};
        double[] covariate = {1, 1, 2, 2, 3, 3, 4, 4};
        int[] folds = {0, 1, 0, 1, 0, 1, 0, 1};
        IndependentHypothesisWeighting.Result first =
            IndependentHypothesisWeighting.adjust(p, covariate, folds,
                new IndependentHypothesisWeighting.Options(2, 0.5, 5));
        double[] changed = p.clone();
        changed[0] = 0.99; changed[2] = 0.99;
        IndependentHypothesisWeighting.Result second =
            IndependentHypothesisWeighting.adjust(changed, covariate, folds,
                new IndependentHypothesisWeighting.Options(2, 0.5, 5));
        for (int index = 0; index < folds.length; index++)
            if (folds[index] == 0)
                assertEquals(first.weights()[index], second.weights()[index], 0.0);
        for (int fold = 0; fold < 2; fold++) {
            double mean = 0.0;
            int count = 0;
            for (int index = 0; index < folds.length; index++)
                if (folds[index] == fold) { mean += first.weights()[index]; count++; }
            assertEquals(1.0, mean / count, 1e-15);
        }
    }

    @Test void hierarchicalBonferroniSpendsACompleteFamilyAndGatesChildren() {
        String[] ids = {"gene", "tissue-a", "tissue-b", "phenotype"};
        String[] parents = {"", "gene", "gene", "tissue-a"};
        double[] p = {0.001, 0.01, 0.9, 0.001};
        HierarchicalTesting.Result result = HierarchicalTesting.adjust(
            ids, parents, p, 0.05);
        assertEquals(1.0, Arrays.stream(result.weights()).sum(), 1e-15);
        assertTrue(result.rejected()[0]);
        assertTrue(result.rejected()[1]);
        assertFalse(result.rejected()[2]);
        assertTrue(result.adjustedPValues()[3]
            >= result.adjustedPValues()[1]);
    }
}
