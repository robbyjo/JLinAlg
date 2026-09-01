/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class GenomicRelationshipMatrixTest {
    @Test
    void acceleratedBuilderFindsDuplicateCrypticPairsAndAuditsVariants() {
        double[][] variants = {
            {0, 0, 2, 2},
            {0, 0, 2, 2},
            {0, 0, 1, 1},
            {Double.NaN, Double.NaN, Double.NaN, 1}
        };
        GenomicRelationshipMatrix grm =
            GenomicRelationshipMatrix.fromVariantDosages(
                variants, List.of("a", "b", "c", "d"),
                new GenomicRelationshipOptions(0, 0.75),
                BackendPolicy.CPU);

        assertEquals(4, grm.variantsConsidered());
        assertEquals(3, grm.variantsUsed());
        assertEquals(1, grm.variantsExcluded());
        assertEquals(grm.relationship("a", "a"),
            grm.relationship("a", "b"), 1e-14);
        assertEquals(0.5 * grm.relationship("a", "b"),
            grm.kinshipCoefficient("a", "b"), 0);
        assertEquals(2, grm.relatedPairs(1).size());
        assertEquals(grm.relatedPairs(1),
            grm.relatedPairsByKinship(0.5));
        assertTrue(grm.computationBackend().isPresent());
        assertArrayEquals(grm.relationshipMatrix(),
            grm.varianceComponent("cryptic").covariance(), 0);
        assertArrayEquals(new double[] {
            grm.relationship("a", "a"), grm.relationship("a", "c"),
            grm.relationship("c", "a"), grm.relationship("c", "c")
        }, grm.varianceComponent(
            "repeated", List.of("a", "c")).covariance(), 0);
    }

    @Test
    void suppliedMatrixRetainsSampleAlignmentAndPairOrdering() {
        GenomicRelationshipMatrix grm = new GenomicRelationshipMatrix(
            List.of("x", "y", "z"), new double[] {
                1, 0.2, 0.5,
                0.2, 1, 0.3,
                0.5, 0.3, 1
            });

        List<RelatednessPair> pairs = grm.relatedPairs(0.25);
        assertEquals(2, pairs.size());
        assertEquals("x", pairs.get(0).firstId());
        assertEquals("z", pairs.get(0).secondId());
        assertEquals(0.5, pairs.get(0).relationship(), 0);
        assertEquals(2, grm.indexOf("z"));
    }
}
