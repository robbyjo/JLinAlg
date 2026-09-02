/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class GammCovariancesTest {
    @Test
    void alignsRepeatedPedigreeIndividualsToObservations() {
        Pedigree pedigree = Pedigree.of(List.of(
            new PedigreeIndividual("sire", null, null),
            new PedigreeIndividual("dam", null, null),
            new PedigreeIndividual("child", "sire", "dam")));
        VarianceComponent component = GammCovariances.pedigree(
            "animal", pedigree, List.of("child", "sire", "child"));

        assertEquals(3, component.dimension());
        assertArrayEquals(new double[] {
            1.0, 0.5, 1.0,
            0.5, 1.0, 0.5,
            1.0, 0.5, 1.0
        }, component.covariance(), 1e-12);
    }
}
