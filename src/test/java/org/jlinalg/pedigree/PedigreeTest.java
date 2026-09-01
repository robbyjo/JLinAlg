/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import org.junit.jupiter.api.Test;

class PedigreeTest {
    @Test
    void buildsRelationshipMatrixInInputOrderAndAccountsForInbreeding() {
        Pedigree pedigree = Pedigree.of(List.of(
            new PedigreeIndividual("C", "S", "D"),
            PedigreeIndividual.founder("S"),
            PedigreeIndividual.founder("D"),
            new PedigreeIndividual("F", "S", "D"),
            new PedigreeIndividual("G", "C", "F")));

        assertEquals(List.of("C", "S", "D", "F", "G"),
            pedigree.individualIds());
        assertArrayEquals(new double[] {
            1.00, 0.50, 0.50, 0.50, 0.75,
            0.50, 1.00, 0.00, 0.50, 0.50,
            0.50, 0.00, 1.00, 0.50, 0.50,
            0.50, 0.50, 0.50, 1.00, 0.75,
            0.75, 0.50, 0.50, 0.75, 1.25
        }, pedigree.relationshipMatrix(), 1e-15);
        assertArrayEquals(new double[] {0.0, 0.0, 0.0, 0.0, 0.25},
            pedigree.inbreedingCoefficients(), 1e-15);
        assertEquals(0.75, pedigree.relationship("G", "C"), 1e-15);
    }

    @Test
    void rejectsDuplicateUnknownAndCyclicAncestry() {
        assertThrows(IllegalArgumentException.class, () -> Pedigree.of(List.of(
            PedigreeIndividual.founder("A"),
            PedigreeIndividual.founder("A"))));
        assertThrows(IllegalArgumentException.class, () -> Pedigree.of(List.of(
            new PedigreeIndividual("A", "missing", null))));
        assertThrows(IllegalArgumentException.class, () -> Pedigree.of(List.of(
            new PedigreeIndividual("A", "B", null),
            new PedigreeIndividual("B", "A", null))));
    }

    @Test
    void relationshipArraysAreDefensiveCopies() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("A")));
        double[] matrix = pedigree.relationshipMatrix();
        matrix[0] = 99.0;
        assertEquals(1.0, pedigree.relationshipMatrix()[0]);
    }

    @Test
    void sparseRelationshipInverseMultipliesRelationshipToIdentity() {
        Pedigree pedigree = Pedigree.of(List.of(
            new PedigreeIndividual("C", "S", "D"),
            PedigreeIndividual.founder("S"),
            PedigreeIndividual.founder("D"),
            new PedigreeIndividual("F", "S", "D"),
            new PedigreeIndividual("G", "C", "F")));
        SparseSymmetricMatrix inverse =
            pedigree.sparseRelationshipMatrixInverse();
        double[] relationship = pedigree.relationshipMatrix();
        double[] inverseDense = inverse.toDense();
        int size = pedigree.size();
        double[] product = new double[size * size];
        for (int row = 0; row < size; row++) {
            for (int column = 0; column < size; column++) {
                for (int shared = 0; shared < size; shared++) {
                    product[row * size + column] +=
                        inverseDense[row * size + shared]
                            * relationship[shared * size + column];
                }
            }
        }

        double[] identity = new double[size * size];
        for (int index = 0; index < size; index++) {
            identity[index * size + index] = 1.0;
        }
        assertArrayEquals(identity, product, 2e-15);
        assertEquals(inverse.get(0, 1), inverse.get(1, 0), 0.0);
    }
}
