/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;

import java.util.List;
import java.util.Objects;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;

/** A pedigree incidence term paired with its sparse additive precision. */
public record PedigreeRandomEffectTerm(
        RandomEffectTerm randomEffect,
        SparsePrecisionMatrix precision) {

    public PedigreeRandomEffectTerm {
        Objects.requireNonNull(randomEffect, "randomEffect");
        Objects.requireNonNull(precision, "precision");
        if (randomEffect.coefficients() != precision.dimension())
            throw new IllegalArgumentException(
                "pedigree random-effect and precision dimensions must match");
    }

    /** Creates a pedigree term, retaining unobserved ancestors as coefficients. */
    public static PedigreeRandomEffectTerm of(
            String name, List<String> observationIndividualIds,
            Pedigree pedigree) {
        if (pedigree == null || observationIndividualIds == null
                || observationIndividualIds.isEmpty())
            throw new IllegalArgumentException(
                "pedigree and observation identifiers are required");
        int rows = observationIndividualIds.size();
        int[] rowStarts = new int[rows + 1];
        int[] columns = new int[rows];
        double[] values = new double[rows];
        for (int row = 0; row < rows; row++) {
            rowStarts[row] = row;
            columns[row] = pedigree.indexOf(observationIndividualIds.get(row));
            values[row] = 1.0;
        }
        rowStarts[rows] = rows;
        RandomEffectTerm term = RandomEffectTerm.ofSparseCsr(name, rows,
            pedigree.size(), rowStarts, columns, values, pedigree.individualIds());
        SparseSymmetricMatrix inverse = pedigree.sparseRelationshipMatrixInverse();
        return new PedigreeRandomEffectTerm(term, new SparsePrecisionMatrix(
            inverse.dimension(), inverse.rowPointers(),
            inverse.columnIndices(), inverse.values()));
    }
}
