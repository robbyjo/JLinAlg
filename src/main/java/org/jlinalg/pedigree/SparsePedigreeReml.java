/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.reml.RemlOptions;

/** Pedigree REML that never materializes dense A or observation covariance. */
public final class SparsePedigreeReml {
    private SparsePedigreeReml() { }

    public static SparsePedigreeRemlResult fit(
            double[] response,
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null)
            throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, fixed, response.length, fixedEffects[0].length,
            observationIndividualIds, pedigree, options, backendPolicy);
    }

    public static SparsePedigreeRemlResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        if (pedigree == null || observationIndividualIds == null
                || observationIndividualIds.size() != rows)
            throw new IllegalArgumentException(
                "one known pedigree individual is required per observation");
        int[] rowStarts = new int[rows + 1];
        int[] individualColumns = new int[rows];
        double[] values = new double[rows];
        for (int row = 0; row < rows; row++) {
            String id = observationIndividualIds.get(row);
            if (id == null)
                throw new IllegalArgumentException(
                    "observation pedigree identifiers must not be null");
            rowStarts[row] = row;
            individualColumns[row] = pedigree.indexOf(id);
            values[row] = 1.0;
        }
        rowStarts[rows] = rows;
        RandomEffectTerm animal = RandomEffectTerm.ofSparseCsr(
            "additive genetic", rows, pedigree.size(), rowStarts,
            individualColumns, values, pedigree.individualIds());
        SparseSymmetricMatrix inverse =
            pedigree.sparseRelationshipMatrixInverse();
        SparsePrecisionMatrix precision = new SparsePrecisionMatrix(
            inverse.dimension(), inverse.rowPointers(),
            inverse.columnIndices(), inverse.values());
        SparseLinearMixedModelResult fitted =
            SparseLinearMixedModel.fitWithPrecision(response, fixedEffects,
                rows, columns, List.of(animal), List.of(precision),
                options, backendPolicy);
        return new SparsePedigreeRemlResult(fitted,
            pedigree.individualIds(),
            fitted.randomEffects("additive genetic").estimates());
    }
}
