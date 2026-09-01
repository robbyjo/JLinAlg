/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.reml.RemlOptions;

/** Sparse REML combining one or more pedigree and ordinary random terms. */
public final class SparsePedigreeMixedModel {
    private SparsePedigreeMixedModel() { }

    public static SparseLinearMixedModelResult fit(
            double[] response,
            double[][] fixedEffects,
            List<PedigreeRandomEffectTerm> pedigreeEffects,
            List<RandomEffectTerm> ordinaryEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null)
            throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, fixed, response.length, fixedEffects[0].length,
            pedigreeEffects, ordinaryEffects, options, backendPolicy);
    }

    public static SparseLinearMixedModelResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<PedigreeRandomEffectTerm> pedigreeEffects,
            List<RandomEffectTerm> ordinaryEffects,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (pedigreeEffects == null || pedigreeEffects.isEmpty()
                || ordinaryEffects == null)
            throw new IllegalArgumentException(
                "at least one pedigree term and an ordinary-term list are required");
        List<RandomEffectTerm> terms = new ArrayList<>();
        List<SparsePrecisionMatrix> precisions = new ArrayList<>();
        for (PedigreeRandomEffectTerm value : pedigreeEffects) {
            if (value == null)
                throw new IllegalArgumentException("pedigree terms must not be null");
            terms.add(value.randomEffect());
            precisions.add(value.precision());
        }
        for (RandomEffectTerm value : ordinaryEffects) {
            if (value == null)
                throw new IllegalArgumentException("ordinary terms must not be null");
            terms.add(value);
            precisions.add(SparsePrecisionMatrix.identity(value.coefficients()));
        }
        return SparseLinearMixedModel.fitWithPrecision(response, fixedEffects,
            rows, columns, terms, precisions, options, backendPolicy);
    }
}
