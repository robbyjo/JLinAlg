/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.pedigree.Pedigree;

/** Cox model with an additive pedigree-correlated Gaussian log frailty. */
public final class CoxPedigreeFrailty {
    private CoxPedigreeFrailty() { }

    public static CoxPedigreeResult fit(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree) {
        return fit(survival, fixedEffects, observationIndividualIds,
            pedigree, null, CoxMixedOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    public static CoxPedigreeResult fit(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            double[] offset,
            CoxMixedOptions options,
            BackendPolicy backendPolicy) {
        if (survival == null || pedigree == null
                || observationIndividualIds == null
                || observationIndividualIds.size() != survival.observations())
            throw new IllegalArgumentException(
                "one pedigree individual is required per survival row");
        int rows = survival.observations();
        int animals = pedigree.size();
        double[][] incidence = new double[rows][animals];
        for (int row = 0; row < rows; row++) {
            String id = observationIndividualIds.get(row);
            if (id == null)
                throw new IllegalArgumentException(
                    "pedigree survival identifiers must not be null");
            incidence[row][pedigree.indexOf(id)] = 1;
        }
        CoxRandomEffectTerm term = new CoxRandomEffectTerm(
            "pedigree", incidence,
            pedigree.sparseRelationshipMatrixInverse().toDense(),
            pedigree.individualIds());
        CoxMixedResult result = CoxMixedModel.fit(survival, fixedEffects,
            List.of(term), offset, options, backendPolicy);
        return new CoxPedigreeResult(result, pedigree.individualIds());
    }
}
