/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.survival;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.genetics.GenomicRelationshipMatrix;

/** Cox Gaussian frailty driven by an empirical kinship or genomic GRM. */
public final class CoxKinshipFrailty {
    private static final double DEFAULT_RELATIVE_REGULARIZATION = 1e-8;

    private CoxKinshipFrailty() { }

    public static CoxMixedResult fit(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<String> observationSampleIds,
            GenomicRelationshipMatrix relationship) {
        return fit(survival, fixedEffects, observationSampleIds,
            relationship, null, CoxMixedOptions.defaults(),
            DEFAULT_RELATIVE_REGULARIZATION, BackendPolicy.PREFERRED);
    }

    public static CoxMixedResult fit(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<String> observationSampleIds,
            GenomicRelationshipMatrix relationship,
            double[] offset,
            CoxMixedOptions options,
            double relativeDiagonalRegularization,
            BackendPolicy backendPolicy) {
        if (survival == null || relationship == null
                || observationSampleIds == null
                || observationSampleIds.size() != survival.observations())
            throw new IllegalArgumentException(
                "one GRM sample ID is required per survival row");
        int rows = survival.observations();
        int samples = relationship.samples();
        double[][] incidence = new double[rows][samples];
        for (int row = 0; row < rows; row++) {
            String id = observationSampleIds.get(row);
            if (id == null)
                throw new IllegalArgumentException(
                    "survival sample IDs must not be null");
            incidence[row][relationship.indexOf(id)] = 1;
        }
        CoxRandomEffectTerm kinship = CoxRandomEffectTerm.fromCovariance(
            "kinship", incidence, relationship.relationshipMatrix(),
            relationship.sampleIds(), relativeDiagonalRegularization,
            backendPolicy);
        return CoxMixedModel.fit(survival, fixedEffects,
            List.of(kinship), offset, options, backendPolicy);
    }
}
