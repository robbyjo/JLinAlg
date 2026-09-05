/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

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

    /**
     * Fits through the sparse pedigree-precision kernel.
     *
     * <p>This path supports one-stratum right-censored data with distinct event
     * times. Repeated observations may map to the same pedigree individual.
     * It uses a diagonal approximation to the profiled random-effect
     * information while retaining the exact penalized score.</p>
     */
    public static CoxPedigreeResult fitSparse(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree) {
        return fitSparse(survival, fixedEffects, observationIndividualIds,
            pedigree, null, CoxMixedOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    /** Fits a sparse pedigree Cox model with explicit controls. */
    public static CoxPedigreeResult fitSparse(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            double[] offset,
            CoxMixedOptions options,
            BackendPolicy backendPolicy) {
        validate(survival, observationIndividualIds, pedigree);
        PedigreeRandomEffectTerm term = PedigreeRandomEffectTerm.of(
            "pedigree", observationIndividualIds, pedigree);
        return fitSparse(survival, fixedEffects, term, offset, options,
            backendPolicy);
    }

    /**
     * Fits a fully sparse pedigree term created with
     * {@link PedigreeRandomEffectTerm#ofSparse} or
     * {@link PedigreeRandomEffectTerm#ofUninbred}.
     */
    public static CoxPedigreeResult fitSparse(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            PedigreeRandomEffectTerm pedigreeEffect) {
        return fitSparse(survival, fixedEffects, pedigreeEffect, null,
            CoxMixedOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a fully sparse pedigree term with explicit controls. */
    public static CoxPedigreeResult fitSparse(
            CoxSurvivalData survival,
            double[][] fixedEffects,
            PedigreeRandomEffectTerm pedigreeEffect,
            double[] offset,
            CoxMixedOptions options,
            BackendPolicy backendPolicy) {
        if (pedigreeEffect == null)
            throw new IllegalArgumentException(
                "sparse pedigree effect is required");
        CoxMixedResult result = SparseCoxMixedModel.fit(survival,
            fixedEffects, pedigreeEffect, List.of(), offset, options,
            backendPolicy);
        return new CoxPedigreeResult(result,
            pedigreeEffect.randomEffect().coefficientNames(),
            pedigreeEffect.randomEffect().name());
    }

    private static void validate(
            CoxSurvivalData survival,
            List<String> observationIndividualIds,
            Pedigree pedigree) {
        if (survival == null || pedigree == null
                || observationIndividualIds == null
                || observationIndividualIds.size() != survival.observations())
            throw new IllegalArgumentException(
                "one pedigree individual is required per survival row");
        for (String id : observationIndividualIds)
            if (id == null)
                throw new IllegalArgumentException(
                    "pedigree survival identifiers must not be null");
    }
}
