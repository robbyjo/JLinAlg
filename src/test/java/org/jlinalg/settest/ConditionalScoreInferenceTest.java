/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

class ConditionalScoreInferenceTest {
    private static final String TARGET = "1:100:A:G";
    private static final String CONDITION = "1:200:C:T";

    @Test void alignsAlleleSwapsConditionsWithinCohortThenPools() {
        var first = study("one", "null-one", model("binomial", List.of()),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        // Both source dosages are reversed. Scores change sign, while the
        // cross-covariance keeps its sign because both axes are reversed.
        var second = study("two", "null-two", model("binomial", List.of()),
            List.of("1:100:G:A", "1:200:T:C"), new double[]{-3,-2},
            new double[]{5,2,2,4});
        var result = ConditionalScoreInference.condition(List.of(first,second),
            List.of(TARGET), List.of(CONDITION));
        assertArrayEquals(new double[]{3.5}, result.state().scores(), 1e-14);
        assertArrayEquals(new double[]{7.5}, result.state().information(), 1e-14);
        assertArrayEquals(new int[]{2}, result.cohortCounts());
        assertArrayEquals(new String[]{"++"}, result.directions());
        assertEquals(List.of("one","two"), result.cohorts());
        var association = SummarySetTests.singleVariant(TARGET,
            result.state().scores()[0], result.state().information()[0]);
        assertEquals(3.5/7.5, association.beta(), 1e-14);
    }

    @Test void modelContractsAndCohortIdentityAreCheckedBeforePooling() {
        var first = study("one", "null-one", model("binomial", List.of()),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        var incompatible = study("two", "null-two", model("poisson", List.of()),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        IllegalArgumentException mismatch = assertThrows(
            IllegalArgumentException.class, () -> ConditionalScoreInference.condition(
                List.of(first,incompatible), List.of(TARGET), List.of(CONDITION)));
        assertTrue(mismatch.getMessage().contains("family"));

        var duplicate = study("two", "null-one", model("binomial", List.of()),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        IllegalArgumentException repeated = assertThrows(
            IllegalArgumentException.class, () -> ConditionalScoreInference.condition(
                List.of(first,duplicate), List.of(TARGET), List.of(CONDITION)));
        assertTrue(repeated.getMessage().contains("double-count"));
    }

    @Test void fittedConditioningSetsAlignSwapsButDifferentSetsReject() {
        var first = study("one", "null-one",
            model("binomial", List.of("2:10:A:C")),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        var swapped = study("two", "null-two",
            model("binomial", List.of("2:10:C:A")),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        assertDoesNotThrow(() -> ConditionalScoreInference.condition(
            List.of(first,swapped), List.of(TARGET), List.of(CONDITION)));
        var exactWithoutIntercept = study("exact-no-intercept", "null-exact-no-intercept",
            model("binomial", List.of("2:10:A:C"), false),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        assertDoesNotThrow(() -> ConditionalScoreInference.condition(
            List.of(exactWithoutIntercept), List.of(TARGET), List.of(CONDITION)));
        var swappedWithoutIntercept = study("swap-no-intercept", "null-swap-no-intercept",
            model("binomial", List.of("2:10:C:A"), false),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        IllegalArgumentException orientation = assertThrows(
            IllegalArgumentException.class, () -> ConditionalScoreInference.condition(
                List.of(exactWithoutIntercept, swappedWithoutIntercept),
                List.of(TARGET), List.of(CONDITION)));
        assertTrue(orientation.getMessage().contains("fitted conditioning set"));
        var other = study("three", "null-three",
            model("binomial", List.of("2:11:A:C")),
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        assertThrows(IllegalArgumentException.class,
            () -> ConditionalScoreInference.condition(List.of(first,other),
                List.of(TARGET), List.of(CONDITION)));
    }

    @Test void missingAllelesAndRankDeficientConditioningInformationReject() {
        var missing = study("one", "null-one", model("binomial", List.of()),
            List.of("1:100:A:C", CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        IllegalArgumentException absent = assertThrows(
            IllegalArgumentException.class, () -> ConditionalScoreInference.condition(
                List.of(missing), List.of(TARGET), List.of(CONDITION)));
        assertTrue(absent.getMessage().contains("lacks requested"));

        var singular = study("one", "null-one", model("binomial", List.of()),
            List.of(TARGET, CONDITION, "1:300:G:T"), new double[]{1,1,1},
            new double[]{2,0,0, 0,1,1, 0,1,1});
        IllegalArgumentException rank = assertThrows(
            IllegalArgumentException.class, () -> ConditionalScoreInference.condition(
                List.of(singular), List.of(TARGET),
                List.of(CONDITION,"1:300:G:T")));
        assertTrue(rank.getMessage().contains("singular or ill-conditioned"));
    }

    @Test void refAltSwapRequiresTranslationInvariantNullScore() {
        var withoutIntercept = model("binomial", List.of(), false);
        var exact = study("exact", "null-exact", withoutIntercept,
            List.of(TARGET, CONDITION), new double[]{2,1},
            new double[]{4,1,1,2});
        assertDoesNotThrow(() -> ConditionalScoreInference.condition(
            List.of(exact), List.of(TARGET), List.of(CONDITION)));
        var reversed = study("one", "null-one", withoutIntercept,
            List.of("1:100:G:A", CONDITION), new double[]{-2,1},
            new double[]{4,-1,-1,2});
        IllegalArgumentException failure = assertThrows(
            IllegalArgumentException.class, () ->
                ConditionalScoreInference.condition(List.of(reversed),
                    List.of(TARGET), List.of(CONDITION)));
        assertTrue(failure.getMessage().contains("no intercept"));
    }

    private static ConditionalScoreStudy study(String name, String nullId,
            ConditionalScoreStudy.ModelContract contract, List<String> keys,
            double[] scores, double[] covariance) {
        return new ConditionalScoreStudy(name, nullId, "none", contract, keys,
            new SetTestScoreState(scores, covariance, scores.length));
    }

    private static ConditionalScoreStudy.ModelContract model(String family,
            List<String> fittedConditions) {
        return model(family, fittedConditions, true);
    }

    private static ConditionalScoreStudy.ModelContract model(String family,
            List<String> fittedConditions, boolean intercept) {
        List<String> covariates = new ArrayList<>();
        if (intercept) covariates.add("(Intercept)");
        covariates.add("age");
        covariates.addAll(fittedConditions);
        return new ConditionalScoreStudy.ModelContract(
            "jlinalg-conditional-score-v1", "GRCh38", "glm", family,
            "trait ~ age + <omics>", covariates,
            "as-specified-in-formula", "1", "0",
            "fixed complete phenotype/covariate cases aligned to genotype source",
            "additive dosage",
            "analysis-sample mean imputation",
            "unstandardized; includes inverse fitted-null dispersion for GLMs",
            "model-based; nuisance covariates projected out", "none", "",
            "normal approximation only; no SPA or finite-sample guarantee",
            fittedConditions);
    }
}
