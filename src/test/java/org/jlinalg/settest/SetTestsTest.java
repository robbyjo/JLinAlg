/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import jdistlib.ChiSquare;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gwas.RemlAssociationScanner;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.pipeline.VariantFilterOptions;
import org.jlinalg.pipeline.VariantRecord;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class SetTestsTest {
    private static final double[][] INTERCEPT = {
        {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
    };
    private static final double[] FIRST =
        {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
    private static final double[] SECOND =
        {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    private static final double[] RESPONSE =
        {1.0, 1.3, 1.1, 1.4, 2.5, 2.9, 2.6, 3.0, 4.1, 4.6, 4.2, 4.7};

    @Test
    void burdenReturnsEffectUncertaintyAndAlleleOrientation() {
        LinearSetTestNullModel nullModel = nullModel();
        VariantSet alternate = new VariantSet("gene",
            List.of(member("v1", FIRST, EffectAllele.ALTERNATE)));
        VariantSet reference = new VariantSet("gene-ref",
            List.of(member("v1", FIRST, EffectAllele.REFERENCE)));

        SetTestResult alt = SetTests.burden(
            alternate, nullModel, options(100));
        SetTestResult ref = SetTests.burden(
            reference, nullModel, options(100));

        assertTrue(alt.beta() > 1);
        assertTrue(alt.standardError() > 0);
        assertTrue(alt.pValue() < 1e-5);
        assertEquals(-alt.beta(), ref.beta(), 1e-12);
        assertEquals(alt.pValue(), ref.pValue(), 1e-12);
    }

    @Test
    void skatAndSkatOAreFiniteAndSeedReproducible() {
        VariantSet set = new VariantSet("gene", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("v2", SECOND, EffectAllele.ALTERNATE)));
        SetTestOptions options = options(250);

        SetTestResult skat = SetTests.skat(set, nullModel(), options);
        SkatOResult first = SetTests.skatO(set, nullModel(), options);
        SkatOResult second = SetTests.skatO(set, nullModel(), options);

        assertTrue(skat.statistic() > 0);
        assertTrue(skat.pValue() > 0 && skat.pValue() <= 1);
        assertTrue(skat.eigenvalues().length > 0);
        assertEquals(3, first.components().size());
        assertEquals(first.adjustedPValue(), second.adjustedPValue(), 0);
        assertEquals(250, first.simulations());
    }

    @Test
    void quadraticFormUsesExactScaledChiSquareForOneEigenvalue() {
        QuadraticFormDistribution.Tail result =
            QuadraticFormDistribution.survival(3.5, new double[] {2});
        assertEquals(ChiSquare.cumulative(3.5 / 2, 1, false, false),
            result.pValue(), 1e-15);
        assertEquals("exact-scaled-chi-square", result.method());
    }

    @Test
    void setMembershipFilteringIsAudited() {
        VariantSet set = new VariantSet("filtered", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("singleton",
                new double[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                EffectAllele.ALTERNATE)));
        SetTestOptions options = new SetTestOptions(
            VariantFilterOptions.builder().minimumMac(2).build(),
            SetTestMissingPolicy.MEAN_IMPUTE,
            new double[] {0, 1}, 20, 7);

        SetTestResult result = SetTests.skat(set, nullModel(), options);

        assertEquals(1, result.includedVariants());
        assertEquals(1, result.excludedVariants().size());
        assertEquals("singleton",
            result.excludedVariants().get(0).variant().id());
    }

    @Test
    void relatedSampleTestsReuseOneRemlProjection() {
        double[] kinship = new double[RESPONSE.length * RESPONSE.length];
        for (int row = 0; row < RESPONSE.length; row++)
            for (int column = 0; column < RESPONSE.length; column++)
                kinship[row * RESPONSE.length + column] =
                    Math.pow(0.2, Math.abs(row - column));
        RemlAssociationScanner scanner = RemlAssociationScanner.prepare(
            RESPONSE, INTERCEPT,
            List.of(
                new VarianceComponent(
                    "cryptic-relatedness", RESPONSE.length, kinship),
                VarianceComponent.identity("residual", RESPONSE.length)),
            RemlOptions.builder().initialVariances(0.2, 1).build(),
            BackendPolicy.CPU);
        RemlSetTestNullModel nullModel = new RemlSetTestNullModel(scanner);
        VariantSet set = new VariantSet("related-gene", List.of(
            member("v1", FIRST, EffectAllele.ALTERNATE),
            member("v2", SECOND, EffectAllele.ALTERNATE)));

        SetTestResult burden = SetTests.burden(set, nullModel, options(100));
        SetTestResult skat = SetTests.skat(set, nullModel, options(100));
        SkatOResult skatO = SetTests.skatO(set, nullModel, options(100));

        assertTrue(Double.isFinite(burden.beta()));
        assertTrue(burden.standardError() > 0);
        assertTrue(burden.pValue() > 0 && burden.pValue() <= 1);
        assertTrue(skat.statistic() > 0);
        assertTrue(skat.pValue() > 0 && skat.pValue() <= 1);
        assertTrue(skatO.adjustedPValue() > 0
            && skatO.adjustedPValue() <= 1);
    }

    private static LinearSetTestNullModel nullModel() {
        return LinearSetTestNullModel.prepare(
            RESPONSE, INTERCEPT, OlsOptions.defaults(), BackendPolicy.CPU);
    }

    private static SetTestOptions options(int simulations) {
        return new SetTestOptions(VariantFilterOptions.defaults(),
            SetTestMissingPolicy.MEAN_IMPUTE,
            new double[] {0, 0.5, 1}, simulations, 42);
    }

    private static WeightedVariant member(
            String id, double[] dosage, EffectAllele effectAllele) {
        return new WeightedVariant(new VariantRecord(
            id, "1", 1, "A", "G", dosage, Double.NaN),
            effectAllele, 1);
    }
}
