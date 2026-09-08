/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import static org.junit.jupiter.api.Assertions.*;
import java.util.List;
import org.junit.jupiter.api.Test;

class ConditionalAssociationTest {
    private static ConditionalAssociation a(String id, double beta, double se) {
        return new ConditionalAssociation(id, beta, se, "locus");
    }

    @Test void twoSnpRegressionUsesScaledEffectsAndVarianceInflation() {
        var model = new ConditionalAssociationModel(List.of(a("lead", .5, .05), a("tag", .8, .1)),
            new double[][] {{1,.8},{.8,1}});
        var tag = model.condition(new int[] {0}, .95).get(1);
        assertEquals(0, tag.effect(), 1e-14);
        assertEquals(.1/.6, tag.standardError(), 1e-14);
        assertEquals(1, tag.pValue(), 1e-14);
        assertEquals(List.of("lead"), tag.conditionedOn());
        assertEquals(.36, tag.residualGenotypeVariance(), 1e-14);
        // A selected SNP is jointly tested against every other selected SNP.
        var joint = model.condition(new int[] {0,1}, .95);
        assertEquals(.5, joint.get(0).effect(), 1e-14);
        assertEquals(0, joint.get(1).effect(), 1e-14);
    }

    @Test void leadChangesRecomputeStatisticsAndRemoveProxyLead() {
        var values = List.of(a("proxy",1.3,.1), a("causal1",1,.1), a("causal2",1,.1));
        double[][] r = {{1,.65,.65},{.65,1,0},{.65,0,1}};
        var selected = SecondarySignalClumper.select(values,r,.05,3,.95);
        assertTrue(selected.converged());
        assertEquals(List.of("add proxy","add causal1","add causal2","remove proxy"), selected.changes());
        assertEquals(List.of("causal1","causal2"),selected.selectedVariantIds());
        assertEquals(0,selected.associations().get(0).effect(),1e-14);
        assertEquals(1,selected.associations().get(0).pValue(),1e-13);
        assertEquals(List.of("causal1","causal2"),selected.associations().get(0).conditionedOn());
    }

    @Test void correlatedTagIsNotSelectedFromItsMarginalSignificance() {
        var selection = SecondarySignalClumper.select(List.of(a("lead",1,.1),a("tag",.8,.1)),
            new double[][]{{1,.8},{.8,1}},5e-8,2,.95);
        assertEquals(List.of("lead"),selection.selectedVariantIds());
        assertEquals(1,selection.associations().get(1).pValue(),1e-13);
    }

    @Test void emptyConditioningAndAlleleChangesAreExplicitAndConsistent() {
        var values = List.of(a("a",.5,.05),a("b",.3,.1));
        var model = new ConditionalAssociationModel(values,new double[][]{{1,.4},{.4,1}});
        assertEquals(.3,model.condition(new int[0],.95).get(1).effect(),1e-14);
        assertTrue(model.condition(new int[0],.95).get(1).conditionedOn().isEmpty());
        var flipped = new ConditionalAssociationModel(List.of(a("a",-.5,.05),a("b",.3,.1)),
            new double[][]{{1,-.4},{-.4,1}});
        assertEquals(model.condition(new int[]{0},.95).get(1), flipped.condition(new int[]{0},.95).get(1));
    }

    @Test void invalidLdAndConditioningAreRejected() {
        var a = List.of(a("a",1,.1),a("b",1,.1));
        for (double[][] r : List.of(new double[][]{{1,1},{1,1}},
                new double[][]{{1,Double.NaN},{.2,1}},new double[][]{{1,.3},{.2,1}},
                new double[][]{{Double.NaN,0},{0,1}}))
            assertThrows(IllegalArgumentException.class,()->new ConditionalAssociationModel(a,r));
        var model = new ConditionalAssociationModel(a,new double[][]{{1,0},{0,1}});
        assertThrows(IllegalArgumentException.class,()->model.condition(new int[]{0,0},.95));
        assertThrows(IllegalArgumentException.class,()->model.condition(new int[]{2},.95));
        assertThrows(IllegalArgumentException.class,()->model.condition(new int[]{0},Double.NaN));
        assertThrows(IllegalArgumentException.class,()->new ConditionalAssociationModel(
            List.of(a("a",1,.1),a("a",2,.1)),new double[][]{{1,0},{0,1}}));
    }
}
