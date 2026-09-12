/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;

class SummaryScoreModelsTest {
    @Test void conditionalMatchesHandSchurComplementAndRejectsSingularBlock() {
        var joint=new SetTestScoreState(new double[]{2,3,4},new double[]{5,1,2,1,6,1,2,1,4},3);
        var adjusted=SummaryScoreModels.condition(joint,new int[]{0,1},new int[]{2});
        assertArrayEquals(new double[]{0,2},adjusted.scores(),1e-13);
        assertArrayEquals(new double[]{4,.5,.5,5.75},adjusted.information(),1e-13);
        assertThrows(IllegalArgumentException.class,()->SummaryScoreModels.condition(joint,new int[]{0},new int[]{0}));
        var singular=new SetTestScoreState(new double[]{1,2,2},new double[]{2,0,0,0,1,1,0,1,1},3);
        assertThrows(IllegalArgumentException.class,()->SummaryScoreModels.condition(singular,new int[]{0},new int[]{1,2}));
    }
    @Test void heterogeneousKernelKeepsOpposingCohortEffects() {
        var a=new SetTestScoreState(new double[]{3,0},new double[]{1,0,0,1},2);
        var b=new SetTestScoreState(new double[]{-3,0},new double[]{1,0,0,1},2);
        var het=SummaryScoreModels.heterogeneous(List.of(a,b));
        var r=SummarySetTests.skat("g",het,new double[]{1,1,1,1});
        assertEquals(18,r.statistic());
        // chi-square(4) survival at 18 = exp(-9)*(1+9).
        assertEquals(10*Math.exp(-9),r.pValue(),1e-15);
        assertEquals(1,SummarySetTests.skat("g",ScoreMetaAnalysis.pool(new double[][]{a.scores(),b.scores()},new double[][]{a.information(),b.information()},1).state(),new double[]{1,1}).pValue());
    }
    @Test void noncentralSeriesMatchesIndependentRChiSquare()throws Exception {
        var rows=Files.readAllLines(Path.of("src/test/resources/raremetal/noncentral-tail.tsv"));
        for(String row:rows.subList(1,rows.size())) {
            String[] f=row.split("\t");double expected=Double.parseDouble(f[2]);
            double actual=SkatODeterministicDistribution.tail(Double.parseDouble(f[0]),new double[]{1},new double[]{Double.parseDouble(f[1])},1,Math.max(1e-300,expected*1e-11));
            assertEquals(expected,actual,Math.max(1e-300,expected*2e-9),row);
        }
    }
    @Test void deterministicSkatoMatchesIndependentPolarIntegrationIncludingRareTails()throws Exception {
        var rows=Files.readAllLines(Path.of("src/test/resources/raremetal/advanced-skato.tsv"));
        var defaults=SetTestOptions.defaults();
        var options=new SetTestOptions(defaults.variantFilter(),defaults.missingPolicy(),defaults.skatORhoGrid(),0,0,SkatOCalibration.DETERMINISTIC);
        for(String row:rows.subList(1,rows.size())) {
            String[] f=row.split("\t");double expected=Double.parseDouble(f[4]);
            double[] scores=numbers(f[1]),weights=new double[scores.length];Arrays.fill(weights,1);
            var r=SummarySetTests.skatO(f[0],new SetTestScoreState(scores,numbers(f[2]),scores.length),weights,options);
            assertEquals(Double.parseDouble(f[3]),r.minimumComponentPValue(),Double.parseDouble(f[3])*1e-7,f[0]+" component");
            assertEquals(expected,r.adjustedPValue(),expected*3e-6,f[0]+" adjusted");
            assertEquals(0,r.simulations());
        }
    }
    @Test void vtMatchesIndependentCorrelatedPolarProbability()throws Exception {
        var state=new SetTestScoreState(new double[]{2,-1},new double[]{2,.3,.3,1},2);
        var r=SummaryScoreModels.variableThreshold("g",state,new double[]{1,1},new double[]{.01,.03},200000,9843);
        double expected=Double.parseDouble(Files.readString(Path.of("src/test/resources/raremetal/advanced-vt.txt")).trim());
        assertEquals(expected,r.adjustedPValue(),5*r.monteCarloStandardError());
        assertEquals(.01,r.selectedMaf());assertEquals(1,r.selectedBurden().beta());
        assertTrue(r.adjustedPValue()>r.selectedBurden().pValue());
        var same=SummaryScoreModels.variableThreshold("g",state,new double[]{1,1},new double[]{.01,.01},100,1);
        assertEquals(0,same.simulations());assertEquals(same.selectedBurden().pValue(),same.adjustedPValue());
        var zero=new SetTestScoreState(new double[]{1,-1},new double[]{1,-1,-1,1},2);
        assertNull(SummaryScoreModels.variableThreshold("zero",zero,new double[]{1,1},new double[]{.01,.01},100,1).selectedBurden());
    }
    @Test void noncentralMixtureMatchesIndependentNormalConvolution()throws Exception {
        var rows=Files.readAllLines(Path.of("src/test/resources/raremetal/noncentral-mixture.tsv"));
        for(String row:rows.subList(1,rows.size())) {
            String[] f=row.split("\t");double expected=Double.parseDouble(f[1]);
            assertEquals(expected,SkatODeterministicDistribution.tail(Double.parseDouble(f[0]),new double[]{.8,1.2},new double[]{.5,2},1,expected*1e-11),expected*1e-7,row);
        }
    }
    private static double[] numbers(String s){return Arrays.stream(s.split(",")).mapToDouble(Double::parseDouble).toArray();}
}
