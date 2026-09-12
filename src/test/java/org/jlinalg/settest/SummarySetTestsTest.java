/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
class SummarySetTestsTest {
    @Test void independentRQuadratureAndSkatoFixtures()throws Exception {
        var rows=java.nio.file.Files.readAllLines(java.nio.file.Path.of("src/test/resources/raremetal/independent-scores.tsv"));
        for(String row:rows.subList(1,rows.size())) {
            String[] f=row.split("\t");
            double[] u=java.util.Arrays.stream(f[1].split(",")).mapToDouble(Double::parseDouble).toArray();
            double[] v=java.util.Arrays.stream(f[2].split(",")).mapToDouble(Double::parseDouble).toArray();
            var state=new SetTestScoreState(u,v,u.length);double[] w={1,1,1};
            var burden=SummarySetTests.burden(f[0],state,w);
            assertEquals(Double.parseDouble(f[3]),burden.beta(),1e-13);
            var skat=SummarySetTests.skat(f[0],state,w);
            double expected=Double.parseDouble(f[11]);
            assertEquals(expected,skat.pValue(),expected*1e-8,f[0]+" spherical tail");
            if(!f[0].equals("tail")) {
                var defaults=SetTestOptions.defaults();
                var options=new SetTestOptions(defaults.variantFilter(),defaults.missingPolicy(),defaults.skatORhoGrid(),200000,1234);
                var result=SummarySetTests.skatO(f[0],state,w,options);
                assertEquals(Double.parseDouble(f[10]),result.adjustedPValue(),.015,f[0]+" R SKAT-O vs score-null simulation");
            }
        }
    }
    @Test void analyticBurdenAndRankOneKernel() {
        var state=new SetTestScoreState(new double[]{2,-3},new double[]{4,1,1,9},2);
        var b=SummarySetTests.burden("gene",state,new double[]{1,2});
        assertEquals(-4.0/44,b.beta(),1e-15);assertEquals(1/Math.sqrt(44),b.standardError(),1e-15);
        var single=SummarySetTests.singleVariant("v",8,4);
        assertEquals(6.334248366623996e-5,single.pValue(),1e-16);
        var skat=SummarySetTests.skat("v",new SetTestScoreState(new double[]{8},new double[]{4},1),new double[]{2});
        assertEquals(single.pValue(),skat.pValue(),1e-15);
        var skato=SummarySetTests.skatO("v",new SetTestScoreState(new double[]{8},new double[]{4},1),new double[]{2},SetTestOptions.defaults());
        assertEquals(single.pValue(),skato.adjustedPValue(),1e-15);
    }
    @Test void invalidCovarianceAndMissingnessFailExplicitly() {
        assertThrows(IllegalArgumentException.class,()->SummarySetTests.burden("g",new SetTestScoreState(new double[]{1,2},new double[]{1,2,2,1},2),new double[]{1,1}));
        assertThrows(IllegalArgumentException.class,()->SummarySetTests.validate(new SetTestScoreState(new double[]{1,-1},new double[]{1,1,1,1},2)));
        var p=ScoreMetaAnalysis.pool(new double[][]{{2,Double.NaN},{3,4}},new double[][]{{4,0,0,0},{9,1,1,16}},2);
        assertArrayEquals(new int[]{0},p.indices());assertArrayEquals(new double[]{5},p.state().scores());assertArrayEquals(new double[]{13},p.state().information());assertEquals("++",p.directions()[0]);
        assertThrows(IllegalArgumentException.class,()->ScoreMetaAnalysis.pool(new double[][]{{1,1},{1,1}},
            new double[][]{{1,2,2,1},{10,0,0,10}},1));
        var cancelled=ScoreMetaAnalysis.pool(new double[][]{{1e16},{1},{-1e16}},new double[][]{{1},{1},{1}},1);
        assertEquals(1,cancelled.state().scores()[0]);
    }
    @Test void zeroBurdenAndRankOneSkatoRemainDefined() {
        var state=new SetTestScoreState(new double[]{1,-1},new double[]{1,-1,-1,1},2);
        var burden=SummarySetTests.burden("g",state,new double[]{1,1});
        assertEquals(1,burden.pValue());assertTrue(Double.isNaN(burden.beta()));
        var skato=SummarySetTests.skatO("g",state,new double[]{1,1},SetTestOptions.defaults());
        assertEquals(.3173105078629141,skato.adjustedPValue(),1e-12);assertEquals(0,skato.simulations());
    }
    @Test void extremeNormalTailRetainsLogProbability() {
        var r=SummarySetTests.singleVariant("v",40,1);assertEquals(0,r.pValue());assertTrue(r.negativeLog10PValue()>349);
    }
}
