/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

import static org.junit.jupiter.api.Assertions.*;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import org.junit.jupiter.api.Test;

class CoreStatisticsAuditTest {
    private static Map<String, double[]> reference() throws Exception {
        Map<String, double[]> values = new HashMap<>();
        try (var input = CoreStatisticsAuditTest.class.getResourceAsStream("/core-audit/reference.csv")) {
            assertNotNull(input);
            try (var reader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8))) {
                reader.readLine();
                for (String line; (line = reader.readLine()) != null;) {
                    String[] fields = line.replace("\"", "").split(",");
                    values.put(fields[0], new double[]{Double.parseDouble(fields[1]),
                        Double.parseDouble(fields[2]), Double.parseDouble(fields[3])});
                }
            }
        }
        return values;
    }

    private static void check(double[] expected, StatisticalTestResult actual) {
        assertEquals(expected[0], actual.statistic(), Math.max(1e-10, Math.abs(expected[0])*1e-12));
        assertEquals(expected[1], actual.pValue(), Math.max(1e-300, Math.abs(expected[1])*2e-10));
        if (!Double.isNaN(expected[2])) assertEquals(expected[2], actual.estimates().values().iterator().next(), 1e-13);
    }

    @Test void pearsonIsScaleInvariantIncludingExtremeExponents() throws Exception {
        double[] expected = reference().get("pearson");
        for (double scale : new double[]{1, 1e100, 1e-100, 1e300, 1e-300}) {
            double[] x={1,2,3,4,5,6}, y={2,1,4,3,6,5};
            for(int i=0;i<x.length;i++){x[i]*=scale;y[i]*=scale;}
            check(expected, StatisticalTests.correlation(x,y));
        }
        assertEquals(1, StatisticalTests.correlation(new double[]{-1e308,0,1e308},
            new double[]{-1e308,0,1e308}).estimates().get("correlation"), 1e-15);
        assertTrue(Double.isNaN(StatisticalTests.correlation(new double[]{1,1,1},
            new double[]{1,2,3}).pValue()));
    }

    @Test void signedZerosAreTiesForBothRankTests() throws Exception {
        double[] x={-0.0,0.0,1,2,3,4}, y={1,2,3,2,5,4};
        check(reference().get("kendall-zero"),StatisticalTests.correlation(x,y,CorrelationMethod.KENDALL,Alternative.TWO_SIDED));
        check(reference().get("spearman-zero"),StatisticalTests.correlation(x,y,CorrelationMethod.SPEARMAN,Alternative.TWO_SIDED));
    }

    @Test void rankTestsMatchRIncludingApproximationCutoffs() throws Exception {
        Map<String,double[]> r=reference();
        for(int n:new int[]{12,30,49,50,180,1300,5000}) {
            double[] x=new double[n],y=new double[n];
            for(int i=0;i<n;i++){x[i]=i;y[i]=(i*37)%(n+1);}
            if(n<=50) check(r.get("kendall-exact-"+n),StatisticalTests.correlation(x,y,
                CorrelationMethod.KENDALL,Alternative.TWO_SIDED,.95,true,false));
            check(r.get("kendall-asym-"+n),StatisticalTests.correlation(x,y,
                CorrelationMethod.KENDALL,Alternative.TWO_SIDED,.95,false,false));
            check(r.get("spearman-"+n),StatisticalTests.correlation(x,y,
                CorrelationMethod.SPEARMAN,Alternative.TWO_SIDED));
            if(n>=50) assertEquals("normal approximation",StatisticalTests.correlation(x,y,
                CorrelationMethod.KENDALL,Alternative.TWO_SIDED).pValueMethod());
        }
    }

    @Test void exactKendallRetainsTinyUpperTailAndAvoidsFactorialOverflow() {
        double[] x=new double[30]; for(int i=0;i<x.length;i++)x[i]=i;
        double expected=2;for(int i=2;i<=30;i++)expected/=i;
        assertEquals(expected,StatisticalTests.correlation(x,x,CorrelationMethod.KENDALL,
            Alternative.TWO_SIDED).pValue(),expected*1e-13);
        double[] a=new double[180],b=new double[180];
        for(int i=0;i<a.length;i++){a[i]=i;b[i]=(i*37)%181;}
        double p=StatisticalTests.correlation(a,b,CorrelationMethod.KENDALL,
            Alternative.TWO_SIDED,.95,true,false).pValue();
        assertTrue(p>0 && p<1);
    }

    @Test void fastKendallPairCountsMatchIndependentQuadraticOracle() {
        Random random=new Random(903);
        for(int trial=0;trial<100;trial++) {
            int n=2+random.nextInt(130);double[]x=new double[n],y=new double[n];
            for(int i=0;i<n;i++){x[i]=random.nextInt(12);y[i]=random.nextInt(8);}
            long concordant=0,discordant=0;
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
                double s=(x[i]-x[j])*(y[i]-y[j]);
                if(s>0)concordant++;if(s<0)discordant++;
            }
            assertArrayEquals(new long[]{concordant,discordant},StatisticsSupport.kendallPairs(x,y));
        }
    }

    @Test void anovaSupportsConstantGroupsAndLargeOrSmallUnits() throws Exception {
        check(reference().get("anova-constant-group"),StatisticalTests.oneWayAnova(
            new double[][]{{1,1,1},{2,3,4},{3,5,4}},true));
        for(boolean equal:new boolean[]{true,false})for(double scale:new double[]{1,1e200,1e-200}) {
            double[][]g={{1,2,1.5,3},{2,3,4,4.1},{3,5,4,6}};
            for(double[]group:g)for(int i=0;i<group.length;i++)group[i]*=scale;
            check(reference().get("anova-"+(equal?"TRUE":"FALSE")),StatisticalTests.oneWayAnova(g,equal));
        }
    }

    @Test void tAndVarianceTestsPreserveInferenceWhenUnitsChange() throws Exception {
        for(double scale:new double[]{1,1e200,1e-200}) {
            double[]x={scale,2*scale,4*scale,3*scale},y={3*scale,5*scale,4*scale,6*scale};
            var one=StatisticalTests.studentT(x,0);
            assertEquals(reference().get("t-one")[1],one.pValue(),1e-14);
            assertEquals(3.8729833462074175,one.statistic(),1e-13);
            assertEquals(.030466291662170977,one.pValue(),1e-14);
            assertEquals(2.5,one.estimates().get("mean")/scale,1e-14);
            for(boolean equal:new boolean[]{true,false}) {
                var two=StatisticalTests.studentT(x,y,0,equal,Alternative.TWO_SIDED,.95);
                assertEquals(reference().get("t-welch")[1],two.pValue(),1e-14);
                assertEquals(-2.1908902300206643,two.statistic(),1e-13);
                assertEquals(.07098765432098764,two.pValue(),1e-14);
                assertEquals(-4.2337146951647044,two.confidenceInterval().orElseThrow().lower()/scale,1e-13);
            }
            assertEquals(1,StatisticalTests.variance(x,y).statistic(),1e-14);
            assertEquals(reference().get("variance")[1],StatisticalTests.variance(x,y).pValue(),1e-14);
        }
    }
    @Test void proportionTrendIsIndependentOfScoreUnits() {
        int[] successes={2,5,10,15},trials={20,20,20,20};
        var expected=StatisticalTests.proportionTrend(successes,trials,new double[]{1,2,3,4});
        for(double scale:new double[]{1e200,1e-200}){
            var actual=StatisticalTests.proportionTrend(successes,trials,new double[]{scale,2*scale,3*scale,4*scale});
            assertEquals(expected.statistic(),actual.statistic(),1e-13);
            assertEquals(expected.pValue(),actual.pValue(),1e-14);
        }
    }
    @Test void pairedSummaryAndLargeNullRetainFiniteEstimableQuantities() {
        var pair=StatisticalTests.pairedT(new double[]{1e308,1e308,1e308},
            new double[]{0,1e307,2e307},0,Alternative.TWO_SIDED,.95);
        assertEquals(9,pair.estimates().get("mean difference")/1e307,1e-14);
        assertEquals(15.5884572681199,pair.statistic(),1e-12);
        var one=StatisticalTests.studentT(new double[]{1,2,3},1e300);
        assertEquals(-Math.sqrt(3),one.statistic()/1e300,1e-14);
        assertEquals(0,one.pValue());
        var two=StatisticalTests.studentT(new double[]{1,2,3},new double[]{1,2,3},
            1e300,false,Alternative.TWO_SIDED,.95);
        assertEquals(-Math.sqrt(1.5),two.statistic()/1e300,1e-14);
    }
}
