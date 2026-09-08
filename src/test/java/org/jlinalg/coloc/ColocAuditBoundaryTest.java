/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.coloc;

import static org.junit.jupiter.api.Assertions.*;
import java.util.List;
import org.junit.jupiter.api.Test;

class ColocAuditBoundaryTest {
    private static final double N = Double.NEGATIVE_INFINITY;
    private ColocSusieResult fit(double[] a, double[] b, double[] w1, double[] w2) {
        var ids = java.util.stream.IntStream.range(0,a.length).mapToObj(i->"v"+i).toList();
        return ColocSusie.analyze(new ColocSusieInput(ids,new double[][]{a}),
            new ColocSusieInput(ids,new double[][]{b}),
            new ColocOptions(1e-4,1e-4,5e-6,0,false,w1,w2));
    }
    @Test void disjointEvidenceIsH3NotNull() {
        var fit = fit(new double[]{100,N},new double[]{N,100},null,null);
        assertEquals(1,fit.signalPairs().get(0).posteriorH3(),1e-15);
        assertEquals(0,fit.signalPairs().get(0).posteriorH4());
        assertArrayEquals(new double[2],fit.sharedVariantPosterior());
    }
    @Test void dominantSharedVariantDoesNotEraseDistinctConfigurations() {
        var fit = fit(new double[]{1000,950},new double[]{1000,950},null,null);
        // H3/H4 = (p1*p2/p12) * 2 exp(-50)/(1+exp(-100)).
        double expected = .004*Math.exp(-50)/(1+Math.exp(-100));
        double actual = fit.signalPairs().get(0).posteriorH3()/fit.signalPairs().get(0).posteriorH4();
        assertEquals(expected,actual,expected*2e-12);
    }
    @Test void weightedH4PosteriorAndWeightRescalingAgree() {
        double[] a={2,-1,.5}, b={1,.2,-.5};
        var fit=fit(a,b,new double[]{1,2,3},new double[]{3,1,2});
        var large=fit(a,b,new double[]{5e307,1e308,1.5e308},new double[]{3e-300,1e-300,2e-300});
        double norm=3*Math.exp(3)+2*Math.exp(-.8)+6;
        assertArrayEquals(new double[]{3*Math.exp(3)/norm,2*Math.exp(-.8)/norm,6/norm},fit.sharedVariantPosterior(),3e-15);
        assertArrayEquals(fit.signalPairs().get(0).hypothesisPosteriors(),large.signalPairs().get(0).hypothesisPosteriors(),2e-14);
        assertArrayEquals(fit.sharedVariantPosterior(),large.sharedVariantPosterior(),3e-14);
    }
    @Test void traitSwapAndVariantPermutationPreserveEvidence() {
        var a = new ColocSusieInput(List.of("a","b","c"),new double[][]{{3,1,2}});
        var b = new ColocSusieInput(List.of("c","a","b"),new double[][]{{5,4,2}});
        var ab=ColocSusie.analyze(a,b); var ba=ColocSusie.analyze(b,a);
        double[] h=ab.signalPairs().get(0).hypothesisPosteriors(), g=ba.signalPairs().get(0).hypothesisPosteriors();
        assertArrayEquals(new double[]{h[0],h[2],h[1],h[3],h[4]},g,2e-15);
        assertArrayEquals(new double[]{ab.sharedVariantPosterior()[2],ab.sharedVariantPosterior()[0],ab.sharedVariantPosterior()[1]},ba.sharedVariantPosterior(),2e-15);
    }
}
