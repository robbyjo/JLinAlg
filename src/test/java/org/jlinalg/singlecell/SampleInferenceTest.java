/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.singlecell;

import java.util.*;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class SampleInferenceTest {
    private final String[] condition={"a","a","a","a","b","b","b","b"};
    private final double[] response={1,2,4,5,4,6,5,9};
    @Test void independentAndPairedUncertaintyMatchHandCalculations() {
        var independent=SampleInference.design(new String[]{"1","2","3","4","5","6","7","8"},condition,"a","b",new double[8][0],List.of(),false);
        var fit=SampleInference.fit(response,independent);
        assertEquals(3,fit.coefficients()[1],1e-12);assertEquals(Math.sqrt(2),fit.standardErrors()[1],1e-12);assertEquals(6,fit.residualDegreesOfFreedom());
        var paired=SampleInference.design(new String[]{"1","2","3","4","1","2","3","4"},condition,"a","b",new double[8][0],List.of(),true);
        var pair=SampleInference.fit(response,paired);
        assertEquals(3,pair.coefficients()[1],1e-12);assertEquals(Math.sqrt(.5),pair.standardErrors()[1],1e-12);assertEquals(3,pair.residualDegreesOfFreedom());
    }
    @Test void pseudoreplicationIncompletePairsAndConfoundingReject() {
        String[] ids={"1","2","3","4","1","2","3","4"};
        assertThrows(IllegalArgumentException.class,()->SampleInference.design(ids,condition,"a","b",new double[8][0],List.of(),false));
        String[] incomplete=ids.clone();incomplete[7]="5";
        assertThrows(IllegalArgumentException.class,()->SampleInference.design(incomplete,condition,"a","b",new double[8][0],List.of(),true));
        double[][] confound=new double[8][1];for(int i=4;i<8;i++)confound[i][0]=1;
        assertThrows(IllegalArgumentException.class,()->SampleInference.design(ids,condition,"a","b",confound,List.of("batch"),true));
    }
    @Test void replicatedGaussianNullHasReasonableFalsePositiveRate() {
        String[] ids={"1","2","3","4","5","6","7","8"};var design=SampleInference.design(ids,condition,"a","b",new double[8][0],List.of(),false);
        Random random=new Random(170);int rejected=0;
        for(int b=0;b<400;b++){double[] y=new double[8];for(int i=0;i<8;i++)y[i]=random.nextGaussian();if(SampleInference.fit(y,design).pValues()[1]<.05)rejected++;}
        assertTrue(rejected>=7&&rejected<=36,"Gaussian null rejections: "+rejected);
    }
}
