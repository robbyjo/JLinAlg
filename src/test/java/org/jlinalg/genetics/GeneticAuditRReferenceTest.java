/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.coloc.*;
import org.jlinalg.susie.*;
import org.junit.jupiter.api.Test;

class GeneticAuditRReferenceTest {
    @Test void matchesIndependentRConfigurationEnumerationAndSingleEffectPosterior() throws Exception {
        var ref=new Properties();
        try(var in=getClass().getResourceAsStream("/r-reference/genetic-audit-reference.properties")){assertNotNull(in);ref.load(in);}
        double[] a=new double[64],b=new double[64],z=new double[32];
        for(int i=0;i<64;i++){a[i]=5*Math.sin((i+1)*.37)+.1*(i+1);b[i]=4*Math.cos((i+1)*.23)+.07*(i+1);}
        for(int i=0;i<32;i++)z[i]=5*Math.sin((i+1)*.37);
        var names=java.util.stream.IntStream.range(0,64).mapToObj(i->"s"+i).toList();
        var coloc=ColocSusie.analyze(new ColocSusieInput(names,new double[][]{a}),new ColocSusieInput(names,new double[][]{b}));
        check(ref,"coloc64",coloc.signalPairs().get(0).hypothesisPosteriors(),0);
        check(ref,"coloc64",coloc.sharedVariantPosterior(),5);
        double[][] ld=new double[32][32];for(int i=0;i<32;i++)ld[i][i]=1;
        var fit=Susie.fitSummary(z,ld,400,null,new SusieOptions(1,100,1e-8,.2,false,.95,0),BackendPolicy.CPU);
        assertTrue(fit.converged());
        check(ref,"ser32",fit.pip(),0); check(ref,"ser32",fit.posteriorMean(),32);check(ref,"ser32",fit.logBayesFactors(),64);
        var grm=GenomicRelationshipMatrix.fromVariantDosages(new double[][]{{0,1,2},{0,Double.NaN,2},{1,2,0}},List.of("a","b","c"),new GenomicRelationshipOptions(.01,.5),BackendPolicy.CPU);
        check(ref,"grm",grm.relationshipMatrix(),0);
    }
    private static void check(Properties ref,String key,double[] values,int offset){
        for(int i=0;i<values.length;i++){
            double expected=Double.parseDouble(ref.getProperty(key+"."+(offset+i)));
            assertEquals(expected,values[i],3e-12*Math.max(1,Math.abs(expected)),key+"["+i+"]");
        }
    }
}
