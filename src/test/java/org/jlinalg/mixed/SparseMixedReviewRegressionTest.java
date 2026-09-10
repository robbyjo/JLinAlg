/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;
import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.*;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

public final class SparseMixedReviewRegressionTest {
    private static Properties reference() throws Exception {
        Properties p=new Properties();
        try(var in=SparseMixedReviewRegressionTest.class.getResourceAsStream("/r-reference/mixed-review.properties")){p.load(in);}
        return p;
    }
    private static double ref(Properties p,String name){return Double.parseDouble(p.getProperty(name));}
    @Test void responseShiftCannotChangeVarianceEstimatesOrLikelihood() throws Exception {
        int n=40; double[] y=new double[n], shifted=new double[n],x=new double[n];
        String[] g=new String[n]; double[][] z=new double[n][1];
        for(int i=0;i<n;i++){y[i]=2+2*Math.sin(i/5.0)+Math.cos(2.1*i);shifted[i]=y[i]+1e8;x[i]=1;g[i]="g"+(i/5);z[i][0]=1;}
        var term=RandomEffectTerm.randomIntercept("g",Arrays.asList(g));
        var block=CorrelatedRandomEffectBlock.of("g",Arrays.asList(g),List.of("intercept"),z);
        var options=RemlOptions.defaults();
        var a=SparseLinearMixedModel.fit(y,x,n,1,List.of(term),options,BackendPolicy.CPU);
        var b=SparseLinearMixedModel.fit(shifted,x,n,1,List.of(term),options,BackendPolicy.CPU);
        var c=SparseUnstructuredCorrelatedModel.fit(y,x,n,1,List.of(block),options,BackendPolicy.CPU).correlatedFit();
        var d=SparseUnstructuredCorrelatedModel.fit(shifted,x,n,1,List.of(block),options,BackendPolicy.CPU).correlatedFit();
        assertArrayEquals(a.varianceComponents(),b.varianceComponents(),2e-6);
        assertEquals(a.logLikelihood(),b.logLikelihood(),1e-6);
        assertEquals(a.beta()[0]+1e8,b.beta()[0],1e-7);
        assertArrayEquals(c.randomEffects().get(0).covariance(),d.randomEffects().get(0).covariance(),2e-6);
        assertEquals(c.logLikelihood(),d.logLikelihood(),1e-6);
        Properties p=reference();
        assertEquals(ref(p,"randomVariance"),a.varianceComponents()[0],1e-5);
        assertEquals(ref(p,"residualVariance"),a.varianceComponents()[1],1e-5);
        assertEquals(ref(p,"logLik"),a.logLikelihood(),1e-7);
        System.out.println("shift-stable variances="+Arrays.toString(b.varianceComponents())+" LL="+b.logLikelihood());
    }
    @Test void varianceBoundsConstrainPhysicalVariancesInBothFitters() throws Exception {
        int n=40;double[] y=new double[n],x=new double[n];String[] g=new String[n];double[][] z=new double[n][1];
        for(int i=0;i<n;i++){y[i]=2+2*Math.sin(i/5.0)+Math.cos(2.1*i);x[i]=1;g[i]="g"+(i/5);z[i][0]=1;}
        var options=RemlOptions.builder().varianceBounds(.01,.1).build();
        var a=SparseLinearMixedModel.fit(y,x,n,1,List.of(RandomEffectTerm.randomIntercept("g",Arrays.asList(g))),options,BackendPolicy.CPU);
        var b=SparseUnstructuredCorrelatedModel.fit(y,x,n,1,List.of(CorrelatedRandomEffectBlock.of("g",Arrays.asList(g),List.of("i"),z)),options,BackendPolicy.CPU).correlatedFit();
        Properties p=reference();
        for(double value:a.varianceComponents())assertTrue(value>=.01-1e-12&&value<=.1+1e-12);
        assertEquals(ref(p,"boundedRandom"),a.varianceComponents()[0],1e-7);
        assertEquals(ref(p,"boundedResidual"),a.varianceComponents()[1],1e-7);
        assertEquals(ref(p,"boundedLogLik"),a.logLikelihood(),1e-6);
        assertEquals(ref(p,"boundedRandom"),b.randomEffects().get(0).covariance()[0],1e-7);
        assertEquals(ref(p,"boundedResidual"),b.residualVariance(),1e-7);
        assertEquals(ref(p,"boundedLogLik"),b.logLikelihood(),1e-6);
        var finiteDf = options.toBuilder().degreesOfFreedomMethod(
            org.jlinalg.inference.DegreesOfFreedomMethod.KENWARD_ROGER).build();
        assertThrows(IllegalArgumentException.class, () -> SparseLinearMixedModel.fit(y,x,n,1,
            List.of(RandomEffectTerm.randomIntercept("g",Arrays.asList(g))),finiteDf,BackendPolicy.CPU));
        assertThrows(IllegalArgumentException.class, () -> SparseUnstructuredCorrelatedModel.fit(y,x,n,1,
            List.of(CorrelatedRandomEffectBlock.of("g",Arrays.asList(g),List.of("i"),z)),finiteDf,BackendPolicy.CPU));
        System.out.println("bounded variances="+Arrays.toString(a.varianceComponents())+" LL="+a.logLikelihood());
    }
    @Test void correlatedBoundsApplyToMarginalVariancesNotCholeskyDiagonals() {
        int n=80;double[] y=new double[n],x=new double[n];String[] g=new String[n];
        for(int i=0;i<n;i++){int j=i%8,group=i/8;x[i]=j-3.5;g[i]="g"+group;
            y[i]=2+.3*x[i]+2*Math.sin(group)*(1+.5*x[i])+((j%4==0||j%4==3)?.3:-.3);}
        var table=ModelTable.builder(n).numeric("y",y).numeric("x",x).categorical("g",g).build();
        var fit=MixedFormula.compile("y~x+(1+x|g)",table).fitSparseUnstructured(
            RemlOptions.builder().varianceBounds(.01,.1).maximumIterations(500).build(),BackendPolicy.CPU).correlatedFit();
        assertTrue(fit.converged());
        double[] covariance=fit.randomEffects().get(0).covariance();
        for(double variance:new double[]{covariance[0],covariance[3],fit.residualVariance()})
            assertTrue(variance>=.01-1e-12&&variance<=.1+1e-12,"physical variance="+variance);
    }
    @Test void boundaryCorrelationsHaveAnEstimableOneSidedProfile() throws Exception {
        Properties p=reference();
        for(int sign:new int[]{1,-1}){
            int n=80;double[] y=new double[n],x=new double[n];String[] g=new String[n];
            for(int i=0;i<n;i++){int j=i%8,group=i/8;x[i]=j-3.5;g[i]="g"+group;
                y[i]=2+.3*x[i]+2*Math.sin(group)*(1+sign*.5*x[i])+((j%4==0||j%4==3)?.3:-.3);}
            var table=ModelTable.builder(n).numeric("y",y).numeric("x",x).categorical("g",g).build();
            var interval=MixedFormula.compile("y~x+(1+x|g)",table).profileCorrelation(0,.95,
                RemlOptions.builder().maximumIterations(500).build(),BackendPolicy.CPU);
            System.out.println("boundary correlation "+interval);
            if(sign==1){assertEquals(1,interval.upper(),1e-10);assertEquals(ref(p,"correlationLower"),interval.lower(),2e-6);assertTrue(interval.lowerFound());}
            else{assertEquals(-1,interval.lower(),1e-10);assertEquals(ref(p,"correlationUpper"),interval.upper(),2e-6);assertTrue(interval.upperFound());}
        }
    }
    @Test void profilesAnyPairInAThreeCoefficientBlock() {
        int n=96; double[] y=new double[n],x=new double[n],z=new double[n];
        String[] g=new String[n];
        for(int i=0;i<n;i++){
            int j=i%8,group=i/8;
            x[i]=j-3.5;z[i]=(j%3)-1;g[i]="g"+group;
            double randomIntercept=1.2*Math.sin(.7*group);
            double randomX=.35*Math.cos(.5*group);
            double randomZ=.25*Math.sin(1.1*group);
            y[i]=2+.2*x[i]-.15*z[i]+randomIntercept
                +randomX*x[i]+randomZ*z[i]+.08*Math.cos(1.7*i);
        }
        var table=ModelTable.builder(n).numeric("y",y).numeric("x",x)
            .numeric("z",z).categorical("g",g).build();
        var formula=MixedFormula.compile("y~x+z+(1+x+z|g)",table);
        var options=RemlOptions.builder().maximumIterations(500).build();
        var fitted=formula.fitCorrelated(options,BackendPolicy.CPU);
        double[] pev=fitted.randomEffects().get(0)
            .predictionErrorCovariances();
        assertEquals(12*9,pev.length);
        for(int group=0;group<12;group++)for(int first=0;first<3;first++){
            assertTrue(pev[group*9+first*3+first]>=0);
            for(int second=0;second<3;second++)
                assertEquals(pev[group*9+first*3+second],
                    pev[group*9+second*3+first],2e-10);
        }
        var interval=formula.profileCorrelation(
            0,1,2,.90,options,BackendPolicy.CPU);
        assertEquals(fitted.randomEffects().get(0).correlation(1,2),
            interval.estimate(),2e-4);
        assertTrue(interval.estimate()>=-1&&interval.estimate()<=1);
        assertTrue(interval.lowerFound()||interval.upperFound()
            || interval.lower()==-1||interval.upper()==1);
    }
    public static void main(String[]args)throws Exception{
        var t=new SparseMixedReviewRegressionTest();t.responseShiftCannotChangeVarianceEstimatesOrLikelihood();
        t.varianceBoundsConstrainPhysicalVariancesInBothFitters();t.boundaryCorrelationsHaveAnEstimableOneSidedProfile();
        t.profilesAnyPairInAThreeCoefficientBlock();
    }
}
