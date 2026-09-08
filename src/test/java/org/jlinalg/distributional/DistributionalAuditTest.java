/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PenalizedPredictor;
import org.junit.jupiter.api.Test;

public class DistributionalAuditTest {
    @Test void hurdleAndZipHighCountDensitiesMatchR() throws Exception {
        List<String> rows=lines("distributional-poisson-extremes.tsv");
        for(String row:rows.subList(1,rows.size())) {
            String[] f=row.split("\t"); double y=Double.parseDouble(f[0]),mu=Double.parseDouble(f[1]);
            double zip=Double.parseDouble(f[2]),hurdle=Double.parseDouble(f[3]);
            assertEquals(hurdle,new HurdlePoissonFamily().logLikelihood(y,new double[]{mu,.2}),2e-14*Math.max(1,Math.abs(hurdle)),row);
            assertEquals(zip,new ZeroInflatedPoissonFamily().logLikelihood(y,new double[]{mu,.2}),2e-14*Math.max(1,Math.abs(zip)),row);
            // The mixed helper accepts eta: account for exp(log(mu)) rounding.
            double actualMean=Math.exp(Math.log(mu));
            double expected=y==0?Math.log(.2+.8*Math.exp(-actualMean)):
                Math.log(.8)+jdistlib.Poisson.density(y,actualMean,true);
            assertEquals(expected,SparseZeroInflatedMixedModel.likelihoodDerivatives(y,Math.log(mu),Math.log(.25),1,false)[0],2e-14*Math.max(1,Math.abs(expected)),row);
        }
    }
    @Test void betaMixedObservedInformationIncludesNoncanonicalScoreTerm() throws Exception {
        var constructor=Class.forName("org.jlinalg.distributional.BetaMixedModel$LaplaceBetaFamily").getDeclaredConstructors()[0];
        constructor.setAccessible(true);
        var family=(org.jlinalg.glm.GlmFamily)constructor.newInstance(3.,.01,1e6);
        var exact=(org.jlinalg.glmm.LaplaceFamilyDerivatives)family;
        for(double eta:new double[]{-4,-2,0,2,4})for(double y:new double[]{.01,.2,.8,.99}) {
            double h=1e-5,weight=.7;
            double score=(family.logLikelihood(y,family.inverseLink(eta+h),weight,1)-family.logLikelihood(y,family.inverseLink(eta-h),weight,1))/(2*h);
            double information=-(exact.linearPredictorScore(y,eta+h,weight)-exact.linearPredictorScore(y,eta-h,weight))/(2*h);
            assertEquals(score,exact.linearPredictorScore(y,eta,weight),2e-8);
            assertEquals(information,exact.linearPredictorInformation(y,eta,weight),2e-8);
        }
        assertTrue(exact.linearPredictorInformation(.8,-2,1)<0);
    }
    private static List<String> lines(String name) throws Exception {
        try(var in=DistributionalAuditTest.class.getResourceAsStream("/r-reference/"+name)) {
            return new String(Objects.requireNonNull(in).readAllBytes(),StandardCharsets.UTF_8).lines().toList();
        }
    }
    private static Properties reference(String name) throws Exception {
        Properties p=new Properties();
        try(var in=DistributionalAuditTest.class.getResourceAsStream("/r-reference/"+name)){p.load(in);}
        return p;
    }
    @Test void betaAndNbLogDensitiesStayAccurateAtExtremePrecisionAndCounts() throws Exception {
        double max=0;
        List<String> densityRows=lines("distributional-audit-densities.tsv");
        for(String line:densityRows.subList(1,densityRows.size())) {
            String[] f=line.split("\t");double y=Double.parseDouble(f[1]),mu=Double.parseDouble(f[2]),p=Double.parseDouble(f[3]);
            DistributionalFamily family=f[0].equals("beta")?new BetaMeanPrecisionFamily():new NegativeBinomialMeanDispersionFamily();
            double actual=family.logLikelihood(y,new double[]{mu,p}),expected=Double.parseDouble(f[4]);
            // R 4.6.1's small-count branch substitutes -mu for
            // -size*log1p(mu/size). JDistlib 0.10.2 fixes that approximation.
            // Keep the historical R fixture, but use the exact finite-product
            // value here; the independent grid below checks this more tightly.
            if(f[0].equals("nb") && y==2 && mu==1e5 && p==1e14)
                expected=-99977.66724625262;
            max=Math.max(max,Math.abs(actual-expected));
            assertEquals(expected,actual,2e-10*Math.max(1,Math.abs(expected)),line);
        }
        System.out.println("distributional extreme-density maximum absolute error="+max);
    }
    @Test void negativeBinomialSmallCountsMatchIndependentFiniteProduct() throws Exception {
        List<String> rows=lines("jdistlib-upgrade-nb.tsv");
        double maximumScaledError=0;
        for(String row:rows.subList(1,rows.size())) {
            String[] f=row.split("\t");
            double y=Double.parseDouble(f[0]),mu=Double.parseDouble(f[1]),size=Double.parseDouble(f[2]);
            double expected=Double.parseDouble(f[4]);
            double actual=new NegativeBinomialMeanDispersionFamily().logLikelihood(y,new double[]{mu,size});
            double scale=Math.max(1,Math.abs(expected));
            maximumScaledError=Math.max(maximumScaledError,Math.abs(actual-expected)/scale);
            // Tight gate for the branch improved in 0.10.2. Outside it, both
            // dependency versions retain up to 1.51e-9 scaled error on this
            // grid; record that pre-existing limitation, not false FP64 parity.
            double tolerance=y < 1e-10*size ? 5e-14 : 2e-9;
            assertEquals(expected,actual,tolerance*scale,row);
        }
        System.out.println("NB finite-product grid maximum scaled error="+maximumScaledError);
    }
    @Test void hurdleInformationHandlesSmallAndLargeMeansWithoutClippingOrOverflow() {
        var family=new HurdlePoissonFamily();
        for(double mu:new double[]{1e-10,1e-6,.1,1,20,1000}) {
            double[] score=new double[2],info=new double[4];
            family.derivatives(1,new double[]{mu,.2},score,info);
            assertTrue(Double.isFinite(info[0])&&info[0]>0);
            if(mu<1e-3) {
                assertEquals(.8*(mu/2+mu*mu/6),info[0],mu*mu*mu*.01+1e-25);
                assertEquals(-mu/2-mu*mu/12,score[0],1e-25+mu*mu*mu*.01);
            }
            if(mu==1000)assertEquals(800,info[0],1e-12);
            if(mu>=.1) {
                double h=1e-5;double plus=family.logLikelihood(1,new double[]{mu*Math.exp(h),.2});
                double minus=family.logLikelihood(1,new double[]{mu*Math.exp(-h),.2});
                assertEquals((plus-minus)/(2*h),score[0],1e-7*Math.max(1,mu));
            }
        }
    }
    @Test void largeMeanHurdleFitMatchesIndependentRMaximum() throws Exception {
        List<String> data=lines("distributional-audit-hurdle.tsv");int n=data.size()-1;
        double[] y=new double[n];double[][] x=new double[n][1];
        for(int r=0;r<n;r++){y[r]=Double.parseDouble(data.get(r+1));x[r][0]=1;}
        var predictor=PenalizedPredictor.linear(x);
        var fit=DistributionalModel.fit(y,List.of(predictor,predictor),new HurdlePoissonFamily(),DistributionalOptions.defaults(),BackendPolicy.CPU);
        var p=reference("distributional-audit-hurdle.properties");
        assertTrue(fit.converged(),fit.convergenceMessage());
        assertEquals(Double.parseDouble(p.getProperty("mu")),fit.parameter("mu").fittedValues()[0],1e-4);
        assertEquals(Double.parseDouble(p.getProperty("zero")),fit.parameter("zeroProbability").fittedValues()[0],1e-10);
        assertEquals(Double.parseDouble(p.getProperty("logLik")),fit.logLikelihood(),1e-7);
    }
    @Test void tinyStepLimitsCannotCertifyNonstationaryFits() {
        double[] y={.1,.15,.2,.35,.4,.5,.6,.65,.8,.9};double[][] x=new double[10][2],p=new double[10][1];
        for(int r=0;r<10;r++){x[r][0]=1;x[r][1]=r;p[r][0]=1;}
        var fit=BetaRegression.fit(y,x,p,new BetaRegressionOptions(10,1e-8,1e-12,BetaMeanLink.LOGIT,BetaPrecisionLink.LOG),BackendPolicy.CPU);
        assertFalse(fit.converged());
        var generic=DistributionalModel.fit(y,List.of(PenalizedPredictor.linear(x),PenalizedPredictor.linear(p)),
            new BetaMeanPrecisionFamily(),new DistributionalOptions(10,1e-8,1e-12,.95),BackendPolicy.CPU);
        assertFalse(generic.converged());
    }
    @Test void nbObservedCovarianceMatchesGlmmTmbRatherThanScoringMetric() throws Exception {
        List<String> data=lines("distributional-audit-nb.tsv");int n=data.size()-1;
        double[] y=new double[n];double[][] x=new double[n][2],z=new double[n][2];
        for(int r=0;r<n;r++){String[] f=data.get(r+1).split("\t");y[r]=Double.parseDouble(f[0]);x[r][0]=z[r][0]=1;x[r][1]=Double.parseDouble(f[1]);z[r][1]=Double.parseDouble(f[2]);}
        var fit=DistributionalModel.fit(y,List.of(PenalizedPredictor.linear(x),PenalizedPredictor.linear(z)),new NegativeBinomialMeanDispersionFamily(),
            new DistributionalOptions(500,1e-9,3,.95),BackendPolicy.CPU);
        var p=reference("distributional-audit-nb.properties");
        assertTrue(fit.converged(),fit.convergenceMessage());
        assertEquals(Double.parseDouble(p.getProperty("logLik")),fit.logLikelihood(),2e-6);
        double[] coefficients={fit.parameter("mu").coefficients()[0],fit.parameter("mu").coefficients()[1],fit.parameter("size").coefficients()[0],fit.parameter("size").coefficients()[1]};
        for(int c=0;c<4;c++)assertEquals(Double.parseDouble(p.getProperty("coefficient"+c)),coefficients[c],2e-4);
        double[] covariance=fit.covariance();double maximum=0;
        System.out.println("NB beta="+Arrays.toString(coefficients)+" LL="+fit.logLikelihood()+" cov="+Arrays.toString(covariance));
        for(int c=0;c<16;c++){double error=Math.abs(Double.parseDouble(p.getProperty("covariance"+c))-covariance[c]);maximum=Math.max(maximum,error);assertTrue(error<2e-5,"covariance entry "+c+" error="+error);}
        System.out.println("NB observed covariance maximum absolute R error="+maximum);
    }
    @Test void zeroInflatedNbRetainsPoissonLimitAndObservedDerivatives() {
        for(double y:new double[]{0,2,100}) {
            double[] nb=SparseZeroInflatedMixedModel.likelihoodDerivatives(y,Math.log(2),-.7,1e14,true);
            double[] poisson=SparseZeroInflatedMixedModel.likelihoodDerivatives(y,Math.log(2),-.7,1,false);
            assertArrayEquals(poisson,nb,2e-8);
            double h=1e-4;
            double[] plus=SparseZeroInflatedMixedModel.likelihoodDerivatives(y,Math.log(2)+h,-.7,1e14,true);
            double[] minus=SparseZeroInflatedMixedModel.likelihoodDerivatives(y,Math.log(2)-h,-.7,1e14,true);
            assertEquals(nb[1],(plus[0]-minus[0])/(2*h),2e-6);
            assertEquals(nb[3],-(plus[1]-minus[1])/(2*h),2e-6);
        }
    }
}
