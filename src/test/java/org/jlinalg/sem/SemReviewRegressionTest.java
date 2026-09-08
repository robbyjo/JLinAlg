/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

/** Regressions for the six independent statistical SEM review findings. */
class SemReviewRegressionTest {
    private static Properties reference() throws Exception {
        Properties p=new Properties();try(var in=SemReviewRegressionTest.class.getResourceAsStream("/r-reference/sem-review.properties")){p.load(in);}return p;
    }
    private static double value(Properties p,String key){return Double.parseDouble(p.getProperty(key));}
    private static int[][] table(int a,int b,int c,int d) {
        int[][] data=new int[a+b+c+d][2];int[][] counts={{a,b},{c,d}};int r=0;
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int h=0;h<counts[i][j];h++)data[r++]=new int[]{i,j};return data;
    }
    private static SemModel fixedCorrelation(double rho){return SemModel.builder("x","y").fixedVariance("x",1).fixedVariance("y",1).fixedCovariance("x","y",rho).build();}
    @Test void acceptedHighCorrelationsMatchAnalyticProbabilityAndDerivative() throws Exception {
        Properties p=reference();
        for(double rho:new double[]{-.99989,-.9998,-.999,-.99,0,.99,.999,.9998,.99989})
            assertEquals(.25+Math.asin(rho)/(2*Math.PI),SemOrdinal.bvn(0,0,rho),2e-14);
        assertEquals(value(p,"orthant"),SemOrdinal.bvn(0,0,.9998),2e-14);
        double h=1e-8,derivative=(SemOrdinal.bvn(0,0,.9998+h)-SemOrdinal.bvn(0,0,.9998-h))/(2*h);
        assertEquals(value(p,"orthantDerivative"),derivative,1e-7);
        SemOrdinal.Result fit=SemOrdinal.fit(table(4968,32,32,4968),new int[]{2,2},fixedCorrelation(.9998));
        assertTrue(fit.converged());assertTrue(fit.informationAvailable());
        assertEquals(value(p,"nearBoundaryLL"),fit.pairwiseLogLikelihood(),2e-9);
    }
    @Test void rareRectangleProbabilitiesMatchIndependentPositiveRIntegration() throws Exception {
        Properties p=reference();
        for(int i=1;i<=5;i++) {
            double[] b=Arrays.stream(p.getProperty("rectangle."+i+".bounds").replace("Inf","Infinity").split(",")).mapToDouble(Double::parseDouble).toArray();
            double expected=value(p,"rectangle."+i+".probability");
            double actual=SemNormalRectangle.rectangle(b[0],b[1],b[2],b[3],b[4]);
            assertTrue(actual>0);assertEquals(1,actual/expected,2e-9,"rectangle "+i);
        }
    }
    @Test void zeroCellsNeverEnterCaseScoresAndObservedTinyCellsRemainUsable() {
        SemOrdinal.Result fit=SemOrdinal.fit(table(1,98,0,1),new int[]{2,2},fixedCorrelation(.95));
        assertTrue(fit.converged());assertTrue(fit.informationAvailable());
        // The absent cell has negligible mass here; the other three masses are .01,.98,.01.
        assertEquals(2*Math.log(.01)+98*Math.log(.98),fit.pairwiseLogLikelihood(),1e-8);
        SemOrdinal.Result observedRare=SemOrdinal.fit(table(1,97,1,1),new int[]{2,2},fixedCorrelation(.95));
        assertTrue(observedRare.converged(),"score="+observedRare.scoreNorm());assertTrue(observedRare.informationAvailable());
        assertTrue(Double.isFinite(observedRare.pairwiseLogLikelihood()));
    }
    @Test void ordinalRankChecksSuppressAllInferenceIncludingThresholds() throws Exception {
        SemModel overparameterized=SemModel.builder("x","y").fixedVariance("x",1).fixedVariance("y",1)
            .regression("b","y","x",.2).covariance("c","x","y",.2).build();
        assertThrows(IllegalArgumentException.class,()->SemOrdinal.fit(table(35,15,15,35),new int[]{2,2},overparameterized));
        // Parameter count alone passes: 5 structural labels vs 6 correlations,
        // but the unfixed factor scale makes the analytic moment Jacobian singular.
        SemModel redundant=SemModel.builder("z1","z2","z3","z4").latent("f")
            .loading("a","z1","f",.8).loading("b","z2","f",.8).loading("c","z3","f",.8).loading("d","z4","f",.8)
            .fixedVariance("z1",1).fixedVariance("z2",1).fixedVariance("z3",1).fixedVariance("z4",1).build();
        int[][] data=Arrays.stream(SemJointTest.data("ordinal")).map(row->Arrays.stream(row).mapToInt(v->(int)v).toArray()).toArray(int[][]::new);
        SemOrdinal.Result fit=SemOrdinal.fit(data,new int[]{3,3,3,3},redundant);
        assertFalse(fit.informationAvailable());assertTrue(Arrays.stream(fit.parameterCovariance()).allMatch(Double::isNaN));
        assertTrue(fit.parameters().stream().allMatch(parameter->Double.isNaN(parameter.standardError()) && Double.isNaN(parameter.pValue())));
    }
    @Test void redundantLatentMeanCannotChangeTheReportedSignificance() {
        SemModel valid=SemModel.builder("x1","x2","x3","x4").latent("f").fixedLoading("x1","f",1)
            .loading("x2","f",.8).loading("x3","f",.8).loading("x4","f",.8).meanStructure().build();
        double[] cov={2,1.3,1,1,1.3,2,1,1,1,1,2,1,1,1,1,2},mean={1,2,3,4};
        SemFitResult identified=Sem.fitMoments(cov,mean,500,valid,SemOptions.defaults(),BackendPolicy.CPU);
        SemFitResult redundant=Sem.fitMoments(cov,mean,500,valid.toBuilder().intercept("f",.2).build(),SemOptions.defaults(),BackendPolicy.CPU);
        assertEquals(identified.logLikelihood(),redundant.logLikelihood(),1e-7);
        assertTrue(identified.fitTestsAvailable());assertEquals(2,identified.degreesOfFreedom());
        assertFalse(redundant.informationAvailable());assertFalse(redundant.fitTestsAvailable());assertEquals(-1,redundant.degreesOfFreedom());
        for(double diagnostic:new double[]{redundant.chiSquare(),redundant.pValue(),redundant.cfi(),redundant.tli(),redundant.rmsea(),redundant.aic(),redundant.bic()})assertTrue(Double.isNaN(diagnostic));
    }
    @Test void meanInclusiveSrmrMatchesLavaanAndAnalyticDefinition() throws Exception {
        SemModel model=SemModel.builder("x","y").fixedVariance("x",1).fixedVariance("y",1).fixedIntercept("x",0).fixedIntercept("y",0).build();
        double[][] data={{0,0},{0,2},{2,0},{2,2}};Properties p=reference();
        SemFitResult fit=Sem.fit(data,model);assertEquals(value(p,"srmr.fixedMean"),fit.srmr(),1e-14);
        assertEquals(Math.sqrt(2./5),fit.srmr(),1e-14);assertEquals(value(p,"srmr.fixedMean.chi"),fit.chiSquare(),1e-14);
        double[][] missing=Arrays.copyOf(data,8);missing[4]=new double[]{0,Double.NaN};missing[5]=new double[]{2,Double.NaN};missing[6]=new double[]{Double.NaN,0};missing[7]=new double[]{Double.NaN,2};
        assertEquals(Math.sqrt(2./5),SemFiml.fit(missing,model).fit().srmr(),1e-10);
    }
    @Test void stableAuxiliaryMomentsAndSingularH1NeverDestroyTargetFit() throws Exception {
        double m=1e8;double[][] data={{m-1,m-1},{m-1,m+1},{m+1,m-1},{m+1,m+1},{m-1,Double.NaN},{m+1,Double.NaN},{Double.NaN,m-1},{Double.NaN,m+1}};
        SemModel model=SemModel.builder("x","y").intercept("x",m).intercept("y",m).build();
        SemFitResult fit=SemFiml.fit(data,model).fit();Properties p=reference();
        assertTrue(fit.converged());assertTrue(fit.informationAvailable());assertTrue(fit.fitTestsAvailable());
        assertEquals(value(p,"largeLocationLL"),fit.logLikelihood(),2e-10);
        assertEquals(1,fit.impliedCovariance()[0],1e-10);assertEquals(1,fit.impliedCovariance()[3],1e-10);
        SemModel constrained=SemModel.builder("x","y").fixedIntercept("x",0).fixedIntercept("y",0).build();
        SemFitResult singular=Sem.fit(new double[][]{{2,-1},{2,1},{2,-1},{2,1}},constrained);
        assertTrue(singular.converged());assertTrue(singular.informationAvailable());
        assertFalse(singular.fitTestsAvailable());assertTrue(Double.isNaN(singular.pValue()));
        assertEquals(value(p,"singularSampleLL"),singular.logLikelihood(),1e-10);
        assertEquals(4,singular.parameter("x~~x").estimate(),1e-7);
    }
}
