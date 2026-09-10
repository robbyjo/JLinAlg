/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class MrAuditBoundaryTest {
    @Test void ivwReconstructionDoesNotOverflowAFiniteEffectOrStandardError() {
        var values=List.of(instrument(0,1e-200,1,1e200,1e200),instrument(1,1e-200,1,1,1));
        var fit=MendelianRandomization.ivw(values,false,.95);
        assertEquals(1,fit.estimate()/1e200,1e-13);
        assertEquals(1,fit.standardError()/1e200,1e-13);
        assertEquals(1,fit.statistic(),1e-13);
    }
    private HarmonizedInstrument instrument(int i,double x,double sx,double y,double sy) {
        return new HarmonizedInstrument("s"+i,"A","C",x,sx,y,sy,.2,.2,false,false);
    }
    @Test void dispersionMomentTruncatesAfterAveragingNotBefore() {
        var values=new ArrayList<HarmonizedInstrument>();
        for(int i=0;i<8;i++)values.add(instrument(i,1,1e-6,i%4==0?Math.sqrt(2):i%4==1?-Math.sqrt(2):0,1));
        var fit=RobustMendelianRandomization.raps(values);
        assertTrue(fit.converged());
        assertEquals(0,fit.estimate().estimate(),1e-15);
        assertEquals(0,fit.overdispersion(),1e-14);
        assertTrue(fit.iterations()<=100);
    }
    @Test void mixtureUsesLogDensitiesBelowOrdinaryProbabilityUnderflow() throws Exception {
        var values=List.of(instrument(0,1,1e-6,-400,1),instrument(1,1,1e-6,0,1),instrument(2,1,1e-6,400,1));
        var fit=ContaminationMixture.fit(values,101);
        var ref=new Properties();
        try(var in=getClass().getResourceAsStream("/r-reference/genetic-audit-reference.properties")){assertNotNull(in);ref.load(in);}
        assertEquals(Double.parseDouble(ref.getProperty("mixture.0")),fit.estimate().estimate(),1e-12);
        assertEquals(Double.parseDouble(ref.getProperty("mixture.1")),fit.logLikelihood(),1e-8);
        assertEquals(Double.parseDouble(ref.getProperty("mixture.2")),fit.estimate().standardError(),1e-10);
    }
    @Test void winnerMatchesIndependentRConditionalScoreIncludingUnderflowingSelection() {
        double[] t={5.46,6,40.01,41};
        double[] expected={.36935412377725574,4.570504684075494,.10718042663600426,40.48105838703463};
        for(int i=0;i<t.length;i++) {
            double threshold=i<2?5.45:40;
            assertEquals(expected[i],WinnerCurseCorrection.correct(t[i],1,threshold),3e-11);
            assertEquals(-expected[i]*1e-100,WinnerCurseCorrection.correct(-t[i]*1e-100,1e-100,threshold),3e-110);
        }
        assertThrows(IllegalArgumentException.class,()->WinnerCurseCorrection.correct(4,1,5));
    }
    @Test void eggerDoesNotLoseSmallExposureDifferencesAtLargeOffsets() {
        var values=new ArrayList<HarmonizedInstrument>();
        for(int i=0;i<6;i++){double x=1e8+i;values.add(instrument(i,x,1,2*x+3,1));}
        var fit=MendelianRandomization.egger(values,.95);
        assertEquals(2,fit.slope().estimate(),1e-15);
        assertEquals(3,fit.intercept(),1e-15);
        assertEquals(0,fit.slope().cochranQ(),1e-15);
        assertEquals(1/Math.sqrt(17.5),fit.slope().standardError(),1e-15);
    }
    @Test void ivwAndSteigerAreInvariantToExtremeUnits() {
        for(double scale:new double[]{1,1e-200,1e200}) {
            var values=new ArrayList<HarmonizedInstrument>();
            for(int i=0;i<4;i++) values.add(instrument(i,(i+1)*scale,scale,(i+1)*2*scale,2*scale));
            var fit=MendelianRandomization.ivw(values,false,.95);
            assertEquals(2,fit.estimate(),1e-14);
            assertEquals(2/Math.sqrt(30),fit.standardError(),1e-14);
            var steiger=SteigerFiltering.analyze(values,new double[]{100,100,100,100},new double[]{100,100,100,100});
            double expected=0;for(int i=1;i<=4;i++)expected+=i*i/(i*i+98.0);
            assertEquals(expected,steiger.exposureVarianceExplained(),1e-15);
            assertEquals(expected,steiger.outcomeVarianceExplained(),1e-15);
        }
    }
    @Test @SuppressWarnings("deprecation")
    void multivariableMarginalStrengthIsNotPresentedAsConditionalStrength() {
        var values=new ArrayList<MultivariableInstrument>();
        for(int i=0;i<8;i++)values.add(new MultivariableInstrument("s"+i,new double[]{i+1,2*(i+1)+.01*Math.sin(i)},new double[]{.1,.1},.5*(i+1),.1));
        var fit=MultivariableMendelianRandomization.fit(values,List.of("a","b"),false,BackendPolicy.CPU);
        assertEquals(2550,fit.marginalFStatistics()[0],1e-10);
        assertThrows(UnsupportedOperationException.class,fit::conditionalFStatistics);
    }

    @Test @SuppressWarnings("deprecation")
    void covarianceAwareFitReportsTrueConditionalStrength() {
        var values=new ArrayList<MultivariableInstrument>();
        var sampling=new ArrayList<double[][]>();
        int n=10;
        double[][] outcomeCovariance=new double[n][n];
        for(int i=0;i<n;i++){
            double x1=.08+.02*i;
            double x2=.18-.012*i+.025*Math.sin(i);
            values.add(new MultivariableInstrument("c"+i,
                new double[]{x1,x2},new double[]{.02,.03},
                .6*x1-.25*x2,.04));
            sampling.add(new double[][]{{.0004,.00012},{.00012,.0009}});
            outcomeCovariance[i][i]=.0016;
        }
        var fit=MultivariableMendelianRandomization.generalizedFit(values,
            List.of("x1","x2"),false,outcomeCovariance,sampling,
            BackendPolicy.CPU);
        assertArrayEquals(new double[]{.6,-.25},fit.beta(),1e-10);
        assertTrue(fit.conditionalStrengthAvailable());
        for(double value:fit.conditionalFStatistics())
            assertTrue(Double.isFinite(value)&&value>0);
        assertFalse(java.util.Arrays.equals(fit.marginalFStatistics(),
            fit.conditionalFStatistics()));
    }
}
