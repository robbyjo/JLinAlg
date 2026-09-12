/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class MetaArrayTest {
    @Test void primitiveRegressionRetainsEstimatesWhenEffectUnitsChange() {
        double[][] mods={{-1},{0},{1},{2}};
        for(var estimator:TauSquaredEstimator.values()) {
            var options=MetaAnalysisOptions.builder().tauSquaredEstimator(estimator).build();
            var reference=MetaRegression.fit(new double[]{.2,.5,.1,.7},new double[]{.1,.2,.15,.25},
                mods,List.of("dose"),true,options,BackendPolicy.CPU);
            for(double unit:new double[]{1e-6,1e6}) {
                var actual=MetaRegression.fit(new double[]{.2*unit,.5*unit,.1*unit,.7*unit},
                    new double[]{.1*unit,.2*unit,.15*unit,.25*unit},mods,List.of("dose"),true,options,BackendPolicy.CPU);
                for(int i=0;i<2;i++) assertEquals(reference.beta()[i],actual.beta()[i]/unit,1e-7);
                assertEquals(reference.tauSquared(),actual.tauSquared()/(unit*unit),1e-7);
            }
        }
    }
    @Test void missingBatchMatchesCompleteCaseScalarModelsAndInference() {
        double na=Double.NaN;
        double[] beta={.2,na,.5,.1,.7},se={.1,na,.2,.15,.25};
        var studies=List.of(new MetaStudy("a",.2,.1),new MetaStudy("b",.5,.2),
            new MetaStudy("c",.1,.15),new MetaStudy("d",.7,.25));
        var prepared=MetaAnalysis.prepareBatch(beta,se,1,5,1);
        assertEquals("+?+++",prepared.direction(0));assertArrayEquals(new int[]{4},prepared.cohortCounts());
        for(var method:MetaAnalysisMethod.values())for(var tau:TauSquaredEstimator.values())for(var inference:MetaInferenceMethod.values()) {
            var options=MetaAnalysisOptions.builder().method(method).tauSquaredEstimator(tau).inferenceMethod(inference).build();
            var fit=prepared.fit(options);var reference=MetaAnalysis.fit(studies,options,BackendPolicy.CPU);
            assertEquals(reference.pooledEffectSize(),fit.pooledEffectSizes()[0],3e-8);
            assertEquals(reference.standardError(),fit.standardErrors()[0],3e-8);
            assertEquals(reference.pValue(),fit.pValues()[0],3e-8);
            assertEquals(reference.tauSquared(),fit.tauSquared()[0],3e-8);
        }
    }
    @Test void oneAndZeroCohortRowsHaveExplicitUnavailableHeterogeneity() {
        double na=Double.NaN;
        var prepared=MetaAnalysis.prepareBatch(new double[]{.2,na,na,na},new double[]{.1,na,na,na},2,2,1);
        var fit=prepared.fit(MetaAnalysisOptions.builder().inferenceMethod(MetaInferenceMethod.HARTUNG_KNAPP).build());
        assertEquals(.2,fit.pooledEffectSizes()[0]);assertEquals(.1,fit.standardErrors()[0]);
        assertEquals(.0455002638963584,fit.pValues()[0],1e-14);
        assertTrue(Double.isNaN(fit.tauSquared()[0]));assertTrue(Double.isNaN(fit.cochranQ()[0]));
        assertTrue(Double.isNaN(fit.pooledEffectSizes()[1]));
        assertTrue(Double.isNaN(MetaAnalysis.prepareBatch(new double[]{.2,na},new double[]{.1,na},1,2,2)
            .fit(MetaAnalysisOptions.fixedEffect()).pooledEffectSizes()[0]));
        assertThrows(IllegalArgumentException.class,()->MetaAnalysis.prepareBatch(new double[]{.2,na},new double[]{.1,.1},1,2,1));
    }
    @Test void parallelMissingRowsUseTheirOwnDegreesOfFreedom() {
        int n=4097;double[] beta=new double[n*3],se=new double[n*3];
        for(int i=0;i<n;i++)for(int j=0;j<3;j++) {
            beta[i*3+j]=j==2&&i%2==0?Double.NaN:.2+j*.1;
            se[i*3+j]=Double.isNaN(beta[i*3+j])?Double.NaN:.1;
        }
        var batch=MetaAnalysis.prepareBatch(beta,se,n,3,1);
        var options=MetaAnalysisOptions.builder().method(MetaAnalysisMethod.FIXED_EFFECT).inferenceMethod(MetaInferenceMethod.STUDENT_T).build();
        assertArrayEquals(batch.fit(options,1).pValues(),batch.fit(options,3).pValues());
    }
    @Test void randomBatchEstimationRespectsEffectUnits() {
        for(var estimator:TauSquaredEstimator.values())for(double unit:new double[]{1e-6,1,1e6}) {
            double[] beta={-2*unit,-unit,0,unit,2*unit},se={.1*unit,.1*unit,.1*unit,.1*unit,.1*unit};
            var fit=MetaAnalysis.prepareBatch(beta,se,1,5,1).fit(MetaAnalysisOptions.builder().tauSquaredEstimator(estimator).build());
            assertEquals(2.49,fit.tauSquared()[0]/(unit*unit),3e-6);
        }
    }
    @Test void primitiveRegressionMatchesStudyApi() {
        double[] effects={.2,.5,.1,.7},se={.1,.2,.15,.25};double[][] mods={{-1},{0},{1},{2}};
        var studies=new ArrayList<MetaStudy>();for(int i=0;i<4;i++)studies.add(new MetaStudy("c"+i,effects[i],se[i]));
        for(var method:MetaAnalysisMethod.values())for(var estimator:TauSquaredEstimator.values()) {
            var options=MetaAnalysisOptions.builder().method(method).tauSquaredEstimator(estimator).build();
            var array=MetaRegression.fit(effects,se,mods,List.of("dose"),true,options,BackendPolicy.CPU);
            var scalar=MetaRegression.fit(studies,mods,List.of("dose"),true,options,BackendPolicy.CPU);
            assertArrayEquals(scalar.beta(),array.beta());assertArrayEquals(scalar.standardErrors(),array.standardErrors());
        }
    }
}
