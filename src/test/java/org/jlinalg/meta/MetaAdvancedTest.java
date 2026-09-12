/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.junit.jupiter.api.Test;

class MetaAdvancedTest {
    @Test void sparseHeterogeneityIsInvariantToEffectUnits() {
        for(var estimator:List.of(TauSquaredEstimator.REML,TauSquaredEstimator.PAULE_MANDEL))
            for(double unit:new double[]{1e-100,1e-6,1,1e6,1e100}) {
                var studies=new java.util.ArrayList<MetaStudy>();
                double[] precision=new double[5];
                for(int i=0;i<5;i++) {
                    studies.add(new MetaStudy("s"+i,(i-2)*unit,.1*unit));
                    precision[i]=100/(unit*unit);
                }
                var result=SparseMetaAnalysis.fit(studies,new SparsePrecisionMatrix(5,
                    new int[]{0,1,2,3,4,5},new int[]{0,1,2,3,4},precision),
                    MetaAnalysisOptions.builder().tauSquaredEstimator(estimator).build(),BackendPolicy.CPU);
                // Equal sampling variances: REML and PM both give s_y^2 - s_e^2.
                assertEquals(2.49,result.tauSquared()/(unit*unit),2e-6,estimator+" unit="+unit);
                assertEquals(Math.sqrt(.5),result.standardError()/unit,3e-7);
                assertEquals(0,result.pooledEffectSize()/unit,2e-14);
            }
    }

    @Test void sparseHeterogeneityRetainsTheZeroBoundaryInDifferentUnits() {
        for(var estimator:List.of(TauSquaredEstimator.REML,TauSquaredEstimator.PAULE_MANDEL))
            for(double unit:new double[]{1e-6,1,1e6})
                for(double spread:new double[]{0,.01}) {
                    var studies=new java.util.ArrayList<MetaStudy>();
                    double[] precision=new double[5];
                    for(int i=0;i<5;i++) {
                        studies.add(new MetaStudy("s"+i,(2+(i-2)*spread)*unit,.1*unit));
                        precision[i]=100/(unit*unit);
                    }
                    var result=SparseMetaAnalysis.fit(studies,new SparsePrecisionMatrix(5,
                        new int[]{0,1,2,3,4,5},new int[]{0,1,2,3,4},precision),
                        MetaAnalysisOptions.builder().tauSquaredEstimator(estimator).build(),BackendPolicy.CPU);
                    assertEquals(0,result.tauSquared());
                    assertEquals(2,result.pooledEffectSize()/unit,2e-14);
                }
    }
    private static final List<MetaStudy> STUDIES = List.of(
        new MetaStudy("a", 0.2, 0.1), new MetaStudy("b", 0.5, 0.2),
        new MetaStudy("c", 0.1, 0.15), new MetaStudy("d", 0.7, 0.25));

    @Test
    void constructsCommonEffectSizes() {
        assertEquals(Math.log(2.25), MetaEffectSizes.logOddsRatio(20, 80, 10, 90).effectSize(), 1e-12);
        assertEquals(1.0 / Math.sqrt(4.0), MetaEffectSizes.fisherZCorrelation(0.5, 7).standardError(), 1e-12);
        assertTrue(MetaEffectSizes.standardizedMeanDifference(2, 1, 1, 1, 20, 20).standardError() > 0);
    }

    @Test
    void correlatedIdentityMatchesFixedEffect() {
        double[][] covariance = new double[4][4];
        for (int i = 0; i < 4; i++) covariance[i][i] = STUDIES.get(i).variance();
        MetaCorrelatedResult result = MetaCorrelatedAnalysis.fit(STUDIES, covariance,
            MetaAnalysisOptions.fixedEffect(), BackendPolicy.CPU);
        assertEquals(MetaAnalysis.fit(STUDIES, MetaAnalysisOptions.fixedEffect(), BackendPolicy.CPU)
            .pooledEffectSize(), result.pooledEffectSize(), 1e-12);
    }

    @Test
    void clusterRobustAndPublicationDiagnosticsAreFinite() {
        MetaClusterRobustResult robust = MetaClusterRobust.fit(STUDIES,
            new String[] {"one", "one", "two", "two"},
            new double[][] {{-1}, {0}, {1}, {2}}, List.of("dose"),
            MetaAnalysisOptions.randomEffects(), BackendPolicy.CPU);
        assertTrue(Double.isFinite(robust.beta()[0]));
        MetaPublicationBiasResult bias = MetaPublicationBias.diagnose(STUDIES);
        assertTrue(Double.isFinite(bias.eggerPValue()));
        assertTrue(Double.isFinite(bias.rankPValue()));
    }

    @Test void sparsePrecisionRemlMatchesDenseCorrelatedAnalysis() {
        int n=STUDIES.size();double variance=.04,rho=.2;
        double[][] covariance=new double[n][n];
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)
            covariance[i][j]=variance*Math.pow(rho,Math.abs(i-j));
        double denominator=variance*(1-rho*rho);
        int[] starts={0,2,5,8,10};
        int[] columns={0,1,0,1,2,1,2,3,2,3};
        double[] values={1/denominator,-rho/denominator,
            -rho/denominator,(1+rho*rho)/denominator,-rho/denominator,
            -rho/denominator,(1+rho*rho)/denominator,-rho/denominator,
            -rho/denominator,1/denominator};
        var precision=new SparsePrecisionMatrix(n,starts,columns,values);
        var options=MetaAnalysisOptions.builder()
            .tauSquaredEstimator(TauSquaredEstimator.REML)
            .maximumIterations(1000).tolerance(1e-11).build();
        var dense=MetaCorrelatedAnalysis.fit(STUDIES,covariance,options,
            BackendPolicy.CPU);
        var sparse=SparseMetaAnalysis.fit(STUDIES,precision,options,
            BackendPolicy.CPU);
        assertEquals(dense.tauSquared(),sparse.tauSquared(),2e-9);
        assertEquals(dense.pooledEffectSize(),sparse.pooledEffectSize(),2e-10);
        assertEquals(dense.covariance()[0],sparse.covariance()[0],2e-10);
        assertEquals(dense.generalizedQ(),sparse.generalizedQ(),2e-10);
    }
}
