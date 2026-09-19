/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

final class MultipleImputationTest {
    @Test void binaryBootstrapRetainsObservedDataUncertainty() {
        double[][] data = new double[1000][2];
        for (int i=0;i<data.length;i++) { data[i][0]=i<50?0:i<100?1:Double.NaN; data[i][1]=7; }
        var completed=MiceImputer.impute(data,new VariableType[]{VariableType.BINARY,VariableType.CONTINUOUS},
            new MiceOptions(500,1,5,271828,1e-6)).datasets();
        var pooled=poolMeans(completed);
        // Observed n=100 Bernoulli(.5) information: variance .25/100.
        assertEquals(.0025,pooled.totalVariance(),.0004);
        assertEquals(.5,pooled.estimate(),.01);
    }

    @Test void multinomialBootstrapPreservesUnequalCategoryProbabilities() {
        double[][] data=new double[300][2];
        for(int i=0;i<data.length;i++)data[i][0]=i<60?0:i<90?1:i<100?2:Double.NaN;
        var completed=MiceImputer.impute(data,new VariableType[]{VariableType.CATEGORICAL,VariableType.CONTINUOUS},
            new MiceOptions(300,1,5,271828,1e-6)).datasets();
        double[] probability=new double[3];
        for(var matrix:completed)for(int i=100;i<data.length;i++)probability[(int)matrix[i][0]]+=1.0/(200*completed.size());
        assertArrayEquals(new double[]{.6,.3,.1},probability,.015);
    }

    @Test void tiedDonorsPreserveMarginalMeanAndContinuousUncertainty() {
        for(boolean reverse:new boolean[]{false,true}) {
            double[][] data=new double[200][2];
            for(int i=0;i<data.length;i++)data[i][0]=i<100?(reverse?100-i:i+1):Double.NaN;
            var completed=MiceImputer.impute(data,new VariableType[]{VariableType.CONTINUOUS,VariableType.CONTINUOUS},
                new MiceOptions(300,1,5,271828,1e-6)).datasets();
            double mean=0;double min=100,max=1;
            for(var matrix:completed)for(int i=100;i<data.length;i++) {
                mean+=matrix[i][0]/(100*completed.size());min=Math.min(min,matrix[i][0]);max=Math.max(max,matrix[i][0]);
            }
            assertEquals(50.5,mean,1.0);assertEquals(1,min);assertEquals(100,max);
            // Sample variance of 1..100 divided by the 100 observed values.
            assertEquals(8.416666666666666,poolMeans(completed).totalVariance(),1.7);
        }
    }

    private static RubinPooling.Estimate poolMeans(java.util.List<double[][]> completed) {
        double[] mean=new double[completed.size()],variance=mean.clone();
        for(int k=0;k<mean.length;k++) {
            var data=completed.get(k);int n=data.length;
            for(var row:data)mean[k]+=row[0]/n;
            for(var row:data)variance[k]+=Math.pow(row[0]-mean[k],2)/(n*(n-1.0));
        }
        return RubinPooling.pool(mean,variance,Double.POSITIVE_INFINITY);
    }
    @Test void mixedChainsAreReproducibleAndPreserveObservedCells() {
        double missing = Double.NaN;
        double[][] data = {
            {20, 0, 1}, {21, 0, 1}, {22, 0, 2}, {missing, 1, 2},
            {30, 1, missing}, {31, 1, 2}, {32, missing, 1}, {33, 1, 2}
        };
        VariableType[] types = {VariableType.CONTINUOUS,
            VariableType.BINARY, VariableType.CATEGORICAL};
        MiceOptions options = new MiceOptions(3, 5, 3, 42, 1e-4);
        MiceResult first = MiceImputer.impute(data, types, options);
        MiceResult second = MiceImputer.impute(data, types, options);
        assertEquals(3, first.datasets().size());
        for (int chain = 0; chain < first.datasets().size(); chain++) {
            double[][] completed = first.datasets().get(chain);
            double[][] repeated = second.datasets().get(chain);
            for (int row = 0; row < completed.length; row++) {
                assertArrayEquals(completed[row], repeated[row], 0.0);
                for (int column = 0; column < completed[row].length; column++) {
                    assertFalse(Double.isNaN(completed[row][column]));
                    if (!Double.isNaN(data[row][column]))
                        assertEquals(data[row][column], completed[row][column], 0.0);
                }
            }
            assertTrue(completed[4][2] == 1.0 || completed[4][2] == 2.0);
        }
        assertEquals(1, first.diagnostics().get(0).missingCount());
    }

    @Test void rubinPoolingMatchesHandCalculationAndFiniteDfCorrection() {
        RubinPooling.Estimate pooled = RubinPooling.pool(
            new double[] {1, 2, 3}, new double[] {4, 4, 4}, 20);
        assertEquals(2.0, pooled.estimate(), 0.0);
        assertEquals(4.0, pooled.withinVariance(), 0.0);
        assertEquals(1.0, pooled.betweenVariance(), 0.0);
        assertEquals(16.0 / 3.0, pooled.totalVariance(), 1e-15);
        assertTrue(pooled.degreesOfFreedom() > 0.0
            && pooled.degreesOfFreedom() < 20.0);
        assertTrue(pooled.fractionMissingInformation() > 0.0
            && pooled.fractionMissingInformation() < 1.0);
    }
}
