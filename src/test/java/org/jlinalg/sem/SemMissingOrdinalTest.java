/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.junit.jupiter.api.Test;

class SemMissingOrdinalTest {
    @Test void availablePairsMatchLavaanWithVaryingCoverage() throws Exception {
        Path root=Path.of("src/test/resources/feasible-extensions");
        int[][] data=Files.readAllLines(root.resolve("ordinal-missing.csv")).stream()
            .map(s->Arrays.stream(s.split(",")).mapToInt(Integer::parseInt).toArray()).toArray(int[][]::new);
        SemModel model=SemJointTest.ordinalModel();
        var fit=SemOrdinal.fitPairwiseMissing(data,new int[]{3,3,3,3},model);
        assertTrue(fit.converged(),"score="+fit.scoreNorm());assertTrue(fit.informationAvailable());
        Properties ref=new Properties();try(var in=Files.newInputStream(root.resolve("ordinal-missing.properties"))){ref.load(in);}
        double[] beta=numbers(ref,"estimates"),covariance=numbers(ref,"covariance");
        assertArrayEquals(beta,fit.parameters().stream().mapToDouble(SemParameterEstimate::estimate).toArray(),4e-5);
        assertArrayEquals(covariance,fit.parameterCovariance(),4e-5);
        int[][] extra=Arrays.copyOf(data,data.length+2);extra[data.length]=new int[]{-1,-1,-1,-1};extra[data.length+1]=new int[]{1,-1,-1,-1};
        var same=SemOrdinal.fitPairwiseMissing(extra,new int[]{3,3,3,3},model);
        assertEquals(data.length,same.observations());assertEquals(fit.pairwiseLogLikelihood(),same.pairwiseLogLikelihood(),1e-10);
        int[][] repeated=new int[2*extra.length][];int[] ids=new int[repeated.length];
        for(int i=0;i<extra.length;i++){repeated[2*i]=repeated[2*i+1]=extra[i];ids[2*i]=ids[2*i+1]=i;}
        var cluster=SemOrdinal.fitPairwiseMissing(repeated,new int[]{3,3,3,3},model,SemOptions.defaults(),ids);
        assertArrayEquals(fit.parameterCovariance(),cluster.parameterCovariance(),1e-8);
        assertThrows(IllegalArgumentException.class,()->SemOrdinal.fit(data,new int[]{3,3,3,3},model));
        var stopped=SemOrdinal.fitPairwiseMissing(data,new int[]{3,3,3,3},model,
            new SemOptions(100,1e-16,org.jlinalg.model.MissingDataPolicy.ERROR),null);
        assertFalse(stopped.converged());assertFalse(stopped.informationAvailable());
        assertTrue(Arrays.stream(stopped.parameterCovariance()).allMatch(Double::isNaN));
        assertTrue(stopped.parameters().stream().allMatch(p->Double.isNaN(p.standardError()) && Double.isNaN(p.pValue())));
    }
    @Test void completeRowsGiveIdenticalResultsThroughBothApis() {
        List<int[]> rows=new ArrayList<>();
        for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int r=0;r<(a==b?15:35);r++)rows.add(new int[]{a,b});
        int[][] data=rows.toArray(int[][]::new);
        SemModel model=SemModel.builder("x","y").regression("b","y","x",.2).fixedVariance("x",1).fixedVariance("y",1).build();
        var complete=SemOrdinal.fit(data,new int[]{2,2},model);
        var available=SemOrdinal.fitPairwiseMissing(data,new int[]{2,2},model);
        assertTrue(available.converged());assertEquals(complete.pairwiseLogLikelihood(),available.pairwiseLogLikelihood(),1e-12);
        assertArrayEquals(complete.parameterCovariance(),available.parameterCovariance(),1e-12);
    }
    @Test void disjointBinaryPairsHaveAnalyticProbitCorrelationAndVariance() {
        List<int[]> data=new ArrayList<>();int[][] pairs={{0,1},{0,2},{1,2}};
        int[] equalCounts={20,60,45},unequalCounts={30,40,55};
        for(int h=0;h<3;h++)for(int a=0;a<2;a++)for(int b=0;b<2;b++)
            for(int r=0;r<(a==b?equalCounts[h]:unequalCounts[h]);r++) {
                int[] row={-1,-1,-1};row[pairs[h][0]]=a;row[pairs[h][1]]=b;data.add(row);
            }
        SemModel model=SemModel.builder("x","y","z").fixedVariance("x",1).fixedVariance("y",1).fixedVariance("z",1)
            .covariance("xy","x","y",0).covariance("xz","x","z",0).covariance("yz","y","z",0).build();
        var fit=SemOrdinal.fitPairwiseMissing(data.toArray(int[][]::new),new int[]{2,2,2},model);
        assertTrue(fit.converged());assertTrue(fit.informationAvailable());
        for(int h=0;h<3;h++) {
            double n=2.0*(equalCounts[h]+unequalCounts[h]),probability=2*equalCounts[h]/n;
            double rho=Math.sin(Math.PI*(probability-.5));
            assertEquals(rho,fit.parameters().get(h).estimate(),2e-7);
            double variance=Math.PI*Math.PI*(1-rho*rho)*probability*(1-probability)/n;
            assertEquals(variance,fit.parameterCovariance()[h*6+h],1e-8);
        }
        for(double[] thresholds:fit.thresholds())assertEquals(0,thresholds[0],1e-8);
        double[] copy=fit.parameterCovariance();copy[0]=-1;assertTrue(fit.parameterCovariance()[0]>0);
    }
    @Test void invalidCoverageMissingCodesAndFailedFitsDoNotProvideInference() {
        SemModel model=SemModel.builder("x","y","z").fixedVariance("x",1).fixedVariance("y",1).fixedVariance("z",1).build();
        assertThrows(IllegalArgumentException.class,()->SemOrdinal.fitPairwiseMissing(new int[][]{{0,0,-1},{1,1,-1},{0,1,-1}},new int[]{2,2,2},model));
        assertThrows(IllegalArgumentException.class,()->SemOrdinal.fitPairwiseMissing(new int[][]{{-2,0,0},{1,1,1},{0,1,0}},new int[]{2,2,2},model));
        assertThrows(IllegalArgumentException.class,()->SemOrdinal.fitPairwiseMissing(new int[][]{{0,0,0},{1,1,1},{0,1,0}},new int[]{2,2,2},model,SemOptions.defaults(),new int[2]));
    }
    private static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
}
