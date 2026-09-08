/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;

/** Exact LP objectives and KKT certificates, including nonunique quantreg optima. */
public class ExactQuantileTest {
    private static final Path ROOT = Path.of("src/test/resources/exact-quantile");
    record Data(double[] y,double[][] x) { }
    static Data data(String name) throws Exception {
        List<String> lines = Files.readAllLines(ROOT.resolve(name+".csv"));
        int n=lines.size(),p=lines.get(0).split(",").length-1; double[] y=new double[n]; double[][] x=new double[n][p];
        for(int i=0;i<n;i++) { String[] row=lines.get(i).split(",");y[i]=Double.parseDouble(row[0]);for(int j=0;j<p;j++)x[i][j]=Double.parseDouble(row[j+1]); }
        return new Data(y,x);
    }
    @Test void frozenRqObjectivesAndUniqueCoefficientsAgree() throws Exception {
        Map<String,Double> coefficients=new HashMap<>();
        for(String line:Files.readAllLines(ROOT.resolve("coefficients.csv")).subList(1,Files.readAllLines(ROOT.resolve("coefficients.csv")).size())) {
            String[] a=line.split(",");coefficients.put(a[0]+":"+Double.parseDouble(a[1])+":"+a[2],Double.parseDouble(a[3]));
        }
        List<String> refs=Files.readAllLines(ROOT.resolve("reference.csv"));
        double largestObjectiveError=0;
        for(String line:refs.subList(1,refs.size())) {
            String[] a=line.split(",");Data d=data(a[0]);double tau=Double.parseDouble(a[1]),expected=Double.parseDouble(a[2]);
            var result=QuantileLinearProgram.solve(d.y,d.x,tau);
            assertTrue(result.fit().converged(),a[0]+" tau="+tau+" "+result.certificate());
            double error=Math.abs(expected-result.fit().objective());largestObjectiveError=Math.max(largestObjectiveError,error);
            assertEquals(expected,result.fit().objective(),2e-7*Math.max(1,expected),a[0]+" objective tau="+tau);
            checkKkt(d,result,tau);
            // rq may choose a different basic solution on a face; only compare coefficients where unique.
            if(a[3].equals("1"))for(int j=0;j<d.x[0].length;j++)
                assertEquals(coefficients.get(a[0]+":"+tau+":"+j),result.fit().beta()[j],2e-6,a[0]+" beta "+j+" tau="+tau);
        }
        System.out.println("quantreg cases="+(refs.size()-1)+" max objective error="+largestObjectiveError);
    }
    @Test void evenMedianAndTiesCertifyTheOptimalFace() {
        var result=QuantileLinearProgram.solve(new double[]{1,2},new double[][]{{1},{1}},.5);
        assertTrue(result.fit().converged());assertTrue(result.fit().beta()[0]>=1&&result.fit().beta()[0]<=2);
        assertEquals(.5,result.fit().objective(),1e-10);
        var tied=QuantileRegression.fitExact(new double[]{0,0,0,1,1,1},new double[][]{{1},{1},{1},{1},{1},{1}},.5);
        assertTrue(tied.converged());assertEquals(1.5,tied.objective(),1e-9);
    }
    @Test void scalingColumnsAndEffectsPreservesTheLp() throws Exception {
        Data d=data("regular");var base=QuantileLinearProgram.solve(d.y,d.x,.5);
        double[][] scaled=Arrays.stream(d.x).map(double[]::clone).toArray(double[][]::new);
        double[] scales={1e12,1e-9,1e5};for(double[] row:scaled)for(int j=0;j<3;j++)row[j]*=scales[j];
        var fit=QuantileLinearProgram.solve(d.y,scaled,.5);assertTrue(fit.fit().converged(),fit.certificate().toString());
        for(int j=0;j<3;j++)assertEquals(base.fit().beta()[j],fit.fit().beta()[j]*scales[j],1e-7);
        assertEquals(base.fit().objective(),fit.fit().objective(),1e-8);
        double[] ys=d.y.clone();for(int i=0;i<ys.length;i++)ys[i]*=1e-9;
        var small=QuantileLinearProgram.solve(ys,d.x,.5);assertTrue(small.fit().converged());
        assertEquals(base.fit().objective(),small.fit().objective()/1e-9,1e-7);
        for(int j=0;j<3;j++)assertEquals(base.fit().beta()[j],small.fit().beta()[j]/1e-9,1e-7);
    }
    @Test void largeInterceptShiftPreservesLossButMustPassOriginalUnitCertificate() throws Exception {
        Data d=data("regular");var base=QuantileLinearProgram.solve(d.y,d.x,.5);
        double[] shifted=d.y.clone();for(int i=0;i<shifted.length;i++)shifted[i]+=1e8;
        var strict=QuantileLinearProgram.solve(shifted,d.x,.5);
        // At this offset one response ulp is 1.49e-8: an accurate internal LP
        // solution need not certify its reconstructed original-unit residuals.
        assertFalse(strict.fit().converged(),strict.certificate().toString());
        assertEquals(base.fit().objective(),strict.fit().objective(),1e-6);
        var relaxed=QuantileLinearProgram.solve(shifted,d.x,.5,new QuantileLinearProgram.Options(200,1e-8));
        assertTrue(relaxed.fit().converged(),relaxed.certificate().toString());
        assertEquals(base.fit().beta()[0],relaxed.fit().beta()[0]-1e8,1e-7);
        for(int j=1;j<3;j++)assertEquals(base.fit().beta()[j],relaxed.fit().beta()[j],1e-7);
        assertEquals(base.fit().objective(),relaxed.fit().objective(),1e-6);
        checkKkt(new Data(shifted,d.x),relaxed,.5);
        System.out.println("large-shift strict="+strict.certificate()+" explicit tolerance="+relaxed.certificate());
    }
    @Test void budgetsAndInvalidSystemsNeverReportSuccess() throws Exception {
        Data d=data("regular");var limited=QuantileLinearProgram.solve(d.y,d.x,.9,new QuantileLinearProgram.Options(1,1e-12));
        assertFalse(limited.fit().converged());assertEquals(1,limited.fit().iterations());
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(new double[]{1,2,3},new double[][]{{1,2},{1,2},{1,2}},.5));
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(new double[]{1,2},new double[][]{{1},{Double.NaN}},.5));
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(new double[]{1,Double.POSITIVE_INFINITY},new double[][]{{1},{1}},.5));
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(new double[]{1,2},new double[][]{{1},{1,2}},.5));
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(new double[]{1,2},new double[][]{{1,0,1},{0,1,1}},.5));
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(d.y,d.x,0));
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(d.y,d.x,Double.NaN));
        assertThrows(IllegalArgumentException.class,()->QuantileRegression.fitExact(new double[]{1,2,3},new double[][]{{1,1},{1,1+1e-14},{1,1-1e-14}},.5));
    }
    @Test void deterministicAndDefensiveResults() throws Exception {
        Data d=data("ties");double[] original=d.y.clone();double[][] x=Arrays.stream(d.x).map(double[]::clone).toArray(double[][]::new);
        var a=QuantileLinearProgram.solve(d.y,d.x,.25);var b=QuantileLinearProgram.solve(d.y,d.x,.25);
        assertArrayEquals(a.fit().beta(),b.fit().beta());assertArrayEquals(a.dualWeights(),b.dualWeights());
        assertArrayEquals(original,d.y);for(int i=0;i<x.length;i++)assertArrayEquals(x[i],d.x[i]);
        double saved=a.dualWeights()[0];double[] exported=a.dualWeights();exported[0]=99;assertEquals(saved,a.dualWeights()[0]);
    }
    @Test void exactWrapperRetainsTheExistingSmoothedEntryPoint() {
        double[] y={1,2,3};double[][] x={{1},{1},{1}};
        var exact=QuantileRegression.fitExact(y,x,.5);
        var smooth=QuantileRegression.fit(y,x,.5);
        assertTrue(exact.converged());assertEquals(1,exact.objective(),1e-9);
        assertTrue(smooth.objective()>exact.objective());
    }
    static void checkKkt(Data d,QuantileLinearProgram.Result result,double tau) {
        double[] z=result.dualWeights(),beta=result.fit().beta();double primal=0,dual=0,maxError=0;
        for(int i=0;i<d.y.length;i++) {
            double fitted=0;for(int j=0;j<beta.length;j++)fitted+=d.x[i][j]*beta[j];
            double residual=d.y[i]-fitted;primal+=residual>=0?tau*residual:(tau-1)*residual;dual+=d.y[i]*z[i];
            assertTrue(z[i]>=tau-1-1e-10&&z[i]<=tau+1e-10);
            maxError=Math.max(maxError,Math.abs(residual)*(residual>=0?tau-z[i]:1-tau+z[i]));
        }
        for(int j=0;j<beta.length;j++) {double sum=0,norm=0;for(int i=0;i<d.y.length;i++){sum+=d.x[i][j]*z[i];norm+=Math.abs(d.x[i][j]);}assertEquals(0,sum,1e-8*Math.max(1,norm));}
        assertEquals(primal,dual,1e-7*Math.max(1,primal));assertTrue(maxError<1e-7*Math.max(1,primal));
        // Independent naive X*beta and the solver's compensated dot can differ
        // by a response ulp per row after a large intercept translation.
        double rounding=d.y.length*Math.ulp(Arrays.stream(d.y).map(Math::abs).max().orElse(1));
        assertEquals(primal,result.fit().objective(),1e-10*Math.max(1,primal)+rounding);
    }
    public static void main(String[] args) throws Exception {
        ExactQuantileTest test=new ExactQuantileTest();int count=0;
        for(var method:ExactQuantileTest.class.getDeclaredMethods())if(method.isAnnotationPresent(Test.class)){System.out.println("Checking "+method.getName());method.invoke(test);count++;}
        System.out.println(count+" exact quantile tests passed");
    }
}
