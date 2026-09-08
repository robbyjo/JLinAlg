/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.penalized.PolyhedralSelectiveInference;
import org.junit.jupiter.api.Test;

final class RegressionRParityTest {
    private static final int N=160;
    private static final Map<String,double[]> REF=reference();
    private static Map<String,double[]> reference() {
        Map<String,List<Double>> values=new HashMap<>();
        try(var stream=RegressionRParityTest.class.getResourceAsStream("/regression/r-accuracy.csv");
            var reader=new BufferedReader(new InputStreamReader(Objects.requireNonNull(stream),StandardCharsets.UTF_8))) {
            reader.readLine();String line;while((line=reader.readLine())!=null){String[] fields=line.split(",");values.computeIfAbsent(fields[0],ignored->new ArrayList<>()).add(Double.valueOf(fields[2]));}
        }catch(IOException e){throw new UncheckedIOException(e);}
        Map<String,double[]> result=new HashMap<>();values.forEach((k,v)->result.put(k,v.stream().mapToDouble(Double::doubleValue).toArray()));return result;
    }
    private record Data(double[] z,double[] x,double[][] design,double[] y,double[][] responses,int[] classes) { }
    private static Data data() {
        double[] z=new double[N],x=new double[N],y=new double[N];double[][] design=new double[N][2],responses=new double[N][2];int[] classes=new int[N];
        for(int i=0;i<N;i++) {
            z[i]=-3+6.*i/(N-1);x[i]=Math.cos(i*1.7)+.2*z[i];design[i]=new double[]{1,x[i]};
            y[i]=1+.4*x[i]+Math.sin(z[i])+.15*Math.cos(i*.37);responses[i]=new double[]{y[i],-.5+.8*x[i]+.1*Math.sin(i*.8)};
            double l1=Math.exp(.3+.5*x[i]),l2=Math.exp(-.2-.3*x[i]),u=((i*37+11)%101+.5)/101;
            classes[i]=u<1/(1+l1+l2)?0:u<(1+l1)/(1+l1+l2)?1:2;
        }
        return new Data(z,x,design,y,responses,classes);
    }
    @Test void multivariateOlsMatchesR() {
        var d=data();var fit=MultivariateRegression.fit(d.responses,d.design,BackendPolicy.CPU);
        assertArrayEquals(REF.get("ols_beta"),fit.coefficients(),1e-12);
        assertArrayEquals(REF.get("ols_cov"),fit.residualCovariance(),1e-12);
    }
    @Test void multinomialMatchesNnet() {
        var d=data();var fit=MultinomialRegression.fit(d.classes,d.design,3);
        assertTrue(fit.converged());assertArrayEquals(REF.get("multinomial_beta"),fit.coefficients(),1e-5);
        assertEquals(REF.get("multinomial_loglik")[0],fit.logLikelihood(),1e-9);
    }
    @Test void quantileMatchesSmoothedObjectiveAndApproximatesRq() {
        var d=data();for(double tau:new double[]{.25,.5,.75}) {
            var fit=QuantileRegression.fit(d.y,d.design,tau);
            assertTrue(fit.converged(),"tau="+tau);
            assertArrayEquals(REF.get("smooth_rq_"+tau),fit.coefficients(),1e-5);
            assertArrayEquals(REF.get("rq_"+tau),fit.coefficients(),.003);
        }
    }
    @Test void kernelAndSupersmootherMatchR() {
        var d=data();assertArrayEquals(REF.get("kernel"),KernelRegression.fit(d.z,d.y,.4).fittedValues(),1e-12);
        assertArrayEquals(REF.get("supsmu"),SuperSmoother.fit(d.z,d.y).fittedValues(),1e-10);
        double[] weights=new double[N],periodic=new double[N];for(int i=0;i<N;i++){weights[i]=1+i%3;periodic[i]=(d.z[i]+3)/6;}
        assertArrayEquals(REF.get("supsmu_weighted"),SuperSmoother.fit(d.z,d.y,weights,0,false,5).fittedValues(),1e-10);
        assertArrayEquals(REF.get("supsmu_fixed"),SuperSmoother.fit(d.z,d.y,null,.3,false,0).fittedValues(),1e-10);
        assertArrayEquals(REF.get("supsmu_periodic"),SuperSmoother.fit(periodic,d.y,null,0,true,0).fittedValues(),1e-10);
    }
    @Test void robinsonAndHc3MatchRMatrixAlgebra() {
        var d=data();var fit=PartiallyLinearRegression.fitWithInference(d.y,d.design,d.z,.4,BackendPolicy.CPU);
        assertArrayEquals(REF.get("robinson_beta"),fit.fit().coefficients(),1e-12);
        assertArrayEquals(REF.get("robinson_slope_cov"),fit.slopeCovariance(),1e-12);
        assertArrayEquals(REF.get("robinson_homoskedastic_cov"),fit.homoskedasticSlopeCovariance(),1e-12);
        assertEquals(REF.get("robinson_noise_df")[0],fit.residualNoiseDegreesOfFreedom(),1e-10);
        assertArrayEquals(REF.get("robinson_smooth"),fit.fit().smoothEffect(),1e-12);
        assertEquals(0,Arrays.stream(fit.fit().smoothEffect()).sum(),1e-10);
    }
    @Test void conditionalLassoMatchesSelectiveInference() {
        var d=data();double[][] x=new double[N][4];double[] y=new double[N];
        for(int i=0;i<N;i++) {
            x[i]=new double[]{Math.cos(i*.7)+.2*d.z[i],Math.sin(i*.23)-.1*d.z[i],Math.cos(i*.33),Math.sin(i*.41)};
            y[i]=.5+.9*x[i][0]-.6*x[i][1]+.5*Math.sin(i*1.13);
        }
        var fit=PolyhedralSelectiveInference.fit(y,x,.12,1,.5,true,.95,BackendPolicy.CPU);
        for(int j=0;j<fit.effects().size();j++) {
            var e=fit.effects().get(j);
            assertEquals(REF.get("selective_beta")[j],e.estimate(),1e-8);
            assertEquals(REF.get("selective_p")[j],e.pValue(),1e-8);
            assertEquals(REF.get("selective_ci")[2*j],e.confidenceLower(),1e-5);
            assertEquals(REF.get("selective_ci")[2*j+1],e.confidenceUpper(),1e-5);
            assertEquals(REF.get("selective_limits")[2*j],e.truncationLower(),1e-6);
            assertEquals(REF.get("selective_limits")[2*j+1],e.truncationUpper(),1e-6);
        }
    }
}
