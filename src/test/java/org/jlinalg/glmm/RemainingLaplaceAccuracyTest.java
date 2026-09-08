/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.*;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.reml.*;
import org.junit.jupiter.api.Test;

class RemainingLaplaceAccuracyTest {
    static final Path ROOT=Path.of("src/test/resources/r-reference/remaining-mixed");
    static final GlmmLaplaceOptions OPTIONS=new GlmmLaplaceOptions(100,100,1e-8,1,1e-8,100,null);
    static Properties reference() throws Exception {var p=new Properties();try(var reader=Files.newBufferedReader(ROOT.resolve("reference.properties"))){p.load(reader);}return p;}
    static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
    record Data(double[] y,double[] x,double[] offset,List<String> groups) { }
    static Data data(String kind) throws Exception {
        var lines=Files.readAllLines(ROOT.resolve(kind+".tsv"));int n=lines.size()-1;
        double[] y=new double[n],x=new double[2*n],offset=new double[n];List<String> groups=new ArrayList<>();
        for(int i=0;i<n;i++){var a=lines.get(i+1).split("\t");y[i]=Double.parseDouble(a[0]);x[2*i]=1;x[2*i+1]=Double.parseDouble(a[1]);groups.add(a[2]);offset[i]=Double.parseDouble(a[3]);}
        return new Data(y,x,offset,groups);
    }
    static GlmmLaplaceResult fit(Data d,String kind,boolean dense,GlmmLaplaceOptions options) {
        var family=kind.equals("rare")?GlmFamilies.binomial():GlmFamilies.poisson();
        return dense?GlmmLaplace.fit(d.y,d.x,d.y.length,2,family,List.of(VarianceComponent.randomIntercept("g",d.groups)),null,d.offset,options,BackendPolicy.CPU)
            :SparseGlmmLaplace.fit(d.y,d.x,d.y.length,2,family,List.of(RandomEffectTerm.randomIntercept("g",d.groups)),null,d.offset,options,BackendPolicy.CPU);
    }
    @Test void actualMarginalOptimumAndJointNuisanceCovarianceMatchIndependentR() throws Exception {
        var p=reference();
        for(String kind:List.of("rare","poisson"))for(boolean dense:new boolean[]{false,true}) {
            var fit=fit(data(kind),kind,dense,OPTIONS);assertTrue(fit.converged(),kind+" dense="+dense);
            assertArrayEquals(numbers(p,kind+".beta"),fit.beta(),2e-6);
            assertArrayEquals(numbers(p,kind+".variance"),fit.varianceComponents(),3e-6);
            assertEquals(numbers(p,kind+".LL")[0],fit.marginalLogLikelihood(),2e-8);
            assertArrayEquals(numbers(p,kind+".covariance"),fit.fixedEffectCovariance(),3e-6);
        }
    }
    @Test void largeDesignUnitsAndTinyVarianceStartsDoNotHideInteriorOptimum() throws Exception {
        Data d=data("rare");double[] scaled=d.x.clone();for(int i=0;i<d.y.length;i++)scaled[2*i+1]*=1e12;
        var control=new GlmmLaplaceOptions(100,100,1e-8,1,1e-10,100,new double[]{1e-10});
        var fit=fit(new Data(d.y,scaled,d.offset,d.groups),"rare",false,control);var p=reference();
        assertTrue(fit.converged());assertEquals(numbers(p,"rare.beta")[1],fit.beta()[1]*1e12,2e-6);
        assertEquals(numbers(p,"rare.LL")[0],fit.marginalLogLikelihood(),2e-8);
    }
    @Test void covarianceBackTransformRetainsRepresentableSubnormalVariance() throws Exception {
        Data d=data("rare");double[] scaled=d.x.clone();
        for(int i=0;i<d.y.length;i++)scaled[2*i+1]*=1e160;
        var fit=fit(new Data(d.y,scaled,d.offset,d.groups),"rare",false,OPTIONS);
        assertTrue(fit.converged());assertTrue(fit.fixedEffectCovariance()[3]>0);
        assertTrue(fit.standardErrors()[1]>0);
        double reference=numbers(reference(),"rare.covariance")[3];
        assertEquals(reference,fit.fixedEffectCovariance()[3]*1e160*1e160,1e-3);
        assertEquals(numbers(reference(),"rare.beta")[1],fit.beta()[1]*1e160,2e-6);
    }
    @Test void covarianceBackTransformAvoidsIntermediateScaleProductUnderflow() {
        assertEquals(1e300,SparseGlmmLaplace.unscaleCovariance(1e-20,1e-160,1e-160),2e284);
        assertEquals(1e20,SparseGlmmLaplace.unscaleCovariance(1e20,1e-300,1e300),2e4);
        assertEquals(1e-20,SparseGlmmLaplace.unscaleCovariance(1e-20,1e300,1e-300),2e-36);
        assertEquals(-1e20,SparseGlmmLaplace.unscaleCovariance(-1e20,1e-300,1e300),2e4);
    }
    @Test void exhaustionAndUnmodeledDispersionDoNotReturnSuccessfulInference() throws Exception {
        Data d=data("rare");var fit=fit(d,"rare",false,new GlmmLaplaceOptions(1,1,1e-8,1,1e-8,100,null));
        assertFalse(fit.converged());for(double v:fit.fixedEffectCovariance())assertTrue(Double.isNaN(v));
        assertThrows(IllegalArgumentException.class,()->SparseGlmmLaplace.fit(d.y,d.x,d.y.length,2,GlmFamilies.gaussian(),
            List.of(RandomEffectTerm.randomIntercept("g",d.groups)),null,d.offset,OPTIONS,BackendPolicy.CPU));
    }
    @Test void pqlMatchesFixedDispersionWorkingRemlRatherThanMarginalML() throws Exception {
        double[] y={0,1,1,1,2,3,3,4,5,7,8,9};double[][] x=new double[y.length][1];List<Integer> groups=new ArrayList<>();
        for(int i=0;i<y.length;i++){x[i][0]=1;groups.add(i/3);}
        var options=GlmmPqlOptions.builder().maximumIterations(100).relativeTolerance(1e-7).remlOptions(
            RemlOptions.builder().initialVariances(.5).scoreTolerance(1e-7).build()).build();
        var fit=GlmmPql.fit(y,x,GlmFamilies.poisson(),List.of(VarianceComponent.randomIntercept("g",groups)),null,null,options,BackendPolicy.CPU);
        var p=reference();assertTrue(fit.converged(),fit.convergenceMessage());
        assertArrayEquals(numbers(p,"pql.beta"),fit.beta(),2e-6);
        assertArrayEquals(numbers(p,"pql.variance"),fit.varianceComponents(),2e-6);
        assertArrayEquals(numbers(p,"pql.random"),fit.randomLinearPredictor(),2e-6);
    }
    @Test void nuisanceGradientAndInformationStencilsRespectNarrowBounds() {
        LaplaceOptimization.Objective narrow=x->{
            assertTrue(x[1]>=0&&x[1]<=1e-6,"nuisance stencil left its bounds");
            return x[0]*x[0]+x[1];
        };
        assertArrayEquals(new double[]{0,1},LaplaceOptimization.gradient(narrow,new double[]{0,0},
            new double[]{Double.NEGATIVE_INFINITY,0},new double[]{Double.POSITIVE_INFINITY,1e-6}),1e-10);
        LaplaceOptimization.Objective interior=x->{
            assertTrue(x[1]>=0&&x[1]<=1e-4,"information stencil left its bounds");
            return .5*x[0]*x[0]+.5e8*(x[1]-5e-5)*(x[1]-5e-5);
        };
        assertArrayEquals(new double[]{1},LaplaceOptimization.fixedCovariance(interior,new double[]{0,5e-5},
            new double[]{Double.NEGATIVE_INFINITY,0},new double[]{Double.POSITIVE_INFINITY,1e-4},1),1e-9);
    }
    @Test void roundedPositiveBinomialTailRetainsLaplaceDensityAndMode() {
        // At negligible, tightly bounded random variance this approaches the
        // independently solved GLM score: 1e22*logistic(-b-50)-3*logistic(b)=0.
        var controls=new GlmmLaplaceOptions(100,100,1e-8,1,1e-12,1.0001e-12,new double[]{1e-12});
        var fit=SparseGlmmLaplace.fit(new double[]{1,0,0,0},new double[]{1,1,1,1},4,1,
            GlmFamilies.binomial(),List.of(RandomEffectTerm.randomIntercept("g",List.of("a","a","a","a"))),
            new double[]{1e22,1,1,1},new double[]{50,0,0,0},controls,BackendPolicy.CPU);
        assertTrue(fit.converged());
        assertEquals(.17001035679438581,fit.beta()[0],2e-7);
        assertEquals(-3.9724844032274467,fit.marginalLogLikelihood(),2e-9);
    }
    @Test void pqlWorkingModelUsesTailResidualAndUnclippedPrecision() throws Exception {
        var method=GlmmPql.class.getDeclaredMethod("workingModel",double[].class,double[].class,
            double[].class,double[].class,GlmFamily.class);
        method.setAccessible(true);
        Object model=method.invoke(null,new double[]{1},new double[]{50},new double[]{1e22},
            new double[]{0},GlmFamilies.binomial());
        var response=model.getClass().getDeclaredMethod("response");response.setAccessible(true);
        var covariance=model.getClass().getDeclaredMethod("residualCovariance");covariance.setAccessible(true);
        assertEquals(51,((double[])response.invoke(model))[0],1e-14);
        assertEquals(Math.exp(50)/1e22,((double[])covariance.invoke(model))[0],2e-15);
        model=method.invoke(null,new double[]{1},new double[]{0},new double[]{1e-20},
            new double[]{0},GlmFamilies.poisson());
        assertEquals(1e20,((double[])covariance.invoke(model))[0],1e5);
    }
}
