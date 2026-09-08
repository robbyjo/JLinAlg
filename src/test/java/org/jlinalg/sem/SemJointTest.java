/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

/** Frozen lavaan comparisons, including every parameter covariance entry. */
public class SemJointTest {
    public static SemModel latent(boolean means) {
        SemModel.Builder b=SemModel.builder("x1","x2","x3","y1","y2","y3").latent("f","g")
            .fixedLoading("x1","f",1).loading("l2","x2","f",.8).loading("l3","x3","f",1.1)
            .fixedLoading("y1","g",1).loading("l5","y2","g",.9).loading("l6","y3","g",.7)
            .regression("b","g","f",.5).variance("vf","f",1).variance("vg","g",.6);
        if(means)b.meanStructure();return b.build();
    }
    public static SemModel ordinalModel() {
        return SemModel.builder("z1","z2","z3","z4").latent("f").fixedVariance("f",1)
            .loading("a","z1","f",.8).loading("b","z2","f",.8)
            .loading("c","z3","f",.8).loading("d","z4","f",.8)
            .fixedVariance("z1",1).fixedVariance("z2",1).fixedVariance("z3",1).fixedVariance("z4",1).build();
    }
    public static double[][] data(String name) throws IOException {
        try(var in=SemJointTest.class.getResourceAsStream("/r-reference/sem-joint/"+name+".csv")) {
            if(in==null)throw new IOException("missing SEM fixture: "+name);
            return new String(in.readAllBytes(),java.nio.charset.StandardCharsets.UTF_8).lines().skip(1)
                .map(line->Arrays.stream(line.split(",")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);
        }
    }
    public static Properties reference(String name) throws IOException {
        Properties p=new Properties();try(InputStream in=SemJointTest.class.getResourceAsStream("/r-reference/sem-joint/"+name+".properties")){p.load(in);}return p;
    }
    static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
    static void compare(Properties r,List<SemParameterEstimate> parameters,double[] covariance,double tolerance) {
        List<String> labels=List.of(r.getProperty("labels").split(","));double[] estimates=numbers(r,"estimates"),expected=numbers(r,"covariance");int k=labels.size();
        assertEquals(k,parameters.size());
        for(int i=0;i<k;i++) {
            int a=labels.indexOf(parameters.get(i).label());assertTrue(a>=0,parameters.get(i).label());
            assertEquals(estimates[a],parameters.get(i).estimate(),tolerance,parameters.get(i).label());
            for(int j=0;j<k;j++) {
                int b=labels.indexOf(parameters.get(j).label());
                assertEquals(expected[a*k+b],covariance[i*k+j],tolerance,parameters.get(i).label()+" covariance "+parameters.get(j).label());
            }
        }
    }
    @Test void latentMeasurementAndStructuralPathsMatchLavaan() throws IOException {
        SemFitResult fit=Sem.fit(data("continuous"),latent(false),SemOptions.defaults(),BackendPolicy.CPU);
        assertTrue(fit.converged(),"score="+fit.scoreNorm());
        compare(reference("latent"),fit.parameters(),fit.parameterCovariance(),2e-6);
        assertEquals(numbers(reference("latent"),"logLikelihood")[0],fit.logLikelihood(),1e-6);
    }
    @Test void structuralMeansAndFullCovarianceMatchLavaan() throws IOException {
        SemFitResult fit=Sem.fit(data("continuous"),latent(true),SemOptions.defaults(),BackendPolicy.CPU);
        assertTrue(fit.converged(),"score="+fit.scoreNorm());
        compare(reference("means"),fit.parameters(),fit.parameterCovariance(),2e-6);
        assertEquals(numbers(reference("means"),"logLikelihood")[0],fit.logLikelihood(),1e-6);
        SemFitResult mediation=Sem.fit(data("mediation"),SemModel.builder("x","m","y")
            .regression("a","m","x",.4).regression("b","y","m",.6).regression("c","y","x",.1).meanStructure().build());
        compare(reference("mediation"),mediation.parameters(),mediation.parameterCovariance(),2e-6);
        double empiricalM=Arrays.stream(data("mediation")).mapToDouble(row->row[1]).average().orElseThrow();
        assertEquals(empiricalM,mediation.impliedMeans()[1],2e-6);
        assertTrue(Math.abs(empiricalM-mediation.parameter("m~1").estimate())>1e-3);
        assertEquals(numbers(reference("mediation"),"indirect.se")[0],SemInference.indirect(mediation,.95,"a","b").standardError(),2e-7);
    }
    @Test void constrainedPatternFimlMatchesLavaan() throws IOException {
        SemFimlResult result=SemFiml.fit(data("missing"),latent(true));SemFitResult fit=result.fit();
        assertTrue(fit.converged(),"score="+fit.scoreNorm());assertTrue(result.patternCount()>5);
        compare(reference("fiml"),fit.parameters(),fit.parameterCovariance(),3e-6);
        assertEquals(numbers(reference("fiml"),"logLikelihood")[0],fit.logLikelihood(),1e-6);
        assertEquals(numbers(reference("fiml"),"chisq")[0],fit.chiSquare(),2e-5);
        assertEquals(numbers(reference("fiml"),"cfi")[0],fit.cfi(),2e-6);
        assertEquals(numbers(reference("fiml"),"tli")[0],fit.tli(),2e-6);
        assertEquals(numbers(reference("fiml"),"rmsea")[0],fit.rmsea(),2e-6);
        assertTrue(Arrays.stream(SemInference.robust(fit).parameterCovariance()).allMatch(Double::isFinite));
    }
    @Test void completeFimlDoesNotInflateCovarianceAndEmptyRowsContributeNothing() {
        double[][] d={{-2,0},{-1,1},{0,-1},{1,1},{2,-1}};
        SemModel model=SemModel.builder("x","y").covariance("x","y",0).build();
        SemFitResult fit=SemFiml.fit(d,model).fit();assertTrue(fit.converged(),"score="+fit.scoreNorm());
        assertEquals(2,fit.impliedCovariance()[0],1e-7);assertEquals(.8,fit.impliedCovariance()[3],1e-7);
        double[][] withEmpty=Arrays.copyOf(d,6);withEmpty[5]=new double[]{Double.NaN,Double.NaN};
        SemFitResult same=SemFiml.fit(withEmpty,model).fit();assertEquals(fit.logLikelihood(),same.logLikelihood(),1e-9);assertEquals(5,same.observations());
    }
    @Test void constrainedMeanChangesJointCovarianceFitAndMatchesLavaan() throws IOException {
        // A fixed observed intercept changes the likelihood and the joint covariance optimum.
        SemModel model=latent(false).toBuilder().fixedIntercept("x1",1).build();
        SemFitResult fit=Sem.fit(data("continuous"),model);
        assertTrue(fit.converged(),"score="+fit.scoreNorm());
        compare(reference("fixedMean"),fit.parameters(),fit.parameterCovariance(),2e-6);
        assertEquals(numbers(reference("fixedMean"),"logLikelihood")[0],fit.logLikelihood(),1e-6);
        assertEquals(1,fit.impliedMeans()[0]);
    }
    @Test void invalidDataAndNegativeCorrelationPcaRegression() {
        SemModel model=SemModel.builder("x","y").build();
        assertThrows(IllegalArgumentException.class,()->SemFiml.fit(new double[][]{{1,Double.NaN},{2,Double.NaN},{3,Double.NaN}},model));
        assertThrows(IllegalArgumentException.class,()->Sem.fitCovariance(new double[]{1,2,0,1},100,model,SemOptions.defaults(),BackendPolicy.CPU));
        LatentMeasurement.Result pca=LatentMeasurement.fit(new double[][]{{-2,2},{-1,1},{0,0},{1,-1},{2,-2}},1);
        assertTrue(pca.loadings()[0]*pca.loadings()[1]<-1);
        assertThrows(IllegalArgumentException.class,()->SemOrdinal.fit(new int[][]{{0,0},{1,1},{0,1}},new int[]{2,2},model));
    }
    @Test void robustClusterAndEfficientModificationIndexMatchLavaan() throws IOException {
        SemFitResult fit=Sem.fit(data("continuous"),latent(true));Properties reference=reference("robust");
        SemInference.RobustResult robust=SemInference.robust(fit);
        compare(reference,fit.parameters(),robust.parameterCovariance(),3e-6);
        assertEquals(numbers(reference,"scalingFactor")[0],robust.scalingFactor(),2e-5);
        int[] clusters=java.util.stream.IntStream.range(0,fit.observations()).map(i->i/5).toArray();
        Properties clustered=reference("means");clustered.setProperty("covariance",clustered.getProperty("clusterCovariance"));
        compare(clustered,fit.parameters(),SemInference.robust(fit,clusters).parameterCovariance(),3e-6);
        assertThrows(IllegalArgumentException.class,()->SemInference.robust(fit,new int[fit.observations()]));
        SemInference.ModificationIndex mi=SemInference.modificationIndices(fit,SemInference.Modification.covariance("x1","x2")).get(0);
        assertEquals(numbers(reference("means"),"mi.x1.x2")[0],mi.chiSquare(),2e-5);
        assertEquals(numbers(reference("means"),"epc.x1.x2")[0],mi.expectedChange(),2e-6);
    }
    @Test void jointOrdinalPmlMatchesLavaan() throws IOException {
        int[][] data=Arrays.stream(data("ordinal")).map(row->Arrays.stream(row).mapToInt(v->(int)v).toArray()).toArray(int[][]::new);
        SemOrdinal.Result fit=SemOrdinal.fit(data,new int[]{3,3,3,3},ordinalModel());
        assertTrue(fit.converged(),"score="+fit.scoreNorm());
        compare(reference("ordinal"),fit.parameters(),fit.parameterCovariance(),3e-5);
        assertEquals(numbers(reference("ordinal"),"pairwiseLogLikelihood")[0],fit.pairwiseLogLikelihood(),2e-5);
        for(double[] t:fit.thresholds())assertTrue(t[0]<t[1]);
        assertEquals(.25+Math.asin(.7)/(2*Math.PI),SemOrdinal.bvn(0,0,.7),1e-13);
    }
    @Test void binaryNegativeAssociationHasExactProbitLikelihoodAndClusterCovariance() {
        int[][] data=new int[100][2];int row=0;
        for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int r=0;r<(a==b?15:35);r++)data[row++]=new int[]{a,b};
        SemModel model=SemModel.builder("x","y").regression("b","y","x",.2).fixedVariance("x",1).fixedVariance("y",1).build();
        SemOrdinal.Result fit=SemOrdinal.fit(data,new int[]{2,2},model);
        assertTrue(fit.converged(),"score="+fit.scoreNorm());
        assertEquals(-Math.tan(.2*Math.PI),fit.parameter("b").estimate(),2e-7);
        assertEquals(0,fit.thresholds()[0][0],2e-7);assertEquals(0,fit.thresholds()[1][0],2e-7);
        assertEquals(30*Math.log(.15)+70*Math.log(.35),fit.pairwiseLogLikelihood(),1e-9);
        int[][] replicated=new int[200][];int[] ids=new int[200];
        for(int r=0;r<100;r++){replicated[2*r]=replicated[2*r+1]=data[r];ids[2*r]=ids[2*r+1]=r;}
        SemOrdinal.Result clustered=SemOrdinal.fit(replicated,new int[]{2,2},model,SemOptions.defaults(),ids);
        assertArrayEquals(fit.parameterCovariance(),clustered.parameterCovariance(),2e-8);
    }
    @Test void analyticRamScoresMatchFiniteDifferences() throws IOException {
        SemModel m=latent(true);double[] x=RamFit.initial(m);List<RamFit.Pattern> patterns=RamFit.patterns(data("missing"));
        double[] score=RamFit.evaluate(m,x,patterns,600).gradient();
        for(int j=0;j<x.length;j++) {
            double[] a=x.clone(),b=x.clone();a[j]+=1e-5;b[j]-=1e-5;
            double derivative=(RamFit.evaluate(m,a,patterns,600).value()-RamFit.evaluate(m,b,patterns,600).value())/2e-5;
            assertEquals(derivative,score[j],2e-8,"score "+j);
        }
    }
    @Test void unidentifiedModelDoesNotFabricateInformation() {
        SemModel model=SemModel.builder("x","y","z","w").latent("f")
            .loading("x","f",.8).loading("y","f",.8).loading("z","f",.8).loading("w","f",.8).build();
        double[] cov={1,.4,.4,.4,.4,1,.4,.4,.4,.4,1,.4,.4,.4,.4,1};
        SemFitResult fit=Sem.fitCovariance(cov,200,model,SemOptions.defaults(),BackendPolicy.CPU);
        assertFalse(fit.informationAvailable());
    }
    /** Isolated package checks without rebuilding concurrently edited packages. */
    public static void main(String[] args) throws Exception {
        SemJointTest joint=new SemJointTest();int count=0;
        for(Class<?> c:List.of(SemTest.class,SemExtensionsTest.class,SemJointTest.class,SemReviewRegressionTest.class)) {
            Object test=c==SemJointTest.class?joint:c.getDeclaredConstructor().newInstance();
            for(var method:c.getDeclaredMethods())if(method.isAnnotationPresent(Test.class)) {
                System.out.println("RUN "+method.getName());
                try{method.invoke(test);count++;}catch(java.lang.reflect.InvocationTargetException e){throw new AssertionError(method.getName(),e.getCause());}
            }
        }
        System.out.println("SEM passed "+count+" tests");
    }
}
