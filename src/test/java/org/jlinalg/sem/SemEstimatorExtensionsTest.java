/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.Test;
final class SemEstimatorExtensionsTest {
    @Test void fimlRobustScalingMatchesIndependentObservedInformationCalculation() throws Exception {
        var fit=SemFiml.fit(SemJointTest.data("missing"),SemJointTest.latent(true)).fit();
        var robust=SemInference.robust(fit);double[] reference=vector("fiml-robust-reference.tsv");
        assertEquals(reference[0],fit.chiSquare(),1e-5);
        assertEquals(reference[1],robust.scalingFactor(),1e-5);
        assertEquals(reference[2],robust.chiSquare(),2e-5);
        assertTrue(Double.isFinite(robust.pValue()));
    }
    static double[][] matrix(String name) throws Exception {
        return Files.readAllLines(Path.of("src/test/resources/estimator-extensions/"+name)).stream()
            .filter(s->!s.isBlank()).map(s->Arrays.stream(s.split("\\s+")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);
    }
    static double[] vector(String name) throws Exception {return Arrays.stream(matrix(name)).flatMapToDouble(Arrays::stream).toArray();}
    @Test void dwlsWlsmvMatchesLavaanIncludingSandwichAndAdjustment() throws Exception {
        var b=SemModel.builder("y1","y2","y3","y4").latent("F").fixedVariance("F",1);
        for(int j=1;j<=4;j++)b.loading("l"+j,"y"+j,"F",.6).fixedVariance("y"+j,1);
        var moments=new SemDwls.Moments(600,new int[]{3,3,3,3},vector("dwls-moments.tsv"),vector("dwls-gamma.tsv"));
        var fit=SemDwls.fit(moments,b.build());double[] expected=vector("dwls-reference.tsv");
        assertTrue(fit.converged());assertEquals(2,fit.degreesOfFreedom());
        for(int j=0;j<4;j++)assertEquals(expected[j],fit.estimates()[j],2e-5);
        assertEquals(expected[4],fit.wlsmvStatistic(),2e-5);assertEquals(expected[5],fit.scalingFactor(),2e-5);assertEquals(expected[6],fit.shift(),2e-5);
        double[] covariance=vector("dwls-covariance.tsv");int k=fit.estimates().length;
        for(int i=0;i<4;i++)for(int j=0;j<4;j++)assertEquals(covariance[i*4+j],fit.parameterCovariance()[i*k+j],2e-6);
    }
    private static SemModel mixedModel(boolean slope) {
        var b=SemModel.builder("x","y").meanStructure().intercept("mu","x",.3).variance("vx","x",1)
            .fixedVariance("y",1).fixedIntercept("y",0);
        if(slope)b.regression("b","y","x",.4);return b.build();
    }
    @Test void mixedMarLikelihoodMatchesIndependentRConditionalDensity() throws Exception {
        var fit=SemMixed.fit(matrix("mixed-data.tsv"),new int[]{0,2},mixedModel(true));
        double[] reference=vector("mixed-reference.tsv");assertTrue(fit.converged());assertTrue(fit.informationAvailable());
        double[] estimates=fit.estimates();
        assertEquals(reference[0],estimates[fit.parameterLabels().indexOf("mu")],2e-5);
        assertEquals(Math.exp(2*reference[1]),estimates[fit.parameterLabels().indexOf("vx")],2e-5);
        assertEquals(reference[2],estimates[fit.parameterLabels().indexOf("b")],2e-5);
        assertEquals(reference[3],estimates[fit.parameterLabels().indexOf("y|t1")],2e-5);
        assertEquals(reference[4],fit.logLikelihood(),1e-7);
        int[] order={fit.parameterLabels().indexOf("mu"),fit.parameterLabels().indexOf("vx"),fit.parameterLabels().indexOf("b"),fit.parameterLabels().indexOf("y|t1")};
        double[] factors={1,2*Math.exp(2*reference[1]),1,1},covariance=vector("mixed-covariance.tsv");
        for(int i=0;i<4;i++)for(int j=0;j<4;j++)assertEquals(covariance[i*4+j]*factors[i]*factors[j],fit.parameterCovariance()[order[i]*4+order[j]],2e-5);
        var nullFit=SemMixed.fit(matrix("mixed-data.tsv"),new int[]{0,2},mixedModel(false));
        var mi=SemMixed.modificationIndices(nullFit,SemInference.Modification.regression("y","x")).get(0);
        assertTrue(mi.chiSquare()>10);assertTrue(mi.expectedChange()>0);
    }
    @Test void independentTrivariateRectangleFactorizesAndCorrelatedOrthantIsAnalytic() {
        double[] lo={Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY},hi={0,0,0};
        assertEquals(.125,SemRectangle.probability(lo,hi,RamFit.identity(3)),1e-12);
        double[] covariance={1,.2,.3,.2,1,.4,.3,.4,1};
        assertEquals(.125+(Math.asin(.2)+Math.asin(.3)+Math.asin(.4))/(4*Math.PI),SemRectangle.probability(lo,hi,covariance),1e-9);
    }
    @Test void fourDimensionalEquicorrelatedOrthantHasOrderStatisticProbability() {
        double[] lo=new double[4],hi=new double[4],cov=new double[16];Arrays.fill(lo,Double.NEGATIVE_INFINITY);Arrays.fill(cov,.5);
        for(int i=0;i<4;i++)cov[i*4+i]=1;
        // P(Z0 exceeds four independent normals) = 1/5.
        assertEquals(.2,SemRectangle.probability(lo,hi,cov),1e-9);
    }
    @Test void partialThresholdInvarianceSharesActualThresholdNotIncrement() {
        java.util.ArrayList<double[]> a=new java.util.ArrayList<>(),b=new java.util.ArrayList<>();
        int[] ca={10,20,30},cb={20,30,10};
        for(int category=0;category<3;category++) {
            for(int i=0;i<ca[category];i++)a.add(new double[]{category,i%2});
            for(int i=0;i<cb[category];i++)b.add(new double[]{category,i%2});
        }
        var model=SemModel.builder("x","y").fixedVariance("x",1).fixedVariance("y",1).build();
        var fit=SemMultigroup.fit(List.of(a.toArray(double[][]::new),b.toArray(double[][]::new)),new int[]{3,2},List.of(model,model),Set.of("x|t2"),SemOptions.defaults());
        assertTrue(fit.informationAvailable());
        assertEquals(jdistlib.Normal.quantile(2.0/3,0,1,true,false),fit.estimates()[fit.parameterLabels().indexOf("x|t2")],1e-6);
        assertEquals(jdistlib.Normal.quantile(2.0/9,0,1,true,false),fit.estimates()[fit.parameterLabels().indexOf("group0:x|t1")],1e-6);
        assertEquals(jdistlib.Normal.quantile(4.0/15,0,1,true,false),fit.estimates()[fit.parameterLabels().indexOf("group1:x|t1")],1e-6);
    }
    @Test void multigroupMeanInvarianceIsPooledGaussianLikelihood() {
        double[][] a={{1,0},{2,1},{3,0},{4,-1},{5,0},{6,1}},b={{2,0},{3,1},{4,0},{5,-1},{6,0},{7,1}};
        SemModel model=SemModel.builder("x","z").meanStructure().intercept("mean","x",3).variance("variance","x",3).fixedVariance("z",1).fixedIntercept("z",0).build();
        var groups=List.of(a,b);var models=List.of(model,model);
        var separate=SemMultigroup.fit(groups,new int[]{0,0},models,Set.of(),SemOptions.defaults());
        var shared=SemMultigroup.fit(groups,new int[]{0,0},models,Set.of("mean","variance"),SemOptions.defaults());
        assertTrue(separate.informationAvailable());assertTrue(shared.informationAvailable());
        assertEquals(4,shared.estimates()[shared.parameterLabels().indexOf("mean")],1e-6);
        assertEquals(19.0/6,shared.estimates()[shared.parameterLabels().indexOf("variance")],1e-6);
        assertEquals(12*Math.log((19.0/6)/(35.0/12)),SemMultigroup.compare(separate,shared).chiSquare(),1e-8);
    }
    @Test void saturatedMixedMomentPreparationRetainsFullGamma() throws Exception {
        var moments=SemDwls.prepare(matrix("mixed-data.tsv"),new int[]{0,2},"x","y");
        assertEquals(4,moments.statistics().length);assertEquals(16,moments.gamma().length);
        assertTrue(Arrays.stream(moments.gamma()).allMatch(Double::isFinite));
        var fit=SemDwls.fit(moments,mixedModel(true));assertTrue(fit.converged());assertEquals(0,fit.degreesOfFreedom());assertTrue(fit.dwlsStatistic()<1e-8);
    }
}
