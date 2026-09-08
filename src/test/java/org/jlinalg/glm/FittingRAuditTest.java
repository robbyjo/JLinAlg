package org.jlinalg.glm;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.*;
import org.jlinalg.gee.*;
import org.jlinalg.penalized.*;
import org.junit.jupiter.api.Test;

/** Frozen R 4.6.1 / glmnet 5.0 values; generator also reconstructs GEE analytically. */
class FittingRAuditTest {
    private final double[][] x=new double[60][2], design=new double[60][3];
    private final double[] y=new double[60], weights=new double[60], offset=new double[60], counts=new double[60], gamma=new double[60];
    private final int[] id=new int[60], waves=new int[60];
    FittingRAuditTest() {
        for(int i=0;i<60;i++) {
            x[i][0]=Math.sin(i*.37);x[i][1]=Math.cos(i*.19);
            design[i][0]=1;design[i][1]=x[i][0];design[i][2]=x[i][1];
            weights[i]=1+i%5;offset[i]=.2*Math.sin(i*.11);
            y[i]=1+.7*x[i][0]-.4*x[i][1]+.3*Math.sin(i*1.7)+offset[i];
            counts[i]=(i*7+3)%11;gamma[i]=Math.exp(.3+.5*x[i][0]+.7*Math.sin(i));
            id[i]=i/3;waves[i]=i%3;
        }
    }
    @Test void weightedOffsetOlsAndGaussian() throws IOException {
        var o=Ols.fit(y,design,weights,offset,OlsOptions.defaults(),BackendPolicy.CPU);
        check("ols.beta",o.coefficients(),1e-12);check("ols.covariance",o.covariance(),1e-12);
        check("ols.p",o.pValues(),1e-12);check("ols.ll",new double[]{o.logLikelihood()},1e-10);
        var g=Glm.fit(y,design,GlmFamilies.gaussian(),weights,offset,GlmOptions.defaults(),BackendPolicy.CPU);
        assertTrue(g.converged());check("gaussian.beta",g.coefficients(),1e-12);
        check("gaussian.p",g.pValues(),1e-12);check("gaussian.ll",new double[]{g.logLikelihood()},1e-10);
    }
    @Test void weightedOffsetPoisson() throws IOException {
        var g=Glm.fit(counts,design,GlmFamilies.poisson(),weights,offset,
            GlmOptions.builder().relativeTolerance(1e-12).maximumIterations(200).build(),BackendPolicy.CPU);
        assertTrue(g.converged(),g.convergenceMessage());check("poisson.beta",g.coefficients(),1e-8);
        check("poisson.covariance",g.covariance(),1e-8);check("poisson.p",g.pValues(),1e-7);
        check("poisson.ll",new double[]{g.logLikelihood()},1e-9);
    }
    @Test void gammaPearsonInferenceAndProfileLikelihood() throws IOException {
        var g=Glm.fit(gamma,design,GlmFamilies.gamma(),null,null,
            GlmOptions.builder().relativeTolerance(1e-12).maximumIterations(200).build(),BackendPolicy.CPU);
        assertTrue(g.converged(),g.convergenceMessage());check("gamma.beta",g.coefficients(),1e-7);
        check("gamma.dispersion",new double[]{g.dispersion()},1e-8);
        check("gamma.profile.ll",new double[]{g.logLikelihood()},1e-9);
    }
    @Test void fixedCorrelationWeightedGeeAgainstIndependentRMatrixEquations() throws IOException {
        var g=Gee.fit(y,design,id,waves,GlmFamilies.gaussian(),weights,offset,
            GeeOptions.builder().correlation(GeeCorrelation.FIXED)
                .fixedAssociation(new double[][]{{1,.2,.2},{.2,1,.2},{.2,.2,1}})
                .fixedDispersion(1).build(),BackendPolicy.CPU);
        assertTrue(g.converged(),g.convergenceMessage());check("gee.beta",g.coefficients(),1e-12);
        check("gee.covariance",g.covariance(),1e-12);
    }
    @Test void inverseGaussianAndRareBinomialLikelihoods() throws IOException {
        var ig=Glm.fit(gamma,design,GlmFamilies.inverseGaussian(),null,null,
            GlmOptions.builder().relativeTolerance(1e-12).maximumIterations(200).build(),BackendPolicy.CPU);
        assertTrue(ig.converged(),ig.convergenceMessage());check("ig.beta",ig.coefficients(),1e-7);
        check("ig.covariance",ig.covariance(),1e-7);check("ig.profile.ll",new double[]{ig.logLikelihood()},1e-9);
        var bin=Glm.fit(new double[]{1e-14,2e-14,3e-14},new double[][]{{1},{1},{1}},GlmFamilies.binomial(),
            new double[]{1e14,1e14,1e14},null,GlmOptions.defaults(),BackendPolicy.CPU);
        assertTrue(bin.converged());check("binomial.rare.ll",new double[]{bin.logLikelihood()},1e-10);
        var tail=Glm.fit(new double[]{1,0,0,0},new double[][]{{1},{1},{1},{1}},GlmFamilies.binomial(),
            new double[]{1e22,1,1,1},new double[]{50,0,0,0},GlmOptions.defaults(),BackendPolicy.CPU);
        assertTrue(tail.converged());check("binomial.tail",new double[]{tail.coefficients()[0],tail.logLikelihood(),tail.standardErrors()[0]},1e-10);
    }
    @Test void weightedLassoAndWeightedFoldRiskAgainstGlmnet() throws IOException {
        var options=ElasticNetOptions.builder().alpha(1).observationWeights(weights).relativeTolerance(1e-14).build();
        var path=PenalizedRegression.path(y,x,new double[]{.2,.08,.02},options);
        double[] beta=new double[9];
        for(int j=0;j<3;j++){var f=path.fit(j);assertTrue(f.converged());beta[3*j]=f.intercept();System.arraycopy(f.coefficients(),0,beta,3*j+1,2);}
        check("lasso.beta",beta,1e-7);
        var cv=PenalizedRegressionCrossValidation.fit(y,x,new double[]{.2,.08,.02},5,17,options);
        check("cv.mse",cv.meanSquaredErrors(),1e-8);check("cv.se",cv.standardErrors(),1e-8);
    }
    @Test void ridgeAndElasticNetWithUnpenalizedPredictorMatchExactRActiveSetEnumeration() throws IOException {
        for(double alpha:new double[]{0,.35}) {
            var options=ElasticNetOptions.builder().alpha(alpha).observationWeights(weights)
                .penaltyFactors(0,2).relativeTolerance(1e-14).build();
            var path=PenalizedRegression.path(y,x,new double[]{.2,.08,.02},options);
            double[] beta=new double[9];
            for(int j=0;j<3;j++){var f=path.fit(j);assertTrue(f.converged());beta[3*j]=f.intercept();System.arraycopy(f.coefficients(),0,beta,3*j+1,2);}
            check(alpha==0?"enet.0":"enet.0.35",beta,1e-7);
        }
    }
    private static void check(String key,double[] values,double tolerance) throws IOException {
        Properties p=new Properties();
        try(InputStream input=FittingRAuditTest.class.getResourceAsStream("/fitting-audit-v030/reference.properties")) {
            assertNotNull(input);p.load(input);
        }
        String encoded=p.getProperty(key);assertNotNull(encoded,key);
        double[] expected=Arrays.stream(encoded.split(",")).mapToDouble(Double::parseDouble).toArray();
        assertArrayEquals(expected,values,tolerance,key);
        double maximumError=0;
        for(int i=0;i<values.length;i++)maximumError=Math.max(maximumError,Math.abs(expected[i]-values[i]));
        System.out.println("fitting R accuracy "+key+" max_abs_error="+maximumError);
    }
}
