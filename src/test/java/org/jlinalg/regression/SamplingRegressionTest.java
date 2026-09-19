/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.jlinalg.glm.*;

class SamplingRegressionTest {
    private static double[][] data() throws IOException {
        try(var reader=new BufferedReader(new InputStreamReader(SamplingRegressionTest.class.getResourceAsStream("/r-reference/regression-inference/survey.tsv")))){
            return reader.lines().skip(1).map(s->Arrays.stream(s.split("\t")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);
        }
    }
    private static double[] reference(String key) throws IOException {
        Properties p=new Properties();try(var in=SamplingRegressionTest.class.getResourceAsStream("/r-reference/regression-inference/reference.properties")){p.load(in);}
        return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();
    }
    private static double[] col(double[][] data,int column){return Arrays.stream(data).mapToDouble(r->r[column]).toArray();}
    private static double[][] design(double[][] data){return Arrays.stream(data).map(r->new double[]{1,r[1],r[2]}).toArray(double[][]::new);}
    @Test void rareEventCorrectionsMatchIndependentEquationsAndRemainSeparate() throws IOException {
        double[][] data=data(),x=design(data);double[] y=col(data,3);
        var ml=RareEventsLogistic.fit(y,x,false,null);
        var bias=RareEventsLogistic.fit(y,x,true,null);
        var prior=RareEventsLogistic.fit(y,x,true,.04);
        assertArrayEquals(reference("rare.ml.beta"),ml.coefficients(),2e-7);
        assertArrayEquals(reference("rare.ml.cov"),ml.covariance(),2e-7);
        assertArrayEquals(reference("rare.bias.beta"),bias.coefficients(),2e-7);
        assertArrayEquals(reference("rare.prior.beta"),prior.coefficients(),2e-7);
        var priorOnly=RareEventsLogistic.fit(y,x,false,.04);
        assertEquals(ml.coefficients()[1],priorOnly.coefficients()[1],1e-12);
        assertEquals(ml.coefficients()[0]+priorOnly.interceptPriorCorrection(),priorOnly.coefficients()[0],1e-12);
        assertThrows(IllegalArgumentException.class,()->RareEventsLogistic.fit(y,x,true,0.0));
        assertThrows(IllegalArgumentException.class,()->RareEventsLogistic.fit(new double[4],new double[][]{{1},{1},{1},{1}},true,null));
    }
    @Test void surveyMatchesPackageAcrossFamiliesAndIsWeightScaleInvariant() throws IOException {
        double[][] data=data(),x=design(data);double[] w=col(data,5);
        String[] strata=Arrays.stream(data).map(r->Double.toString(r[6])).toArray(String[]::new),psu=Arrays.stream(data).map(r->Double.toString(r[7])).toArray(String[]::new);
        var design=new SurveyRegression.Design(w,strata,psu);
        String[] names={"gaussian","logit","probit","poisson"};GlmFamily[] families={GlmFamilies.gaussian(),GlmFamilies.binomial(),GlmFamilies.probit(),GlmFamilies.poisson()};
        for(int i=0;i<families.length;i++){
            double[] y=col(data,i==0?0:i==3?4:3);var fit=SurveyRegression.fit(y,x,families[i],design);
            assertArrayEquals(reference("survey."+names[i]+".beta"),fit.coefficients(),2e-6,names[i]);
            assertArrayEquals(reference("survey."+names[i]+".cov"),fit.covariance(),2e-6,names[i]);
            assertEquals(55,fit.degreesOfFreedom());
            double[] scaled=w.clone();for(int j=0;j<scaled.length;j++)scaled[j]*=1e12;
            var again=SurveyRegression.fit(y,x,families[i],new SurveyRegression.Design(scaled,strata,psu));
            assertArrayEquals(fit.covariance(),again.covariance(),1e-10);
        }
        String[] singleton=psu.clone();Arrays.fill(singleton,"one");
        assertThrows(IllegalArgumentException.class,()->SurveyRegression.fit(col(data,0),x,GlmFamilies.gaussian(),new SurveyRegression.Design(w,strata,singleton)));
    }
}
