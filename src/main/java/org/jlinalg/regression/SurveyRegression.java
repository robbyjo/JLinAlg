/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.*;

/** With-replacement, stratified PSU Taylor-linearized GLM covariance.
 * Supports a one-stage design (or ultimate-PSU approximation), no FPC.
 * Single-PSU strata are rejected; they are not silently treated as certainty units. */
public final class SurveyRegression {
    private SurveyRegression() { }
    public record Design(double[] weights,String[] strata,String[] psu) {
        public Design {
            if(weights==null||strata==null||psu==null||weights.length==0||strata.length!=weights.length||psu.length!=weights.length)
                throw new IllegalArgumentException("matching survey weights, strata and PSU identifiers required");
            weights=weights.clone();strata=strata.clone();psu=psu.clone();
            for(int i=0;i<weights.length;i++)if(!(weights[i]>0)||!Double.isFinite(weights[i])||strata[i]==null||strata[i].isBlank()||psu[i]==null||psu[i].isBlank())
                throw new IllegalArgumentException("positive finite survey weights and nonblank design labels required");
        }
        @Override public double[] weights(){return weights.clone();}
        @Override public String[] strata(){return strata.clone();}
        @Override public String[] psu(){return psu.clone();}
    }
    public static Result fit(double[] y,double[][] x,GlmFamily family,Design design) {
        if(design==null||y==null||y.length!=design.weights.length||x==null||x.length!=y.length)
            throw new IllegalArgumentException("response/design rows must match the sampling design");
        if(family!=GlmFamilies.gaussian()&&family!=GlmFamilies.binomial()&&family!=GlmFamilies.probit()&&family!=GlmFamilies.poisson())
            throw new IllegalArgumentException("survey regression supports Gaussian, logit, probit and Poisson families");
        if(family==GlmFamilies.binomial()||family==GlmFamilies.probit())for(double v:y)if(v!=0&&v!=1)
            throw new IllegalArgumentException("survey binomial response must be individual binary outcomes");
        int n=y.length;double max=Arrays.stream(design.weights).max().orElseThrow(),mean=0;
        for(double w:design.weights)mean+=w/max/n;
        double[] weights=new double[n];for(int i=0;i<n;i++)weights[i]=(design.weights[i]/max)/mean;
        var fit=Glm.fit(y,x,family,weights,null,GlmOptions.defaults(),BackendPolicy.CPU);
        if(!fit.converged()||fit.rank()!=fit.parameters())throw new IllegalArgumentException("survey mean model failed convergence or rank check");
        int p=fit.parameters();double[] bread=fit.covariance(),eta=fit.linearPredictor(),mu=fit.fittedMeans();
        if(fit.estimatedDispersion()){
            if(!(fit.dispersion()>0))throw new IllegalArgumentException("positive residual dispersion required");
            for(int j=0;j<bread.length;j++)bread[j]/=fit.dispersion();
        }
        Map<String,Map<String,double[]>> strata=new LinkedHashMap<>();
        for(int i=0;i<n;i++){
            var units=strata.computeIfAbsent(design.strata[i],v->new LinkedHashMap<>());
            double[] total=units.computeIfAbsent(design.psu[i],v->new double[p]);
            double scalar=weights[i]*family.residualAtPredictor(y[i],eta[i],mu[i])*family.meanDerivative(eta[i])/family.varianceAtPredictor(eta[i],mu[i]);
            for(int j=0;j<p;j++)for(int k=0;k<p;k++)total[j]+=bread[j*p+k]*x[i][k]*scalar;
        }
        int psus=0;double[] covariance=new double[p*p];
        for(var units:strata.values()) {
            int m=units.size();if(m<2)throw new IllegalArgumentException("each stratum must contain at least two PSUs");psus+=m;
            double[] average=new double[p];for(double[] unit:units.values())for(int j=0;j<p;j++)average[j]+=unit[j]/m;
            for(double[] unit:units.values())for(int j=0;j<p;j++)for(int k=0;k<p;k++)
                covariance[j*p+k]+=(double)m/(m-1)*(unit[j]-average[j])*(unit[k]-average[k]);
        }
        for(double value:covariance)if(!Double.isFinite(value))throw new IllegalArgumentException("nonfinite survey covariance");
        return new Result(fit.coefficients(),covariance,psus-strata.size(),psus,strata.size(),family.name());
    }
    /** Tests use the conservative design df = PSUs - strata, explicitly rather than residual model df. */
    public record Result(double[] coefficients,double[] covariance,int degreesOfFreedom,int psus,int strata,String family) {
        public Result{coefficients=coefficients.clone();covariance=covariance.clone();}
        @Override public double[] coefficients(){return coefficients.clone();}
        @Override public double[] covariance(){return covariance.clone();}
    }
}
