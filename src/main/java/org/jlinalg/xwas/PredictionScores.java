/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.xwas;

import java.util.Arrays;
import org.jlinalg.ols.Ols;
import org.jlinalg.penalized.*;

/** Training-only Gaussian penalized scores and independent-cohort evaluation. */
public final class PredictionScores {
    private PredictionScores() { }
    public record Model(double intercept,double[] weights,double lambda,double alpha,double trainingMean) {
        public Model {SummaryMath.finite(weights);weights=weights.clone();
            if(!Double.isFinite(intercept)||!Double.isFinite(trainingMean))throw new IllegalArgumentException("nonfinite model");}
        @Override public double[] weights(){return weights.clone();}
        public double[] predict(double[][] features) {
            double[] result=new double[features.length];
            for(int i=0;i<features.length;i++){SummaryMath.finite(features[i]);result[i]=intercept+SummaryMath.dot(features[i],weights);}
            SummaryMath.finite(result);return result;
        }
    }
    public record Evaluation(int observations,double rmse,double predictiveR2,double calibrationIntercept,
            double calibrationSlope,double auc) { }
    public static Model train(double[] y,double[][] x,double[] lambdas,double alpha,int folds,long seed) {
        var options=ElasticNetOptions.builder().alpha(alpha).relativeTolerance(1e-10).build();
        PenalizedRegressionResult fit;
        if(lambdas.length==1)fit=PenalizedRegression.fit(y,x,lambdas[0],options);
        else fit=PenalizedRegressionCrossValidation.fit(y,x,lambdas,folds,seed,options).minimumErrorFit();
        if(!fit.converged())throw new IllegalArgumentException("score training failed: "+fit.convergenceMessage());
        return new Model(fit.intercept(),fit.coefficients(),fit.lambda(),alpha,Arrays.stream(y).average().orElseThrow());
    }
    /** Prediction R2 uses a frozen baseline prediction, usually the training outcome mean. */
    public static Evaluation evaluate(double[] y,double[] prediction,double[] baseline,boolean binary) {
        SummaryMath.finite(y);SummaryMath.finite(prediction);SummaryMath.finite(baseline);
        int n=y.length;
        if(n<4||prediction.length!=n||baseline.length!=n)throw new IllegalArgumentException("evaluation requires at least four aligned observations");
        double ss=0,base=0;
        for(int i=0;i<n;i++){ss+=Math.pow(y[i]-prediction[i],2);base+=Math.pow(y[i]-baseline[i],2);}
        if(!(base>0))throw new IllegalArgumentException("baseline error is zero; predictive R2 is undefined");
        double intercept=Double.NaN,slope=Double.NaN;
        double min=Arrays.stream(prediction).min().orElseThrow(),max=Arrays.stream(prediction).max().orElseThrow();
        if(max>min) {
            double[][] design=new double[n][2];for(int i=0;i<n;i++){design[i][0]=1;design[i][1]=prediction[i];}
            double[] coefficients=Ols.fit(y,design).coefficients();intercept=coefficients[0];slope=coefficients[1];
        }
        double auc=Double.NaN;
        if(binary) {
            Integer[] order=new Integer[n];int cases=0;
            for(int i=0;i<n;i++){if(y[i]!=0&&y[i]!=1)throw new IllegalArgumentException("binary outcomes must be 0/1");cases+=(int)y[i];order[i]=i;}
            if(cases==0||cases==n)throw new IllegalArgumentException("AUC requires both cases and controls");
            Arrays.sort(order,(a,b)->Double.compare(prediction[a],prediction[b]));
            double ranks=0;
            for(int start=0;start<n;) {
                int end=start+1;while(end<n&&prediction[order[end]]==prediction[order[start]])end++;
                double rank=(start+1+end)/2.0;
                for(int i=start;i<end;i++)if(y[order[i]]==1)ranks+=rank;
                start=end;
            }
            auc=(ranks-cases*(cases+1.0)/2)/(cases*(double)(n-cases));
        }
        return new Evaluation(n,Math.sqrt(ss/n),1-ss/base,intercept,slope,auc);
    }
}
