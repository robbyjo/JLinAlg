/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import org.jlinalg.compute.*;

/** Observed information from analytic-score differences; no ridge or pseudoinverse. */
final class LikelihoodInference {
    private LikelihoodInference() { }
    static double[] covariance(RegressionOptimizer.Objective objective,double[] point,int n) {
        int p=point.length;double[] h=hessian(objective,point,1e-4),fine=hessian(objective,point,5e-5);
        for(int j=0;j<h.length;j++)if(!Double.isFinite(h[j])||Math.abs(h[j]-fine[j])>2e-5*(1+Math.abs(fine[j])))
            throw new IllegalArgumentException("observed-information step refinement failed");
        double[] identity=new double[p*p];for(int j=0;j<p;j++)identity[j*p+j]=1;
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            double[] cov=context.backend().dpotrf(fine,p).solve(identity,p);
            double[] score=objective.evaluate(point).gradient();
            for(int j=0;j<p;j++) {
                double step=0;for(int l=0;l<p;l++)step+=cov[j*p+l]*score[l];
                if(!Double.isFinite(step)||Math.abs(step)>1e-3)
                    throw new IllegalArgumentException("unresolved Newton step; possible separation or boundary likelihood");
            }
            for(int j=0;j<cov.length;j++)cov[j]/=n;
            // Verify the inverse in optimizer coordinates before transforming it.
            for(int j=0;j<p;j++)for(int k=0;k<p;k++) {
                double product=0;for(int l=0;l<p;l++)product+=fine[j*p+l]*cov[l*p+k]*n;
                if(!Double.isFinite(product)||Math.abs(product-(j==k?1:0))>1e-7)
                    throw new IllegalArgumentException("singular or unresolved observed information");
            }
            return cov;
        } catch(IllegalStateException e){throw new IllegalArgumentException("positive definite observed information required",e);}
    }
    private static double[] hessian(RegressionOptimizer.Objective objective,double[] point,double delta) {
        int p=point.length;double[] h=new double[p*p];
        for(int j=0;j<p;j++) {
            double step=delta*(1+Math.abs(point[j]));double[] a=point.clone(),b=point.clone();a[j]+=step;b[j]-=step;
            double[] hi=objective.evaluate(a).gradient(),lo=objective.evaluate(b).gradient();
            for(int k=0;k<p;k++)h[k*p+j]=(hi[k]-lo[k])/(2*step);
        }
        for(int j=0;j<p;j++)for(int k=0;k<j;k++)h[j*p+k]=h[k*p+j]=(h[j*p+k]+h[k*p+j])/2;
        return h;
    }
}
