/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 * Adapted from R stats/src/ppr.f: SUPSMU, Copyright (C) 1984 Jerome H.
 * Friedman; modifications by B. D. Ripley and R Core, GPL-2.0-or-later.
 * https://svn.r-project.org/R/trunk/src/library/stats/src/ppr.f
 */
package org.jlinalg.regression;

import java.util.Arrays;
import java.util.Comparator;

/** Friedman's running-lines supersmoother with smoothed cross-validation spans. */
public final class SuperSmoother {
    private SuperSmoother() { }
    private static final double[] SPANS={.05,.2,.5};
    public static Result fit(double[] x,double[] y) {return fit(x,y,SPANS);}
    public static Result fit(double[] x,double[] y,double[] spans) {return fit(x,y,null,spans,false,0);}
    /** R-style controls: span zero selects CV; positive span fixes the window. */
    public static Result fit(double[] x,double[] y,double[] weights,double span,boolean periodic,double bass) {
        if(!Double.isFinite(span)||span<0||span>1)throw new IllegalArgumentException("span must be in [0,1]");
        return fit(x,y,weights,span==0?SPANS:new double[]{span},periodic,bass);
    }
    private static Result fit(double[] predictor,double[] response,double[] weights,double[] spans,boolean periodic,double bass) {
        if(predictor==null||response==null||predictor.length!=response.length||predictor.length<4
                ||spans==null||spans.length==0||!Double.isFinite(bass)||bass<0||bass>10
                ||(weights!=null&&weights.length!=predictor.length))throw new IllegalArgumentException("invalid supersmoother inputs");
        int n=predictor.length;
        if(periodic&&n<5)throw new IllegalArgumentException("periodic supersmoothing requires at least five observations");
        double[] orderedSpans=spans.clone();Arrays.sort(orderedSpans);
        for(int i=0;i<orderedSpans.length;i++)if(!(orderedSpans[i]>0)||orderedSpans[i]>1
                ||(i>0&&orderedSpans[i]==orderedSpans[i-1]))throw new IllegalArgumentException("spans must be distinct and in (0,1]");
        Integer[] order=new Integer[n];
        for(int i=0;i<n;i++) {
            order[i]=i;
            if(!Double.isFinite(predictor[i])||!Double.isFinite(response[i])
                    ||(weights!=null&&(!(weights[i]>0)||!Double.isFinite(weights[i]))))throw new IllegalArgumentException("finite observations and positive weights required");
            if(periodic&&(predictor[i]<0||predictor[i]>1))throw new IllegalArgumentException("periodic predictors must be in [0,1]");
        }
        Arrays.sort(order,Comparator.comparingDouble((Integer i)->predictor[i]).thenComparingDouble(i->response[i]));
        double[] x=new double[n],y=new double[n],w=new double[n];
        for(int i=0;i<n;i++){x[i]=predictor[order[i]];y[i]=response[order[i]];w[i]=weights==null?1:weights[order[i]];}
        double[] fitted=new double[n],selected=new double[n];
        if(x[n-1]==x[0]) {
            double sum=0,weight=0;for(int i=0;i<n;i++){sum+=w[i]*y[i];weight+=w[i];}
            Arrays.fill(fitted,sum/weight);Arrays.fill(selected,1);
        } else {
            int lo=Math.max(0,n/4-1),hi=Math.max(lo+1,3*(n/4)-1);
            while(x[hi]<=x[lo]){if(hi<n-1)hi++;if(lo>0)lo--;}
            double small=Math.pow(.001*(x[hi]-x[lo]),2);
            if(orderedSpans.length==1){fitted=smooth(x,y,w,orderedSpans[0],periodic,small).fitted;Arrays.fill(selected,orderedSpans[0]);}
            else {
                int k=orderedSpans.length;double middle=orderedSpans[k/2];
                double[][] candidates=new double[k][],errors=new double[k][];
                for(int j=0;j<k;j++) {
                    Running fit=smooth(x,y,w,orderedSpans[j],periodic,small);
                    candidates[j]=fit.fitted;
                    errors[j]=smooth(x,fit.cv,w,middle,periodic,small).fitted;
                }
                for(int i=0;i<n;i++) {
                    int best=0;for(int j=1;j<k;j++)if(errors[j][i]<errors[best][i])best=j;
                    selected[i]=orderedSpans[best];double min=errors[best][i];
                    if(bass>0&&min>0&&min<errors[k-1][i])selected[i]+=(orderedSpans[k-1]-selected[i])*Math.pow(Math.max(1e-7,min/errors[k-1][i]),10-bass);
                }
                selected=smooth(x,selected,w,middle,periodic,small).fitted;
                double[] blend=new double[n];
                for(int i=0;i<n;i++) {
                    selected[i]=Math.max(orderedSpans[0],Math.min(orderedSpans[k-1],selected[i]));
                    int right=1;while(right<k-1&&orderedSpans[right]<selected[i])right++;
                    double fraction=(selected[i]-orderedSpans[right-1])/(orderedSpans[right]-orderedSpans[right-1]);
                    blend[i]=(1-fraction)*candidates[right-1][i]+fraction*candidates[right][i];
                }
                fitted=smooth(x,blend,w,orderedSpans[0],periodic,small).fitted;
            }
        }
        double[] originalFits=new double[n],originalSpans=new double[n];
        for(int i=0;i<n;i++){originalFits[order[i]]=fitted[i];originalSpans[order[i]]=selected[i];}
        return new Result(predictor,response,originalFits,originalSpans);
    }
    private record Running(double[] fitted,double[] cv) { }
    private static Running smooth(double[] x,double[] y,double[] w,double span,boolean periodic,double small) {
        int n=x.length,half=Math.max(2,(int)(.5*span*n+.5)),count=Math.min(n,2*half+1);
        Moments m=new Moments();
        for(int i=0;i<count;i++){int index=periodic?i-half-1:i;int wrapped=Math.floorMod(index,n);m.add(x[wrapped]+Math.floorDiv(index,n),y[wrapped],w[wrapped]);}
        double[] fitted=new double[n],cv=new double[n];
        for(int j=0;j<n;j++) {
            int out=j-half-1,in=j+half;
            if(periodic||(out>=0&&in<n)) {
                int a=Math.floorMod(out,n),b=Math.floorMod(in,n);
                m.remove(x[a]+Math.floorDiv(out,n),y[a],w[a]);m.add(x[b]+Math.floorDiv(in,n),y[b],w[b]);
            }
            double dx=x[j]-m.x,slope=m.variance>small?m.covariance/m.variance:0;
            fitted[j]=m.y+slope*dx;
            double leverage=1/m.weight+(m.variance>small?dx*dx/m.variance:0),denominator=1-w[j]*leverage;
            cv[j]=denominator>0?Math.abs(y[j]-fitted[j])/denominator:(j>0?cv[j-1]:0);
        }
        for(int start=0;start<n;) {
            int end=start+1;double sum=w[start]*fitted[start],weight=w[start];
            while(end<n&&x[end]==x[start]){sum+=w[end]*fitted[end];weight+=w[end++];}
            if(end>start+1)Arrays.fill(fitted,start,end,sum/weight);
            start=end;
        }
        return new Running(fitted,cv);
    }
    private static final class Moments {
        double weight,x,y,variance,covariance;
        void add(double xi,double yi,double wi) {
            double old=weight;weight+=wi;x=(old*x+wi*xi)/weight;y=(old*y+wi*yi)/weight;
            double factor=old>0?weight*wi*(xi-x)/old:0;variance+=factor*(xi-x);covariance+=factor*(yi-y);
        }
        void remove(double xi,double yi,double wi) {
            double old=weight;weight-=wi;double factor=weight>0?old*wi*(xi-x)/weight:0;
            variance-=factor*(xi-x);covariance-=factor*(yi-y);
            if(weight>0){x=(old*x-wi*xi)/weight;y=(old*y-wi*yi)/weight;}
        }
    }
    /** Fitted values retain input order; predictions linearly interpolate the fitted
     * curve, with constant extrapolation outside the training range. */
    public record Result(double[] predictor,double[] response,double[] fittedValues,double[] selectedSpans) {
        public Result {predictor=predictor.clone();response=response.clone();fittedValues=fittedValues.clone();selectedSpans=selectedSpans.clone();}
        public double[] predictor(){return predictor.clone();}
        public double[] response(){return response.clone();}
        public double[] fittedValues(){return fittedValues.clone();}
        public double[] selectedSpans(){return selectedSpans.clone();}
        public double[] predict(double[] query) {
            if(query==null)throw new IllegalArgumentException("queries are required");
            int n=predictor.length;Integer[] order=new Integer[n];for(int i=0;i<n;i++)order[i]=i;
            Arrays.sort(order,Comparator.comparingDouble(i->predictor[i]));
            double[] result=new double[query.length];
            for(int i=0;i<query.length;i++) {
                double q=query[i];if(!Double.isFinite(q))throw new IllegalArgumentException("queries must be finite");
                if(q<=predictor[order[0]]){result[i]=fittedValues[order[0]];continue;}
                if(q>=predictor[order[n-1]]){result[i]=fittedValues[order[n-1]];continue;}
                int left=0,right=n-1;
                while(right-left>1){int mid=(left+right)/2;if(predictor[order[mid]]<=q)left=mid;else right=mid;}
                int a=order[left],b=order[right];double fraction=(q-predictor[a])/(predictor[b]-predictor[a]);
                result[i]=(1-fraction)*fittedValues[a]+fraction*fittedValues[b];
            }
            return result;
        }
    }
}
