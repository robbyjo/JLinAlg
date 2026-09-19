/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.network;

import java.util.*;
import org.jlinalg.penalized.ElasticNetOptions;
import org.jlinalg.penalized.PenalizedRegression;

/** Bounded dense-matrix network methods. Edges describe association, not causality. */
public final class NetworkAnalysis {
    private NetworkAnalysis() { }
    /** Standardized Gaussian neighborhood LASSO coefficients, [target][predictor]. */
    public static double[][] neighborhoodLasso(double[][] samples,double lambda,int maxIterations,double tolerance) {
        if(!Double.isFinite(lambda)||lambda<=0 ||maxIterations<1||!Double.isFinite(tolerance)||tolerance<=0)
            throw new IllegalArgumentException("Require positive lambda, iterations and tolerance");
        double[][] x=standardize(samples);int n=x.length,p=x[0].length;
        double[][] coefficients=new double[p][p];
        for(int target=0;target<p;target++) {
            double[] y=new double[n];double[][] predictors=new double[n][p-1];
            for(int i=0;i<n;i++){y[i]=x[i][target];int col=0;for(int j=0;j<p;j++)if(j!=target)predictors[i][col++]=x[i][j];}
            // The underlying solver compares squared coordinate contributions.
            var options=ElasticNetOptions.builder().fitIntercept(false).standardize(false).maximumIterations(maxIterations).relativeTolerance(tolerance*tolerance).build();
            var fit=PenalizedRegression.lasso(y,predictors,lambda,options);
            if(!fit.converged())throw new IllegalArgumentException("Neighborhood LASSO did not converge for feature "+target);
            double[] beta=fit.coefficients();int col=0;for(int j=0;j<p;j++)if(j!=target)coefficients[target][j]=beta[col++];
        }
        return coefficients;
    }
    /** Center and scale each feature to population variance one. Never omit missing rows. */
    public static double[][] standardize(double[][] samples) {
        if(samples==null||samples.length<4||samples[0]==null||samples[0].length<2||samples[0].length>2000)
            throw new IllegalArgumentException("Require >=4 independent samples and 2..2000 features");
        int n=samples.length,p=samples[0].length;double[][] x=new double[n][p];double[] mean=new double[p],sd=new double[p];
        for(double[] row:samples) {
            if(row==null||row.length!=p)throw new IllegalArgumentException("Ragged matrix");
            for(int j=0;j<p;j++) {if(!Double.isFinite(row[j]))throw new IllegalArgumentException("Matrix must be finite");mean[j]+=row[j]/n;}
        }
        for(double[] row:samples)for(int j=0;j<p;j++){double d=row[j]-mean[j];sd[j]+=d*d/n;}
        for(int j=0;j<p;j++){sd[j]=Math.sqrt(sd[j]);if(!(sd[j]>0)||!Double.isFinite(sd[j]))throw new IllegalArgumentException("Constant or numerically invalid feature: "+j);}
        for(int i=0;i<n;i++)for(int j=0;j<p;j++)x[i][j]=(samples[i][j]-mean[j])/sd[j];
        return x;
    }
    /** Pearson correlation matrix. */
    public static double[][] correlation(double[][] samples) {
        double[][] x=standardize(samples);int n=x.length,p=x[0].length;double[][] result=new double[p][p];
        for(int a=0;a<p;a++)for(int b=a;b<p;b++){double v=0;for(int i=0;i<n;i++)v+=x[i][a]*x[i][b]/n;result[a][b]=result[b][a]=Math.max(-1,Math.min(1,v));}
        return result;
    }
    /** Benjamini-Hochberg correction over the complete supplied family. */
    public static double[] bh(double[] p) {
        Integer[] order=new Integer[p.length];for(int i=0;i<p.length;i++){if(!Double.isFinite(p[i])||p[i]<0||p[i]>1)throw new IllegalArgumentException("Invalid p value");order[i]=i;}
        Arrays.sort(order,Comparator.comparingDouble(i->p[i]));double[] q=new double[p.length];double running=1;
        for(int j=p.length-1;j>=0;j--){running=Math.min(running,p[order[j]]*p.length/(j+1));q[order[j]]=running;}return q;
    }
    /** Weighted random walk with restart; dangling mass returns to the seed distribution. */
    public static double[] propagate(List<Map<Integer,Double>> graph,double[] seeds,double restart,double tolerance,int maxIterations) {
        int n=graph.size();if(seeds.length!=n||n==0||!(restart>0&&restart<=1)||!(tolerance>0)||maxIterations<1)throw new IllegalArgumentException("Invalid propagation controls");
        double sum=0;for(double v:seeds){if(!Double.isFinite(v)||v<0)throw new IllegalArgumentException("Invalid seed weight");sum+=v;}
        if(!(sum>0)||!Double.isFinite(sum))throw new IllegalArgumentException("Require positive finite seed mass");
        double[] base=seeds.clone(),degree=new double[n];for(int i=0;i<n;i++){base[i]/=sum;for(var e:graph.get(i).entrySet()) {if(e.getKey()<0||e.getKey()>=n||!Double.isFinite(e.getValue())||e.getValue()<0)throw new IllegalArgumentException("Invalid graph weight");degree[i]+=e.getValue();}if(!Double.isFinite(degree[i]))throw new IllegalArgumentException("Degree overflow");}
        double[] score=base.clone();
        for(int iteration=0;iteration<maxIterations;iteration++) {
            double[] next=new double[n];double dangling=0;
            for(int i=0;i<n;i++)if(degree[i]==0)dangling+=score[i];else for(var e:graph.get(i).entrySet())next[e.getKey()]+=(1-restart)*score[i]*e.getValue()/degree[i];
            double change=0;for(int i=0;i<n;i++){next[i]+=(restart+(1-restart)*dangling)*base[i];change+=Math.abs(next[i]-score[i]);}score=next;
            if(change<tolerance)return score;
        }
        throw new IllegalArgumentException("Network propagation did not converge");
    }
}
