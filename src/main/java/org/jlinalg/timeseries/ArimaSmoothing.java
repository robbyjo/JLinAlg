/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.Arrays;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Exact Gaussian historical state conditioning, including symbolic diffuse levels.
 * Uses dense observation conditioning: O(n^3 + n^2 r^2), bounded to 1024 dates.
 * Uncertainty conditions on the supplied dynamics and innovation variance. */
public final class ArimaSmoothing {
    private ArimaSmoothing() { }
    public static Result smooth(double[] series,double[] ar,double[] ma,
            double[] differencingPolynomial,double innovationVariance) {
        if(series==null||series.length==0||series.length>1024||ar==null||ma==null
                ||differencingPolynomial==null||differencingPolynomial.length==0
                ||Math.max(ar.length,ma.length+1)+differencingPolynomial.length-1>128
                ||differencingPolynomial[0]!=1||!(innovationVariance>0)||!Double.isFinite(innovationVariance))
            throw new IllegalArgumentException("finite variance and 1..1024 dates required");
        for(double[] a:new double[][]{ar,ma,differencingPolynomial})
            for(double v:a)if(!Double.isFinite(v))throw new IllegalArgumentException("nonfinite dynamics");
        for(double v:series)if(Double.isInfinite(v))throw new IllegalArgumentException("use NaN for missing dates");
        ArimaStateSpace model=new ArimaStateSpace(ar,ma,differencingPolynomial);
        int n=series.length,r=model.stateSize(),d=r-model.stationarySize(),m=0;
        for(double v:series)if(Double.isFinite(v))m++;
        if(m<=d)throw new IllegalArgumentException("insufficient observations to identify diffuse states");
        int[] dates=new int[m];double[] y=new double[m];
        for(int t=0,i=0;t<n;t++)if(Double.isFinite(series[t])){dates[i]=t;y[i++]=series[t];}
        double[] transition=model.transitionMatrix(),z=model.observationVector();
        double[][] powers=new double[n][],p=new double[n][],diffuse=new double[n][];
        powers[0]=MatrixOps.identity(r);p[0]=model.initialCovariance();diffuse[0]=new double[r*d];
        for(int j=0;j<d;j++)diffuse[0][(model.stationarySize()+j)*d+j]=1;
        for(int t=1;t<n;t++) {
            powers[t]=mul(transition,r,r,powers[t-1],r);
            p[t]=mul(mul(transition,r,r,p[t-1],r),r,r,transpose(transition,r,r),r);
            double[] q=model.noiseMatrix();for(int j=0;j<q.length;j++)p[t][j]+=q[j];
            diffuse[t]=mul(transition,r,r,diffuse[t-1],d);
        }
        double[] v=new double[m*m],h=new double[m*d];
        for(int i=0;i<m;i++) {
            for(int j=0;j<d;j++)for(int k=0;k<r;k++)h[i*d+j]+=z[k]*diffuse[dates[i]][k*d+j];
            for(int j=0;j<=i;j++) {
                double[] c=cross(dates[i],dates[j],r,p,powers);double value=0;
                for(int a=0;a<r;a++)for(int b=0;b<r;b++)value+=z[a]*c[a*r+b]*z[b];
                v[i*m+j]=v[j*m+i]=value;
            }
        }
        double[][] means=new double[n][r],covariances=new double[n][r*r];
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            var factor=context.backend().dpotrf(v,m);
            double[] vy=factor.solve(y),vh=d==0?new double[0]:factor.solve(h,d);
            double[] infoInverse=new double[d*d],level=new double[d];
            if(d>0) {
                double[] information=mul(transpose(h,m,d),d,m,vh,d);
                var df=context.backend().dpotrf(information,d);
                infoInverse=df.solve(MatrixOps.identity(d),d);
                level=df.solve(mul(transpose(h,m,d),d,m,vy,1));
                for(int i=0;i<m;i++)for(int j=0;j<d;j++)vy[i]-=vh[i*d+j]*level[j];
            }
            for(int t=0;t<n;t++) {
                double[] c=new double[r*m];
                for(int j=0;j<m;j++) {
                    double[] block=cross(t,dates[j],r,p,powers);
                    for(int a=0;a<r;a++)for(int b=0;b<r;b++)c[a*m+j]+=block[a*r+b]*z[b];
                }
                double[] mean=mul(c,r,m,vy,1),dc=mul(c,r,m,vh,d);
                for(int a=0;a<r;a++)for(int j=0;j<d;j++) {
                    mean[a]+=diffuse[t][a*d+j]*level[j];
                    dc[a*d+j]=diffuse[t][a*d+j]-dc[a*d+j];
                }
                double[] reduction=mul(c,r,m,factor.solve(transpose(c,r,m),r),r);
                double[] addition=mul(mul(dc,r,d,infoInverse,d),r,d,transpose(dc,r,d),r);
                for(int j=0;j<r*r;j++)covariances[t][j]=innovationVariance*(p[t][j]-reduction[j]+addition[j]);
                for(int j=0;j<r;j++) {
                    double variance=covariances[t][j*r+j];
                    if(variance < -1e-8*innovationVariance*(1+Math.abs(p[t][j*r+j])))
                        throw new IllegalArgumentException("smoothing covariance lost precision");
                    covariances[t][j*r+j]=Math.max(0,variance);
                }
                means[t]=mean;
            }
        }
        double[] signal=new double[n],variance=new double[n];
        for(int t=0;t<n;t++)for(int a=0;a<r;a++) {
            signal[t]+=z[a]*means[t][a];
            for(int b=0;b<r;b++)variance[t]+=z[a]*covariances[t][a*r+b]*z[b];
        }
        for(int t=0;t<n;t++)variance[t]=Math.max(0,variance[t]);
        return new Result(means,covariances,signal,variance);
    }
    private static double[] cross(int t,int s,int r,double[][] p,double[][] powers) {
        if(t>=s)return mul(powers[t-s],r,r,p[s],r);
        return transpose(mul(powers[s-t],r,r,p[t],r),r,r);
    }
    private static double[] mul(double[] a,int rows,int shared,double[] b,int cols) {
        double[] c=new double[rows*cols];
        for(int i=0;i<rows;i++)for(int k=0;k<shared;k++)for(int j=0;j<cols;j++)c[i*cols+j]+=a[i*shared+k]*b[k*cols+j];
        return c;
    }
    private static double[] transpose(double[] a,int rows,int cols) {
        double[] b=new double[a.length];for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)b[j*rows+i]=a[i*cols+j];return b;
    }
    public record Result(double[][] states,double[][] stateCovariances,double[] signal,double[] signalVariance) {
        public Result {states=copy(states);stateCovariances=copy(stateCovariances);signal=signal.clone();signalVariance=signalVariance.clone();}
        public double[][] states(){return copy(states);}
        public double[][] stateCovariances(){return copy(stateCovariances);}
        public double[] signal(){return signal.clone();}
        public double[] signalVariance(){return signalVariance.clone();}
        private static double[][] copy(double[][] a){return Arrays.stream(a).map(double[]::clone).toArray(double[][]::new);}
    }
}
