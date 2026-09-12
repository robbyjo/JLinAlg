/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import java.util.Arrays;
import jdistlib.Normal;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.MatrixOps;

/** Product-kernel local linear regression for continuous, unordered and ordered
 * predictors. Categorical values are finite numeric codes (ordered codes have
 * meaningful integer distances). Uses Gaussian, Aitchison-Aitken and geometric
 * kernels. Local linear terms are included for continuous coordinates only. */
public final class ProductKernelRegression {
    private ProductKernelRegression() { }
    public enum Type { CONTINUOUS, UNORDERED, ORDERED }
    public static double predict(double[][] x,double[] y,double[] query,Type[] types,double[] bandwidths) {
        validate(x,y,query,types,bandwidths);
        return dot(weights(x,query,types,bandwidths),y);
    }
    /** Automatic pointwise inference for independent rows with smooth conditional
     * means/variances, interior continuous query, and positive design density.
     * Categorical strata are matched exactly. Continuous bandwidths are a robust
     * scale times n_stratum^(-1/(d+2)); this undersmooths a second-order local
     * linear estimator: n h^d -> infinity and n h^(d+4) -> 0. The rule depends
     * on predictors only. No finite-sample or boundary coverage guarantee. */
    public static Inference infer(double[][] x,double[] y,double[] query,Type[] types,double level) {
        if(types==null||query==null||types.length!=query.length||!(level>0&&level<1))
            throw new IllegalArgumentException("types, query and confidence level required");
        double[] h=new double[types.length];
        for(int j=0;j<h.length;j++)h[j]=types[j]==Type.CONTINUOUS?1:0;
        validate(x,y,query,types,h);
        int d=(int)Arrays.stream(types).filter(t->t==Type.CONTINUOUS).count();
        int[] stratum=java.util.stream.IntStream.range(0,y.length).filter(i->{
            for(int j=0;j<types.length;j++)if(types[j]!=Type.CONTINUOUS&&x[i][j]!=query[j])return false;
            return true;
        }).toArray();
        if(stratum.length<=d+3)throw new IllegalArgumentException("insufficient rows in query stratum");
        for(int j=0;j<h.length;j++)if(types[j]==Type.CONTINUOUS) {
            double[] column=new double[stratum.length];for(int i=0;i<column.length;i++)column[i]=x[stratum[i]][j];
            Arrays.sort(column);double mean=Arrays.stream(column).average().orElseThrow(),ss=0;
            for(double value:column)ss+=Math.pow(value-mean,2);
            double sd=Math.sqrt(ss/(column.length-1)),iqr=(quantile(column,.75)-quantile(column,.25))/1.349;
            double scale=iqr>0?Math.min(sd,iqr):sd;
            if(!(scale>0)||!Double.isFinite(scale)||!(query[j]>column[0]&&query[j]<column[column.length-1]))
                throw new IllegalArgumentException("inference needs varying continuous predictors and interior queries");
            h[j]=scale*Math.pow(stratum.length,-1.0/(d+2));
        }
        double[] influence=weights(x,query,types,h);double estimate=dot(influence,y),variance=0,sumSquares=dot(influence,influence);
        for(int i:stratum) {
            double[] row=weights(x,x[i],types,h);double residual=y[i]-dot(row,y),remaining=1-row[i];
            if(!(remaining>1e-8))throw new IllegalArgumentException("insufficient local support for HC3 inference");
            variance+=Math.pow(influence[i]*residual/remaining,2);
        }
        if(!(variance>0)||!Double.isFinite(variance))throw new IllegalArgumentException("unresolved local variance");
        double se=Math.sqrt(variance),critical=Normal.quantile(.5+level/2,0,1,true,false);
        return new Inference(estimate,se,estimate-critical*se,estimate+critical*se,h,1/sumSquares,stratum.length,level);
    }
    private static double[] weights(double[][] x,double[] query,Type[] types,double[] h) {
        int n=x.length,d=(int)Arrays.stream(types).filter(t->t==Type.CONTINUOUS).count(),p=d+1;
        double[] log=new double[n];double maximum=Double.NEGATIVE_INFINITY;
        int[] levels=new int[types.length];
        for(int j=0;j<types.length;j++)if(types[j]==Type.UNORDERED) {
            final int column=j;levels[j]=(int)Arrays.stream(x).mapToDouble(row->row[column]).distinct().count();
        }
        for(int i=0;i<n;i++) {
            for(int j=0;j<types.length;j++) {
                double delta=x[i][j]-query[j];
                log[i]+=switch(types[j]) {
                    case CONTINUOUS -> -.5*Math.pow(delta/h[j],2);
                    case UNORDERED -> Math.log(delta==0?1-h[j]:(levels[j]>1?h[j]/(levels[j]-1):0));
                    case ORDERED -> delta==0?0:Math.abs(delta)*Math.log(h[j]);
                };
            }
            maximum=Math.max(maximum,log[i]);
        }
        if(!Double.isFinite(maximum))throw new IllegalArgumentException("query has no kernel support");
        double[] design=new double[n*p],w=new double[n];
        for(int i=0;i<n;i++) {
            w[i]=Math.exp(.5*(log[i]-maximum));design[i*p]=w[i];
            for(int j=0,k=1;j<types.length;j++)if(types[j]==Type.CONTINUOUS)design[i*p+k++]=w[i]*(x[i][j]-query[j])/h[j];
        }
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            double[] inverse=LeastSquaresSolver.solve(design,new double[n],n,p,false,context.backend()).unscaledCovariance();
            double[] result=new double[n];
            for(int i=0;i<n;i++)for(int j=0;j<p;j++)result[i]+=inverse[j]*design[i*p+j]*w[i];
            return result;
        }
    }
    private static void validate(double[][] x,double[] y,double[] query,Type[] types,double[] h) {
        if(y==null||query==null||types==null||h==null||y.length<2||query.length==0
                ||types.length!=query.length||h.length!=query.length)throw new IllegalArgumentException("invalid product kernel dimensions");
        MatrixOps.rowMajor(x,y.length);
        if(x[0].length!=query.length)throw new IllegalArgumentException("query dimension differs");
        for(double v:y)if(!Double.isFinite(v))throw new IllegalArgumentException("finite response required");
        for(int j=0;j<h.length;j++) {
            if(types[j]==null||!Double.isFinite(query[j])||!Double.isFinite(h[j])
                    ||(types[j]==Type.CONTINUOUS?!(h[j]>0):!(h[j]>=0&&h[j]<1)))
                throw new IllegalArgumentException("invalid type, query or bandwidth");
            if(types[j]!=Type.CONTINUOUS) {
                boolean known=false;for(double[] row:x){known|=row[j]==query[j];
                    if(types[j]==Type.ORDERED&&row[j]!=Math.rint(row[j]))throw new IllegalArgumentException("ordered codes must be integers");}
                if(!known)throw new IllegalArgumentException("query category is absent from training data");
            }
        }
    }
    private static double quantile(double[] a,double p){double x=p*(a.length-1);int i=(int)x;return a[i]+(x-i)*(a[Math.min(i+1,a.length-1)]-a[i]);}
    private static double dot(double[] a,double[] b){double result=0;for(int i=0;i<a.length;i++)result+=a[i]*b[i];return result;}
    public record Inference(double estimate,double standardError,double lower,double upper,
            double[] bandwidths,double effectiveSampleSize,int stratumSize,double confidenceLevel) {
        public Inference{bandwidths=bandwidths.clone();}
        public double[] bandwidths(){return bandwidths.clone();}
    }
}
