/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Box-constrained BFGS for the complete profiled Laplace likelihood. */
final class LaplaceOptimization {
    private LaplaceOptimization() { }
    interface Objective { double value(double[] point); }
    record Fit(double[] point, double value, int iterations, boolean converged) { }
    static Fit fit(Objective objective, double[] initial, double[] lower, double[] upper,
            int limit, double tolerance, int fixed, double initialVarianceStep) {
        int n=initial.length;double[] x=initial.clone();
        for(int i=0;i<n;i++)x[i]=Math.max(lower[i],Math.min(upper[i],x[i]));
        double f=objective.value(x);double[][] inverse=initialMetric(n,fixed,initialVarianceStep);
        double[] g=gradient(objective,x,lower,upper);int iteration=0;
        for(;iteration<limit;iteration++) {
            if(norm(g,x,lower,upper)<=tolerance)break;
            double[] d=new double[n];
            for(int i=0;i<n;i++)for(int j=0;j<n;j++)d[i]-=inverse[i][j]*g[j];
            projectDirection(d,x,lower,upper);
            if(!(dot(d,g)<0)) {inverse=identity(n);for(int i=0;i<n;i++)d[i]=-g[i];projectDirection(d,x,lower,upper);}
            double max=0;for(double v:d)max=Math.max(max,Math.abs(v));
            if(max>3)for(int i=0;i<n;i++)d[i]*=3/max;
            double[] next=x.clone();double nextF=Double.POSITIVE_INFINITY;boolean accepted=false;
            for(double alpha=1;alpha>=1e-10;alpha*=.5) {
                double slope=0;for(int i=0;i<n;i++){next[i]=Math.max(lower[i],Math.min(upper[i],x[i]+alpha*d[i]));slope+=g[i]*(next[i]-x[i]);}
                nextF=objective.value(next);
                if(slope<0 && Double.isFinite(nextF) && nextF<=f+1e-4*slope){accepted=true;break;}
            }
            if(!accepted)break;
            double[] ng=gradient(objective,next,lower,upper),s=new double[n],y=new double[n],hy=new double[n];
            for(int i=0;i<n;i++){s[i]=next[i]-x[i];y[i]=ng[i]-g[i];}
            double sy=dot(s,y);
            if(sy>1e-12*Math.sqrt(dot(s,s)*dot(y,y))) {
                for(int i=0;i<n;i++)for(int j=0;j<n;j++)hy[i]+=inverse[i][j]*y[j];
                double a=(sy+dot(y,hy))/(sy*sy);
                for(int i=0;i<n;i++)for(int j=0;j<n;j++)inverse[i][j]+=a*s[i]*s[j]-(hy[i]*s[j]+s[i]*hy[j])/sy;
            }else inverse=identity(n);
            x=next.clone();g=ng;f=nextF;
        }
        return new Fit(x,f,iteration,Double.isFinite(f)&&norm(g,x,lower,upper)<=tolerance);
    }
    static double[] gradient(Objective objective,double[] x,double[] lower,double[] upper) {
        double[] g=new double[x.length];
        for(int i=0;i<x.length;i++) {
            double h=Math.min(1e-5*Math.max(1,Math.abs(x[i])),(upper[i]-lower[i])/4);
            double[] a=x.clone(),b=x.clone();
            if(!(h>0)||x[i]+h==x[i]||x[i]-h==x[i]){g[i]=Double.NaN;continue;}
            if(x[i]-h<lower[i]) {a[i]+=h;b[i]+=2*h;g[i]=(-3*objective.value(x)+4*objective.value(a)-objective.value(b))/(2*h);}
            else if(x[i]+h>upper[i]) {a[i]-=h;b[i]-=2*h;g[i]=(3*objective.value(x)-4*objective.value(a)+objective.value(b))/(2*h);}
            else {a[i]+=h;b[i]-=h;g[i]=(objective.value(a)-objective.value(b))/(2*h);}
        }
        return g;
    }
    static double[] fixedCovariance(Objective objective,double[] x,double[] lower,double[] upper,int fixed) {
        // At an active variance bound report covariance conditional on that bound.
        int[] active=new int[x.length];int count=0;
        for(int i=0;i<x.length;i++)if(i<fixed||(x[i]>lower[i]+1e-5&&x[i]<upper[i]-1e-5))active[count++]=i;
        double[][] hessian=new double[count][count];double center=objective.value(x);
        double[] steps=new double[count];
        for(int i=0;i<count;i++) {
            int a=active[i];steps[i]=Math.min(2e-4*Math.max(1,Math.abs(x[a])),
                .25*Math.min(x[a]-lower[a],upper[a]-x[a]));
        }
        for(int i=0;i<count;i++) {
            int a=active[i];double ha=steps[i];
            double[] plus=x.clone(),minus=x.clone();plus[a]+=ha;minus[a]-=ha;
            hessian[i][i]=(objective.value(plus)-2*center+objective.value(minus))/(ha*ha);
            for(int j=0;j<i;j++) {
                int b=active[j];double hb=steps[j];
                double[] pp=x.clone(),pm=x.clone(),mp=x.clone(),mm=x.clone();
                pp[a]+=ha;pp[b]+=hb;pm[a]+=ha;pm[b]-=hb;mp[a]-=ha;mp[b]+=hb;mm[a]-=ha;mm[b]-=hb;
                hessian[i][j]=(objective.value(pp)-objective.value(pm)-objective.value(mp)+objective.value(mm))/(4*ha*hb);hessian[j][i]=hessian[i][j];
            }
        }
        double[][] inverse=QuadratureOptimizer.inversePositiveDefinite(hessian);
        double[] covariance=new double[fixed*fixed];
        if(inverse==null){java.util.Arrays.fill(covariance,Double.NaN);return covariance;}
        for(int i=0;i<fixed;i++)for(int j=0;j<fixed;j++)covariance[i*fixed+j]=inverse[i][j];
        return covariance;
    }
    private static void projectDirection(double[] d,double[] x,double[] lo,double[] hi) {
        for(int i=0;i<d.length;i++)if((x[i]<=lo[i]&&d[i]<0)||(x[i]>=hi[i]&&d[i]>0))d[i]=0;
    }
    private static double norm(double[] g,double[] x,double[] lo,double[] hi) {
        double max=0;for(int i=0;i<g.length;i++)if(!((x[i]<=lo[i]&&g[i]>0)||(x[i]>=hi[i]&&g[i]<0)))max=Math.max(max,Math.abs(g[i]));return max;
    }
    private static double dot(double[] a,double[] b){double sum=0;for(int i=0;i<a.length;i++)sum+=a[i]*b[i];return sum;}
    private static double[][] identity(int n){double[][] a=new double[n][n];for(int i=0;i<n;i++)a[i][i]=1;return a;}
    private static double[][] initialMetric(int n,int fixed,double step){double[][] a=identity(n);for(int i=fixed;i<n;i++)a[i][i]=step*step;return a;}
}
