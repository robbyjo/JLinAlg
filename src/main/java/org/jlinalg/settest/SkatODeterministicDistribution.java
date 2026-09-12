/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import java.util.*;
import java.util.function.DoubleUnaryOperator;
import jdistlib.ChiSquare;
import jdistlib.Normal;
import jdistlib.accelerator.ComputeBackend;

/** Conditional Gaussian integration, without moment matching or Monte Carlo.
 * Positive-series remainder <= 1e-10 times min-p; scaled quadrature tolerance
 * 1e-6 times min-p. Unresolved spectra/integrals fail, never approximate silently.
 * Quadrature errors are estimates, not interval-arithmetic guarantees. */
final class SkatODeterministicDistribution {
    private SkatODeterministicDistribution() { }
    static double adjustedPValue(double minimum, double[] rho,
            List<SkatOResult.Component> components, SetTestScoreState state, ComputeBackend backend) {
        if(minimum==1)return 1;
        if(!(minimum>=1e-250))throw failure("min-p below supported 1e-250 tail");
        int n=state.variants();double[] v=state.information(),a=new double[n];double b=0;
        for(int i=0;i<n;i++){for(int j=0;j<n;j++)a[i]+=v[i*n+j];b+=a[i];}
        if(!(b>0))return minimum; // all rho kernels are proportional when burden is null
        double[] critical=new double[rho.length];
        double upper=Normal.quantile(minimum*1e-10/2,0,1,false,false);
        for(int h=0;h<rho.length;h++) {
            double[] eigen=components.get(h).result().eigenvalues();
            if(components.get(h).result().pValueMethod().contains("saddlepoint"))throw failure("component tail used a saddlepoint approximation");
            critical[h]=strictCritical(eigen,minimum);
            if(rho[h]==1)upper=Math.min(upper,Math.sqrt(critical[h]/b));
        }
        double[] c=v.clone();for(int i=0;i<n;i++)a[i]/=Math.sqrt(b);
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)c[i*n+j]-=a[i]*a[j];
        var eig=backend.dsyev(c,n);double[] d=eig.eigenvalues(),e=eig.eigenvectors();
        double scale=Arrays.stream(v).map(Math::abs).max().orElseThrow(),zero=0;
        List<Double> lambda=new ArrayList<>(),shift=new ArrayList<>();
        for(int j=0;j<n;j++) {
            double projection=0;for(int i=0;i<n;i++)projection+=e[i*n+j]*a[i];
            if(d[j]<-scale*1e-10)throw failure("conditional covariance is indefinite");
            if(d[j]<=scale*1e-12)zero+=projection*projection;
            else {lambda.add(d[j]);shift.add(projection*projection/d[j]);}
        }
        double[] eigen=lambda.stream().mapToDouble(Double::doubleValue).toArray();
        double[] noncentral=shift.stream().mapToDouble(Double::doubleValue).toArray();
        final double burden=b,offset=zero;
        DoubleUnaryOperator f=t->{
            double bound=Double.POSITIVE_INFINITY;
            for(int h=0;h<rho.length;h++)if(rho[h]<1)
                bound=Math.min(bound,(critical[h]-rho[h]*burden*t*t)/(1-rho[h]));
            double p=tail(bound-offset*t*t,eigen,noncentral,t*t,minimum*1e-10);
            return Math.sqrt(2/Math.PI)*Math.exp(-t*t/2)*(p/minimum);
        };
        // Split where the active quadratic threshold changes, and use short
        // panels so rare-tail peaks cannot hide between quadrature nodes.
        TreeSet<Double> knots=new TreeSet<>();knots.add(0.0);knots.add(upper);
        for(double t=.125;t<upper;t+=.125)knots.add(t);
        for(int i=0;i<rho.length;i++)for(int j=0;j<i;j++)if(rho[i]<1&&rho[j]<1) {
            double denominator=burden*(rho[i]/(1-rho[i])-rho[j]/(1-rho[j]));
            double sq=(critical[i]/(1-rho[i])-critical[j]/(1-rho[j]))/denominator;
            if(sq>0&&sq<upper*upper)knots.add(Math.sqrt(sq));
        }
        double integral=0,left=0;int panels=knots.size()-1;
        for(double right:knots)if(right>left) {
            double mid=(left+right)/2,fa=f.applyAsDouble(left),fm=f.applyAsDouble(mid),fb=f.applyAsDouble(right);
            integral+=integrate(f,left,right,fa,fm,fb,(right-left)*(fa+4*fm+fb)/6,1e-6/panels,20);
            left=right;
        }
        double result=minimum*integral+2*Normal.cumulative(upper,0,1,false,false);
        if(!Double.isFinite(result)||result<minimum*(1-1e-5)||result>minimum*rho.length*(1+1e-5))
            throw failure("union-probability bounds violated");
        return Math.min(1,Math.max(minimum,result));
    }
    private static double strictCritical(double[] eigen,double p) {
        double lo=0,hi=Arrays.stream(eigen).max().orElseThrow()*ChiSquare.quantile(p,eigen.length,false,false);
        for(int i=0;i<100;i++) {
            double mid=(lo+hi)/2;
            var tail=QuadraticFormDistribution.survival(mid,eigen);
            if(tail.method().contains("saddlepoint"))throw failure("unresolved component quantile");
            if(tail.pValue()>p)lo=mid;else hi=mid;
        }
        return (lo+hi)/2;
    }
    /** Positive gamma mixture for sum lambda*(Z+delta)^2. The generating
     * function supplies a Chernoff bound on omitted mixture weight. */
    static double tail(double q,double[] lambda,double[] shifts,double t2,double tolerance) {
        if(q<=0)return 1;
        if(lambda.length==0)return 0;
        double beta=Arrays.stream(lambda).min().orElseThrow();
        double[] r=new double[lambda.length],powers=new double[lambda.length];double logC=0,maxR=0;
        for(int i=0;i<r.length;i++) {
            r[i]=1-beta/lambda[i];maxR=Math.max(maxR,r[i]);powers[i]=1;
            logC+=.5*Math.log(beta/lambda[i])-.5*shifts[i]*t2;
        }
        double z=maxR==0?1.5:1+.5*(1/maxR-1),logG=logC;
        for(int i=0;i<r.length;i++)logG+=-.5*Math.log1p(-r[i]*z)+.5*shifts[i]*t2*(1-r[i])*z/(1-r[i]*z);
        if(!Double.isFinite(logG)||Math.exp(logC)==0)throw failure("noncentral series dynamic range");
        int limit=8192;double[] weights=new double[limit+1],coefficients=new double[limit+1];weights[0]=Math.exp(logC);double sum=0;
        for(int k=0;k<=limit;k++) {
            if(k>0) {
                for(int i=0;i<r.length;i++) {
                    coefficients[k]+=.5*powers[i]*(r[i]+shifts[i]*t2*(1-r[i])*k);
                    powers[i]*=r[i];
                }
                double value=0;for(int j=1;j<=k;j++)value+=coefficients[j]*weights[k-j];weights[k]=value/k;
            }
            sum+=weights[k]*ChiSquare.cumulative(q/beta,lambda.length+2.0*k,false,false);
            if(logG-(k+1)*Math.log(z)<Math.log(tolerance))return Math.min(1,sum);
        }
        throw failure("noncentral series exceeded 8192 terms");
    }
    private static double integrate(DoubleUnaryOperator f,double a,double b,double fa,double fm,double fb,double whole,double tol,int depth) {
        double m=(a+b)/2,fl=f.applyAsDouble((a+m)/2),fr=f.applyAsDouble((m+b)/2);
        double l=(m-a)*(fa+4*fl+fm)/6,r=(b-m)*(fm+4*fr+fb)/6,delta=l+r-whole;
        if(Math.abs(delta)<=15*tol)return l+r+delta/15;
        if(depth==0)throw failure("quadrature did not meet relative tolerance");
        return integrate(f,a,m,fa,fl,fm,l,tol/2,depth-1)+integrate(f,m,b,fm,fr,fb,r,tol/2,depth-1);
    }
    private static IllegalArgumentException failure(String reason) {
        return new IllegalArgumentException("deterministic SKAT-O unresolved: "+reason);
    }
}
