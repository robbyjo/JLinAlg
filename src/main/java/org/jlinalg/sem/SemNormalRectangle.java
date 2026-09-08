/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import jdistlib.Normal;

/** Error-controlled bivariate normal probabilities, with positive-integrand
 * rectangle integration when subtraction of CDFs is ill-conditioned. */
final class SemNormalRectangle {
    private SemNormalRectangle() { }
    private static final double[][] Q16=rule(16), Q32=rule(32);
    static double interval(double lo,double hi) {
        if(lo>=hi)return 0;
        if(lo>=0)return Normal.cumulative(lo,0,1,false,false)-Normal.cumulative(hi,0,1,false,false);
        return Normal.cumulative(hi,0,1,true,false)-Normal.cumulative(lo,0,1,true,false);
    }
    static double cdf(double h,double k,double rho) {
        if(h==Double.NEGATIVE_INFINITY||k==Double.NEGATIVE_INFINITY)return 0;
        if(h==Double.POSITIVE_INFINITY)return Normal.cumulative(k,0,1,true,false);
        if(k==Double.POSITIVE_INFINITY)return Normal.cumulative(h,0,1,true,false);
        double bound=Math.asin(rho);
        // r=sin(t) cancels the 1/sqrt(1-r*r) singularity of Plackett's formula.
        DoubleUnaryOperator f=t->{double r=Math.sin(t),c=Math.cos(t),z=(h-r*k)/c;return Math.exp(-.5*(z*z+k*k))/(2*Math.PI);};
        double integral=bound>=0?integrate(f,0,bound,2e-15,1e-12,0):-integrate(f,bound,0,2e-15,1e-12,0);
        double value=Normal.cumulative(h,0,1,true,false)*Normal.cumulative(k,0,1,true,false)+integral;
        if(value<1e-8)return positiveRectangle(Double.NEGATIVE_INFINITY,h,Double.NEGATIVE_INFINITY,k,rho);
        return value;
    }
    static double rectangle(double a,double b,double c,double d,double rho) {
        double p1=cdf(b,d,rho),p2=cdf(a,d,rho),p3=cdf(b,c,rho),p4=cdf(a,c,rho);
        double probability=(p1-p2)-(p3-p4);
        if(probability>1e-5 && probability>1e-8*(p1+p2+p3+p4))return probability;
        return positiveRectangle(a,b,c,d,rho);
    }
    private static double positiveRectangle(double a,double b,double c,double d,double rho) {
        if(interval(a,b)>interval(c,d)){double t=a;a=c;c=t;t=b;b=d;d=t;}
        final double lo=c,hi=d,sigma=Math.sqrt((1-rho)*(1+rho));
        double lower=Math.max(-40,a),upper=Math.min(40,b);if(lower>=upper)return 0;
        DoubleUnaryOperator f=x->Math.exp(-x*x/2)/Math.sqrt(2*Math.PI)*interval((lo-rho*x)/sigma,(hi-rho*x)/sigma);
        List<Double> cuts=new ArrayList<>();cuts.add(lower);cuts.add(upper);
        // Resolve narrow conditional transitions even if an initial quadrature
        // rule would put all its nodes outside their support.
        if(Math.abs(rho)>1e-8)for(double endpoint:new double[]{lo,hi})if(Double.isFinite(endpoint))
            for(double shift:new double[]{-8,0,8}){double x=(endpoint+shift*sigma)/rho;if(x>lower&&x<upper)cuts.add(x);}
        for(double x=Math.ceil(lower/2)*2;x<upper;x+=2)if(x>lower)cuts.add(x);
        cuts.sort(Double::compare);double total=0;
        for(int i=1;i<cuts.size();i++)total+=integrate(f,cuts.get(i-1),cuts.get(i),1e-300,2e-11,0);
        return total;
    }
    private static double integrate(DoubleUnaryOperator f,double lo,double hi,double abs,double rel,int depth) {
        if(lo==hi)return 0;
        double small=quadrature(f,lo,hi,Q16),large=quadrature(f,lo,hi,Q32);
        if(Math.abs(large-small)<=Math.max(abs,rel*Math.abs(large)))return large;
        if(depth>=20)throw new IllegalArgumentException("bivariate normal integration did not meet its error tolerance");
        double mid=(lo+hi)/2;
        return integrate(f,lo,mid,abs/2,rel,depth+1)+integrate(f,mid,hi,abs/2,rel,depth+1);
    }
    private static double quadrature(DoubleUnaryOperator f,double lo,double hi,double[][] q) {
        double sum=0,half=(hi-lo)/2,mid=(hi+lo)/2;
        for(int i=0;i<q[0].length;i++)sum+=q[1][i]*f.applyAsDouble(mid+half*q[0][i]);return sum*half;
    }
    private static double[][] rule(int n) {
        double[][] q=new double[2][n];
        for(int i=0;i<(n+1)/2;i++) {
            double z=Math.cos(Math.PI*(i+.75)/(n+.5)),derivative=0;
            for(int it=0;it<30;it++) {
                double a=1,b=0;for(int j=1;j<=n;j++){double c=b;b=a;a=((2*j-1)*z*b-(j-1)*c)/j;}
                derivative=n*(z*a-b)/(z*z-1);double next=z-a/derivative;
                if(Math.abs(next-z)<1e-15){z=next;break;}z=next;
            }
            q[0][i]=-z;q[0][n-1-i]=z;q[1][i]=q[1][n-1-i]=2/((1-z*z)*derivative*derivative);
        }
        return q;
    }
}
