/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

/** Deterministic conditional-normal rectangle integration (up to four variables).
 * Refinement is checked; unresolved probabilities are rejected. */
final class SemRectangle {
    private SemRectangle() { }
    static double probability(double[] lo,double[] hi,double[] covariance) {
        int d=lo.length;if(d==0)return 1;
        if(d>4)throw new IllegalArgumentException("joint ordinal likelihood supports at most four observed ordinal responses per row");
        double[] scale=new double[d],a=lo.clone(),b=hi.clone(),cor=covariance.clone();
        for(int i=0;i<d;i++) {
            scale[i]=Math.sqrt(covariance[i*d+i]);
            if(!(scale[i]>0))throw new IllegalArgumentException("nonpositive conditional variance");
            a[i]/=scale[i];b[i]/=scale[i];
        }
        for(int i=0;i<d;i++)for(int j=0;j<d;j++)cor[i*d+j]/=scale[i]*scale[j];
        RamFit.chol(cor,d);
        if(d==1)return SemNormalRectangle.interval(a[0],b[0]);
        if(d==2)return SemNormalRectangle.rectangle(a[0],b[0],a[1],b[1],cor[1]);
        double previous=integrate(a,b,cor,8);
        for(int order=16;order<=256;order*=2) {
            double value=integrate(a,b,cor,order);
            if(Math.abs(value-previous)<=1e-10*Math.max(value,1e-280))return value;
            previous=value;
        }
        throw new IllegalArgumentException("ordinal rectangle refinement tolerance not met");
    }
    private static double integrate(double[] lo,double[] hi,double[] cor,int order) {
        int d=lo.length;
        if(d==2)return SemNormalRectangle.rectangle(lo[0],hi[0],lo[1],hi[1],cor[1]);
        // Integrate in the normal coordinate. Inverse-CDF integration has
        // endpoint derivative singularities and converges poorly for orthants.
        // Beyond +/-40 the normal density is below double precision range.
        double left=Math.max(-40,lo[0]),right=Math.min(40,hi[0]),width=right-left;
        if(!(width>0))throw new IllegalArgumentException("ordinal tail is beyond representable integration range");
        double[][] rule=rule(order);double total=0;
        double[] conditional=new double[(d-1)*(d-1)],sd=new double[d-1];
        for(int j=1;j<d;j++)sd[j-1]=Math.sqrt(1-cor[j*d]*cor[j*d]);
        for(int i=1;i<d;i++)for(int j=1;j<d;j++)conditional[(i-1)*(d-1)+j-1]=(cor[i*d+j]-cor[i*d]*cor[j*d])/(sd[i-1]*sd[j-1]);
        for(int k=0;k<order;k++) {
            double u=(rule[0][k]+1)/2;
            double x=left+u*width;
            double[] a=new double[d-1],b=new double[d-1];
            for(int j=1;j<d;j++){a[j-1]=(lo[j]-cor[j*d]*x)/sd[j-1];b[j-1]=(hi[j]-cor[j*d]*x)/sd[j-1];}
            total+=rule[1][k]*Math.exp(-.5*x*x)/Math.sqrt(2*Math.PI)*integrate(a,b,conditional,order)/2;
        }
        return width*total;
    }
    private static final java.util.Map<Integer,double[][]> RULES=new java.util.concurrent.ConcurrentHashMap<>();
    private static double[][] rule(int n) {return RULES.computeIfAbsent(n,key->{
        double[][] result=new double[2][n];
        for(int i=0;i<(n+1)/2;i++) {
            double x=Math.cos(Math.PI*(i+.75)/(n+.5)),derivative=0;
            for(int iteration=0;iteration<30;iteration++) {
                double p=1,previous=0;for(int j=1;j<=n;j++){double old=p;p=((2*j-1)*x*p-(j-1)*previous)/j;previous=old;}
                derivative=n*(x*p-previous)/(x*x-1);double step=p/derivative;x-=step;if(Math.abs(step)<1e-15)break;
            }
            result[0][i]=-x;result[0][n-1-i]=x;
            result[1][i]=result[1][n-1-i]=2/((1-x*x)*derivative*derivative);
        }return result;
    });}
}
