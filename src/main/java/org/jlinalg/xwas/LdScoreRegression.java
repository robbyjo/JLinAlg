/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.xwas;

import java.util.*;

/** One-step, unpartitioned LD score regression with two weight updates and common block deletes. */
public final class LdScoreRegression {
    private LdScoreRegression() { }
    /** Pairs use upper triangle row order (0,0), (0,1), ..., (1,1), ... . */
    public record Result(int traits, double[] geneticCovariance, double[] samplingCovariance,
            double[] intercepts, double[] interceptStandardErrors, double[][] deleteValues,
            double[] geneticCorrelations, double[] correlationStandardErrors) {
        public Result {
            geneticCovariance=geneticCovariance.clone(); samplingCovariance=samplingCovariance.clone();
            intercepts=intercepts.clone(); interceptStandardErrors=interceptStandardErrors.clone();
            deleteValues=Arrays.stream(deleteValues).map(double[]::clone).toArray(double[][]::new);
            geneticCorrelations=geneticCorrelations.clone(); correlationStandardErrors=correlationStandardErrors.clone();
        }
        @Override public double[] geneticCovariance(){return geneticCovariance.clone();}
        @Override public double[] samplingCovariance(){return samplingCovariance.clone();}
        @Override public double[] intercepts(){return intercepts.clone();}
        @Override public double[] interceptStandardErrors(){return interceptStandardErrors.clone();}
        @Override public double[][] deleteValues(){return Arrays.stream(deleteValues).map(double[]::clone).toArray(double[][]::new);}
        @Override public double[] geneticCorrelations(){return geneticCorrelations.clone();}
        @Override public double[] correlationStandardErrors(){return correlationStandardErrors.clone();}
    }
    /** z and sampleSizes are variant by trait, already aligned to one effect allele. */
    public static Result fit(double[] ld, double[] weightLd, double[][] z,
            double[][] sampleSizes, int[] block, double referenceVariants) {
        SummaryMath.finite(ld); SummaryMath.finite(weightLd);
        int n=ld.length;
        if(n<6 || weightLd.length!=n || z.length!=n || sampleSizes.length!=n || block.length!=n
                || !(referenceVariants>0) || !Double.isFinite(referenceVariants))
            throw new IllegalArgumentException("invalid LDSC dimensions or reference variant count");
        int t=z[0].length, pairs=t*(t+1)/2;
        if(t<1 || t>32) throw new IllegalArgumentException("LDSC supports 1 to 32 traits");
        int b=0;
        for(int i=0;i<n;i++) {
            SummaryMath.finite(z[i]); SummaryMath.finite(sampleSizes[i]);
            if(z[i].length!=t || sampleSizes[i].length!=t || ld[i]<0 || weightLd[i]<0 || block[i]<0)
                throw new IllegalArgumentException("invalid LDSC row");
            for(double count:sampleSizes[i]) if(!(count>0)) throw new IllegalArgumentException("N must be positive");
            if(i>0 && block[i]!=block[i-1] && block[i]!=block[i-1]+1)
                throw new IllegalArgumentException("blocks must be contiguous and numbered from zero");
            b=Math.max(b,block[i]+1);
        }
        if(block[0]!=0 || b<3 || b>n) throw new IllegalArgumentException("at least three blocks are required");
        int[] sizes=new int[b]; for(int v:block)sizes[v]++;
        if(Arrays.stream(sizes).max().orElseThrow()-Arrays.stream(sizes).min().orElseThrow()>1)
            throw new IllegalArgumentException("jackknife blocks must have equal row counts within one row");
        double[] s=new double[t*t], intercept=new double[t*t], ise=new double[t*t];
        double[][] deletes=new double[b][pairs];
        Fit[][] fits=new Fit[t][t];
        for(int i=0;i<t;i++) fits[i][i]=regress(ld,weightLd,z,sampleSizes,block,b,referenceVariants,i,i,null,null);
        for(int i=0;i<t;i++)for(int j=i+1;j<t;j++)
            fits[i][j]=regress(ld,weightLd,z,sampleSizes,block,b,referenceVariants,i,j,fits[i][i],fits[j][j]);
        int k=0;
        for(int i=0;i<t;i++)for(int j=i;j<t;j++,k++) {
            Fit f=fits[i][j];
            s[i*t+j]=s[j*t+i]=f.estimate[0]; intercept[i*t+j]=intercept[j*t+i]=f.estimate[1];
            ise[i*t+j]=ise[j*t+i]=Math.sqrt(SummaryMath.jackknifeCovariance(f.deletes)[3]);
            for(int r=0;r<b;r++) deletes[r][k]=f.deletes[r][0];
        }
        double[] v=SummaryMath.jackknifeCovariance(deletes), rg=new double[t*t], rgse=new double[t*t];
        Arrays.fill(rg,Double.NaN);Arrays.fill(rgse,Double.NaN);
        for(int i=0;i<t;i++)for(int j=i;j<t;j++) {
            if(s[i*t+i]<=0 || s[j*t+j]<=0)continue;
            rg[i*t+j]=rg[j*t+i]=s[i*t+j]/Math.sqrt(s[i*t+i]*s[j*t+j]);
            double[][] dr=new double[b][1]; boolean valid=true;
            for(int r=0;r<b;r++) {
                double h1=fits[i][i].deletes[r][0],h2=fits[j][j].deletes[r][0];
                if(h1<=0 || h2<=0){valid=false;break;}
                dr[r][0]=fits[i][j].deletes[r][0]/Math.sqrt(h1*h2);
            }
            if(valid)rgse[i*t+j]=rgse[j*t+i]=Math.sqrt(SummaryMath.jackknifeCovariance(dr)[0]);
        }
        return new Result(t,s,v,intercept,ise,deletes,rg,rgse);
    }
    private record Fit(double[] estimate,double[][] deletes) { }
    private static double clamp(double x,double lo,double hi){return Math.max(lo,Math.min(hi,x));}
    private static Fit regress(double[] ld,double[] wld,double[][] z,double[][] ns,int[] blocks,
            int b,double m,int a,int c,Fit h1,Fit h2) {
        int n=ld.length; double[] x=new double[n],y=new double[n],w=new double[n];
        double sy=0,sx=0;
        for(int r=0;r<n;r++){x[r]=Math.sqrt(ns[r][a])*Math.sqrt(ns[r][c])*ld[r]/m;y[r]=z[r][a]*z[r][c];sy+=y[r];sx+=x[r];}
        double[] estimate={(sy-(a==c?n:0))/sx,a==c?1:0};
        for(int update=0;update<=2;update++) {
            for(int r=0;r<n;r++) {
                double l=Math.max(1,ld[r]),variance;
                if(a==c){double e=estimate[1]+clamp(estimate[0],0,1)*ns[r][a]*l/m;variance=2*e*e;}
                else {
                    double e1=h1.estimate[1]+clamp(h1.estimate[0],0,1)*ns[r][a]*l/m;
                    double e2=h2.estimate[1]+clamp(h2.estimate[0],0,1)*ns[r][c]*l/m;
                    double ec=estimate[1]+clamp(estimate[0],-1,1)*Math.sqrt(ns[r][a])*Math.sqrt(ns[r][c])*l/m;
                    variance=e1*e2+ec*ec;
                }
                w[r]=1/(Math.max(1,wld[r])*variance);
                if(!(w[r]>0)||!Double.isFinite(w[r]))throw new IllegalArgumentException("invalid LDSC regression weight");
            }
            estimate=wls(x,y,w,blocks,-1);
        }
        Moments[] pieces=new Moments[b];for(int r=0;r<b;r++)pieces[r]=new Moments();
        Moments total=new Moments();
        for(int r=0;r<n;r++){total.add(x[r],y[r],w[r]);pieces[blocks[r]].add(x[r],y[r],w[r]);}
        double[][] deleted=new double[b][];
        for(int r=0;r<b;r++) {
            Moments part=pieces[r];double remain=total.weight-part.weight;
            double mx=total.mx+(total.mx-part.mx)*part.weight/remain;
            double my=total.my+(total.my-part.my)*part.weight/remain;
            double factor=total.weight*part.weight/remain;
            double xx=total.xx-part.xx-factor*Math.pow(total.mx-part.mx,2);
            double xy=total.xy-part.xy-factor*(total.mx-part.mx)*(total.my-part.my);
            if(n-part.count<3||!(xx>1e-14*remain*Math.max(1,mx*mx)))throw new IllegalArgumentException("delete-block LDSC design is rank deficient");
            double slope=xy/xx;deleted[r]=new double[]{slope,my-slope*mx};SummaryMath.finite(deleted[r]);
        }
        return new Fit(estimate,deleted);
    }
    private static final class Moments {
        double weight,mx,my,xx,xy;int count;
        void add(double x,double y,double w) {
            double next=weight+w,dx=x-mx,dy=y-my;
            xx+=w*weight/next*dx*dx;xy+=w*weight/next*dx*dy;
            mx+=w/next*dx;my+=w/next*dy;weight=next;count++;
        }
    }
    // Centering avoids subtracting raw second moments in the slope normal equation.
    private static double[] wls(double[] x,double[] y,double[] w,int[] block,int omit) {
        double sw=0,mx=0,my=0;int count=0;
        for(int r=0;r<x.length;r++)if(block[r]!=omit){sw+=w[r];mx+=w[r]*x[r];my+=w[r]*y[r];count++;}
        mx/=sw;my/=sw;double xx=0,xy=0;
        for(int r=0;r<x.length;r++)if(block[r]!=omit){xx+=w[r]*(x[r]-mx)*(x[r]-mx);xy+=w[r]*(x[r]-mx)*(y[r]-my);}
        if(count<3 || !(xx>1e-14*sw*Math.max(1,mx*mx)))throw new IllegalArgumentException("LDSC design has insufficient LD-score variation");
        double slope=xy/xx;double[] result={slope,my-slope*mx};SummaryMath.finite(result);return result;
    }
}
