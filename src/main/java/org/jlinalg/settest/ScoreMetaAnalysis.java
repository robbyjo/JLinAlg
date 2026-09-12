/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import java.util.Arrays;

/** Fixed-effect aggregation of aligned, independent cohort score summaries.
 * Arrays are cohort by variant (scores) and cohort by row-major variant pair
 * (covariances). NaN scores denote unavailable variants; corresponding covariance
 * entries are ignored. Inputs must already share allele and phenotype units.
 */
public final class ScoreMetaAnalysis {
    private ScoreMetaAnalysis() { }

    /** Pool variants meeting the minimum number of informative cohorts.
     * A cohort is informative for a variant only when its variance is positive.
     * Returns null when no variants meet the threshold.
     */
    public static Pooled pool(double[][] scores, double[][] covariances, int minimumCohorts) {
        if(scores==null || covariances==null || scores.length==0 || scores.length!=covariances.length
                || minimumCohorts<1 || scores[0]==null || scores[0].length==0)
            throw new IllegalArgumentException("invalid cohort score dimensions or minimum cohorts");
        int n=scores[0].length,k=scores.length;
        int[] counts=new int[n];
        for(int c=0;c<k;c++) {
            if(scores[c]==null || scores[c].length!=n || covariances[c]==null
                    || covariances[c].length!=Math.multiplyExact(n,n))
                throw new IllegalArgumentException("cohort score dimensions differ");
            for(int i=0;i<n;i++) if(!Double.isNaN(scores[c][i])) {
                double variance=covariances[c][i*n+i];
                if(!Double.isFinite(scores[c][i]) || !Double.isFinite(variance) || variance<0
                        || variance==0 && scores[c][i]!=0)
                    throw new IllegalArgumentException("invalid cohort score or diagonal variance");
                if(variance>0) counts[i]++;
            }
            int present=0;for(int i=0;i<n;i++)if(informative(scores[c],covariances[c],n,i))present++;
            if(present>0) {
                int[] index=new int[present];for(int i=0,j=0;i<n;i++)if(informative(scores[c],covariances[c],n,i))index[j++]=i;
                double[] cu=new double[present],cv=new double[Math.multiplyExact(present,present)];
                for(int i=0;i<present;i++){cu[i]=scores[c][index[i]];for(int j=0;j<present;j++)cv[i*present+j]=covariances[c][index[i]*n+index[j]];}
                SummarySetTests.validate(new SetTestScoreState(cu,cv,present));
            }
        }
        int[] indices=new int[n]; int m=0;
        for(int i=0;i<n;i++) if(counts[i]>=minimumCohorts) indices[m++]=i;
        if(m==0) return null;
        indices=Arrays.copyOf(indices,m);
        double[] u=new double[m],v=new double[Math.multiplyExact(m,m)],uc=new double[m],vc=new double[v.length];
        StringBuilder[] directions=new StringBuilder[m];
        for(int i=0;i<m;i++) directions[i]=new StringBuilder(k);
        for(int c=0;c<k;c++) {
            for(int i=0;i<m;i++) {
                int a=indices[i]; boolean present=informative(scores[c],covariances[c],n,a);
                directions[i].append(!present?'?':scores[c][a]>0?'+':scores[c][a]<0?'-':'0');
                if(!present) continue;
                add(u,uc,i,scores[c][a]);
                for(int j=0;j<m;j++) {
                    int b=indices[j];
                    if(!informative(scores[c],covariances[c],n,b)) continue;
                    double value=covariances[c][a*n+b];
                    if(!Double.isFinite(value)) throw new IllegalArgumentException("missing covariance between informative variants");
                    add(v,vc,i*m+j,value);
                }
            }
        }
        for(int i=0;i<u.length;i++)u[i]+=uc[i];
        for(int i=0;i<v.length;i++)v[i]+=vc[i];
        int[] keptCounts=new int[m]; String[] strings=new String[m];
        for(int i=0;i<m;i++) { keptCounts[i]=counts[indices[i]]; strings[i]=directions[i].toString(); }
        SetTestScoreState state=new SetTestScoreState(u,v,m);
        SummarySetTests.validate(state);
        return new Pooled(state,indices,keptCounts,strings);
    }
    private static boolean informative(double[] u,double[] v,int n,int i) {
        return !Double.isNaN(u[i]) && v[i*n+i]>0;
    }
    private static void add(double[] sums,double[] corrections,int i,double value) {
        double old=sums[i],t=old+value;
        corrections[i]+=Math.abs(old)>=Math.abs(value)?(old-t)+value:(value-t)+old;sums[i]=t;
    }
    /** Immutable pooled state and input-column mapping. */
    public record Pooled(SetTestScoreState state,int[] indices,int[] cohortCounts,String[] directions) {
        public Pooled { indices=indices.clone(); cohortCounts=cohortCounts.clone(); directions=directions.clone(); }
        @Override public int[] indices(){return indices.clone();}
        @Override public int[] cohortCounts(){return cohortCounts.clone();}
        @Override public String[] directions(){return directions.clone();}
    }
}
