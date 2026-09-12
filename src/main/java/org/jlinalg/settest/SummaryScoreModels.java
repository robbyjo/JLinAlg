/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import java.util.*;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Conditional, heterogeneous and adaptive models of Gaussian score summaries. */
public final class SummaryScoreModels {
    private SummaryScoreModels() { }

    /** Extract an ordered principal substate; indices must be distinct. */
    public static SetTestScoreState subset(SetTestScoreState state, int[] indices) {
        int n=state.variants(), m=indices.length;
        if(m==0) throw new IllegalArgumentException("empty score subset");
        boolean[] used=new boolean[n]; double[] u=new double[m],v=new double[m*m];
        double[] sourceU=state.scores(),sourceV=state.information();
        for(int i=0;i<m;i++) {
            int a=indices[i];
            if(a<0||a>=n||used[a])throw new IllegalArgumentException("invalid or duplicate subset index");
            used[a]=true;u[i]=sourceU[a];
            for(int j=0;j<m;j++) {
                if(indices[j]<0||indices[j]>=n)throw new IllegalArgumentException("invalid subset index");
                v[i*m+j]=sourceV[a*n+indices[j]];
            }
        }
        return new SetTestScoreState(u,v,m);
    }

    /** Schur-complement adjustment in each cohort, BEFORE pooling. The joint
     * state must contain all target/condition cross-covariances. A singular or
     * ill-conditioned conditioning block is rejected (relative cutoff 1e-10).
     * This uses null-model information and is a one-step conditional analysis. */
    public static SetTestScoreState condition(SetTestScoreState joint, int[] targets, int[] conditioning) {
        SummarySetTests.validate(joint);
        if(conditioning.length==0)return subset(joint,targets);
        SetTestScoreState tt=subset(joint,targets), cc=subset(joint,conditioning);
        for(int a:targets)for(int b:conditioning)if(a==b)
            throw new IllegalArgumentException("target and conditioning variants overlap");
        int n=joint.variants(),t=targets.length,c=conditioning.length;
        double[] u=tt.scores(),v=tt.information(),full=joint.information(),cu=cc.scores();
        try(BackendContext context=BackendContext.select(BackendPolicy.CPU)) {
            var eig=context.backend().dsyev(cc.information(),c);
            double[] d=eig.eigenvalues(),e=eig.eigenvectors();
            double max=Arrays.stream(d).max().orElseThrow();
            for(double value:d)if(!(value>max*1e-10))
                throw new IllegalArgumentException("conditioning information is singular or ill-conditioned");
            for(int j=0;j<c;j++) {
                double score=0; double[] cross=new double[t];
                for(int a=0;a<c;a++) {
                    score+=e[a*c+j]*cu[a];
                    for(int i=0;i<t;i++)cross[i]+=full[targets[i]*n+conditioning[a]]*e[a*c+j];
                }
                for(int i=0;i<t;i++) {
                    u[i]-=cross[i]*score/d[j];
                    for(int h=0;h<t;h++)v[i*t+h]-=cross[i]*cross[h]/d[j];
                }
            }
        }
        SetTestScoreState result=new SetTestScoreState(u,v,t);
        SummarySetTests.validate(result);
        return result;
    }

    /** Independent cohort-by-variant effects: stack U and block-diagonal V.
     * Use repeated, fixed variant weights with SKAT or SKAT-O on this state.
     * SKAT-O combines the heterogeneous kernel with the pooled burden. */
    public static SetTestScoreState heterogeneous(List<SetTestScoreState> cohorts) {
        if(cohorts.isEmpty())throw new IllegalArgumentException("cohort states required");
        int n=0;for(var s:cohorts){SummarySetTests.validate(s);n=Math.addExact(n,s.variants());}
        double[] u=new double[n],v=new double[Math.multiplyExact(n,n)];int offset=0;
        for(var s:cohorts) {
            int m=s.variants();double[] su=s.scores(),sv=s.information();
            System.arraycopy(su,0,u,offset,m);
            for(int i=0;i<m;i++)System.arraycopy(sv,i*m,v,(i+offset)*n+offset,m);
            offset+=m;
        }
        return new SetTestScoreState(u,v,n);
    }

    /** VT reports the selected burden separately from its correlated-search p.
     * The selected burden is null when no threshold has positive information.
     * Monte Carlo SE describes simulation uncertainty, not biological uncertainty. */
    public record VariableThresholdResult(double selectedMaf, SetTestResult selectedBurden,
            double adjustedPValue, double monteCarloStandardError, int thresholds,
            int simulations, long seed) { }

    /** Search all distinct supplied MAFs using fixed linear weights. Calibrate
     * max |Z| under the joint Gaussian score null, preserving LD among thresholds. */
    public static VariableThresholdResult variableThreshold(String id, SetTestScoreState state,
            double[] weights, double[] maf, int simulations, long seed) {
        SummarySetTests.validate(state);int n=state.variants();
        if(weights.length!=n||maf.length!=n||simulations<1)
            throw new IllegalArgumentException("VT weights/MAFs/budget invalid");
        for(int i=0;i<n;i++)if(!(weights[i]>0)||!Double.isFinite(weights[i])||!(maf[i]>0&&maf[i]<=.5))
            throw new IllegalArgumentException("VT requires positive weights and MAFs in (0,.5]");
        Integer[] order=new Integer[n];for(int i=0;i<n;i++)order[i]=i;
        Arrays.sort(order,Comparator.comparingDouble(i->maf[i]));
        List<Integer> sizes=new ArrayList<>();List<Double> standardErrors=new ArrayList<>();
        SetTestResult selected=null;double chosen=Double.NaN,maxZ=-1,score=0,variance=0;
        double[] v=state.information(),u=state.scores();
        // Each threshold extends the same fixed weighted burden. Prefix updates
        // avoid a covariance decomposition for every nested threshold.
        for(int end=0;end<n;) {
            double threshold=maf[order[end]];
            do {
                int i=order[end];score+=weights[i]*u[i];variance+=weights[i]*weights[i]*v[i*n+i];
                for(int j=0;j<end;j++)variance+=2*weights[i]*weights[order[j]]*v[i*n+order[j]];
                end++;
            }while(end<n&&maf[order[end]]==threshold);
            if(!(variance>0))continue;
            var one=SummarySetTests.singleVariant(id,score,variance);
            var fit=new SetTestResult(id,"vt-selected-burden",n,end,one.statistic(),one.beta(),one.standardError(),
                Double.NaN,one.pValue(),one.log10PValue(),one.pValueMethod(),null,List.of());
            sizes.add(end);standardErrors.add(fit.standardError());
            if(Math.abs(fit.statistic())>maxZ){maxZ=Math.abs(fit.statistic());selected=fit;chosen=threshold;}
        }
        if(selected==null)return new VariableThresholdResult(Double.NaN,null,Double.NaN,Double.NaN,0,0,0);
        if(sizes.size()==1)return new VariableThresholdResult(chosen,selected,selected.pValue(),0,1,0,0);
        int extreme=0;
        try(BackendContext context=BackendContext.select(BackendPolicy.CPU)) {
            var eig=context.backend().dsyev(v,n);double[] d=eig.eigenvalues(),e=eig.eigenvectors();
            double[][] loading=new double[sizes.size()][n];double[] prefix=new double[n];int end=0;
            for(int h=0;h<loading.length;h++) {
                while(end<sizes.get(h)) {
                    int i=order[end++];for(int j=0;j<n;j++)prefix[j]+=weights[i]*e[i*n+j]*Math.sqrt(Math.max(0,d[j]));
                }
                for(int j=0;j<n;j++)loading[h][j]=prefix[j]*standardErrors.get(h);
            }
            Random random=new Random(seed);double[] z=new double[n];
            for(int s=0;s<simulations;s++) {
                for(int j=0;j<n;j++)z[j]=random.nextGaussian();
                for(double[] row:loading) {
                    double value=0;for(int j=0;j<n;j++)value+=row[j]*z[j];
                    if(Math.abs(value)>=maxZ){extreme++;break;}
                }
            }
        }
        double p=(extreme+1.0)/(simulations+1.0);
        return new VariableThresholdResult(chosen,selected,p,Math.sqrt(p*(1-p)/(simulations+1.0)),sizes.size(),simulations,seed);
    }
}
