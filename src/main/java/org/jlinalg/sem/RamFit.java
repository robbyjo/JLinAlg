/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.ChiSquare;
import jdistlib.Normal;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.model.MissingDataPolicy;

/** Shared RAM distribution, pattern sufficient statistics, score and information. */
final class RamFit {
    record Distribution(double[] mean, double[] covariance, double[][] dMean, double[][] dCovariance) { }
    record Pattern(int[] indices, int count, double[] mean, double[] covariance) { }
    record State(SemModel model, double[] point, List<Pattern> patterns, double[][] data,
                 int observations, boolean missing) { }

    static double[] initial(SemModel model) {
        double[] x = new double[model.freeParameterCount()]; boolean[] set = new boolean[x.length];
        for (SemModel.Element e : model.elements()) if (!e.fixed()) {
            int j=model.freeIndex(e.label());
            if (!set[j]) { x[j]=model.varianceParameter(j)?Math.log(e.start()):e.start(); set[j]=true; }
        }
        return x;
    }

    static Distribution distribution(SemModel model, double[] x) {
        int q=model.allVariables().size(), p=model.variables().size(), k=x.length;
        double[] b=identity(q), s=new double[q*q], intercept=new double[q];
        for (SemModel.Element e:model.elements()) {
            double v=e.fixed()?e.start():model.varianceParameter(model.freeIndex(e.label()))
                ?Math.exp(x[model.freeIndex(e.label())]):x[model.freeIndex(e.label())];
            if (!Double.isFinite(v)) throw new IllegalArgumentException("nonfinite RAM parameter");
            switch(e.kind()) {
                case REGRESSION -> b[e.first()*q+e.second()]-=v;
                case VARIANCE -> s[e.first()*q+e.first()]=v;
                case COVARIANCE -> { s[e.first()*q+e.second()]=v; s[e.second()*q+e.first()]=v; }
                case INTERCEPT -> intercept[e.first()]=v;
            }
        }
        checkPositiveSemidefinite(s,q);
        double[] t=inverse(b,q), full=mm(mm(t,s,q),transpose(t,q),q), mu=mv(t,intercept);
        double[] mean=Arrays.copyOf(mu,p), cov=new double[p*p];
        for(int i=0;i<p;i++) for(int j=0;j<p;j++) cov[i*p+j]=full[i*q+j];
        chol(cov,p);
        double[][] dm=new double[k][p], dc=new double[k][p*p];
        for(SemModel.Element e:model.elements()) if(!e.fixed()) {
            int h=model.freeIndex(e.label()), a=e.first(), c=e.second();
            for(int i=0;i<p;i++) {
                if(e.kind()==SemModel.Kind.INTERCEPT) dm[h][i]+=t[i*q+a];
                if(e.kind()==SemModel.Kind.REGRESSION) dm[h][i]+=t[i*q+a]*mu[c];
                for(int j=0;j<p;j++) dc[h][i*p+j]+=switch(e.kind()) {
                    case REGRESSION -> t[i*q+a]*full[c*q+j]+full[i*q+c]*t[j*q+a];
                    case VARIANCE -> Math.exp(x[h])*t[i*q+a]*t[j*q+a];
                    case COVARIANCE -> t[i*q+a]*t[j*q+c]+t[i*q+c]*t[j*q+a];
                    case INTERCEPT -> 0;
                };
            }
        }
        return new Distribution(mean,cov,dm,dc);
    }

    static SemFitResult fit(double[][] input, SemModel model, SemOptions options,
            BackendPolicy policy, boolean missing) {
        if(input==null || model==null || options==null || policy==null || input.length<2)
            throw new IllegalArgumentException("SEM requires data, model, options and backend");
        SemBackendPolicy.requireSupported(policy);
        if(missing && !model.hasMeanStructure()) model=model.toBuilder().meanStructure().build();
        int p=model.variables().size(); List<double[]> rows=new ArrayList<>(); int[] counts=new int[p];
        for(double[] row:input) {
            if(row==null || row.length!=p) throw new IllegalArgumentException("SEM data dimensions do not match observed variables");
            int seen=0;
            for(double v:row) { if(Double.isInfinite(v)) throw new IllegalArgumentException("infinite SEM observation"); if(!Double.isNaN(v))seen++; }
            if(!missing && seen!=p) {
                if(options.missingDataPolicy()==MissingDataPolicy.ERROR) throw new IllegalArgumentException("missing SEM observation");
                continue;
            }
            if(seen==0)continue;
            rows.add(row.clone()); for(int j=0;j<p;j++) if(!Double.isNaN(row[j]))counts[j]++;
        }
        if(rows.size()<2)throw new IllegalArgumentException("too few informative SEM rows");
        for(int c:counts)if(c<2)throw new IllegalArgumentException("each observed variable needs at least two observations");
        double[][] data=rows.toArray(double[][]::new);
        if(!model.hasMeanStructure()) {
            double[] mean=new double[p]; for(double[] r:data)for(int j=0;j<p;j++)mean[j]+=r[j]/data.length;
            for(double[] r:data)for(int j=0;j<p;j++)r[j]-=mean[j];
        }
        List<Pattern> patterns=patterns(data);
        return finish(model,patterns,data,options,policy,missing);
    }

    static List<Pattern> patterns(double[][] data) {
        Map<String,List<double[]>> groups=new LinkedHashMap<>();
        for(double[] r:data) {
            StringBuilder key=new StringBuilder();for(double v:r)key.append(Double.isNaN(v)?'0':'1');
            groups.computeIfAbsent(key.toString(),unused->new ArrayList<>()).add(r);
        }
        List<Pattern> result=new ArrayList<>();
        for(List<double[]> group:groups.values()) {
            int[] idx=java.util.stream.IntStream.range(0,data[0].length).filter(j->!Double.isNaN(group.get(0)[j])).toArray();
            int p=idx.length,n=group.size();if(p==0)continue;
            double[] mean=new double[p],cov=new double[p*p];
            for(double[] r:group)for(int i=0;i<p;i++)mean[i]+=r[idx[i]]/n;
            for(double[] r:group)for(int i=0;i<p;i++)for(int j=0;j<p;j++)cov[i*p+j]+=(r[idx[i]]-mean[i])*(r[idx[j]]-mean[j])/n;
            result.add(new Pattern(idx,n,mean,cov));
        }
        return List.copyOf(result);
    }

    static SemFitResult fitMoments(double[] cov,double[] mean,int n,SemModel model,SemOptions options,BackendPolicy policy) {
        if(model==null || options==null || policy==null || cov==null || mean==null || n<2)
            throw new IllegalArgumentException("invalid SEM moments");
        SemBackendPolicy.requireSupported(policy);
        int p=model.variables().size();
        if(cov.length!=p*p || mean.length!=p)throw new IllegalArgumentException("invalid SEM moment dimensions");
        chol(cov,p);for(double v:mean)if(!Double.isFinite(v))throw new IllegalArgumentException("nonfinite mean");
        return finish(model,List.of(new Pattern(java.util.stream.IntStream.range(0,p).toArray(),n,
            model.hasMeanStructure()?mean.clone():new double[p],cov.clone())),null,options,policy,false);
    }

    static SemOptimizer.Value evaluate(SemModel model,double[] x,List<Pattern> patterns,int n) {
        try { return evaluate(distribution(model,x),patterns,n); }
        catch(IllegalArgumentException ex) { return new SemOptimizer.Value(Double.POSITIVE_INFINITY,null); }
    }
    static SemOptimizer.Value evaluate(Distribution d,List<Pattern> patterns,int n) {
        int p=d.mean.length,k=d.dMean.length; double value=0;double[] gradient=new double[k];
        for(Pattern g:patterns) {
            int m=g.indices.length; double[] cov=sub(d.covariance,g.indices,p), l=chol(cov,m), w=cholInverse(l,m);
            double[] delta=new double[m];for(int i=0;i<m;i++)delta[i]=g.mean[i]-d.mean[g.indices[i]];
            double[] wd=mv(w,delta), scatter=g.covariance.clone();
            for(int i=0;i<m;i++)for(int j=0;j<m;j++)scatter[i*m+j]+=delta[i]*delta[j];
            double trace=0,logdet=0;for(int i=0;i<m;i++) { logdet+=2*Math.log(l[i*m+i]);for(int j=0;j<m;j++)trace+=w[i*m+j]*scatter[j*m+i]; }
            double weight=(double)g.count/n;
            value+=.5*weight*(m*Math.log(2*Math.PI)+logdet+trace);
            double[] ws=mm(mm(w,scatter,m),w,m);
            for(int h=0;h<k;h++) {
                double score=0;
                for(int i=0;i<m;i++) {
                    score-=d.dMean[h][g.indices[i]]*wd[i];
                    for(int j=0;j<m;j++)score+=.5*(w[i*m+j]-ws[i*m+j])*d.dCovariance[h][g.indices[i]*p+g.indices[j]];
                }
                gradient[h]+=weight*score;
            }
        }
        return new SemOptimizer.Value(value,gradient);
    }

    static double[] expectedInformation(Distribution d,List<Pattern> patterns) {
        int k=d.dMean.length,p=d.mean.length;double[] info=new double[k*k];
        for(Pattern g:patterns) {
            int m=g.indices.length;double[] w=cholInverse(chol(sub(d.covariance,g.indices,p),m),m);
            double[][] v=new double[k][],mean=new double[k][];
            for(int h=0;h<k;h++) {
                v[h]=mm(w,sub(d.dCovariance[h],g.indices,p),m);mean[h]=new double[m];
                for(int i=0;i<m;i++)mean[h][i]=d.dMean[h][g.indices[i]];
            }
            for(int a=0;a<k;a++)for(int b=0;b<=a;b++) {
                double s=SemOptimizer.dot(mean[a],mv(w,mean[b]));
                for(int i=0;i<m;i++)for(int j=0;j<m;j++)s+=.5*v[a][i*m+j]*v[b][j*m+i];
                info[a*k+b]+=g.count*s;if(a!=b)info[b*k+a]+=g.count*s;
            }
        }
        return info;
    }
    static double[] hessian(SemOptimizer.Objective objective,double[] point,int n) {
        int k=point.length;double[] h=new double[k*k];
        for(int j=0;j<k;j++) {
            double step=1e-4*(1+Math.abs(point[j]));double[] a=point.clone(),b=point.clone();a[j]+=step;b[j]-=step;
            double[] ga=objective.at(a).gradient(),gb=objective.at(b).gradient();
            if(ga==null||gb==null)throw new IllegalArgumentException("information unavailable at parameter boundary");
            for(int i=0;i<k;i++)h[i*k+j]=n*(ga[i]-gb[i])/(2*step);
        }
        for(int i=0;i<k;i++)for(int j=0;j<i;j++)h[i*k+j]=h[j*k+i]=.5*(h[i*k+j]+h[j*k+i]);
        return h;
    }
    static double[][] scores(Distribution d,double[][] data) {
        double[][] scores=new double[data.length][];
        // Cache pattern precision once; case scores then require only quadratic products.
        int p=d.mean.length,k=d.dMean.length;Map<String,double[]> precisions=new LinkedHashMap<>();
        for(int r=0;r<data.length;r++) {
            double[] row=data[r];int[] idx=java.util.stream.IntStream.range(0,p).filter(j->!Double.isNaN(row[j])).toArray();
            int m=idx.length;String key=Arrays.toString(idx);
            double[] w=precisions.computeIfAbsent(key,z->cholInverse(chol(sub(d.covariance,idx,p),m),m));
            double[] delta=new double[m];for(int i=0;i<m;i++)delta[i]=row[idx[i]]-d.mean[idx[i]];
            double[] wd=mv(w,delta);scores[r]=new double[k];
            for(int h=0;h<k;h++)for(int i=0;i<m;i++) {
                scores[r][h]+=d.dMean[h][idx[i]]*wd[i];
                for(int j=0;j<m;j++)scores[r][h]+=.5*(wd[i]*wd[j]-w[i*m+j])*d.dCovariance[h][idx[i]*p+idx[j]];
            }
        }
        return scores;
    }

    private static SemFitResult finish(SemModel model,List<Pattern> patterns,double[][] data,
            SemOptions options,BackendPolicy policy,boolean missing) {
        return finish(model,patterns,data,options,policy,missing,true);
    }
    private static SemFitResult finish(SemModel model,List<Pattern> patterns,double[][] data,
            SemOptions options,BackendPolicy policy,boolean missing,boolean fitIndices) {
        int p=model.variables().size(),k=model.freeParameterCount(),n=patterns.stream().mapToInt(Pattern::count).sum();
        int df=p*(p+1)/2+(model.hasMeanStructure()?p:0)-k;
        if(df<0)throw new IllegalArgumentException("SEM has more free parameters than observed moments");
        double[] start=initial(model);
        SemOptimizer.Objective objective=x->evaluate(model,x,patterns,n);
        SemOptimizer.Optimum optimum=SemOptimizer.minimize(objective,start,options.maximumEvaluations(),options.tolerance());
        Distribution d=distribution(model,optimum.point());
        double[] rawCovariance=SemInformation.unavailable(k);
        if(identified(d,patterns))try {
            double[] information=missing?hessian(objective,optimum.point(),n):expectedInformation(d,patterns);
            rawCovariance=informationInverse(information,k);
        } catch(IllegalArgumentException unavailable) { /* Likelihood fitting does not require valid Wald inference. */ }
        double[] covariance=naturalCovariance(rawCovariance,model,optimum.point());
        boolean informationAvailable=SemInformation.finiteCovariance(covariance,k);
        if(!informationAvailable)df=-1; // Unknown inferential model dimension, not a label-count chi-square df.
        double ll=-n*objective.at(optimum.point()).value();
        double chi=Double.NaN,cfi=Double.NaN,tli=Double.NaN,rmsea=Double.NaN,srmr=Double.NaN;
        try {
        if(informationAvailable && patterns.size()==1 && patterns.get(0).indices.length==p) {
            Pattern g=patterns.get(0);double[] l=chol(g.covariance,p);double ld=0,bd=0;
            for(int i=0;i<p;i++){ld+=2*Math.log(l[i*p+i]);bd+=Math.log(g.covariance[i*p+i]);}
            double sat=-.5*n*(p*(Math.log(2*Math.PI)+1)+ld);
            chi=Math.max(0,2*(sat-ll));double baseline=n*(bd-ld),baselineDf=p*(p-1)/2.;
            double excess=Math.max(0,chi-df),den=Math.max(excess,baseline-baselineDf);
            cfi=den>0?1-excess/den:1;
            tli=df>0 && baseline>0?1-(chi/df-1)/(baseline/baselineDf-1):Double.NaN;
            rmsea=df>0?Math.sqrt(Math.max(0,(chi-df)/(df*(double)n))):0;
            srmr=srmr(g.covariance,g.mean,d,model.hasMeanStructure());
        }
        else if(informationAvailable && missing && fitIndices) {
            String[] names=model.variables().toArray(String[]::new);
            SemModel.Builder saturated=SemModel.builder(names).meanStructure(),baseline=SemModel.builder(names).meanStructure();
            for(int i=0;i<p;i++) {
                double mean=0,centered=0,count=0;
                // Chan's merge of within- and between-pattern centered moments.
                for(Pattern g:patterns)for(int j=0;j<g.indices.length;j++)if(g.indices[j]==i) {
                    double total=count+g.count,delta=g.mean[j]-mean;
                    centered+=g.count*g.covariance[j*g.indices.length+j]+delta*delta*count*g.count/total;
                    mean+=delta*g.count/total;count=total;
                }
                double variance=centered/count;
                saturated.intercept(names[i],mean).variance(names[i],variance);
                baseline.intercept(names[i],mean).variance(names[i],variance);
                for(int j=0;j<i;j++)saturated.covariance(names[i],names[j],0);
            }
            SemFitResult h1=finish(saturated.build(),patterns,null,options,policy,true,false);
            SemFitResult h0=finish(baseline.build(),patterns,null,options,policy,true,false);
            if(h1.converged() && h0.converged() && h1.informationAvailable()) {
                chi=Math.max(0,2*(h1.logLikelihood()-ll));double baselineChi=Math.max(0,2*(h1.logLikelihood()-h0.logLikelihood())),baselineDf=p*(p-1)/2.;
                double excess=Math.max(0,chi-df),den=Math.max(excess,baselineChi-baselineDf);
                cfi=den>0?1-excess/den:1;tli=df>0 && baselineChi>0?1-(chi/df-1)/(baselineChi/baselineDf-1):Double.NaN;
                rmsea=df>0?Math.sqrt(Math.max(0,(chi-df)/(df*(double)n))):0;
                srmr=srmr(h1.impliedCovariance(),h1.impliedMeans(),d,model.hasMeanStructure());
            }
        }
        } catch(IllegalArgumentException | IllegalStateException auxiliaryFailure) {
            // Singular H1 moments, invalid auxiliary starts or failed auxiliary
            // likelihoods must never discard the already fitted target model.
            chi=cfi=tli=rmsea=srmr=Double.NaN;
        }
        {
            SemFitResult result=new SemFitResult(estimates(model,optimum.point(),covariance),d.covariance,ll,chi,df,
                df>0?ChiSquare.cumulative(chi,df,false,false):Double.NaN,cfi,tli,rmsea,srmr,
                informationAvailable?-2*ll+2*k:Double.NaN,informationAvailable?-2*ll+Math.log(n)*k:Double.NaN,n,optimum.evaluations(),optimum.converged(),
                new BackendProvenance(policy,"CPU","Portable Java RAM pattern kernels",false,false));
            result.attach(d.mean,covariance,new State(model,optimum.point(),patterns,data,n,missing),optimum.scoreNorm());
            return result;
        }
    }
    private static double srmr(double[] sample,double[] mean,Distribution d,boolean means) {
        int p=mean.length;double sum=0;
        for(int i=0;i<p;i++) {
            for(int j=0;j<=i;j++){double residual=(sample[i*p+j]-d.covariance[i*p+j])/Math.sqrt(sample[i*p+i]*sample[j*p+j]);sum+=residual*residual;}
            if(means){double residual=(mean[i]-d.mean[i])/Math.sqrt(sample[i*p+i]);sum+=residual*residual;}
        }
        return Math.sqrt(sum/(p*(p+1)/2.+(means?p:0)));
    }
    private static boolean identified(Distribution d,List<Pattern> patterns) {
        int k=d.dMean.length,p=d.mean.length,rows=0;
        for(Pattern g:patterns){int m=g.indices.length;rows+=m+m*(m+1)/2;}
        double[][] columns=new double[k][rows];int row=0;
        for(Pattern g:patterns)for(int i:g.indices) {
            for(int h=0;h<k;h++)columns[h][row]=d.dMean[h][i]/Math.sqrt(d.covariance[i*p+i]);row++;
            for(int j:g.indices)if(j<=i) {
                double scale=Math.sqrt(d.covariance[i*p+i]*d.covariance[j*p+j]);
                for(int h=0;h<k;h++)columns[h][row]=d.dCovariance[h][i*p+j]/scale;row++;
            }
        }
        return SemInformation.fullColumnRank(columns);
    }
    static List<SemParameterEstimate> estimates(SemModel model,double[] x,double[] covariance) {
        List<SemParameterEstimate> result=new ArrayList<>();int k=x.length;
        for(int h=0;h<k;h++) {
            double estimate=model.varianceParameter(h)?Math.exp(x[h]):x[h],se=Math.sqrt(covariance[h*k+h]),z=estimate/se;
            result.add(new SemParameterEstimate(model.freeParameterLabels().get(h),estimate,se,z,
                2*Normal.cumulative(-Math.abs(z),0,1,true,false)));
        }
        return result;
    }
    static double[] naturalCovariance(double[] raw,SemModel model,double[] x) {
        double[] c=raw.clone();int k=x.length;
        for(int i=0;i<k;i++)for(int j=0;j<k;j++)c[i*k+j]*=(model.varianceParameter(i)?Math.exp(x[i]):1)*(model.varianceParameter(j)?Math.exp(x[j]):1);
        return c;
    }
    static double[] informationInverse(double[] info,int k) {
        if(k==0)return new double[0];
        try {
            double[] normalized=new double[k*k],scale=new double[k];
            for(int i=0;i<k;i++){scale[i]=Math.sqrt(info[i*k+i]);if(!(scale[i]>0))return SemInformation.unavailable(k);}
            for(int i=0;i<k;i++)for(int j=0;j<k;j++)normalized[i*k+j]=info[i*k+j]/(scale[i]*scale[j]);
            double[] l=chol(normalized,k);
            for(int i=0;i<k;i++)if(l[i*k+i]*l[i*k+i]<1e-10)return SemInformation.unavailable(k);
            double[] inverse=cholInverse(l,k);
            for(int i=0;i<k;i++)for(int j=0;j<k;j++)inverse[i*k+j]/=scale[i]*scale[j];
            return inverse;
        }
        catch(IllegalArgumentException e) { return SemInformation.unavailable(k); }
    }
    static double[] identity(int n) {double[] a=new double[n*n];for(int i=0;i<n;i++)a[i*n+i]=1;return a;}
    static double[] transpose(double[] a,int n){double[] b=new double[a.length];for(int i=0;i<n;i++)for(int j=0;j<n;j++)b[j*n+i]=a[i*n+j];return b;}
    static double[] mm(double[] a,double[] b,int n){double[] c=new double[n*n];for(int i=0;i<n;i++)for(int h=0;h<n;h++)for(int j=0;j<n;j++)c[i*n+j]+=a[i*n+h]*b[h*n+j];return c;}
    static double[] mv(double[] a,double[] b){int n=b.length;double[] c=new double[n];for(int i=0;i<n;i++)for(int j=0;j<n;j++)c[i]+=a[i*n+j]*b[j];return c;}
    static double[] sub(double[] a,int[] idx,int p){int n=idx.length;double[] b=new double[n*n];for(int i=0;i<n;i++)for(int j=0;j<n;j++)b[i*n+j]=a[idx[i]*p+idx[j]];return b;}
    static double[] chol(double[] a,int n) {
        double[] l=new double[n*n];
        for(int i=0;i<n;i++)for(int j=0;j<=i;j++) {
            if(!Double.isFinite(a[i*n+j]) || Math.abs(a[i*n+j]-a[j*n+i])>1e-9*(1+Math.abs(a[i*n+j])))throw new IllegalArgumentException("matrix must be finite and symmetric");
            double v=a[i*n+j];for(int h=0;h<j;h++)v-=l[i*n+h]*l[j*n+h];
            if(i==j){if(!(v>1e-13*Math.abs(a[i*n+i])))throw new IllegalArgumentException("matrix is not positive definite");l[i*n+j]=Math.sqrt(v);}
            else l[i*n+j]=v/l[j*n+j];
        }
        return l;
    }
    static double[] cholInverse(double[] l,int n) {
        double[] inv=new double[n*n];
        for(int c=0;c<n;c++) {
            double[] x=new double[n];
            for(int i=0;i<n;i++){double v=i==c?1:0;for(int j=0;j<i;j++)v-=l[i*n+j]*x[j];x[i]=v/l[i*n+i];}
            for(int i=n-1;i>=0;i--){double v=x[i];for(int j=i+1;j<n;j++)v-=l[j*n+i]*inv[j*n+c];inv[i*n+c]=v/l[i*n+i];}
        }
        return inv;
    }
    static double[] inverse(double[] a,int n) {
        double[] b=a.clone(),inv=identity(n);
        for(int j=0;j<n;j++) {
            int pivot=j;for(int i=j+1;i<n;i++)if(Math.abs(b[i*n+j])>Math.abs(b[pivot*n+j]))pivot=i;
            if(Math.abs(b[pivot*n+j])<1e-12)throw new IllegalArgumentException("singular RAM path matrix");
            for(int h=0;h<n;h++){double v=b[j*n+h];b[j*n+h]=b[pivot*n+h];b[pivot*n+h]=v;v=inv[j*n+h];inv[j*n+h]=inv[pivot*n+h];inv[pivot*n+h]=v;}
            double v=b[j*n+j];for(int h=0;h<n;h++){b[j*n+h]/=v;inv[j*n+h]/=v;}
            for(int i=0;i<n;i++)if(i!=j){v=b[i*n+j];for(int h=0;h<n;h++){b[i*n+h]-=v*b[j*n+h];inv[i*n+h]-=v*inv[j*n+h];}}
        }
        return inv;
    }
    private static void checkPositiveSemidefinite(double[] s,int n) {
        double[] a=s.clone();double scale=0;for(int i=0;i<n;i++)scale=Math.max(scale,Math.abs(s[i*n+i]));double tol=1e-12*Math.max(1,scale);
        for(int j=0;j<n;j++) {
            double d=a[j*n+j];if(d < -tol)throw new IllegalArgumentException("RAM disturbance covariance is indefinite");
            if(d<=0) {for(int i=j+1;i<n;i++)if(Math.abs(a[i*n+j])>tol)throw new IllegalArgumentException("invalid zero disturbance variance");continue;}
            for(int i=j+1;i<n;i++)for(int h=j+1;h<n;h++)a[i*n+h]-=a[i*n+j]*a[j*n+h]/d;
        }
    }
}
