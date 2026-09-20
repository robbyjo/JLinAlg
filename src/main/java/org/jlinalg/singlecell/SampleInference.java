/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.singlecell;

import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.*;

/** Biological-sample inference for independent donors or complete two-condition pairs. */
public final class SampleInference {
    private SampleInference() { }

    /** Design rows follow input order; the tested minus reference coefficient is column one. */
    public record Design(double[][] matrix, List<String> terms, int donors) { }

    /** Bound dense design allocation and QR work before allocating a sample design. */
    public static void checkWork(int rows,int columns) {
        if(rows<1 || columns<1 || columns>256 || (long)rows*columns>1_000_000
                || (long)rows*columns*columns>100_000_000)
            throw new IllegalArgumentException("Design exceeds 256 columns, one million entries or 100 million QR work units");
    }

    /** Numeric covariates must be explicitly encoded; paired designs use donor fixed effects. */
    public static Design design(String[] donors, String[] conditions, String reference,
            String tested, double[][] covariates, List<String> names, boolean paired) {
        int n=donors.length;
        if(n<4 || conditions.length!=n || covariates.length!=n || reference.equals(tested))
            throw new IllegalArgumentException("Require >=4 samples and distinct reference/tested conditions");
        Map<String,List<Integer>> groups=new LinkedHashMap<>();
        int controls=0,cases=0;
        for(int i=0;i<n;i++) {
            if(donors[i]==null || donors[i].isBlank())throw new IllegalArgumentException("Blank donor");
            if(conditions[i].equals(reference))controls++;
            else if(conditions[i].equals(tested))cases++;
            else throw new IllegalArgumentException("Select exactly the declared two conditions before inference");
            if(covariates[i].length!=names.size())throw new IllegalArgumentException("Covariate width mismatch");
            for(double value:covariates[i])if(!Double.isFinite(value))throw new IllegalArgumentException("Nonfinite covariate");
            groups.computeIfAbsent(donors[i],k->new ArrayList<>()).add(i);
        }
        if(controls<2 || cases<2)throw new IllegalArgumentException("Require >=2 independent donors per condition or >=2 complete pairs");
        for(List<Integer> rows:groups.values()) {
            if(paired) {
                if(rows.size()!=2 || conditions[rows.get(0)].equals(conditions[rows.get(1)]))
                    throw new IllegalArgumentException("Paired design requires exactly one sample per donor per condition");
            } else if(rows.size()!=1)throw new IllegalArgumentException("Repeated donor: use complete paired design; cells/sections are not replicates");
        }
        List<String> terms=new ArrayList<>(List.of("intercept","tested_minus_reference"));
        terms.addAll(names);
        List<String> ids=new ArrayList<>(groups.keySet());
        if(paired)for(int i=1;i<ids.size();i++)terms.add("donor:"+ids.get(i));
        if(n<=terms.size())throw new IllegalArgumentException("No residual degrees of freedom");
        checkWork(n,terms.size());
        if(names.stream().anyMatch(s->s==null||s.isBlank()) || new HashSet<>(terms).size()!=terms.size())
            throw new IllegalArgumentException("Design term names must be nonblank and unique, including generated terms");
        double[][] x=new double[n][terms.size()];
        for(int i=0;i<n;i++) {
            x[i][0]=1;x[i][1]=conditions[i].equals(tested)?1:0;
            System.arraycopy(covariates[i],0,x[i],2,names.size());
            if(paired)for(int j=1;j<ids.size();j++)x[i][1+names.size()+j]=donors[i].equals(ids.get(j))?1:0;
        }
        // Validate identifiability even when all supplied outcomes are constant.
        Ols.fit(new double[n],x,OlsOptions.defaults(),BackendPolicy.CPU);
        return new Design(x,List.copyOf(terms),groups.size());
    }

    /** Sample-level Gaussian model, with a finite-sample t test under iid Gaussian errors. */
    public static OlsResult fit(double[] response, Design design) {
        return Ols.fit(response,design.matrix(),OlsOptions.defaults(),BackendPolicy.CPU);
    }
}
