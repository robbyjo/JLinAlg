/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.util.*;
import java.nio.file.*;
import java.io.IOException;
import org.jlinalg.singlecell.SampleInference;
import static org.jlinalg.cli.FollowupSupport.*;

/** Common audited biological-sample design. */
final class CellDesign {
    private CellDesign() { }
    static long fitWork(SampleInference.Design design,long responses) {
        long work=responses*design.matrix().length*design.terms().size()*design.terms().size();
        if(work>100_000_000)throw new IllegalArgumentException("Inference exceeds 100 million response-design work units; select fewer hypotheses or simplify the design");
        return work;
    }
    static List<String> validate(CellData data,Map<String,String> o) {
        if(required(o,"reference").equals(required(o,"tested")))throw new IllegalArgumentException("Reference and tested conditions must differ");
        bool(o,"paired",false);
        List<String> cov=o.containsKey("covariates")?List.of(o.get("covariates").split(",",-1)):List.of();
        if(new HashSet<>(cov).size()!=cov.size() || cov.stream().anyMatch(s->s.isBlank()||Set.of("condition","donor_id","sample_id","batch","intercept","tested_minus_reference").contains(s)||s.startsWith("donor:")||s.startsWith("batch:")))
            throw new IllegalArgumentException("Invalid, reserved or duplicate covariate name");
        for(String name:cov){int column=data.samples.column(name);for(String[] row:data.samples.rows())finite(row[column],"numeric covariate "+name);}
        return cov;
    }
    static SampleInference.Design create(CellData data,List<String> ids,Map<String,String> o) {
        boolean paired=bool(o,"paired",false);
        List<String> cov=o.containsKey("covariates")?List.of(o.get("covariates").split(",",-1)):List.of();
        Set<String> sources=new HashSet<>(),assays=new HashSet<>();List<String> batches=new ArrayList<>();
        Map<String,Set<String>> donorBatches=new HashMap<>();
        for(String id:ids) {
            sources.add(data.sample(id,"source"));assays.add(data.sample(id,"assay"));String batch=data.sample(id,"batch");
            if(!batches.contains(batch))batches.add(batch);
            donorBatches.computeIfAbsent(data.sample(id,"donor_id"),key->new HashSet<>()).add(batch);
        }
        if(sources.size()!=1 || assays.size()!=1)throw new IllegalArgumentException("One source and assay per inference run required");
        // Select a batch basis modulo donor effects only. Never residualize against the
        // tested condition: that would silently remove genuine condition/batch confounding.
        List<String> batchTerms=new ArrayList<>();
        if(!paired)batchTerms.addAll(batches.subList(1,batches.size()));
        else if(donorBatches.values().stream().anyMatch(b->b.size()>1)) {
            SampleInference.checkWork(ids.size(),batches.size());
            Map<String,List<Integer>> donorRows=new LinkedHashMap<>();
            for(int i=0;i<ids.size();i++)donorRows.computeIfAbsent(data.sample(ids.get(i),"donor_id"),s->new ArrayList<>()).add(i);
            List<double[]> basis=new ArrayList<>();
            for(String batch:batches.subList(1,batches.size())) {
                double[] v=new double[ids.size()];for(int i=0;i<v.length;i++)v[i]=data.sample(ids.get(i),"batch").equals(batch)?1:0;
                for(var rows:donorRows.values()){double mean=0;for(int i:rows)mean+=v[i]/rows.size();for(int i:rows)v[i]-=mean;}
                double original=0;for(double a:v)original+=a*a;
                for(int pass=0;pass<2;pass++)for(double[] q:basis){double dot=0;for(int i=0;i<v.length;i++)dot+=q[i]*v[i];for(int i=0;i<v.length;i++)v[i]-=dot*q[i];}
                double norm=0;for(double a:v)norm+=a*a;
                if(norm>1e-20*Math.max(1,original)){norm=Math.sqrt(norm);for(int i=0;i<v.length;i++)v[i]/=norm;basis.add(v);batchTerms.add(batch);}
            }
        }
        List<String> terms=new ArrayList<>(cov);
        for(String batch:batchTerms)terms.add("batch:"+batch);
        SampleInference.checkWork(ids.size(),2+terms.size()+(paired?donorBatches.size()-1:0));
        double[][] x=new double[ids.size()][terms.size()];String[] donor=new String[ids.size()],condition=new String[ids.size()];
        for(int i=0;i<ids.size();i++) {
            String id=ids.get(i);donor[i]=data.sample(id,"donor_id");condition[i]=data.sample(id,"condition");
            for(int j=0;j<cov.size();j++)x[i][j]=finite(data.sample(id,cov.get(j)),"numeric covariate "+cov.get(j));
            for(int b=0;b<batchTerms.size();b++)x[i][cov.size()+b]=data.sample(id,"batch").equals(batchTerms.get(b))?1:0;
        }
        return SampleInference.design(donor,condition,required(o,"reference"),required(o,"tested"),x,terms,paired);
    }
    static void write(SampleInference.Design design,List<String> ids,Path file) throws IOException {
        List<String> header=new ArrayList<>(List.of("sample_id"));header.addAll(design.terms());List<String[]> rows=new ArrayList<>();
        for(int i=0;i<ids.size();i++){String[] row=new String[header.size()];row[0]=ids.get(i);for(int j=1;j<row.length;j++)row[j]=Double.toString(design.matrix()[i][j-1]);rows.add(row);}
        table(file,header,rows);
    }
}
