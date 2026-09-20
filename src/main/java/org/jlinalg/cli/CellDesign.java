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
    static SampleInference.Design create(CellData data,List<String> ids,Map<String,String> o) {
        boolean paired=bool(o,"paired",false);
        List<String> cov=o.containsKey("covariates")?List.of(o.get("covariates").split(",",-1)):List.of();
        if(new HashSet<>(cov).size()!=cov.size() || cov.contains("condition") || cov.contains("donor_id"))throw new IllegalArgumentException("Invalid/duplicate covariates");
        Set<String> sources=new HashSet<>(),assays=new HashSet<>();List<String> batches=new ArrayList<>();
        Map<String,Set<String>> donorBatches=new HashMap<>();
        for(String id:ids) {
            sources.add(data.sample(id,"source"));assays.add(data.sample(id,"assay"));String batch=data.sample(id,"batch");
            if(!batches.contains(batch))batches.add(batch);
            donorBatches.computeIfAbsent(data.sample(id,"donor_id"),key->new HashSet<>()).add(batch);
        }
        if(sources.size()!=1 || assays.size()!=1)throw new IllegalArgumentException("One source and assay per inference run required");
        // Donor effects already absorb batch effects constant within each complete pair.
        boolean addBatch=!paired || donorBatches.values().stream().anyMatch(b->b.size()>1);
        List<String> terms=new ArrayList<>(cov);
        if(addBatch)for(int b=1;b<batches.size();b++)terms.add("batch:"+batches.get(b));
        double[][] x=new double[ids.size()][terms.size()];String[] donor=new String[ids.size()],condition=new String[ids.size()];
        for(int i=0;i<ids.size();i++) {
            String id=ids.get(i);donor[i]=data.sample(id,"donor_id");condition[i]=data.sample(id,"condition");
            for(int j=0;j<cov.size();j++)x[i][j]=finite(data.sample(id,cov.get(j)),"numeric covariate "+cov.get(j));
            if(addBatch)for(int b=1;b<batches.size();b++)x[i][cov.size()+b-1]=data.sample(id,"batch").equals(batches.get(b))?1:0;
        }
        return SampleInference.design(donor,condition,required(o,"reference"),required(o,"tested"),x,terms,paired);
    }
    static void write(SampleInference.Design design,List<String> ids,Path file) throws IOException {
        List<String> header=new ArrayList<>(List.of("sample_id"));header.addAll(design.terms());List<String[]> rows=new ArrayList<>();
        for(int i=0;i<ids.size();i++){String[] row=new String[header.size()];row[0]=ids.get(i);for(int j=1;j<row.length;j++)row[j]=Double.toString(design.matrix()[i][j-1]);rows.add(row);}
        table(file,header,rows);
    }
}
