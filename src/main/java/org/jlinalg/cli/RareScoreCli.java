/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.pipeline.VcfVariantSource;
import org.jlinalg.raremetal.GaussianScoreWriter;
import org.jlinalg.settest.LinearSetTestNullModel;

/** Unrelated-sample quantitative-trait cohort summary export. */
final class RareScoreCli {
    private RareScoreCli(){ }
    static int run(String[] args,PrintStream console,PrintStream error) {
        try {
            if(Arrays.asList(args).contains("--help")){console.println(help());return 0;}
            Map<String,String> o=new HashMap<>();
            for(int i=0;i<args.length;i++) {
                String key=args[i];
                if(!Set.of("--vcf","--pheno","--id","--response","--covariates","--out","--genome-build","--cov-window","--max-variants","--grm","--trait-id","--trait-units").contains(key)
                    ||i+1==args.length||o.put(key,args[++i])!=null)throw new IllegalArgumentException("unknown, duplicate, or incomplete option: "+key);
            }
            for(String key:List.of("--vcf","--pheno","--id","--response","--out","--genome-build"))if(!o.containsKey(key))throw new IllegalArgumentException("required: "+key);
            Path prefix=Path.of(o.get("--out")).toAbsolutePath(),score=Path.of(prefix+".score.txt.gz"),cov=Path.of(prefix+".cov.txt.gz"),log=Path.of(prefix+".log");
            PipelinePaths.requireFreshOutputs(score,cov,Path.of(score+".tbi"),Path.of(cov+".tbi"),log);
            try(RunLog journal=RunLog.openPlain(log,false)) {
                VcfVariantSource source=VcfVariantSource.open(Path.of(o.get("--vcf")));
                DelimitedData table=DelimitedData.read(Path.of(o.get("--pheno")));
                int id=table.column(o.get("--id")),response=table.column(o.get("--response"));
                int[] covariates=Arrays.stream(o.getOrDefault("--covariates","").split(",")).filter(s->!s.isBlank()).mapToInt(table::column).toArray();
                Set<Integer> distinct=new HashSet<>();for(int c:covariates)if(c==response||!distinct.add(c))throw new IllegalArgumentException("covariates must be distinct from each other and response");
                Map<String,String[]> people=new HashMap<>();for(String[] row:table.rows())if(row[id].isBlank()||people.put(row[id],row)!=null)throw new IllegalArgumentException("phenotype IDs must be unique and nonblank");
                List<Integer> indices=new ArrayList<>();List<Double> y=new ArrayList<>();List<double[]> x=new ArrayList<>();
                List<String> sampleIds=source.metadata().sampleIds();
                for(int i=0;i<sampleIds.size();i++) {
                    String[] row=people.get(sampleIds.get(i));if(row==null)continue;
                    double outcome=number(row[response]);double[] design=new double[covariates.length+1];design[0]=1;boolean complete=Double.isFinite(outcome);
                    for(int j=0;j<covariates.length;j++){design[j+1]=number(row[covariates[j]]);complete &= Double.isFinite(design[j+1]);}
                    if(complete){indices.add(i);y.add(outcome);x.add(design);}
                }
                if(y.size()<=covariates.length+1)throw new IllegalArgumentException("insufficient complete-case samples");
                var model=o.containsKey("--grm")?null:LinearSetTestNullModel.prepare(y.stream().mapToDouble(Double::doubleValue).toArray(),x.toArray(double[][]::new),OlsOptions.defaults(),BackendPolicy.CPU);
                long window=Long.parseLong(o.getOrDefault("--cov-window","1000000"));int maximum=Integer.parseInt(o.getOrDefault("--max-variants","2000"));
                Path scratch=Files.createTempDirectory(prefix.getParent(),".rare-score-");
                try {
                    Path s=scratch.resolve("score.gz"),c=scratch.resolve("cov.gz");
                    if(o.containsKey("--grm")) {
                        var grm=GrmReader.read(Path.of(o.get("--grm"))).matrix();
                        List<String> aligned=indices.stream().map(sampleIds::get).toList();
                        var related=org.jlinalg.gwas.RemlAssociationScanner.prepare(
                            y.stream().mapToDouble(Double::doubleValue).toArray(),x.toArray(double[][]::new),
                            List.of(grm.varianceComponent("genetic",aligned),org.jlinalg.reml.VarianceComponent.identity("residual",y.size())),
                            org.jlinalg.reml.RemlOptions.defaults(),BackendPolicy.CPU);
                        GaussianScoreWriter.write(source,indices.stream().mapToInt(Integer::intValue).toArray(),related,window,maximum,o.get("--genome-build"),o.getOrDefault("--trait-id",o.get("--response")),o.getOrDefault("--trait-units","original"),s,c);
                        journal.info("related_Gaussian_REML_variance_components="+Arrays.toString(related.nullModel().varianceComponents()));
                    }else GaussianScoreWriter.write(source,indices.stream().mapToInt(Integer::intValue).toArray(),model,window,maximum,o.get("--genome-build"),o.getOrDefault("--trait-id",o.get("--response")),o.getOrDefault("--trait-units","original"),s,c);
                    journal.metadata("command=rare-score\nmodel=Gaussian; intercept and numeric covariates; optional GRM REML\ntrait_units=original\nvariance_estimator=RSS/(N-p) unrelated; REML with GRM\n"
                        +"analyzed_samples="+y.size()+"\ncovariance_window="+window+"\nmissing_dosages=mean imputed\nhwe=exact hard calls; fractional dosages NA\n"
                        +"model_metadata=version-1 quantitative; asymptotic-normal\ntrait_id="+o.getOrDefault("--trait-id",o.get("--response"))+"\ndeclared_trait_units="+o.getOrDefault("--trait-units","original")+"\n"
                        +"genome_build="+o.get("--genome-build")+"\nvcf="+o.get("--vcf")+"\npheno="+o.get("--pheno")+"\nresponse="+o.get("--response")+"\ncovariates="+o.getOrDefault("--covariates","")+"\ngrm="+o.getOrDefault("--grm","none")+"\n");
                    Files.move(s,score);Files.move(Path.of(s+".tbi"),Path.of(score+".tbi"));Files.move(c,cov);Files.move(Path.of(c+".tbi"),Path.of(cov+".tbi"));
                }finally {
                    try(var files=Files.list(scratch)){for(Path p:files.toList())Files.deleteIfExists(p);}Files.deleteIfExists(scratch);
                }
                journal.complete("complete");console.println("Cohort scores: "+score+"; covariance: "+cov);
            }
            return 0;
        }catch(IOException|RuntimeException e){error.println("jlinalg: "+e.getMessage());return 2;}
    }
    private static double number(String text){return text.isBlank()||Set.of("NA","NaN",".").contains(text)?Double.NaN:Double.parseDouble(text);}
    static String help(){return """
        Usage: jlinalg rare-score --vcf cohort.vcf.gz --pheno phenotype.tsv --id sample
          --response trait [--covariates age,sex,pc1] --genome-build GRCh38 --out PREFIX
          [--cov-window 1000000] [--max-variants 2000] [--grm grm.tsv]
          [--trait-id ID] [--trait-units UNITS]
        Fits one Gaussian null model with an intercept and numeric covariates.
        Scope: quantitative trait, diploid additive biallelic dosages; optional GRM REML for related samples.
        Produces RAREMETAL-compatible BGZF score/covariance files and tabix indices.
        Trait ID defaults to the response column; units default to original (no conversion).
        Exports version-1 quantitative model metadata for rare-meta --model-metadata strict.
        Scores use original phenotype units; variance is RSS/(N-p) or GRM-based REML.
        Complete phenotype/covariate cases are aligned by ID; missing genotypes are mean imputed.
        Covariance is saved for forward pairs up to cov-window base pairs apart.
        HWE is exact for hard calls and NA for fractional dosages; related-sample HWE is descriptive.
        """;}
}
