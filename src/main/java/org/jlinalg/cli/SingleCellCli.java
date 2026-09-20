/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.singlecell.SampleInference;
import org.jlinalg.network.NetworkAnalysis;
import org.jlinalg.ols.OlsResult;
import static org.jlinalg.cli.FollowupSupport.*;

/** Quantified, annotated RNA count workflows with biological-replicate inference. */
final class SingleCellCli {
    private SingleCellCli() { }
    static final Set<String> INPUT=Set.of("method","out","counts","cells","samples","features","min-features","min-counts","max-mito");
    static final Set<String> DESIGN=Set.of("reference","tested","paired","covariates");
    static int run(String[] args,PrintStream out,PrintStream err) {
        if(Arrays.asList(args).contains("--help")){out.println(help());return 0;}
        try {
            Set<String> allowed=new HashSet<>(INPUT);allowed.addAll(DESIGN);
            allowed.addAll(Set.of("min-cells","group-column","rscript","r-library","timeout","min-gene-count","min-gene-samples","gene-sets","seed","variable-features","components","clusters"));
            Map<String,String> o=options(args,allowed);String method=required(o,"method");
            Set<String> specific=new HashSet<>();
            switch(method) {
                case "qc" -> { }
                case "pseudobulk" -> specific.addAll(Set.of("min-cells","group-column"));
                case "state" -> {specific.addAll(DESIGN);specific.addAll(Set.of("min-cells","group-column","rscript","r-library","timeout","min-gene-count","min-gene-samples"));}
                case "abundance" -> {specific.addAll(DESIGN);specific.add("group-column");}
                case "pathways" -> {specific.addAll(DESIGN);specific.addAll(Set.of("min-cells","group-column","gene-sets"));}
                case "representation" -> specific.addAll(Set.of("rscript","r-library","timeout","seed","variable-features","components","clusters"));
                default -> throw new IllegalArgumentException("Unknown single-cell method: "+method);
            }
            for(String key:o.keySet())if(!INPUT.contains(key) && !specific.contains(key))throw new IllegalArgumentException("--"+key+" is not applicable to "+method);
            transaction(o,"single-cell-"+method,(dir,manifest)->execute(method,o,dir,manifest));
            out.println("single-cell "+method+" completed: "+o.get("out"));return 0;
        }catch(Exception failure){err.println("jlinalg: "+failure.getMessage());return 2;}
    }
    static String help(){return """
        single-cell --method qc|pseudobulk|state|abundance|pathways|representation
          --counts sparse.tsv --cells cells.tsv --samples samples.tsv --features features.tsv --out NEW_DIR
        Counts: obs_id,feature_id,count (TSV; absent entries are measured zeros).
        Cells: obs_id,sample_id,cell_type [exclude=true|false; other annotations retained].
        Samples: sample_id,donor_id,condition,batch,source,assay,visit [numeric covariates].
        Features: feature_id,mitochondrial (true|false); quantified raw RNA counts only.
        QC: --min-counts 1 --min-features 1 --max-mito 1 (fraction). No automatic doublet/ambient correction.
        Aggregation: --min-cells 10 --group-column cell_type (also supports supplied spatial domains).
        Inference: --reference CONTROL --tested CASE [--paired true --covariates age,sex_numeric].
        Exactly one sample/donor, or one sample/condition/donor for complete pairs; batch adjusted.
        State: limma voom adapter; --rscript PATH [--r-library DIR --timeout 3600]
          --min-gene-count 10 --min-gene-samples 2. Requires installed limma; no implicit installation.
        Abundance: sample-level arcsine-square-root proportions, equal donor weight; relative composition.
        Pathways: --gene-sets TSV (gene_set,feature_id), mean log2 CPM per supplied measured gene set.
        Representation: log1p(10000 counts/library), variance-ranked features, PCA and k-means.
          --rscript PATH [--seed 1 --variable-features 2000 --components 10 --clusters 5].
        Raw counts and exploratory representations remain separate. See single-cell vignette for limits.
        """;}
    private static void execute(String method,Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        CellData data=CellData.read(o);
        if(Set.of("state","abundance","pathways").contains(method))CellDesign.validate(data,o);
        if(method.equals("state")) {
            required(o,"rscript");integer(o,"timeout",3600,1,86400);integer(o,"min-gene-count",10,1,1000000000);integer(o,"min-gene-samples",2,1,10000);
            if(data.sampleIndex.containsKey("feature_id"))throw new IllegalArgumentException("sample_id feature_id conflicts with the wide count-matrix header");
        }
        Pathways pathways=method.equals("pathways")?readPathways(data,o):null;
        boolean[] keep=data.qc(o,dir,manifest);
        manifest.put("measurement","quantified raw RNA counts; supplied labels; omitted sparse entries are observed zeros");
        manifest.put("multiplicity","BH over all returned hypotheses, including failed/filtered rows as p=1 internally; such rows retain NA results. Populations with no retained observations are absent and listed in the input/QC audit.");
        if(method.equals("qc"))return;
        if(method.equals("representation")){representation(data,keep,o,dir,manifest);return;}
        String group=o.getOrDefault("group-column","cell_type");
        data.validateGroup(group);
        if(method.equals("abundance")){abundance(data,keep,group,o,dir,manifest);return;}
        int min=integer(o,"min-cells",10,1,100000);List<CellData.Aggregate> all=data.aggregate(keep,group);data.writeAggregates(all,dir,min);
        manifest.put("aggregation",Map.of("group_column",group,"min_cells",min,"unit","biological sample; all sections within a sample summed"));
        if(method.equals("pseudobulk"))return;
        LinkedHashSet<String> groups=new LinkedHashSet<>();all.forEach(a->groups.add(a.type()));
        CellData.checkOutputSize((long)groups.size()*(method.equals("state")?data.genes.size():pathways.groups.size()));
        if(pathways!=null) {
            CellData.checkOutputSize((long)all.size()*pathways.groups.size());
            long members=pathways.groups.values().stream().mapToLong(Set::size).sum();
            if(members*all.size()>100_000_000)throw new IllegalArgumentException("Pathway scoring exceeds 100 million aggregate-member visits");
        }
        List<String[]> statuses=new ArrayList<>(),results=new ArrayList<>();int counter=0;long fittingWork=0;
        for(String type:groups) {
            List<CellData.Aggregate> selected=all.stream().filter(a->a.type().equals(type)&&a.cells()>=min).toList();
            List<String> ids=selected.stream().map(CellData.Aggregate::sample).toList();
            SampleInference.Design design;
            try{design=CellDesign.create(data,ids,o);}catch(IllegalArgumentException unsupported){
                statuses.add(new String[]{type,"not_tested",unsupported.getMessage(),"NA"});
                List<String> hypotheses=method.equals("state")?data.genes:new ArrayList<>(pathways.groups.keySet());
                for(String feature:hypotheses)results.add(new String[]{type,feature,"NA","NA","NA","NA","NA","NA","NA","NA","unsupported_design_or_insufficient_samples"});
                continue;
            }
            fittingWork+=CellDesign.fitWork(design,method.equals("state")?data.genes.size():pathways.groups.size());
            if(fittingWork>100_000_000)throw new IllegalArgumentException("Inference exceeds 100 million response-design work units across populations");
            Path local=Files.createDirectory(dir.resolve("group-"+(++counter)));CellDesign.write(design,ids,local.resolve("design.tsv"));
            int start=results.size();
            if(method.equals("state"))state(data,selected,type,design,o,local,results);
            else pathways(selected,type,design,pathways,local,results);
            boolean tested=results.subList(start,results.size()).stream().anyMatch(r->r[10].equals("tested"));
            statuses.add(new String[]{type,tested?"tested":"not_tested",Integer.toString(design.donors())+" donors"+(tested?"":"; no eligible feature tests"),local.getFileName().toString()});
        }
        table(dir.resolve("population-status.tsv"),List.of("population","status","detail","output_directory"),statuses);
        adjustWrite(dir.resolve("results.tsv"),results);
        manifest.put("engine",method.equals("state")?"R limma voom with library-size normalization, empirical Bayes t inference":"Native sample-level OLS of prespecified mean log2 CPM scores");
        manifest.put("supported_design","independent donors or complete two-condition pairs; batch and supplied numeric covariates; no arbitrary longitudinal/random-effects model");
        manifest.put("reference",required(o,"reference"));manifest.put("tested",required(o,"tested"));
    }
    private static void state(CellData data,List<CellData.Aggregate> selected,String type,SampleInference.Design design,
            Map<String,String> o,Path dir,List<String[]> results) throws Exception {
        int count=integer(o,"min-gene-count",10,1,1000000000),samples=integer(o,"min-gene-samples",2,1,10000);
        try(var w=Files.newBufferedWriter(dir.resolve("counts.tsv"))) {
            List<String> h=new ArrayList<>(List.of("feature_id"));selected.forEach(a->h.add(a.sample()));row(w,h.toArray(String[]::new));
            for(int j=0;j<data.genes.size();j++){String[] r=new String[selected.size()+1];r[0]=data.genes.get(j);for(int i=0;i<selected.size();i++)r[i+1]=Double.toString(selected.get(i).values()[j]);row(w,r);}
        }
        adapter("state",o,dir,List.of(Integer.toString(count),Integer.toString(samples)));
        DelimitedData output=DelimitedData.read(dir.resolve("fit.tsv"));
        if(!output.header().equals(List.of("feature_id","effect","se","lower","upper","df","statistic","p","status"))
                || !output.rows().stream().map(r->r[0]).toList().equals(data.genes))throw new IOException("Adapter changed the result schema or feature identities/order");
        for(String[] r:output.rows())results.add(new String[]{type,r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],"NA",r[8]});
        if(!Files.isRegularFile(dir.resolve("session.txt")))throw new IOException("Adapter did not record its session");
    }
    static void adapter(String mode,Map<String,String> o,Path dir,List<String> args) throws Exception {
        Path script=dir.resolve("singlecell.R");
        try(InputStream stream=SingleCellCli.class.getResourceAsStream("/org/jlinalg/singlecell/singlecell.R")) {
            if(stream==null)throw new IOException("Missing bundled R adapter");Files.copy(stream,script);
        }
        List<String> command=new ArrayList<>(List.of(executable(required(o,"rscript")),"--vanilla",script.toAbsolutePath().toString(),mode,
            o.containsKey("r-library")?Path.of(o.get("r-library")).toAbsolutePath().toString():"-"));command.addAll(args);
        process(command,dir,integer(o,"timeout",3600,1,86400));
    }
    static void adjustWrite(Path path,List<String[]> results) throws IOException {
        double[] p=new double[results.size()];for(int i=0;i<p.length;i++)p[i]=results.get(i)[8].equals("NA")?1:finite(results.get(i)[8],"p value");
        double[] q=NetworkAnalysis.bh(p);for(int i=0;i<p.length;i++)results.get(i)[9]=results.get(i)[8].equals("NA")?"NA":Double.toString(q[i]);
        table(path,List.of("population","feature_id","effect","se","lower","upper","df","statistic","p","bh","status"),results);
    }
    static String[] fitRow(String group,String feature,double[] y,SampleInference.Design design) {
        OlsResult fit=SampleInference.fit(y,design);double se=fit.standardErrors()[1];
        double scale=0;for(double value:y)scale=Math.max(scale,Math.abs(value));
        double rounding=32*fit.parameters()*Math.ulp(scale);
        if(!(se>0) || !Double.isFinite(se) || Math.sqrt(fit.residualSumOfSquares()/y.length)<=rounding)
            return new String[]{group,feature,Double.toString(fit.coefficients()[1]),"NA","NA","NA",Integer.toString(fit.residualDegreesOfFreedom()),"NA","NA","NA","zero_or_numerically_zero_residual_variance"};
        return new String[]{group,feature,Double.toString(fit.coefficients()[1]),Double.toString(se),Double.toString(fit.confidenceLower()[1]),Double.toString(fit.confidenceUpper()[1]),Integer.toString(fit.residualDegreesOfFreedom()),Double.toString(fit.tStatistics()[1]),Double.toString(fit.pValues()[1]),"NA","tested"};
    }
    private static void abundance(CellData data,boolean[] keep,String group,Map<String,String> o,Path dir,Map<String,Object> manifest) throws IOException {
        Map<String,Map<String,Integer>> counts=new LinkedHashMap<>();Set<String> types=new TreeSet<>();
        for(int i=0;i<keep.length;i++)if(keep[i]){String t=data.cell(i,group);types.add(t);counts.computeIfAbsent(data.cell(i,"sample_id"),s->new HashMap<>()).merge(t,1,Integer::sum);}
        List<String> ids=new ArrayList<>(counts.keySet());SampleInference.Design design=CellDesign.create(data,ids,o);CellDesign.write(design,ids,dir.resolve("design.tsv"));
        CellData.checkOutputSize((long)types.size()*ids.size());
        CellDesign.fitWork(design,types.size());
        Map<String,Integer> totals=new HashMap<>();for(var e:counts.entrySet())totals.put(e.getKey(),e.getValue().values().stream().mapToInt(Integer::intValue).sum());
        List<String[]> rows=new ArrayList<>(),results=new ArrayList<>();
        for(String t:types) {
            double[] y=new double[ids.size()];
            for(int i=0;i<ids.size();i++) {
                int n=totals.get(ids.get(i)),k=counts.get(ids.get(i)).getOrDefault(t,0);
                double proportion=(double)k/n;y[i]=Math.asin(Math.sqrt(proportion));
                rows.add(new String[]{ids.get(i),t,Integer.toString(k),Integer.toString(n),Double.toString(proportion),Double.toString(y[i])});
            }
            results.add(fitRow(t,"relative_abundance",y,design));
        }
        table(dir.resolve("proportions.tsv"),List.of("sample_id","population","population_cells","total_retained_cells","proportion","asin_sqrt"),rows);adjustWrite(dir.resolve("results.tsv"),results);
        manifest.put("estimand","Tested minus reference arcsine-square-root relative cell proportion; all retained cell labels including unknown define the denominator; samples equally weighted, not pooled cells. Gaussian sample-level errors assumed.");
    }
    private record Pathways(Map<String,Set<Integer>> groups,List<String[]> mapping) { }
    private static Pathways readPathways(CellData data,Map<String,String> o) throws IOException {
        DelimitedData sets=CellData.bounded(Path.of(required(o,"gene-sets")));int sc=sets.column("gene_set"),fc=sets.column("feature_id");
        Map<String,Set<Integer>> groups=new LinkedHashMap<>();List<String[]> mapping=new ArrayList<>();Set<List<String>> seen=new HashSet<>();
        for(String[] r:sets.rows()) {
            if(r[sc].isBlank() || r[fc].isBlank() || r[sc].chars().anyMatch(Character::isISOControl) || r[fc].chars().anyMatch(Character::isISOControl) || !seen.add(List.of(r[sc],r[fc])))throw new IllegalArgumentException("Invalid/duplicate gene-set member");
            int j=data.featureIndex.getOrDefault(r[fc],-1);groups.computeIfAbsent(r[sc],s->new LinkedHashSet<>());if(j>=0)groups.get(r[sc]).add(j);
            mapping.add(new String[]{r[sc],r[fc],j>=0?"measured":"outside_panel"});
        }
        return new Pathways(groups,mapping);
    }
    private static void pathways(List<CellData.Aggregate> selected,String type,SampleInference.Design design,
            Pathways input,Path dir,List<String[]> results) throws IOException {
        List<String[]> scores=new ArrayList<>();
        for(var e:input.groups.entrySet()) {
            if(e.getValue().size()<2){results.add(new String[]{type,e.getKey(),"NA","NA","NA","NA","NA","NA","NA","NA","fewer_than_two_measured_members"});continue;}
            double[] y=new double[selected.size()];
            for(int i=0;i<y.length;i++) {
                CellData.Aggregate a=selected.get(i);
                for(int j:e.getValue())y[i]+=Math.log((a.values()[j]+.5)*1e6/(a.library()+1))/Math.log(2)/e.getValue().size();
                scores.add(new String[]{a.sample(),e.getKey(),Double.toString(y[i]),Integer.toString(e.getValue().size())});
            }
            results.add(fitRow(type,e.getKey(),y,design));
        }
        table(dir.resolve("gene-set-members.tsv"),List.of("gene_set","feature_id","status"),input.mapping);
        table(dir.resolve("scores.tsv"),List.of("sample_id","gene_set","mean_log2_cpm","measured_members"),scores);
    }
    private static void representation(CellData data,boolean[] keep,Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        int n=0;for(boolean b:keep)if(b)n++;
        if(n<3 || (long)n*data.genes.size()>2_000_000)throw new IllegalArgumentException("Exploratory representation requires >=3 cells and <=2 million dense entries");
        if(data.featureIndex.containsKey("obs_id"))throw new IllegalArgumentException("feature_id obs_id conflicts with the wide representation header");
        try(var w=Files.newBufferedWriter(dir.resolve("normalized.tsv"))) {
            List<String> header=new ArrayList<>(List.of("obs_id"));header.addAll(data.genes);row(w,header.toArray(String[]::new));
            for(int i=0;i<keep.length;i++)if(keep[i]){String[] r=new String[data.genes.size()+1];r[0]=data.ids.get(i);for(int j=0;j<data.genes.size();j++)r[j+1]=Double.toString(data.normalized(i,j));row(w,r);}
        }
        adapter("representation",o,dir,List.of(Integer.toString(integer(o,"seed",1,0,Integer.MAX_VALUE)),Integer.toString(integer(o,"variable-features",2000,2,50000)),Integer.toString(integer(o,"components",10,1,100)),Integer.toString(integer(o,"clusters",5,2,n-1))));
        manifest.put("representation","Exploratory log1p library-normalized counts; unscaled centered variance-selected PCA; seeded k-means. No batch integration, UMAP, reference annotation, or population inference from clusters.");
    }
}
