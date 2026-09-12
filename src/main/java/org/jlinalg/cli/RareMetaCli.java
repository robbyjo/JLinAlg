/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.raremetal.RareMetalStudy;
import org.jlinalg.raremetal.RareMetalStudy.Score;
import org.jlinalg.settest.*;

/** Cohort score-summary meta-analysis; no genome-wide dense covariance matrix. */
final class RareMetaCli {
    private RareMetaCli() { }
    static int run(String[] args,PrintStream console,PrintStream error) {
        try {
            if(Arrays.asList(args).contains("--help")){console.println(help());return 0;}
            Options o=new Options(args);
            execute(o,console);return 0;
        }catch(IOException|RuntimeException e){error.println("jlinalg: "+e.getMessage());return 2;}
    }
    private static void execute(Options o,PrintStream console)throws IOException {
        List<RareMetalStudy> studies=new ArrayList<>();
        Map<String,BufferedWriter> writers=new LinkedHashMap<>();
        Map<String,Path> temporary=new LinkedHashMap<>(),outputs=new LinkedHashMap<>();
        Path log=Path.of(o.output+".log");
        for(String test:o.tests)outputs.put(test,Path.of(o.output+"."+test+".tsv"));
        List<Path> destinations=new ArrayList<>(outputs.values());destinations.add(log);
        PipelinePaths.requireFreshOutputs(destinations.toArray(Path[]::new));
        try(RunLog journal=RunLog.openPlain(log,false)) {
            try {
                for(int i=0;i<o.scoreFiles.size();i++) {
                    for(int j=0;j<i;j++) if(Files.isSameFile(o.scoreFiles.get(i),o.scoreFiles.get(j)))
                        throw new IllegalArgumentException("one score file cannot represent two independent cohorts");
                    studies.add(new RareMetalStudy(o.scoreFiles.get(i),o.covFiles.get(i),o.callRate,o.hwe));
                    String declared=studies.get(i).genomeBuild();
                    if(declared!=null&&!declared.equals(o.build))throw new IOException("cohort genome build differs: "+o.names.get(i));
                }
                journal.metadata("command=rare-meta\ngenome_build="+o.build+"\ncohort_order="+String.join(",",o.names)
                    +"\ntests="+String.join(",",o.tests)+"\nweights="+o.weights+"\nmaf="+o.maf+"\naf_policy="+o.afPolicy
                    +"\nmin_cohorts_per_variant="+o.minimum+"\ncohort_results="+o.cohortResults
                    +"\nskato_calibration="+o.calibration+"\nsimulations="+o.simulations+"\nseed="+o.seed
                    +"\nmodel=independent cohorts; quantitative-trait null scores; compatible trait units required\n");
                for(int i=0;i<studies.size();i++)journal.info("cohort="+o.names.get(i)+" scores="+o.scoreFiles.get(i)
                    +" covariance="+o.covFiles.get(i)+" samples="+studies.get(i).samples()+" indexed="+studies.get(i).indexed());
                for(String test:o.tests) {
                    Path tmp=Files.createTempFile(o.output.getParent(),".rare-meta-",".tsv");temporary.put(test,tmp);
                    BufferedWriter w=Files.newBufferedWriter(tmp);writers.put(test,w);w.write(header(test));w.newLine();
                }
                if(o.tests.contains("single")) singles(o,studies,writers.get("single"));
                if(o.tests.stream().anyMatch(t->!t.equals("single"))) {
                    if(o.groups!=null) {
                        try(BufferedReader r=RareMetalStudy.open(o.groups)) {
                            Set<String> ids=new HashSet<>();
                            for(String line;(line=r.readLine())!=null;) {
                                if(line.isBlank()||line.startsWith("#"))continue;
                                String[] f=line.trim().split("\\s+");
                                if(f.length<2||!ids.add(f[0]))throw new IOException("group IDs must be unique and have variants");
                                List<Variant> variants=new ArrayList<>();Set<String> unique=new HashSet<>();
                                for(int i=1;i<f.length;i++){Variant v=Variant.parse(f[i]);if(!unique.add(v.chromosome+":"+v.position))throw new IOException("duplicate/multiallelic group position");variants.add(v);}
                                group(o,studies,writers,f[0],variants);
                            }
                        }
                    }else windows(o,studies,writers);
                }
                for(BufferedWriter w:writers.values())w.close();writers.clear();
                for(String test:o.tests)Files.move(temporary.get(test),outputs.get(test));
                journal.complete("complete");console.println("Rare-variant meta-analysis: "+outputs.values());
            }finally {
                for(BufferedWriter w:writers.values())w.close();
                for(RareMetalStudy s:studies)s.close();
                for(Path p:temporary.values())Files.deleteIfExists(p);
            }
        }
    }
    private static String header(String test) {
        String common="scope\tfeature_id\tn_variants\tn_cohorts\tstatus";
        return common+(test.equals("single")?"\tn_samples\talt_af\tdirection\tscore\tvariance\tbeta\tse\tz\tp_value\tneg_log10_p":
            test.equals("burden")?"\tweights\tdirection\tbeta\tse\tz\tp_value\tneg_log10_p\tcalibration":
            test.equals("skat")?"\tweights\tq\tp_value\tneg_log10_p\tcalibration":
            "\tweights\tminimum_component_p\tadjusted_p\tneg_log10_p\tcalibration\tsimulations\tseed\tcomponent_rho_p");
    }
    private static void singles(Options o,List<RareMetalStudy> studies,BufferedWriter writer)throws IOException {
        try(Union union=new Union(studies)) {
            for(Score[] row;(row=union.next())!=null;) {
                Score anchor=Arrays.stream(row).filter(s->s!=null&&!s.alternate().isEmpty()).findFirst()
                    .orElse(Arrays.stream(row).filter(Objects::nonNull).findFirst().orElseThrow());
                double u=0,uc=0,v=0,n=0,alleles=0,afN=0;int count=0;StringBuilder direction=new StringBuilder();
                for(int c=0;c<row.length;c++) {
                    Score s=row[c];int sign=s==null||s.status().equals("invalid_alleles")?1:s.orientation(anchor.reference(),anchor.alternate());
                    if(sign==0)throw new IOException("incompatible alleles at "+anchor.id()+" in "+o.names.get(c));
                    if(s!=null && !s.status().equals("qc_excluded") && Double.isFinite(s.frequency())){afN+=s.samples();alleles+=s.samples()*(sign==1?s.frequency():1-s.frequency());}
                    else if(s==null&&o.afPolicy.equals("raremetal"))afN+=studies.get(c).samples();
                    if(s!=null&&s.informative()){double term=sign*s.score(),t=u+term;uc+=Math.abs(u)>=Math.abs(term)?(u-t)+term:(term-t)+u;u=t;v+=s.variance();n+=s.samples();count++;direction.append(sign*s.score()>0?'+':sign*s.score()<0?'-':'0');}
                    else direction.append('?');
                    if(o.cohortResults)single(writer,o.names.get(c),anchor.id(),s!=null&&s.informative()?1:0,s==null?"missing":s.status(),
                        s==null?0:s.samples(),s==null?Double.NaN:sign==1?s.frequency():1-s.frequency(),""+direction.charAt(c),s==null?Double.NaN:sign*s.score(),s==null?Double.NaN:s.variance());
                }
                single(writer,"meta",anchor.id(),count,count<o.minimum?"below_min_cohorts":count==1?"single_cohort":"ok",n,
                    afN>0?alleles/afN:Double.NaN,direction.toString(),count<o.minimum?Double.NaN:u+uc,count<o.minimum?Double.NaN:v);
            }
        }
    }
    private static void single(BufferedWriter w,String scope,String id,int count,String status,double n,double af,String direction,double u,double v)throws IOException {
        SetTestResult fit=Double.isFinite(u)&&v>0?SummarySetTests.singleVariant(id,u,v):null;
        line(w,scope,id,1,count,status,n,af,direction,u,v,fit==null?Double.NaN:fit.beta(),fit==null?Double.NaN:fit.standardError(),
            fit==null?Double.NaN:fit.statistic(),fit==null?Double.NaN:fit.pValue(),fit==null?Double.NaN:fit.negativeLog10PValue());
    }
    private static void group(Options o,List<RareMetalStudy> studies,Map<String,BufferedWriter> writers,String id,List<Variant> variants)throws IOException {
        if(variants.size()>o.maxVariants)throw new IOException("group exceeds --max-variants: "+id);
        String chr=variants.get(0).chromosome;long low=Long.MAX_VALUE,high=0;
        for(Variant v:variants){if(!chr.equals(v.chromosome))throw new IOException("groups must lie on one chromosome");low=Math.min(low,v.position);high=Math.max(high,v.position);}
        int k=studies.size(),n=variants.size(); Score[][] values=new Score[k][n];int[][] signs=new int[k][n];
        double[] af=new double[n]; int[] counts=new int[n];boolean[] retained=new boolean[n];
        for(int c=0;c<k;c++) {
            Map<Long,Score> region=studies.get(c).region(chr,low,high);
            for(int i=0;i<n;i++) {
                Variant target=variants.get(i);Score s=region.get(target.position);values[c][i]=s;
                signs[c][i]=s==null||s.status().equals("invalid_alleles")?1:s.orientation(target.reference,target.alternate);
                if(signs[c][i]==0)throw new IOException("incompatible alleles for "+target.id()+" in "+o.names.get(c));
            }
        }
        for(int i=0;i<n;i++) {
            double a=0,denom=0;
            for(int c=0;c<k;c++) {
                Score s=values[c][i];
                if(s!=null&&s.informative())counts[i]++;
                if(s!=null&&!s.status().equals("qc_excluded")&&Double.isFinite(s.frequency())){denom+=s.samples();a+=s.samples()*(signs[c][i]==1?s.frequency():1-s.frequency());}
                else if(s==null&&o.afPolicy.equals("raremetal"))denom+=studies.get(c).samples();
            }
            af[i]=denom>0?a/denom:Double.NaN;
            if(af[i]>.5){af[i]=1-af[i];for(int c=0;c<k;c++)signs[c][i]*=-1;}
            retained[i]=counts[i]>=o.minimum&&af[i]>0&&af[i]<=o.maf;
        }
        int m=0;for(boolean b:retained)if(b)m++;
        if(m==0){for(String test:o.tests)if(!test.equals("single"))empty(writers.get(test),test,"meta",id,0,0,"no_variants_after_filters",o.weight(test));return;}
        double[][] u=new double[k][m],cov=new double[k][Math.multiplyExact(m,m)];double[] mafs=new double[m];
        for(int c=0;c<k;c++) {
            Arrays.fill(u[c],Double.NaN);List<Score> present=new ArrayList<>();List<Integer> positions=new ArrayList<>();List<Integer> directions=new ArrayList<>();
            for(int i=0,j=0;i<n;i++)if(retained[i]) {
                mafs[j]=af[i];Score s=values[c][i];
                if(s!=null&&s.informative()){u[c][j]=signs[c][i]*s.score();present.add(s);positions.add(j);directions.add(signs[c][i]);}j++;
            }
            double[] block=studies.get(c).covariance(present);
            for(int i=0;i<present.size();i++)for(int j=0;j<present.size();j++)
                cov[c][positions.get(i)*m+positions.get(j)]=directions.get(i)*directions.get(j)*block[i*present.size()+j];
            if(!present.isEmpty()) {
                double[] cu=new double[present.size()],cv=new double[block.length];
                for(int i=0;i<cu.length;i++){cu[i]=u[c][positions.get(i)];for(int j=0;j<cu.length;j++)cv[i*cu.length+j]=cov[c][positions.get(i)*m+positions.get(j)];}
                SummarySetTests.validate(new SetTestScoreState(cu,cv,cu.length));
            }
        }
        var pooled=ScoreMetaAnalysis.pool(u,cov,1);int totalSamples=studies.stream().mapToInt(RareMetalStudy::samples).sum();
        for(String test:o.tests)if(!test.equals("single")) {
            double[] weights=weights(o.weight(test),mafs,totalSamples);int count=0;StringBuilder direction=new StringBuilder();
            for(int c=0;c<k;c++) {
                boolean any=false;double score=0;int cm=0;
                for(int i=0;i<m;i++)if(Double.isFinite(u[c][i])){any=true;score+=weights[i]*u[c][i];cm++;}
                if(any)count++;direction.append(!any?'?':score>0?'+':score<0?'-':'0');
                if(o.cohortResults) {
                    if(!any)empty(writers.get(test),test,o.names.get(c),id,0,0,"no_informative_variants",o.weight(test));
                    else {
                        var cp=ScoreMetaAnalysis.pool(new double[][]{u[c]},new double[][]{cov[c]},1);
                        int[] ix=cp.indices();double[] cw=new double[ix.length];for(int i=0;i<ix.length;i++)cw[i]=weights[ix[i]];
                        result(writers.get(test),test,o.names.get(c),id,cp.state(),cw,1,"ok",""+direction.charAt(c),o);
                    }
                }
            }
            result(writers.get(test),test,"meta",id,pooled.state(),weights,count,count==1?"single_cohort":"ok",direction.toString(),o);
        }
    }
    private static double[] weights(String method,double[] maf,int samples) {
        double[] weights=new double[maf.length];
        for(int i=0;i<weights.length;i++)weights[i]=switch(method){
            case "equal"->1;case "mb"->1/Math.sqrt(maf[i]*(1-maf[i]));
            case "beta"->25*Math.pow(1-maf[i],24);
            case "raremetal-beta"->25*Math.pow(1-Math.min(1-2.0/samples,Math.max(2.0/samples,maf[i])),24);
            default->throw new IllegalArgumentException("unknown weights: "+method);};
        return weights;
    }
    private static void result(BufferedWriter w,String test,String scope,String id,SetTestScoreState state,double[] weights,int cohorts,String status,String direction,Options o)throws IOException {
        if(test.equals("skat-o")) {
            SetTestOptions defaults=SetTestOptions.defaults();
            var options=new SetTestOptions(defaults.variantFilter(),defaults.missingPolicy(),defaults.skatORhoGrid(),o.simulations,o.seed,o.calibration);
            SkatOResult r=SummarySetTests.skatO(id,state,weights,options);
            String components=String.join(";",r.components().stream().map(c->c.rho()+":"+c.result().pValue()).toList());
            line(w,scope,id,state.variants(),cohorts,status,o.weight(test),r.minimumComponentPValue(),r.adjustedPValue(),r.negativeLog10AdjustedPValue(),
                r.components().isEmpty()?"rank_one_chi_square":o.calibration==SkatOCalibration.ANALYTIC?"analytic_moment_approximation":"gaussian_score_simulation",r.simulations(),r.randomSeed(),components);
        }else {
            SetTestResult r=test.equals("burden")?SummarySetTests.burden(id,state,weights):SummarySetTests.skat(id,state,weights);
            if(test.equals("burden"))line(w,scope,id,state.variants(),cohorts,Double.isNaN(r.beta())?"no_burden_information":status,o.weight(test),direction,r.beta(),r.standardError(),r.statistic(),r.pValue(),r.negativeLog10PValue(),r.pValueMethod());
            else line(w,scope,id,state.variants(),cohorts,status,o.weight(test),r.statistic(),r.pValue(),r.negativeLog10PValue(),r.pValueMethod());
        }
    }
    private static void empty(BufferedWriter w,String test,String scope,String id,int variants,int cohorts,String status,String weight)throws IOException {
        int n=header(test).split("\t").length;Object[] values=new Object[n];Arrays.fill(values,"NA");
        values[0]=scope;values[1]=id;values[2]=variants;values[3]=cohorts;values[4]=status;values[5]=weight;line(w,values);
    }
    private static void line(BufferedWriter w,Object... values)throws IOException {
        for(int i=0;i<values.length;i++){if(i>0)w.write('\t');Object v=values[i];w.write(v instanceof Double d && Double.isNaN(d)?"NA":v.toString());}w.newLine();
    }
    private static void windows(Options o,List<RareMetalStudy> studies,Map<String,BufferedWriter> writers)throws IOException {
        String chromosome="";long last=-1;
        try(Union union=new Union(studies)) {
            for(Score[] row;(row=union.next())!=null;) {
                Score s=Arrays.stream(row).filter(Objects::nonNull).findFirst().orElseThrow();
                if(!chromosome.equals(s.chromosome())){chromosome=s.chromosome();last=-1;}
                long first=Math.max(0,Math.floorDiv(s.position()-o.windowSize-1,o.windowStep)+1),end=(s.position()-1)/o.windowStep;
                for(long index=Math.max(first,last+1);index<=end;index++) {
                    long start=index*o.windowStep+1,stop=start+o.windowSize-1;
                    Map<Long,Variant> members=new TreeMap<>();
                    for(var study:studies)for(Score v:study.region(chromosome,start,stop).values()) {
                        // A monomorphic first cohort may not name the second allele.
                        // Choose a complete pair from any cohort before harmonizing.
                        if(!v.alternate().isEmpty())members.putIfAbsent(v.position(),new Variant(v.chromosome(),v.position(),v.reference(),v.alternate()));
                    }
                    if(!members.isEmpty())group(o,studies,writers,chromosome+":"+start+"-"+stop,new ArrayList<>(members.values()));
                    last=index;
                }
            }
        }
    }
    private record Variant(String chromosome,long position,String reference,String alternate) {
        static Variant parse(String id){String[] f=id.split(":");if(f.length!=4)throw new IllegalArgumentException("variant ID must be CHROM:POS:REF:ALT");long p=Long.parseLong(f[1]);if(p<1)throw new IllegalArgumentException("variant position must be positive");return new Variant(RareMetalStudy.chromosome(f[0]),p,f[2],f[3]);}
        String id(){return chromosome+":"+position+":"+reference+":"+alternate;}
    }
    private static final class Union implements AutoCloseable {
        private final List<RareMetalStudy.Cursor> cursors=new ArrayList<>();private final Score[] current;
        Union(List<RareMetalStudy> studies)throws IOException {
            current=new Score[studies.size()];try{for(int i=0;i<studies.size();i++){cursors.add(studies.get(i).cursor());current[i]=cursors.get(i).next();}}catch(IOException|RuntimeException e){close();throw e;}
        }
        Score[] next()throws IOException {
            Score smallest=null;for(Score s:current)if(s!=null&&(smallest==null||RareMetalStudy.comparePosition(s,smallest)<0))smallest=s;
            if(smallest==null)return null;Score[] row=new Score[current.length];
            for(int i=0;i<current.length;i++)if(current[i]!=null&&RareMetalStudy.comparePosition(current[i],smallest)==0){row[i]=current[i];current[i]=cursors.get(i).next();}return row;
        }
        public void close()throws IOException{for(var c:cursors)c.close();}
    }
    private static final class Options {
        final List<String> names=new ArrayList<>(),tests=new ArrayList<>();final List<Path> scoreFiles=new ArrayList<>(),covFiles=new ArrayList<>();
        Path output,groups;String build,weights="default",afPolicy="observed";int minimum=1,maxVariants=2000,simulations=10000;long seed=20260901L,windowSize,windowStep;double maf=.05,callRate=0,hwe=0;boolean cohortResults;SkatOCalibration calibration=SkatOCalibration.PARAMETRIC_SIMULATION;
        Options(String[] args)throws IOException {
            Map<String,String> map=new HashMap<>();
            for(int i=0;i<args.length;i++){String key=args[i];if(key.equals("--cohort-results")){cohortResults=true;continue;}if(!Set.of("--cohorts","--out","--test","--groups","--genome-build","--weights","--af-policy","--min-cohorts","--max-variants","--simulations","--seed","--window-size","--window-step","--maf","--call-rate","--hwe","--skato-calibration").contains(key)||i+1==args.length||map.put(key,args[++i])!=null)throw new IllegalArgumentException("unknown, duplicate, or incomplete option: "+key);}
            for(String key:List.of("--cohorts","--out","--test","--genome-build"))if(!map.containsKey(key))throw new IllegalArgumentException("required option: "+key);
            output=Path.of(map.get("--out")).toAbsolutePath();build=map.get("--genome-build");
            for(String t:map.get("--test").split(","))if(!Set.of("single","burden","skat","skat-o").contains(t)||tests.contains(t))throw new IllegalArgumentException("invalid or repeated test: "+t);else tests.add(t);
            if(map.containsKey("--groups"))groups=Path.of(map.get("--groups"));
            weights=map.getOrDefault("--weights",weights);if(!Set.of("default","equal","mb","beta","raremetal-beta").contains(weights))throw new IllegalArgumentException("invalid weights");
            afPolicy=map.getOrDefault("--af-policy",afPolicy);if(!Set.of("observed","raremetal").contains(afPolicy))throw new IllegalArgumentException("invalid AF policy");
            minimum=Integer.parseInt(map.getOrDefault("--min-cohorts","1"));maxVariants=Integer.parseInt(map.getOrDefault("--max-variants","2000"));simulations=Integer.parseInt(map.getOrDefault("--simulations","10000"));seed=Long.parseLong(map.getOrDefault("--seed","20260901"));
            windowSize=Long.parseLong(map.getOrDefault("--window-size","0"));windowStep=Long.parseLong(map.getOrDefault("--window-step",Long.toString(windowSize)));
            maf=Double.parseDouble(map.getOrDefault("--maf","0.05"));callRate=Double.parseDouble(map.getOrDefault("--call-rate","0"));hwe=Double.parseDouble(map.getOrDefault("--hwe","0"));
            String c=map.getOrDefault("--skato-calibration","simulation");calibration=switch(c){case "simulation"->SkatOCalibration.PARAMETRIC_SIMULATION;case "analytic"->SkatOCalibration.ANALYTIC;default->throw new IllegalArgumentException("invalid SKAT-O calibration");};
            if(minimum<1||maxVariants<1||simulations<1||!Double.isFinite(maf)||maf<=0||maf>.5)throw new IllegalArgumentException("invalid count or MAF option");
            if(tests.stream().anyMatch(t->!t.equals("single")) && (groups==null&&(windowSize<1||windowStep<1)||groups!=null&&windowSize!=0))throw new IllegalArgumentException("group tests require --groups or positive --window-size/--window-step");
            Path manifest=Path.of(map.get("--cohorts")).toAbsolutePath();DelimitedData table=DelimitedData.read(manifest);
            int nc=table.column("cohort"),sc=table.column("scores"),cc=table.header().contains("covariance")?table.column("covariance"):-1;
            for(String[] row:table.rows()){String name=row[nc];if(name.isBlank()||name.equals("meta")||names.contains(name))throw new IllegalArgumentException("cohort names must be unique, nonblank, and not meta");names.add(name);scoreFiles.add(manifest.getParent().resolve(row[sc]).normalize());covFiles.add(cc<0||row[cc].equals("NA")||row[cc].isBlank()?null:manifest.getParent().resolve(row[cc]).normalize());}
            if(names.isEmpty())throw new IllegalArgumentException("empty cohort manifest");
        }
        String weight(String test){return weights.equals("default")?test.equals("burden")?"equal":"raremetal-beta":weights;}
    }
    static String help(){return """
        Usage: jlinalg rare-meta --cohorts manifest.tsv --genome-build GRCh38
          --test single|burden|skat|skat-o[,TEST...] --out PREFIX
          [--groups groups.txt | --window-size BP --window-step BP]
          [--weights equal|mb|beta|raremetal-beta] [--maf 0.05]
          [--min-cohorts 1] [--af-policy observed|raremetal] [--cohort-results]
          [--call-rate 0] [--hwe 0] [--max-variants 2000]
          [--skato-calibration simulation|analytic] [--simulations 10000] [--seed 20260901]
        Manifest columns: cohort, scores, covariance (tab separated; paths relative to manifest).
        Score/covariance: RAREMETALWORKER or rvtests; .gz/.tbi recommended.
        Groups: GROUP_ID CHROM:POS:REF:ALT ...; biallelic, one chromosome per group.
        Outputs: PREFIX.TEST.tsv and PREFIX.log; existing files are never overwritten.
        Burden reports beta/SE/direction; SKAT and SKAT-O do not have signed effects.
        SKAT-O simulation has Monte Carlo resolution 1/(simulations+1).
        Analytic SKAT-O is a labeled moment approximation, not an exact tail.
        """;}
}
