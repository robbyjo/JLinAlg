/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import org.jlinalg.meta.*;
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
        outputs.put("qc",Path.of(o.output+".qc.tsv"));
        List<Path> destinations=new ArrayList<>(outputs.values());destinations.add(log);
        PipelinePaths.requireFreshOutputs(destinations.toArray(Path[]::new));
        try(RunLog journal=RunLog.openPlain(log,false)) {
            try {
                for(int i=0;i<o.scoreFiles.size();i++) {
                    for(int j=0;j<i;j++) if(Files.isSameFile(o.scoreFiles.get(i),o.scoreFiles.get(j)))
                        throw new IllegalArgumentException("one score file cannot represent two independent cohorts");
                    studies.add(new RareMetalStudy(o.scoreFiles.get(i),o.covFiles.get(i),o.callRate,o.hwe).cacheBytes(o.cacheBytes));
                    String declared=studies.get(i).genomeBuild();
                    if(declared!=null&&!declared.equals(o.build))throw new IOException("cohort genome build differs: "+o.names.get(i));
                }
                journal.metadata("command=rare-meta\ngenome_build="+o.build+"\ncohort_order="+String.join(",",o.names)
                    +"\ntests="+String.join(",",o.tests)+"\nweights="+o.weights+"\nmaf="+o.maf+"\naf_policy="+o.afPolicy
                    +"\nmin_cohorts_per_variant="+o.minimum+"\ncohort_results="+o.cohortResults
                    +"\nskato_calibration="+o.calibration+"\nsimulations="+o.simulations+"\nseed="+o.seed
                    +"\nthreads="+o.threads+"\ncache_bytes_per_cohort="+o.cacheBytes+"\ncondition_file="+o.conditionFile+"\ncondition_missing="+o.conditionMissing+"\nleave_variant_out="+o.leaveVariant+"\nleave_cohort_out="+o.leaveCohort+"\nmodel=independent cohorts; quantitative-trait null scores; compatible trait units required\n");
                for(int i=0;i<studies.size();i++)journal.info("cohort="+o.names.get(i)+" scores="+o.scoreFiles.get(i)
                    +" covariance="+o.covFiles.get(i)+" samples="+studies.get(i).samples()+" indexed="+studies.get(i).indexed());
                for(String test:outputs.keySet()) {
                    Path tmp=Files.createTempFile(o.output.getParent(),".rare-meta-",".tsv");temporary.put(test,tmp);
                    BufferedWriter w=Files.newBufferedWriter(tmp);writers.put(test,w);w.write(header(test));w.newLine();
                }
                if(o.tests.contains("single")) singles(o,studies,writers.get("single"));
                try(GroupDispatcher dispatcher=new GroupDispatcher(o,studies,writers)) {
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
                                dispatcher.submit(f[0],variants);
                            }
                        }
                    }else windows(o,studies,dispatcher);
                }
                }
                for(BufferedWriter w:writers.values())w.close();writers.clear();
                for(String test:outputs.keySet())Files.move(temporary.get(test),outputs.get(test));
                journal.complete("complete");console.println("Rare-variant meta-analysis: "+outputs.values());
            }finally {
                for(BufferedWriter w:writers.values())w.close();
                for(RareMetalStudy s:studies)s.close();
                for(Path p:temporary.values())Files.deleteIfExists(p);
            }
        }
    }
    private static String header(String test) {
        if(test.equals("qc"))return "feature_id\tcohort\tvariant_id\tretained\treason\tn_samples\tpooled_maf\tcohort_information\teffect_orientation\tcondition_status\tcall_rate\thwe_pvalue\trole";
        if(test.startsWith("het-"))return header(test.substring(4))+"\tn_effect_dimensions";
        String common="scope\tfeature_id\tn_variants\tn_cohorts\tstatus";
        if(test.equals("vt"))return common+"\tweights\tselected_maf\tselected_beta\tselected_se\tselected_p\tadjusted_p\tmc_se\tn_thresholds\tsimulations\tseed\tcalibration";
        if(test.equals("burden-fixed")||test.equals("burden-random"))return common+"\tweights\tdirection\tbeta\tse\tp_value\ttau_squared\tcochran_q\tq_p\ti_squared\tcalibration";
        return common+(test.equals("single")?"\tn_samples\talt_af\tdirection\tscore\tvariance\tbeta\tse\tz\tp_value\tneg_log10_p":
            test.equals("burden")?"\tweights\tdirection\tbeta\tse\tz\tp_value\tneg_log10_p\tcalibration":
            (test.equals("skat")||test.equals("het-skat"))?"\tweights\tq\tp_value\tneg_log10_p\tcalibration":
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
        List<Variant> conditions=o.conditions.getOrDefault(id,List.of());
        for(Variant v:conditions){if(!chr.equals(v.chromosome))throw new IOException("conditioning variants must share the group chromosome");
            for(Variant target:variants)if(target.position==v.position)throw new IOException("conditioning and target variants overlap: "+v.id());
            low=Math.min(low,v.position);high=Math.max(high,v.position);}
        if(variants.size()+conditions.size()>o.maxVariants)throw new IOException("target plus conditioning variants exceed --max-variants: "+id);
        int k=studies.size(),n=variants.size();
        Score[][] conditioning=new Score[k][conditions.size()];int[][] conditionSigns=new int[k][conditions.size()];boolean[] excludedCohort=new boolean[k]; Score[][] values=new Score[k][n];int[][] signs=new int[k][n];
        double[] af=new double[n]; int[] counts=new int[n];boolean[] retained=new boolean[n];
        for(int c=0;c<k;c++) {
            Map<Long,Score> region=studies.get(c).region(chr,low,high);
            for(int i=0;i<conditions.size();i++) {
                Variant target=conditions.get(i);Score cs=region.get(target.position);conditioning[c][i]=cs;
                if(cs==null||!cs.informative()) {
                    if(o.conditionMissing.equals("error"))throw new IOException("missing/noninformative conditioning variant "+target.id()+" in "+o.names.get(c));
                    excludedCohort[c]=true;
                }else {
                    conditionSigns[c][i]=cs.orientation(target.reference,target.alternate);
                    if(conditionSigns[c][i]==0)throw new IOException("incompatible conditioning alleles: "+target.id());
                }
            }
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
                if(!excludedCohort[c]&&s!=null&&s.informative())counts[i]++;
                if(s!=null&&!s.status().equals("qc_excluded")&&Double.isFinite(s.frequency())){denom+=s.samples();a+=s.samples()*(signs[c][i]==1?s.frequency():1-s.frequency());}
                else if(s==null&&o.afPolicy.equals("raremetal"))denom+=studies.get(c).samples();
            }
            af[i]=denom>0?a/denom:Double.NaN;
            if(af[i]>.5){af[i]=1-af[i];for(int c=0;c<k;c++)signs[c][i]*=-1;}
            retained[i]=counts[i]>=o.minimum&&af[i]>0&&af[i]<=o.maf;
        }
        for(int c=0;c<k;c++)for(int i=0;i<n;i++) {
            Score sc=values[c][i];String reason=excludedCohort[c]?"missing_condition":sc==null?"missing":!sc.informative()?sc.status():counts[i]<o.minimum?"below_min_cohorts":!(af[i]>0&&af[i]<=o.maf)?"maf_filter":"included";
            line(writers.get("qc"),id,o.names.get(c),variants.get(i).id(),retained[i]&&!excludedCohort[c]&&sc!=null&&sc.informative(),reason,sc==null?0:sc.samples(),af[i],sc==null?Double.NaN:sc.variance(),signs[c][i],conditions.isEmpty()?"unconditional":excludedCohort[c]?"cohort_excluded":"conditioned",sc==null?Double.NaN:sc.callRate(),sc==null?Double.NaN:sc.hwePValue(),"target");
        }
        for(int c=0;c<k;c++)for(int i=0;i<conditions.size();i++) {
            Score sc=conditioning[c][i];line(writers.get("qc"),id,o.names.get(c),conditions.get(i).id(),!excludedCohort[c],sc==null?"missing":sc.status(),sc==null?0:sc.samples(),Double.NaN,sc==null?Double.NaN:sc.variance(),conditionSigns[c][i],excludedCohort[c]?"cohort_excluded":"conditioned",sc==null?Double.NaN:sc.callRate(),sc==null?Double.NaN:sc.hwePValue(),"conditioning");
        }
        int m=0;for(boolean b:retained)if(b)m++;
        if(m==0){for(String test:o.tests)if(!test.equals("single"))empty(writers.get(test),test,"meta",id,0,0,"no_variants_after_filters",o.weight(test));return;}
        double[][] u=new double[k][m],cov=new double[k][Math.multiplyExact(m,m)];double[] mafs=new double[m];
        for(int c=0;c<k;c++) {
            Arrays.fill(u[c],Double.NaN);List<Score> present=new ArrayList<>();List<Integer> positions=new ArrayList<>();List<Integer> directions=new ArrayList<>();
            for(int i=0,j=0;i<n;i++)if(retained[i]) {
                mafs[j]=af[i];Score s=values[c][i];
                if(!excludedCohort[c]&&s!=null&&s.informative()){u[c][j]=signs[c][i]*s.score();present.add(s);positions.add(j);directions.add(signs[c][i]);}j++;
            }
            double[] block;
            if(!conditions.isEmpty()&&!present.isEmpty()) {
                List<Score> joint=new ArrayList<>(present);joint.addAll(Arrays.asList(conditioning[c]));
                int size=joint.size();double[] ju=new double[size],jv=studies.get(c).covariance(joint);int[] js=new int[size];
                for(int i=0;i<size;i++){js[i]=i<present.size()?directions.get(i):conditionSigns[c][i-present.size()];ju[i]=js[i]*joint.get(i).score();}
                for(int i=0;i<size;i++)for(int j=0;j<size;j++)jv[i*size+j]*=js[i]*js[j];
                int[] targets=java.util.stream.IntStream.range(0,present.size()).toArray(),ci=java.util.stream.IntStream.range(present.size(),size).toArray();
                var adjusted=SummaryScoreModels.condition(new SetTestScoreState(ju,jv,size),targets,ci);
                block=adjusted.information();double[] au=adjusted.scores();
                for(int i=0;i<present.size();i++){u[c][positions.get(i)]=au[i];directions.set(i,1);}
            }else block=studies.get(c).covariance(present);
            for(int i=0;i<present.size();i++)for(int j=0;j<present.size();j++)
                cov[c][positions.get(i)*m+positions.get(j)]=directions.get(i)*directions.get(j)*block[i*present.size()+j];
            if(!present.isEmpty()) {
                double[] cu=new double[present.size()],cv=new double[block.length];
                for(int i=0;i<cu.length;i++){cu[i]=u[c][positions.get(i)];for(int j=0;j<cu.length;j++)cv[i*cu.length+j]=cov[c][positions.get(i)*m+positions.get(j)];}
                SummarySetTests.validate(new SetTestScoreState(cu,cv,cu.length));
            }
        }
        List<String> ids=new ArrayList<>();for(int i=0;i<n;i++)if(retained[i])ids.add(variants.get(i).id());
        int totalSamples=studies.stream().mapToInt(RareMetalStudy::samples).sum();
        analyze(o,writers,"meta",id,u,cov,mafs,totalSamples,o.cohortResults);
        if(o.leaveCohort)for(int c=0;c<k;c++) {
            double[][] lu=u.clone(),lv=cov.clone();lu[c]=new double[m];Arrays.fill(lu[c],Double.NaN);lv[c]=new double[m*m];
            analyze(o,writers,"leave_cohort:"+o.names.get(c),id,lu,lv,mafs,totalSamples,false);
        }
        if(o.leaveVariant)for(int omitted=0;omitted<m;omitted++) {
            int size=m-1;double[][] lu=new double[k][size],lv=new double[k][size*size];double[] lm=new double[size];
            for(int i=0,ii=0;i<m;i++)if(i!=omitted) {
                lm[ii]=mafs[i];for(int c=0;c<k;c++) {
                    lu[c][ii]=u[c][i];for(int j=0,jj=0;j<m;j++)if(j!=omitted)lv[c][ii*size+jj++]=cov[c][i*m+j];
                }ii++;
            }
            analyze(o,writers,"leave_variant:"+ids.get(omitted),id,lu,lv,lm,totalSamples,false);
        }
    }
    private static void analyze(Options o,Map<String,BufferedWriter> writers,String scope,String id,
            double[][] u,double[][] cov,double[] mafs,int totalSamples,boolean cohortRows)throws IOException {
        int k=u.length,m=mafs.length;
        var pooled=m==0?null:ScoreMetaAnalysis.pool(u,cov,1);
        for(String test:o.tests)if(!test.equals("single")) {
            if(pooled==null){empty(writers.get(test),test,scope,id,0,0,"no_informative_variants",o.weight(test));continue;}
            double[] weights=weights(o.weight(test),mafs,totalSamples);int count=0;StringBuilder direction=new StringBuilder();
            List<SetTestScoreState> states=new ArrayList<>();List<Double> hetWeights=new ArrayList<>();List<MetaStudy> burdens=new ArrayList<>();StringBuilder burdenDirection=new StringBuilder();
            for(int c=0;c<k;c++) {
                int cm=0;double score=0;
                for(int i=0;i<m;i++)if(Double.isFinite(u[c][i])){score+=weights[i]*u[c][i];cm++;}
                if(cm>0)count++;direction.append(cm==0?'?':score>0?'+':score<0?'-':'0');
                if(cm==0){burdenDirection.append('?');if(cohortRows)empty(writers.get(test),test,o.names.get(c),id,0,0,"no_informative_variants",o.weight(test));continue;}
                if(!cohortRows&&!test.startsWith("het-")&&!test.startsWith("burden-"))continue;
                var cp=ScoreMetaAnalysis.pool(new double[][]{u[c]},new double[][]{cov[c]},1);
                double[] cw=select(weights,cp.indices()),cf=select(mafs,cp.indices());
                states.add(cp.state());for(double w:cw)hetWeights.add(w);
                if(cm==m&&test.startsWith("burden-")) {
                    var burden=SummarySetTests.burden(id,cp.state(),cw);
                    if(Double.isFinite(burden.standardError())){burdens.add(new MetaStudy(o.names.get(c),burden.beta(),burden.standardError()));burdenDirection.append(burden.beta()>0?'+':burden.beta()<0?'-':'0');}else burdenDirection.append('?');
                    if(cohortRows)line(writers.get(test),o.names.get(c),id,m,1,Double.isFinite(burden.standardError())?"cohort_estimate":"no_burden_information",o.weight(test),""+burdenDirection.charAt(c),burden.beta(),burden.standardError(),burden.pValue(),Double.NaN,Double.NaN,Double.NaN,Double.NaN,"cohort_normal");
                }
                if(cm<m){burdenDirection.append('?');if(cohortRows&&test.startsWith("burden-"))empty(writers.get(test),test,o.names.get(c),id,cm,0,"incomplete_burden_mask",o.weight(test));}
                if(cohortRows&&!test.startsWith("burden-"))result(writers.get(test),test,o.names.get(c),id,cp.state(),cw,cf,1,"ok",""+direction.charAt(c),o);
            }
            if(test.startsWith("burden-")) {
                if(burdens.isEmpty()){empty(writers.get(test),test,scope,id,m,0,"no_complete_mask_cohorts",o.weight(test));continue;}
                if(burdens.size()==1) {
                    var only=burdens.get(0);double z=only.effectSize()/only.standardError();
                    double p=2*jdistlib.Normal.cumulative(-Math.abs(z),0,1,true,false);
                    line(writers.get(test),scope,id,m,1,"single_cohort",o.weight(test),burdenDirection,only.effectSize(),only.standardError(),p,Double.NaN,Double.NaN,Double.NaN,Double.NaN,"single_cohort_normal");continue;
                }
                var fit=MetaAnalysis.fit(burdens,test.equals("burden-fixed")?MetaAnalysisOptions.fixedEffect():MetaAnalysisOptions.randomEffects(),org.jlinalg.compute.BackendPolicy.CPU);
                line(writers.get(test),scope,id,m,burdens.size(),burdens.size()==1?"single_cohort":burdens.size()<count?"incomplete_mask_cohorts_excluded":"ok",o.weight(test),burdenDirection,fit.pooledEffectSize(),fit.standardError(),fit.pValue(),fit.tauSquared(),fit.cochranQ(),fit.cochranQPValue(),fit.iSquared(),test.equals("burden-fixed")?"inverse_variance_normal":"REML_normal");
                continue;
            }
            SetTestScoreState state=pooled.state();double[] w=select(weights,pooled.indices()),f=select(mafs,pooled.indices());
            if(test.startsWith("het-")) {
                int dimension=states.stream().mapToInt(SetTestScoreState::variants).sum();
                if(dimension>o.maxVariants)throw new IOException("heterogeneous cohort-by-variant dimension exceeds --max-variants: "+id);
                state=SummaryScoreModels.heterogeneous(states);w=hetWeights.stream().mapToDouble(Double::doubleValue).toArray();
            }
            result(writers.get(test),test,scope,id,state,w,f,count,count==1?"single_cohort":"ok",direction.toString(),o);
        }
    }
    private static double[] select(double[] values,int[] indices){double[] result=new double[indices.length];for(int i=0;i<indices.length;i++)result[i]=values[indices[i]];return result;}
    private static double[] weights(String method,double[] maf,int samples) {
        double[] weights=new double[maf.length];
        for(int i=0;i<weights.length;i++)weights[i]=switch(method){
            case "equal"->1;case "mb"->1/Math.sqrt(maf[i]*(1-maf[i]));
            case "beta"->25*Math.pow(1-maf[i],24);
            case "raremetal-beta"->25*Math.pow(1-Math.min(1-2.0/samples,Math.max(2.0/samples,maf[i])),24);
            default->throw new IllegalArgumentException("unknown weights: "+method);};
        return weights;
    }
    private static void result(BufferedWriter w,String test,String scope,String id,SetTestScoreState state,double[] weights,double[] mafs,int cohorts,String status,String direction,Options o)throws IOException {
        int reportedVariants=test.startsWith("het-")?mafs.length:state.variants();
        if(test.equals("vt")) {
            var r=SummaryScoreModels.variableThreshold(id,state,weights,mafs,o.simulations,o.seed);
            if(r.selectedBurden()==null){empty(w,test,scope,id,state.variants(),cohorts,"no_burden_information",o.weight(test));return;}
            var b=r.selectedBurden();line(w,scope,id,state.variants(),cohorts,status,o.weight(test),r.selectedMaf(),b.beta(),b.standardError(),b.pValue(),r.adjustedPValue(),r.monteCarloStandardError(),r.thresholds(),r.simulations(),r.seed(),r.simulations()==0?"single_threshold_normal":"correlated_gaussian_score_simulation");
        }else if(test.equals("skat-o")||test.equals("het-skat-o")) {
            SetTestOptions defaults=SetTestOptions.defaults();
            var options=new SetTestOptions(defaults.variantFilter(),defaults.missingPolicy(),defaults.skatORhoGrid(),o.simulations,o.seed,o.calibration);
            SkatOResult r=SummarySetTests.skatO(id,state,weights,options);
            String components=String.join(";",r.components().stream().map(c->c.rho()+":"+c.result().pValue()).toList());
            List<Object> fields=new ArrayList<>(List.of(scope,id,reportedVariants,cohorts,status,o.weight(test),r.minimumComponentPValue(),r.adjustedPValue(),r.negativeLog10AdjustedPValue(),
                r.components().isEmpty()?"rank_one_chi_square":o.calibration==SkatOCalibration.DETERMINISTIC?"conditional_gaussian_quadrature_rel_tol_1e-6":o.calibration==SkatOCalibration.ANALYTIC?"analytic_moment_approximation":"gaussian_score_simulation",r.simulations(),r.randomSeed(),components));
            if(test.startsWith("het-"))fields.add(state.variants());line(w,fields.toArray());
        }else {
            SetTestResult r=test.equals("burden")?SummarySetTests.burden(id,state,weights):SummarySetTests.skat(id,state,weights);
            if(test.equals("burden"))line(w,scope,id,state.variants(),cohorts,Double.isNaN(r.beta())?"no_burden_information":status,o.weight(test),direction,r.beta(),r.standardError(),r.statistic(),r.pValue(),r.negativeLog10PValue(),r.pValueMethod());
            else {List<Object> fields=new ArrayList<>(List.of(scope,id,reportedVariants,cohorts,status,o.weight(test),r.statistic(),r.pValue(),r.negativeLog10PValue(),r.pValueMethod()));if(test.startsWith("het-"))fields.add(state.variants());line(w,fields.toArray());}
        }
    }
    private static void empty(BufferedWriter w,String test,String scope,String id,int variants,int cohorts,String status,String weight)throws IOException {
        int n=header(test).split("\t").length;Object[] values=new Object[n];Arrays.fill(values,"NA");
        values[0]=scope;values[1]=id;values[2]=variants;values[3]=cohorts;values[4]=status;values[5]=weight;line(w,values);
    }
    private static void line(BufferedWriter w,Object... values)throws IOException {
        for(int i=0;i<values.length;i++){if(i>0)w.write('\t');Object v=values[i];w.write(v instanceof Double d && Double.isNaN(d)?"NA":v.toString());}w.newLine();
    }
    private static void windows(Options o,List<RareMetalStudy> studies,GroupDispatcher dispatcher)throws IOException {
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
                    if(!members.isEmpty())dispatcher.submit(chromosome+":"+start+"-"+stop,new ArrayList<>(members.values()));
                    last=index;
                }
            }
        }
    }
    private static final class GroupDispatcher implements AutoCloseable {
        final Options options;final List<RareMetalStudy> studies;final Map<String,BufferedWriter> writers;
        final ExecutorService executor;final ArrayDeque<Future<Map<String,String>>> pending=new ArrayDeque<>();final Set<String> conditionGroups=new HashSet<>();
        GroupDispatcher(Options o,List<RareMetalStudy> s,Map<String,BufferedWriter> w){options=o;studies=s;writers=w;executor=Executors.newFixedThreadPool(o.threads);}
        void submit(String id,List<Variant> variants)throws IOException {
            if(options.conditions.containsKey(id))conditionGroups.add(id);
            pending.add(executor.submit(()->{
                Map<String,StringWriter> text=new LinkedHashMap<>();Map<String,BufferedWriter> output=new LinkedHashMap<>();
                for(String key:writers.keySet()){var sink=new StringWriter();text.put(key,sink);output.put(key,new BufferedWriter(sink));}
                group(options,studies,output,id,variants);for(var w:output.values())w.flush();
                Map<String,String> result=new LinkedHashMap<>();text.forEach((key,value)->result.put(key,value.toString()));return result;
            }));
            if(pending.size()>=options.threads)drain();
        }
        void drain()throws IOException {
            try{for(var e:pending.remove().get().entrySet())writers.get(e.getKey()).write(e.getValue());}
            catch(InterruptedException e){Thread.currentThread().interrupt();throw new IOException("group execution interrupted",e);}
            catch(ExecutionException e){throw new IOException("group analysis failed: "+e.getCause().getMessage(),e.getCause());}
        }
        public void close()throws IOException {
            try{while(!pending.isEmpty())drain();if(!conditionGroups.containsAll(options.conditions.keySet()))throw new IOException("condition file names groups that were not analyzed");}finally{executor.shutdownNow();
                try{if(!executor.awaitTermination(60,TimeUnit.SECONDS))throw new IOException("group workers did not terminate");}
                catch(InterruptedException e){Thread.currentThread().interrupt();throw new IOException("group shutdown interrupted",e);}}
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
        final Map<String,List<Variant>> conditions=new HashMap<>();
        Path conditionFile;String conditionMissing="error";boolean leaveVariant,leaveCohort;int threads=1;long cacheBytes=16L*1024*1024;
        Path output,groups;String build,weights="default",afPolicy="observed";int minimum=1,maxVariants=2000,simulations=10000;long seed=20260901L,windowSize,windowStep;double maf=.05,callRate=0,hwe=0;boolean cohortResults;SkatOCalibration calibration=SkatOCalibration.PARAMETRIC_SIMULATION;
        Options(String[] args)throws IOException {
            Map<String,String> map=new HashMap<>();
            for(int i=0;i<args.length;i++){String key=args[i];if(key.equals("--leave-variant-out")){leaveVariant=true;continue;}if(key.equals("--leave-cohort-out")){leaveCohort=true;continue;}if(key.equals("--cohort-results")){cohortResults=true;continue;}if(!Set.of("--cohorts","--out","--test","--groups","--genome-build","--weights","--af-policy","--min-cohorts","--max-variants","--simulations","--seed","--window-size","--window-step","--maf","--call-rate","--hwe","--skato-calibration","--condition","--condition-missing","--threads","--cache-mb").contains(key)||i+1==args.length||map.put(key,args[++i])!=null)throw new IllegalArgumentException("unknown, duplicate, or incomplete option: "+key);}
            for(String key:List.of("--cohorts","--out","--test","--genome-build"))if(!map.containsKey(key))throw new IllegalArgumentException("required option: "+key);
            output=Path.of(map.get("--out")).toAbsolutePath();build=map.get("--genome-build");
            for(String t:map.get("--test").split(","))if(!Set.of("single","burden","skat","skat-o","vt","het-skat","het-skat-o","burden-fixed","burden-random").contains(t)||tests.contains(t))throw new IllegalArgumentException("invalid or repeated test: "+t);else tests.add(t);
            if(map.containsKey("--groups"))groups=Path.of(map.get("--groups"));
            weights=map.getOrDefault("--weights",weights);if(!Set.of("default","equal","mb","beta","raremetal-beta").contains(weights))throw new IllegalArgumentException("invalid weights");
            afPolicy=map.getOrDefault("--af-policy",afPolicy);if(!Set.of("observed","raremetal").contains(afPolicy))throw new IllegalArgumentException("invalid AF policy");
            minimum=Integer.parseInt(map.getOrDefault("--min-cohorts","1"));maxVariants=Integer.parseInt(map.getOrDefault("--max-variants","2000"));simulations=Integer.parseInt(map.getOrDefault("--simulations","10000"));seed=Long.parseLong(map.getOrDefault("--seed","20260901"));
            windowSize=Long.parseLong(map.getOrDefault("--window-size","0"));windowStep=Long.parseLong(map.getOrDefault("--window-step",Long.toString(windowSize)));
            maf=Double.parseDouble(map.getOrDefault("--maf","0.05"));callRate=Double.parseDouble(map.getOrDefault("--call-rate","0"));hwe=Double.parseDouble(map.getOrDefault("--hwe","0"));
            String c=map.getOrDefault("--skato-calibration","simulation");calibration=switch(c){case "simulation"->SkatOCalibration.PARAMETRIC_SIMULATION;case "analytic"->SkatOCalibration.ANALYTIC;case "deterministic"->SkatOCalibration.DETERMINISTIC;default->throw new IllegalArgumentException("invalid SKAT-O calibration");};
            if(minimum<1||maxVariants<1||simulations<1||!Double.isFinite(maf)||maf<=0||maf>.5)throw new IllegalArgumentException("invalid count or MAF option");
            if(tests.stream().anyMatch(t->!t.equals("single")) && (groups==null&&(windowSize<1||windowStep<1)||groups!=null&&windowSize!=0))throw new IllegalArgumentException("group tests require --groups or positive --window-size/--window-step");
            threads=Integer.parseInt(map.getOrDefault("--threads","1"));cacheBytes=Math.multiplyExact(Long.parseLong(map.getOrDefault("--cache-mb","16")),1024L*1024);
            if(threads<1||threads>256||cacheBytes<0)throw new IllegalArgumentException("threads must be 1..256 and cache-mb nonnegative");
            conditionMissing=map.getOrDefault("--condition-missing","error");
            if(!Set.of("error","exclude").contains(conditionMissing))throw new IllegalArgumentException("condition-missing must be error or exclude");
            if(map.containsKey("--condition")) {
                if(tests.contains("single"))throw new IllegalArgumentException("--condition applies to group tests; omit single from --test");
                conditionFile=Path.of(map.get("--condition"));
                try(BufferedReader r=RareMetalStudy.open(conditionFile)) {
                    for(String line;(line=r.readLine())!=null;)if(!line.isBlank()&&!line.startsWith("#")) {
                        String[] f=line.trim().split("\\s+");List<Variant> list=new ArrayList<>();Set<String> seen=new HashSet<>();
                        if(f.length<2||conditions.containsKey(f[0]))throw new IOException("condition file requires unique group IDs and variants");
                        for(int i=1;i<f.length;i++){Variant v=Variant.parse(f[i]);if(!seen.add(v.chromosome+":"+v.position))throw new IOException("duplicate conditioning position");list.add(v);}conditions.put(f[0],List.copyOf(list));
                    }
                }
            }
            Path manifest=Path.of(map.get("--cohorts")).toAbsolutePath();DelimitedData table=DelimitedData.read(manifest);
            int nc=table.column("cohort"),sc=table.column("scores"),cc=table.header().contains("covariance")?table.column("covariance"):-1;
            for(String[] row:table.rows()){String name=row[nc];if(name.isBlank()||name.equals("meta")||names.contains(name))throw new IllegalArgumentException("cohort names must be unique, nonblank, and not meta");names.add(name);scoreFiles.add(manifest.getParent().resolve(row[sc]).normalize());covFiles.add(cc<0||row[cc].equals("NA")||row[cc].isBlank()?null:manifest.getParent().resolve(row[cc]).normalize());}
            if(names.isEmpty())throw new IllegalArgumentException("empty cohort manifest");
            if((leaveVariant||leaveCohort)&&tests.stream().allMatch(t->t.equals("single")))throw new IllegalArgumentException("leave-out diagnostics require a group test");
        }
        String weight(String test){return weights.equals("default")?(test.startsWith("burden")||test.equals("vt"))?"equal":"raremetal-beta":weights;}
    }
    static String help(){return """
        Usage: jlinalg rare-meta --cohorts manifest.tsv --genome-build GRCh38
          --test single|burden|skat|skat-o|vt|het-skat|het-skat-o|burden-fixed|burden-random[,TEST...] --out PREFIX
          [--groups groups.txt | --window-size BP --window-step BP]
          [--weights equal|mb|beta|raremetal-beta] [--maf 0.05]
          [--min-cohorts 1] [--af-policy observed|raremetal] [--cohort-results]
          [--call-rate 0] [--hwe 0] [--max-variants 2000]
          [--skato-calibration simulation|analytic|deterministic] [--simulations 10000] [--seed 20260901]
          [--leave-variant-out] [--leave-cohort-out] [--threads 1] [--cache-mb 16]
          [--condition group-conditions.txt] [--condition-missing error|exclude]
        Conditions: GROUP_ID CHROM:POS:REF:ALT ...; disjoint from targets; complete cross-covariance required.
        Diagnostics retain original weights/MAFs and report omitted identity in scope.
        VT calibrates the correlated threshold search; selected beta/SE are descriptive after selection.
        Heterogeneous kernels stack independent cohort effects; burden-fixed/random require complete masks.
        Deterministic SKAT-O uses conditional Gaussian quadrature (relative tolerance 1e-6); unresolved tails fail.
        Cache budget is per cohort; group results are emitted in input order with bounded concurrency.
        Manifest columns: cohort, scores, covariance (tab separated; paths relative to manifest).
        Score/covariance: RAREMETALWORKER or rvtests; .gz/.tbi recommended.
        Groups: GROUP_ID CHROM:POS:REF:ALT ...; biallelic, one chromosome per group.
        Outputs: PREFIX.TEST.tsv, PREFIX.qc.tsv and PREFIX.log; existing files are never overwritten.
        Burden reports beta/SE/direction; SKAT and SKAT-O do not have signed effects.
        SKAT-O simulation has Monte Carlo resolution 1/(simulations+1).
        Analytic SKAT-O is a labeled moment approximation, not an exact tail.
        """;}
}
