/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import jdistlib.Normal;
import jdistlib.T;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.meta.*;

/** Streaming cross-cohort meta-analysis and numeric moderator regression. */
final class MetaCli {
    private MetaCli() { }
    private static final String HEADER="feature_id\tterm\tmodel\tn_cohorts\tdirection\tstatus"
        +"\tbeta\tse\tstatistic\tp_value\tneg_log10_p\tci_lower\tci_upper"
        +"\ttau_squared\tq\tq_df\tq_p_value\ti_squared\tmoderator_q\tmoderator_q_p_value"
        +"\tprediction_lower\tprediction_upper";

    static int run(String command,String[] arguments,PrintStream console,PrintStream error) {
        try {
            Options options=Options.parse(command,arguments);
            if(options.help) { console.println(help(command)); return 0; }
            execute(options,console);
            return 0;
        } catch(IOException | RuntimeException failure) {
            error.println("jlinalg: "+failure.getMessage()); return 2;
        }
    }

    private static void execute(Options o,PrintStream console) throws IOException {
        Path output=o.output.toAbsolutePath().normalize(),log=Path.of(output+".log");
        PipelinePaths.requireFreshOutputs(output,log);
        List<Path> cohortPaths=new ArrayList<>(o.cohorts.values());
        for(int i=0;i<cohortPaths.size();i++) for(int j=0;j<i;j++)
            if(Files.isSameFile(cohortPaths.get(i),cohortPaths.get(j)))
                throw new IllegalArgumentException("the same input file cannot represent two independent cohorts");
        try (RunLog journal=RunLog.openPlain(log,false)) {
            execute(o,console,output,log,journal);
            journal.complete("complete");
        }
    }

    private static void execute(Options o,PrintStream console,Path output,Path log,RunLog journal) throws IOException {
        double[][] moderators=o.regression?moderators(o):null;
        Path scratch=Files.createTempDirectory(output.getParent(),".jlinalg-meta-");
        List<MetaCohortSource> sources=new ArrayList<>();
        long[] totals=new long[3]; // input union, fitted, excluded
        try {
            for(var cohort:o.cohorts.entrySet()) {
                console.println("Preparing cohort "+cohort.getKey()+": "+cohort.getValue());
                sources.add(new MetaCohortSource(cohort.getValue(),o.idColumn,o.effectColumn,o.seColumn,scratch,o.sortRows));
            }
            Path result=scratch.resolve("result");
            char delimiter=MetaCohortSource.delimiter(output);
            try(BufferedWriter writer=Files.newBufferedWriter(result)) {
                writer.write(HEADER.replace('\t',delimiter)); writer.newLine();
                int k=sources.size();
                while(true) {
                    List<String> ids=new ArrayList<>(o.blockRows);
                    double[] effects=new double[Math.multiplyExact(o.blockRows,k)],errors=new double[effects.length];
                    Arrays.fill(effects,Double.NaN); Arrays.fill(errors,Double.NaN);
                    for(int row=0;row<o.blockRows;row++) {
                        String next=null;
                        for(var source:sources) if(source.current()!=null
                                && (next==null || source.current().id().compareTo(next)<0)) next=source.current().id();
                        if(next==null) break;
                        ids.add(next);
                        for(int j=0;j<k;j++) {
                            var source=sources.get(j); var value=source.current();
                            if(value!=null && value.id().equals(next)) {
                                if(moderators==null || moderators[j]!=null) {
                                    effects[row*k+j]=value.effect(); errors[row*k+j]=value.se();
                                }
                                source.advance();
                            }
                        }
                    }
                    if(ids.isEmpty()) break;
                    effects=Arrays.copyOf(effects,ids.size()*k); errors=Arrays.copyOf(errors,ids.size()*k);
                    var prepared=MetaAnalysis.prepareBatch(effects,errors,ids.size(),k,o.minimum);
                    int[] counts=prepared.cohortCounts();
                    if(o.regression) regression(writer,delimiter,o,ids,effects,errors,moderators,prepared,counts,totals);
                    else analysis(writer,delimiter,o,ids,prepared,counts,totals);
                }
            }
            StringBuilder metadata=new StringBuilder("command="+(o.regression?"meta-regression":"meta-analysis")
                +"\nmodel="+o.options.method()+"\ntau_estimator="+o.options.tauSquaredEstimator()
                +"\ninference="+o.options.inferenceMethod()+"\nconfidence="+o.options.confidenceLevel()
                +"\nmin_cohorts="+o.minimum+"\nfeatures="+totals[0]+"\nfitted_or_passed_through="+totals[1]
                +"\nexcluded="+totals[2]+"\ncohort_order="+String.join(",",o.cohorts.keySet())
                +"\ndirection_symbols=positive:+ negative:- zero:0 missing:?\n"
                +"single_cohort=normal pass-through; heterogeneity unavailable\n"
                +"effect_orientation=as supplied; no automatic allele harmonization\n");
            int index=0;
            for(var cohort:o.cohorts.entrySet()) metadata.append("cohort.").append(++index).append('=')
                .append(cohort.getKey()).append('\t').append(cohort.getValue().toAbsolutePath()).append('\n');
            if(o.regression) metadata.append("moderators=").append(o.moderatorFile.toAbsolutePath())
                .append("\nmoderator_columns=").append(String.join(",",o.moderatorNames))
                .append("\nintercept=").append(o.intercept).append('\n');
            journal.metadata(metadata.toString());
            Files.move(result,output); // no overwrite; publish only after every row succeeds
            console.println("Meta-analysis features: "+totals[0]+"; fitted/pass-through: "+totals[1]+"; excluded: "+totals[2]);
            console.println("Output: "+output+"; cohort order and settings: "+log);
        } finally {
            for(var source:sources) source.close();
            // Only the directory created by this invocation is removed.
            try(var paths=Files.walk(scratch)) {
                for(Path path:paths.sorted(Comparator.reverseOrder()).toList()) Files.deleteIfExists(path);
            }
        }
    }

    private static void analysis(BufferedWriter writer,char delimiter,Options o,List<String> ids,
            PreparedMetaAnalysisBatch prepared,int[] counts,long[] totals) throws IOException {
        MetaAnalysisBatchResult fit=prepared.fit(o.options,o.threads);
        // Take each column once: result accessors defensively copy arrays.
        double[][] columns={fit.pooledEffectSizes(),fit.standardErrors(),fit.statistics(),fit.pValues(),
            fit.negativeLog10PValues(),fit.confidenceLower(),fit.confidenceUpper(),fit.tauSquared(),
            fit.cochranQ(),null,fit.cochranQPValues(),fit.iSquared(),null,null,fit.predictionLower(),fit.predictionUpper()};
        for(int row=0;row<ids.size();row++) {
            String status=counts[row]<o.minimum?"below_min_cohorts":counts[row]==1?"single_cohort":"ok";
            double[] values=new double[16]; Arrays.fill(values,Double.NaN);
            for(int column=0;column<values.length;column++) if(columns[column]!=null) values[column]=columns[column][row];
            if(counts[row]>=Math.max(2,o.minimum)) values[9]=counts[row]-1;
            write(writer,delimiter,ids.get(row),"pooled",o,counts[row],prepared.direction(row),status,values);
            totals[0]++; if(counts[row]>=o.minimum) totals[1]++; else totals[2]++;
        }
    }

    private static void regression(BufferedWriter writer,char delimiter,Options o,List<String> ids,
            double[] effects,double[] errors,double[][] moderators,PreparedMetaAnalysisBatch prepared,
            int[] counts,long[] totals) throws IOException {
        int k=o.cohorts.size(),p=o.moderatorNames.size()+(o.intercept?1:0);
        for(int row=0;row<ids.size();row++) {
            totals[0]++;
            String status=counts[row]<o.minimum?"below_min_cohorts":counts[row]<=p?"insufficient_cohorts":null;
            if(status!=null) { excluded(writer,delimiter,ids.get(row),o,counts[row],prepared.direction(row),status); totals[2]++; continue; }
            double[] y=new double[counts[row]],se=new double[counts[row]];
            double[][] mods=new double[counts[row]][];
            for(int j=0,n=0;j<k;j++) if(!Double.isNaN(effects[row*k+j])) {
                y[n]=effects[row*k+j]; se[n]=errors[row*k+j]; mods[n++]=moderators[j];
            }
            MetaRegressionResult fit;
            try { fit=MetaRegression.fit(y,se,mods,o.moderatorNames,o.intercept,o.options,BackendPolicy.CPU); }
            catch(IllegalArgumentException failure) {
                if(failure.getMessage()!=null && failure.getMessage().startsWith("design matrix is rank deficient")) {
                    excluded(writer,delimiter,ids.get(row),o,counts[row],prepared.direction(row),"rank_deficient"); totals[2]++; continue;
                }
                throw new IllegalArgumentException("feature "+ids.get(row)+": "+failure.getMessage(),failure);
            }
            double[] beta=fit.beta(),ses=fit.standardErrors(),stats=fit.statistics(),ps=fit.pValues(),logps=fit.negativeLog10PValues();
            double critical=o.options.inferenceMethod()==MetaInferenceMethod.NORMAL
                ?Normal.quantile(.5+o.options.confidenceLevel()/2,0,1,true,false)
                :T.quantile(.5+o.options.confidenceLevel()/2,counts[row]-p,true,false);
            for(int term=0;term<p;term++) write(writer,delimiter,ids.get(row),fit.coefficientNames().get(term),o,
                counts[row],prepared.direction(row),"ok",new double[]{beta[term],ses[term],stats[term],ps[term],logps[term],
                    beta[term]-critical*ses[term],beta[term]+critical*ses[term],fit.tauSquared(),fit.residualQ(),
                    fit.residualQDegreesOfFreedom(),fit.residualQPValue(),fit.residualISquared(),
                    fit.moderatorQ(),fit.moderatorQPValue(),Double.NaN,Double.NaN});
            totals[1]++;
        }
    }
    private static void excluded(BufferedWriter writer,char delimiter,String id,Options o,int count,String direction,String status) throws IOException {
        double[] values=new double[16]; Arrays.fill(values,Double.NaN);
        write(writer,delimiter,id,"NA",o,count,direction,status,values);
    }
    private static void write(BufferedWriter writer,char delimiter,String id,String term,Options o,int count,
            String direction,String status,double[] values) throws IOException {
        writer.write(escape(id,delimiter)+delimiter+escape(term,delimiter)+delimiter
            +(o.options.method()==MetaAnalysisMethod.FIXED_EFFECT?"fixed":"random")+delimiter+count+delimiter+direction+delimiter+status);
        for(double value:values) { writer.write(delimiter); writer.write(Double.isNaN(value)?"NA":Double.toString(value)); }
        writer.newLine();
    }
    private static String escape(String value,char delimiter) {
        return value.indexOf(delimiter)>=0 || value.indexOf('"')>=0 ? '"'+value.replace("\"","\"\"")+'"' : value;
    }
    private static double[][] moderators(Options o) throws IOException {
        DelimitedData table=DelimitedData.read(o.moderatorFile);
        int cohort=table.column("cohort");
        int[] columns=o.moderatorNames.stream().mapToInt(table::column).toArray();
        Map<String,double[]> byName=new HashMap<>();
        for(String[] fields:table.rows()) {
            String name=fields[cohort].trim();
            if(name.isEmpty() || byName.containsKey(name)) throw new IllegalArgumentException("duplicate or blank moderator cohort: "+name);
            double[] values=new double[columns.length]; boolean complete=true;
            for(int j=0;j<values.length;j++) { values[j]=MetaCohortSource.number(fields[columns[j]],o.moderatorFile,0); complete&=Double.isFinite(values[j]); }
            byName.put(name,complete?values:null);
        }
        double[][] result=new double[o.cohorts.size()][]; int row=0;
        for(String name:o.cohorts.keySet()) {
            if(!byName.containsKey(name)) throw new IllegalArgumentException("moderator table is missing cohort: "+name);
            result[row++]=byName.get(name);
        }
        return result;
    }

    static String help(String command) {
        return """
            Usage: java -jar jlinalg-<version>.jar COMMAND
              --cohort NAME=FILE [--cohort NAME=FILE ...] --out FILE
              [--model fixed|random] [--tau-estimator reml|dl|pm]
              [--min-cohorts 1] [--id-column feature_id] [--effect-column beta]
              [--se-column se] [--inference normal|t|hk|modified-hk] [--confidence .95]
              [--threads 1] [--block-rows 4096] [--sort-chunk-rows 100000]
              [--max-iterations 200] [--tolerance 1e-10]
            Meta-regression additionally requires:
              --moderator-file FILE --moderators COLUMN[,COLUMN...] [--no-intercept]
            Each cohort: TSV/CSV (optionally .gz), unique feature ID, beta, positive SE.
            Unsorted files are externally sorted; IDs are matched exactly across files.
            Missing pairs: absent ID, NA, NaN, dot, or empty beta/SE. Zero beta is valid.
            Output direction follows --cohort order (+ positive, - negative, 0 zero, ? missing).
            Moderator file: cohort column and numeric moderator columns. Missing moderator
            values exclude that cohort from regression; encode categorical contrasts first.
            Default: random-effects REML, normal inference, minimum one cohort.
            One-cohort analysis rows pass through with normal inference, NA heterogeneity.
            Regression requires more complete cohorts than coefficients and full column rank.
            Below-minimum/unestimable rows are retained with status and NA estimates.
            Effect units and allele orientations must already agree across cohort files.
            Output is TSV (CSV when --out ends .csv), with a .log sidecar. Existing outputs
            are never overwritten. Scratch disk is required alongside the output.
            --threads controls ordinary batch fitting; meta-regression is serial CPU.
            """.replace("COMMAND",command);
    }
    private static final class Options {
        final LinkedHashMap<String,Path> cohorts=new LinkedHashMap<>();
        boolean help,regression,intercept=true;
        Path output,moderatorFile;
        List<String> moderatorNames=List.of();
        String idColumn="feature_id",effectColumn="beta",seColumn="se";
        int minimum=1,threads=1,blockRows=4096,sortRows=100000;
        MetaAnalysisOptions options;
        static Options parse(String command,String[] args) {
            Options o=new Options(); o.regression=command.equals("meta-regression");
            var builder=MetaAnalysisOptions.builder(); Set<String> seen=new HashSet<>();
            for(int i=0;i<args.length;i++) {
                String key=args[i];
                if(key.equals("--help") || key.equals("-h")) { o.help=true; return o; }
                if(!key.equals("--cohort") && !seen.add(key)) throw new IllegalArgumentException("duplicate option: "+key);
                if(key.equals("--no-intercept")) { o.intercept=false; continue; }
                if(i+1==args.length || args[i+1].startsWith("--")) throw new IllegalArgumentException("missing value for "+key);
                String value=args[++i];
                switch(key) {
                    case "--cohort" -> {
                        int equals=value.indexOf('=');
                        if(equals<1 || equals==value.length()-1) throw new IllegalArgumentException("--cohort requires NAME=FILE");
                        String name=value.substring(0,equals);
                        if(!name.matches("[A-Za-z0-9_.-]+") || o.cohorts.putIfAbsent(name,Path.of(value.substring(equals+1)))!=null)
                            throw new IllegalArgumentException("cohort names must be unique and use letters, digits, dot, dash or underscore");
                    }
                    case "--out", "--output" -> o.output=Path.of(value);
                    case "--id-column" -> o.idColumn=value;
                    case "--effect-column" -> o.effectColumn=value;
                    case "--se-column" -> o.seColumn=value;
                    case "--model" -> builder.method(switch(value) {
                        case "fixed" -> MetaAnalysisMethod.FIXED_EFFECT;
                        case "random" -> MetaAnalysisMethod.RANDOM_EFFECT;
                        default -> throw new IllegalArgumentException("--model must be fixed or random"); });
                    case "--tau-estimator" -> builder.tauSquaredEstimator(switch(value.toLowerCase(Locale.ROOT)) {
                        case "reml" -> TauSquaredEstimator.REML;
                        case "dl" -> TauSquaredEstimator.DERSIMONIAN_LAIRD;
                        case "pm" -> TauSquaredEstimator.PAULE_MANDEL;
                        default -> throw new IllegalArgumentException("--tau-estimator must be reml, dl or pm"); });
                    case "--inference" -> builder.inferenceMethod(switch(value) {
                        case "normal" -> MetaInferenceMethod.NORMAL;
                        case "t" -> MetaInferenceMethod.STUDENT_T;
                        case "hk" -> MetaInferenceMethod.HARTUNG_KNAPP;
                        case "modified-hk" -> MetaInferenceMethod.MODIFIED_HARTUNG_KNAPP;
                        default -> throw new IllegalArgumentException("invalid --inference"); });
                    case "--confidence" -> builder.confidenceLevel(Double.parseDouble(value));
                    case "--max-iterations" -> builder.maximumIterations(Integer.parseInt(value));
                    case "--tolerance" -> builder.tolerance(Double.parseDouble(value));
                    case "--min-cohorts" -> o.minimum=positive(value);
                    case "--threads" -> o.threads=positive(value);
                    case "--block-rows" -> o.blockRows=positive(value);
                    case "--sort-chunk-rows" -> o.sortRows=positive(value);
                    case "--moderator-file" -> o.moderatorFile=Path.of(value);
                    case "--moderators" -> o.moderatorNames=Arrays.asList(value.split(",",-1));
                    default -> throw new IllegalArgumentException("unknown meta option: "+key);
                }
            }
            if(o.cohorts.isEmpty() || o.output==null) throw new IllegalArgumentException("--cohort and --out are required");
            if(o.output.toString().toLowerCase(Locale.ROOT).endsWith(".gz"))
                throw new IllegalArgumentException("compressed output is not supported; use .tsv or .csv");
            if(o.regression && (o.moderatorFile==null || o.moderatorNames.isEmpty()))
                throw new IllegalArgumentException("meta-regression requires --moderator-file and --moderators");
            if(!o.regression && (o.moderatorFile!=null || !o.moderatorNames.isEmpty() || !o.intercept))
                throw new IllegalArgumentException("moderator options require meta-regression");
            if(o.moderatorNames.contains("") || new HashSet<>(o.moderatorNames).size()!=o.moderatorNames.size())
                throw new IllegalArgumentException("moderator columns must be unique and nonblank");
            if(o.regression && o.threads!=1) throw new IllegalArgumentException("meta-regression currently requires --threads 1");
            if(new HashSet<>(List.of(o.idColumn,o.effectColumn,o.seColumn)).size()!=3)
                throw new IllegalArgumentException("ID, effect and SE columns must be distinct");
            o.options=builder.build(); return o;
        }
        private static int positive(String value) {
            int result=Integer.parseInt(value);
            if(result<1) throw new IllegalArgumentException("counts must be positive"); return result;
        }
    }
}
