/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.xwas.*;

/** File-based molecular/polygenic/multi-omics score training, application and held-out evaluation. */
final class ScoreCli {
    private ScoreCli() { }
    static int run(String command,String[] args,PrintStream out,PrintStream error) {
        try {
            if(Arrays.asList(args).contains("--help")){out.println(help(command));return 0;}
            if(command.equals("score-train"))train(args);else apply(args);
            out.println(command+" complete");return 0;
        } catch(IOException|RuntimeException ex){error.println("jlinalg: "+ex.getMessage());return 2;}
    }
    private static void train(String[] args) throws IOException {
        var o=XwasFiles.options(args,"--input","--id","--features","--outcome","--lambdas","--alpha","--folds","--seed","--alleles","--name","--out");
        var table=DelimitedData.read(XwasFiles.path(o,"--input"));
        String id=o.getOrDefault("--id","sample"),outcome=XwasFiles.required(o,"--outcome"),name=XwasFiles.id(o.getOrDefault("--name","score"));
        XwasFiles.indexed(table,id);List<String> features=XwasFiles.names(XwasFiles.required(o,"--features"));
        if(features.contains(id)||features.contains(outcome))throw new IllegalArgumentException("features cannot contain ID or outcome columns");
        int n=table.rows().size(),p=features.size();double[] y=new double[n];double[][] x=new double[n][p];
        for(int r=0;r<n;r++){String[] row=table.rows().get(r);y[r]=XwasFiles.number(table,row,outcome);for(int c=0;c<p;c++)x[r][c]=XwasFiles.number(table,row,features.get(c));}
        double[] lambdas=Arrays.stream(XwasFiles.required(o,"--lambdas").split(",")).mapToDouble(XwasFiles::number).toArray();
        var fit=PredictionScores.train(y,x,lambdas,XwasFiles.number(o.getOrDefault("--alpha","0")),Integer.parseInt(o.getOrDefault("--folds","5")),Long.parseLong(o.getOrDefault("--seed","42")));
        Map<String,String[]> alleles=null;DelimitedData alleleTable=null;
        if(o.containsKey("--alleles")){alleleTable=DelimitedData.read(XwasFiles.path(o,"--alleles"));alleles=XwasFiles.indexed(alleleTable,"variant");}
        StringBuilder model=new StringBuilder("model\tvariant\tea\toa\tweight\tintercept\tlambda\talpha\ttraining_mean\n");
        double[] weights=fit.weights();
        for(int c=0;c<p;c++) {
            String ea=".",oa=".";
            if(alleles!=null) {
                String[] row=alleles.get(features.get(c));if(row==null)throw new IllegalArgumentException("missing training feature alleles: "+features.get(c));
                ea=row[alleleTable.column("ea")];oa=row[alleleTable.column("oa")];PredictedOmics.alleleSign(ea,oa,ea,oa);
                for(int r=0;r<n;r++)if(x[r][c]<0||x[r][c]>2)throw new IllegalArgumentException("genotype features require diploid dosages in [0,2]");
            }
            model.append(name).append('\t').append(features.get(c)).append('\t').append(ea).append('\t').append(oa).append('\t').append(weights[c]).append('\t').append(fit.intercept()).append('\t').append(fit.lambda()).append('\t').append(fit.alpha()).append('\t').append(fit.trainingMean()).append('\n');
        }
        StringBuilder ids=new StringBuilder("sample\n");for(String[] row:table.rows())ids.append(row[table.column(id)]).append('\n');
        Path output=XwasFiles.path(o,"--out");Map<Path,String> files=new LinkedHashMap<>();files.put(output,model.toString());files.put(Path.of(output+".training-ids.tsv"),ids.toString());
        files.put(Path.of(output+".metadata.tsv"),"key\tvalue\nmethod\tGaussian elastic net; raw-scale coefficients\nselection\t"+(lambdas.length==1?"prespecified lambda":"minimum training-only CV squared error")+"\nselected_lambda\t"+fit.lambda()+"\nalpha\t"+fit.alpha()+"\nfolds\t"+o.getOrDefault("--folds","5")+"\nseed\t"+o.getOrDefault("--seed","42")+"\ntraining_observations\t"+n+"\ntraining_features\t"+p+"\nmissing\terror\n");
        XwasFiles.publish(files);
    }
    private static void apply(String[] args) throws IOException {
        var o=XwasFiles.options(args,"--input","--weights","--id","--alleles","--outcome","--family","--baseline","--training-ids","--out");
        String family=o.getOrDefault("--family","quantitative");
        if(!Set.of("quantitative","binary").contains(family))throw new IllegalArgumentException("--family must be quantitative or binary");
        if(!o.containsKey("--outcome")&&(o.containsKey("--family")||o.containsKey("--baseline")))throw new IllegalArgumentException("evaluation options require --outcome");
        Path weightPath=XwasFiles.path(o,"--weights");var weights=DelimitedData.read(weightPath);XwasFiles.indexed(weights,"variant");
        var input=DelimitedData.read(XwasFiles.path(o,"--input"));String id=o.getOrDefault("--id","sample");var samples=XwasFiles.indexed(input,id);
        int n=input.rows().size(),p=weights.rows().size();double[] beta=new double[p];double[][] x=new double[n][p];
        double intercept=optional(weights,weights.rows().get(0),"intercept",0),mean=optional(weights,weights.rows().get(0),"training_mean",Double.NaN);
        Map<String,String[]> alleleRows=null;DelimitedData alleles=null;
        if(o.containsKey("--alleles")){alleles=DelimitedData.read(XwasFiles.path(o,"--alleles"));alleleRows=XwasFiles.indexed(alleles,"variant");}
        Set<String> modelNames=new HashSet<>();
        for(int c=0;c<p;c++) {
            String[] row=weights.rows().get(c);String feature=XwasFiles.id(row[weights.column("variant")]);
            if(weights.header().contains("model"))modelNames.add(row[weights.column("model")]);
            if(optional(weights,row,"intercept",0)!=intercept||Double.compare(optional(weights,row,"training_mean",Double.NaN),mean)!=0)
                throw new IllegalArgumentException("model intercept/training mean must be consistent across rows");
            beta[c]=XwasFiles.number(weights,row,"weight");int sign=1;
            String ea=row[weights.column("ea")],oa=row[weights.column("oa")];boolean genetic=!(ea.equals(".")&&oa.equals("."));
            if(genetic) {
                if(alleleRows==null||!alleleRows.containsKey(feature))throw new IllegalArgumentException("genetic weights require target --alleles for every variant");
                String[] a=alleleRows.get(feature);sign=PredictedOmics.alleleSign(a[alleles.column("ea")],a[alleles.column("oa")],ea,oa);
            }
            for(int r=0;r<n;r++) {
                double value=XwasFiles.number(input,input.rows().get(r),feature);
                if(genetic&&(value<0||value>2))throw new IllegalArgumentException("diploid dosage must lie in [0,2]");
                x[r][c]=sign==1?value:2-value;
            }
        }
        if(modelNames.size()>1)throw new IllegalArgumentException("score-apply accepts one model at a time");
        // Model constructor requires a finite mean; it does not affect predictions.
        var model=new PredictionScores.Model(intercept,beta,0,0,Double.isFinite(mean)?mean:0);
        double[] predictions=model.predict(x);StringBuilder text=new StringBuilder("sample\tscore\n");
        for(int r=0;r<n;r++)text.append(input.rows().get(r)[input.column(id)]).append('\t').append(predictions[r]).append('\n');
        Path output=XwasFiles.path(o,"--out");Map<Path,String> files=new LinkedHashMap<>();files.put(output,text.toString());
        if(o.containsKey("--outcome")) {
            Path training=o.containsKey("--training-ids")?XwasFiles.path(o,"--training-ids"):Path.of(weightPath+".training-ids.tsv");
            if(Files.exists(training)) {
                var trainingTable=DelimitedData.read(training);var trainingIds=XwasFiles.indexed(trainingTable,"sample");
                for(String sample:samples.keySet())if(trainingIds.containsKey(sample))throw new IllegalArgumentException("evaluation sample overlaps training IDs: "+sample);
            } else if(o.containsKey("--training-ids")||weights.header().contains("training_mean"))
                throw new IllegalArgumentException("training ID sidecar is required for trained-model evaluation");
            boolean binary=family.equals("binary");
            if(!binary&&!Double.isFinite(mean)&&!o.containsKey("--baseline"))throw new IllegalArgumentException("external quantitative scores require a frozen --baseline prediction column");
            double[] y=new double[n],base=new double[n];
            for(int r=0;r<n;r++){String[] row=input.rows().get(r);y[r]=XwasFiles.number(input,row,o.get("--outcome"));base[r]=o.containsKey("--baseline")?XwasFiles.number(input,row,o.get("--baseline")):binary?.5:mean;}
            var eval=PredictionScores.evaluate(y,predictions,base,binary);
            files.put(Path.of(output+".evaluation.tsv"),"n\trmse\tpredictive_r2\tcalibration_intercept\tcalibration_slope\tauc\n"+n+"\t"+(binary?Double.NaN:eval.rmse())+"\t"+(binary?Double.NaN:eval.predictiveR2())+"\t"+(binary?Double.NaN:eval.calibrationIntercept())+"\t"+(binary?Double.NaN:eval.calibrationSlope())+"\t"+eval.auc()+"\n");
        }
        files.put(Path.of(output+".metadata.tsv"),"key\tvalue\nweights\t"+weightPath+"\nmodel_count\t1\nfeatures\t"+p+"\nmissing\terror\nevaluation\t"+(o.containsKey("--outcome")?"frozen predictions; caller must ensure independent cohorts and preprocessing":"not requested")+"\n");
        XwasFiles.publish(files);
    }
    private static double optional(DelimitedData t,String[] row,String column,double fallback){return t.header().contains(column)?XwasFiles.number(t,row,column):fallback;}
    static String help(String command) {
        return command.equals("score-train")?"""
            Usage: jlinalg score-train --input FILE --features X1,X2 --outcome Y --lambdas 1,0.1,0.01 --out FILE
              [--id sample] [--alpha 0] [--folds 5] [--seed 42] [--name score] [--alleles FILE]
            Gaussian ridge/elastic-net training on a sample-by-feature numeric table; finite complete values required.
            More than one descending lambda invokes training-only CV. Coefficients/intercept are on the input scale.
            Alleles table: variant ea oa. Use for diploid genetic dosages; numeric omics models use dot allele placeholders.
            Output weights also work with twas/pwas when genetic alleles were supplied. Training ID sidecar retained.
            """:"""
            Usage: jlinalg score-apply --input FILE --weights FILE --out FILE [--id sample] [--alleles FILE]
              [--outcome Y --family quantitative|binary --baseline BASELINE_COLUMN --training-ids FILE]
            Input is sample-by-feature; weights: variant ea oa weight, optionally model/intercept/training_mean.
            One model per file. Genetic weights require target allele table; swaps use 2-dosage. Missing features error.
            Quantitative evaluation reports frozen-baseline predictive R2, RMSE and descriptive calibration.
            Binary evaluation reports rank AUC, not absolute-risk calibration. Training overlap is rejected when IDs exist.
            """;
    }
}
