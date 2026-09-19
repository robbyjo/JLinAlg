/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import jdistlib.Normal;
import jdistlib.T;
import org.jlinalg.glm.*;
import org.jlinalg.regression.CensoredRegression;
import org.jlinalg.regression.OrdinalRegression;
import org.jlinalg.regression.RareEventsLogistic;
import org.jlinalg.regression.SurveyRegression;

/** Censored, ordinal, rare-event and survey regression commands. */
final class InferenceCli {
    private InferenceCli() { }
    static int run(String command,String[] args,PrintStream out,PrintStream error) {
        try {
            if (Arrays.asList(args).contains("--help") || Arrays.asList(args).contains("-h")) {
                out.println(help(command)); return 0;
            }
            switch(command) {
                case "censored-regression" -> censored(args);
                case "ordinal-regression" -> ordinal(args);
                case "rare-events-logit" -> rareEvents(args);
                case "survey-regression" -> survey(args);
                default -> throw new IllegalArgumentException("unknown regression command: "+command);
            }
            out.println(command+" complete"); return 0;
        } catch (IOException | RuntimeException ex) { error.println("jlinalg: "+ex.getMessage()); return 2; }
    }

    private static void censored(String[] args) throws IOException {
        var o=XwasFiles.options(args,"--input","--response","--predictors","--censor","--distribution","--predict","--lower","--upper","--time","--out");
        var t=DelimitedData.read(XwasFiles.path(o,"--input"));var predictors=names(o,"--predictors");
        CensoredRegression.Distribution distribution;
        try{distribution=CensoredRegression.Distribution.valueOf(XwasFiles.required(o,"--distribution").toUpperCase(Locale.ROOT));}
        catch(IllegalArgumentException e){throw new IllegalArgumentException("--distribution must be gaussian, lognormal, weibull, or exponential");}
        boolean gaussian=distribution==CensoredRegression.Distribution.GAUSSIAN;
        if(gaussian&&o.containsKey("--time")||!gaussian&&(o.containsKey("--lower")||o.containsKey("--upper")))
            throw new IllegalArgumentException("--lower/--upper apply to Gaussian predictions; --time applies to AFT predictions");
        if(!o.containsKey("--predict")&&(o.containsKey("--lower")||o.containsKey("--upper")||o.containsKey("--time")))
            throw new IllegalArgumentException("prediction controls require --predict");
        rejectResponsePredictor(o,predictors);
        var fit=CensoredRegression.fit(column(t,XwasFiles.required(o,"--response")),design(t,predictors,true),
            integers(column(t,XwasFiles.required(o,"--censor"))),distribution);
        var terms=withIntercept(predictors);double[] beta=fit.coefficients();
        if(distribution!=CensoredRegression.Distribution.EXPONENTIAL){terms.add("log(scale)");beta=Arrays.copyOf(beta,beta.length+1);beta[beta.length-1]=Math.log(fit.scale());}
        Path target=XwasFiles.path(o,"--out");var files=new LinkedHashMap<Path,String>();
        files.put(target,coefficients(terms,beta,fit.covariance(),Double.POSITIVE_INFINITY));
        files.put(Path.of(target+".metadata.tsv"),"key\tvalue\ndistribution\t"+distribution+"\nscale\t"+fit.scale()+"\nlog_likelihood\t"+fit.logLikelihood()+"\nconverged\ttrue\ncensor_codes\t-1 left; 0 observed; 1 right\n");
        if(o.containsKey("--predict")){
            double[][] x=design(DelimitedData.read(XwasFiles.path(o,"--predict")),predictors,true);
            StringBuilder pred=new StringBuilder(gaussian?"row\tlatent_mean\tobserved_mean\tleft_censor_probability\tright_censor_probability\n":"row\tmedian_time\tmedian_lower\tmedian_upper\tsurvival\n");
            double lower=o.containsKey("--lower")?XwasFiles.number(o.get("--lower")):Double.NEGATIVE_INFINITY;
            double upper=o.containsKey("--upper")?XwasFiles.number(o.get("--upper")):Double.POSITIVE_INFINITY;
            double time=gaussian?1:XwasFiles.number(XwasFiles.required(o,"--time"));
            for(int i=0;i<x.length;i++) {
                pred.append(i+1).append('\t');
                if(gaussian)pred.append(fit.location(x[i])).append('\t').append(fit.observedMean(x[i],lower,upper)).append('\t')
                    .append(Double.isFinite(lower)?fit.leftCensoringProbability(x[i],lower):0).append('\t')
                    .append(Double.isFinite(upper)?Normal.cumulative(upper,fit.location(x[i]),fit.scale(),false,false):0);
                else {double[] interval=fit.timeQuantileInterval(x[i],.5,.95);pred.append(interval[0]).append('\t').append(interval[1]).append('\t').append(interval[2]).append('\t').append(fit.survival(x[i],time));}
                pred.append('\n');
            }
            files.put(Path.of(target+".predictions.tsv"),pred.toString());
        }
        XwasFiles.publish(files);
    }

    private static void ordinal(String[] args) throws IOException {
        var o=XwasFiles.options(args,"--input","--response","--predictors","--categories","--link","--predict","--out");
        var t=DelimitedData.read(XwasFiles.path(o,"--input"));var predictors=names(o,"--predictors");rejectResponsePredictor(o,predictors);
        var link=OrdinalRegression.Link.valueOf(o.getOrDefault("--link","logit").toUpperCase(Locale.ROOT));
        int categories=Integer.parseInt(XwasFiles.required(o,"--categories"));
        var fit=OrdinalRegression.fit(integers(column(t,XwasFiles.required(o,"--response"))),design(t,predictors,false),categories,link);
        double[] beta=Arrays.copyOf(fit.coefficients(),predictors.size()+categories-1);System.arraycopy(fit.thresholds(),0,beta,predictors.size(),categories-1);
        var terms=new ArrayList<>(predictors);for(int j=0;j<categories-1;j++)terms.add("threshold_"+j+"_"+(j+1));
        Path target=XwasFiles.path(o,"--out");var files=new LinkedHashMap<Path,String>();
        files.put(target,coefficients(terms,beta,fit.covariance(),Double.POSITIVE_INFINITY));
        files.put(Path.of(target+".metadata.tsv"),"key\tvalue\nlink\t"+link+"\nlog_likelihood\t"+fit.logLikelihood()+"\nconverged\ttrue\nparameterization\tP(Y <= k) = F(threshold[k] - x beta); no predictor intercept\n");
        if(o.containsKey("--predict")){
            double[][] x=design(DelimitedData.read(XwasFiles.path(o,"--predict")),predictors,false);StringBuilder pred=new StringBuilder("row");
            for(int j=0;j<categories;j++)pred.append("\tcategory_").append(j);pred.append('\n');
            for(int i=0;i<x.length;i++){pred.append(i+1);for(double prob:fit.probabilities(x[i]))pred.append('\t').append(prob);pred.append('\n');}
            files.put(Path.of(target+".predictions.tsv"),pred.toString());
        }
        XwasFiles.publish(files);
    }

    private static void rareEvents(String[] args) throws IOException {
        var o=XwasFiles.options(args,"--input","--response","--predictors","--bias-correction","--prevalence","--out");
        var t=DelimitedData.read(XwasFiles.path(o,"--input"));var predictors=names(o,"--predictors");rejectResponsePredictor(o,predictors);
        String bias=o.getOrDefault("--bias-correction","true");if(!Set.of("true","false").contains(bias))throw new IllegalArgumentException("--bias-correction must be true or false");
        var fit=RareEventsLogistic.fit(column(t,o.get("--response")),design(t,predictors,true),Boolean.parseBoolean(bias),o.containsKey("--prevalence")?XwasFiles.number(o.get("--prevalence")):null);
        Path target=XwasFiles.path(o,"--out");var files=new LinkedHashMap<Path,String>();
        files.put(target,coefficients(withIntercept(predictors),fit.coefficients(),fit.covariance(),Double.POSITIVE_INFINITY));
        files.put(Path.of(target+".metadata.tsv"),"key\tvalue\nbias_correction\t"+bias+"\nsample_prevalence\t"+fit.samplePrevalence()+"\nprior_intercept_shift\t"+fit.interceptPriorCorrection()
            +"\nuncertainty\tunderlying sample ML Fisher covariance; first-order; supplied prevalence treated as known\n");
        XwasFiles.publish(files);
    }

    private static void survey(String[] args) throws IOException {
        var o=XwasFiles.options(args,"--input","--response","--predictors","--family","--weights","--strata","--psu","--out");
        var t=DelimitedData.read(XwasFiles.path(o,"--input"));var predictors=names(o,"--predictors");rejectResponsePredictor(o,predictors);
        GlmFamily family=switch(o.getOrDefault("--family","gaussian")){
            case "gaussian"->GlmFamilies.gaussian();case "logit"->GlmFamilies.binomial();case "probit"->GlmFamilies.probit();case "poisson"->GlmFamilies.poisson();
            default->throw new IllegalArgumentException("--family must be gaussian, logit, probit, or poisson");};
        var sampling=new SurveyRegression.Design(column(t,XwasFiles.required(o,"--weights")),labels(t,XwasFiles.required(o,"--strata")),labels(t,XwasFiles.required(o,"--psu")));
        var fit=SurveyRegression.fit(column(t,o.get("--response")),design(t,predictors,true),family,sampling);
        Path target=XwasFiles.path(o,"--out");var files=new LinkedHashMap<Path,String>();
        files.put(target,coefficients(withIntercept(predictors),fit.coefficients(),fit.covariance(),fit.degreesOfFreedom()));
        files.put(Path.of(target+".metadata.tsv"),"key\tvalue\nfamily\t"+family.name()+"\nvariance\twith-replacement stratified PSU Taylor linearization; no FPC\npsus\t"+fit.psus()+"\nstrata\t"+fit.strata()+"\ndegrees_of_freedom\t"+fit.degreesOfFreedom()+"\n");
        XwasFiles.publish(files);
    }
    private static String[] labels(DelimitedData t,String name){int c=t.column(name);return t.rows().stream().map(r->r[c]).toArray(String[]::new);}
    private static int[] integers(double[] values){int[] out=new int[values.length];for(int i=0;i<out.length;i++){if(values[i]!=Math.rint(values[i])||values[i]<Integer.MIN_VALUE||values[i]>Integer.MAX_VALUE)throw new IllegalArgumentException("integer category/censor codes required");out[i]=(int)values[i];}return out;}
    private static void rejectResponsePredictor(Map<String,String> o,List<String> predictors){if(predictors.contains(XwasFiles.required(o,"--response")))throw new IllegalArgumentException("response cannot be a predictor");}

    private static String coefficients(List<String> terms,double[] beta,double[] cov,double df) {
        StringBuilder text=new StringBuilder("term\testimate\tstandard_error\tstatistic\tp_value\tlower\tupper\tdf\n");
        double critical=Double.isInfinite(df)?Normal.quantile(.975,0,1,true,false):T.quantile(.975,df,true,false);
        for(int j=0;j<beta.length;j++) {
            double se=Math.sqrt(cov[j*beta.length+j]);
            if (!(se>0) || !Double.isFinite(se)) throw new IllegalArgumentException("positive finite coefficient uncertainty required");
            double statistic=beta[j]/se;
            double p=Double.isInfinite(df)?2*Normal.cumulative(Math.abs(statistic),0,1,false,false):2*T.cumulative(Math.abs(statistic),df,false,false);
            text.append(terms.get(j)).append('\t').append(beta[j]).append('\t').append(se).append('\t').append(statistic)
                .append('\t').append(p).append('\t').append(beta[j]-critical*se).append('\t').append(beta[j]+critical*se).append('\t').append(df).append('\n');
        }
        return text.toString();
    }
    private static List<String> names(Map<String,String> o,String key) { return o.containsKey(key)?XwasFiles.names(o.get(key)):List.of(); }
    private static ArrayList<String> withIntercept(List<String> names) { var result=new ArrayList<String>();result.add("(Intercept)");result.addAll(names);return result; }
    private static double[] column(DelimitedData t,String name) { int c=t.column(name);return t.rows().stream().mapToDouble(r->XwasFiles.number(r[c])).toArray(); }
    private static double[][] design(DelimitedData t,List<String> names,boolean intercept) {
        int add=intercept?1:0;double[][] result=new double[t.rows().size()][names.size()+add];
        for(int i=0;i<result.length;i++){if(intercept)result[i][0]=1;for(int j=0;j<names.size();j++)result[i][j+add]=XwasFiles.number(t,t.rows().get(i),names.get(j));}
        return result;
    }
    private static String help(String command) {
        return switch(command) {
            case "censored-regression" -> """
                censored-regression --input FILE --response COLUMN --censor COLUMN
                  --distribution gaussian|lognormal|weibull|exponential [--predictors COLUMNS]
                  [--predict FILE] [--lower LIMIT] [--upper LIMIT] [--time TIME] --out FILE.tsv
                Censor codes: -1 left, 0 observed, 1 right. Response holds the value/limit.
                An intercept is included. AFT responses must be positive times.
                Gaussian prediction limits use --lower/--upper; AFT predictions require --time.
                No interval censoring, truncation, competing risks, or frailty.
                """;
            case "ordinal-regression" -> """
                ordinal-regression --input FILE --response COLUMN --categories K
                  [--predictors COLUMNS] [--link logit|probit] [--predict FILE] --out FILE.tsv
                Response codes are 0,...,K-1, all observed. Numeric predictors exclude an
                intercept: ordered thresholds supply the location. Full observed-information covariance.
                """;
            case "rare-events-logit" -> """
                rare-events-logit --input FILE --response COLUMN [--predictors COLUMNS]
                  [--bias-correction true|false] [--prevalence VALUE] --out FILE.tsv
                Binary 0/1 responses. King-Zeng coefficient correction and separate prior
                prevalence correction. Underlying logistic ML must converge without separation.
                Covariance is first-order sample ML Fisher covariance; prevalence is treated as known.
                """;
            case "survey-regression" -> """
                survey-regression --input FILE --response COLUMN [--predictors COLUMNS]
                  --weights COLUMN --strata COLUMN --psu COLUMN
                  [--family gaussian|logit|probit|poisson] --out FILE.tsv
                With-replacement stratified PSU Taylor covariance; t tests use PSUs minus strata df.
                Positive sampling weights and at least two PSUs per stratum required. No FPC.
                """;
            default -> throw new IllegalArgumentException("unknown regression command: "+command);
        };
    }
}
