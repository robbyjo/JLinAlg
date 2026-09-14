/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.Locale;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.prediction.AverageMarginalEffect;
import org.jlinalg.prediction.Averaging;
import org.jlinalg.prediction.ExpectedResponse;
import org.jlinalg.prediction.ModelPredictions;
import org.jlinalg.prediction.RiskRatio;
import org.jlinalg.prediction.ScenarioDifference;

/** Fit a GLM and report response-scale predictions or standardized contrasts. */
final class PredictionCli {
    private PredictionCli() { }
    static int run(String[] args,PrintStream output,PrintStream error){try{Options o=Options.parse(args);if(o.help){output.println(help());return 0;}CliNumericTable table=CliNumericTable.read(o.input);double[] response=table.column(o.response);double[][] design=table.design(o.predictors,o.intercept);GlmFamily family=family(o.family);GlmResult fit=Glm.fit(response,design,family);if(!fit.converged())throw new IllegalArgumentException("GLM did not converge: "+fit.convergenceMessage());String result=switch(o.estimand){case "expected"->expected(fit,family,design,o.level);case "difference"->difference(fit,family,scenarios(table,design,o,true),scenarios(table,design,o,false),o);case "risk-ratio"->ratio(fit,family,scenarios(table,design,o,true),scenarios(table,design,o,false),o);case "ame"->ame(fit,family,design,table.designColumn(o.predictors,o.ameColumn,o.intercept),o);default->throw new IllegalArgumentException("--estimand must be expected, difference, risk-ratio, or ame");};if(o.output==null)output.print(result);else{if(Files.exists(o.output)&&!o.overwrite)throw new IOException("output exists; use --overwrite: "+o.output);Files.writeString(o.output,result,StandardCharsets.UTF_8,o.overwrite?new StandardOpenOption[]{StandardOpenOption.CREATE,StandardOpenOption.TRUNCATE_EXISTING}:new StandardOpenOption[]{StandardOpenOption.CREATE_NEW});}return 0;}catch(IOException|IllegalArgumentException failure){error.println("jlinalg: "+failure.getMessage());return 2;}}
    private static double[][] scenarios(CliNumericTable table,double[][] design,Options o,boolean first){if(o.scenarioColumn==null)throw new IllegalArgumentException("scenario estimands require --scenario-column, --first, and --second");int column=table.designColumn(o.predictors,o.scenarioColumn,o.intercept);return table.withColumnSet(design,column,first?o.first:o.second);}
    private static String expected(GlmResult fit,GlmFamily family,double[][] design,double level){ExpectedResponse[] rows=ModelPredictions.expectedResponses(fit,family,design,null,level);StringBuilder text=new StringBuilder("row\tlinear_predictor\tlink_se\tmean\tmean_se\tlower\tupper\n");for(int i=0;i<rows.length;i++){var r=rows[i];text.append(i+1).append('\t').append(r.linearPredictor()).append('\t').append(r.linkStandardError()).append('\t').append(r.mean().estimate()).append('\t').append(r.mean().standardError()).append('\t').append(r.mean().confidenceLower()).append('\t').append(r.mean().confidenceUpper()).append('\n');}return text.toString();}
    private static String difference(GlmResult fit,GlmFamily family,double[][] first,double[][] second,Options o){ScenarioDifference r=ModelPredictions.scenarioDifference(fit,family,first,second,null,null,o.averaging,o.level);return "estimand\tfirst_mean\tsecond_mean\testimate\tstandard_error\tlower\tupper\nscenario_difference\t"+r.firstScenario().estimate()+"\t"+r.secondScenario().estimate()+"\t"+r.difference()+"\t"+r.standardError()+"\t"+r.confidenceLower()+"\t"+r.confidenceUpper()+"\n";}
    private static String ratio(GlmResult fit,GlmFamily family,double[][] first,double[][] second,Options o){RiskRatio r=ModelPredictions.riskRatio(fit,family,first,second,null,null,o.averaging,o.level);return "estimand\tfirst_mean\tsecond_mean\testimate\tstandard_error\tlog_standard_error\tlower\tupper\nrisk_ratio\t"+r.firstScenario().estimate()+"\t"+r.secondScenario().estimate()+"\t"+r.ratio()+"\t"+r.standardError()+"\t"+r.logStandardError()+"\t"+r.confidenceLower()+"\t"+r.confidenceUpper()+"\n";}
    private static String ame(GlmResult fit,GlmFamily family,double[][] design,int column,Options o){AverageMarginalEffect r=ModelPredictions.averageMarginalEffect(fit,family,design,null,column,o.averaging,o.level);return "estimand\tpredictor_column\testimate\tstandard_error\tlower\tupper\naverage_marginal_effect\t"+column+"\t"+r.estimate()+"\t"+r.standardError()+"\t"+r.confidenceLower()+"\t"+r.confidenceUpper()+"\n";}
    private static GlmFamily family(String value){return switch(value){case "gaussian"->GlmFamilies.gaussian();case "binomial"->GlmFamilies.binomial();case "probit","binomial-probit"->GlmFamilies.probit();case "poisson"->GlmFamilies.poisson();case "gamma"->GlmFamilies.gamma();case "inverse-gaussian"->GlmFamilies.inverseGaussian();case "quasi-binomial"->GlmFamilies.quasiBinomial();case "quasi-poisson"->GlmFamilies.quasiPoisson();default->throw new IllegalArgumentException("unsupported --family: "+value);};}
    private static String help(){return """
        Usage: java -jar jlinalg.jar glm-predict --input TABLE --response Y
          --predictors x1,x2 --family gaussian|binomial|probit|poisson|gamma
          --estimand expected|difference|risk-ratio|ame [--out FILE]
        Scenario estimands require --scenario-column X --first A --second B.
        AME requires --ame-column X. Defaults: intercept, population averaging, 95% CI.
        """;}
    private static final class Options{boolean help,overwrite,intercept=true;Path input,output;String response,family="gaussian",estimand="expected",scenarioColumn,ameColumn;List<String> predictors=List.of();double first=Double.NaN,second=Double.NaN,level=.95;Averaging averaging=Averaging.POPULATION_AVERAGE;static Options parse(String[] args){Options o=new Options();for(int i=0;i<args.length;i++){String a=args[i];switch(a){case "--help","-h"->o.help=true;case "--overwrite"->o.overwrite=true;case "--no-intercept"->o.intercept=false;case "--input"->o.input=Path.of(value(args,++i,a));case "--out"->o.output=Path.of(value(args,++i,a));case "--response"->o.response=value(args,++i,a);case "--predictors"->o.predictors=CliNumericTable.columns(value(args,++i,a));case "--family"->o.family=value(args,++i,a).toLowerCase(Locale.ROOT);case "--estimand"->o.estimand=value(args,++i,a).toLowerCase(Locale.ROOT);case "--scenario-column"->o.scenarioColumn=value(args,++i,a);case "--ame-column"->o.ameColumn=value(args,++i,a);case "--first"->o.first=Double.parseDouble(value(args,++i,a));case "--second"->o.second=Double.parseDouble(value(args,++i,a));case "--confidence-level"->o.level=Double.parseDouble(value(args,++i,a));case "--averaging"->{String v=value(args,++i,a).toLowerCase(Locale.ROOT);o.averaging=switch(v){case "at-average"->Averaging.AT_AVERAGE_COVARIATES;case "population"->Averaging.POPULATION_AVERAGE;default->throw new IllegalArgumentException("--averaging must be population or at-average");};}default->throw new IllegalArgumentException("unknown option: "+a);}}if(o.help)return o;if(o.input==null||o.response==null||o.predictors.isEmpty())throw new IllegalArgumentException("--input, --response, and --predictors are required");if(!(o.level>0&&o.level<1))throw new IllegalArgumentException("confidence level must be in (0,1)");if(List.of("difference","risk-ratio").contains(o.estimand)&&(o.scenarioColumn==null||!Double.isFinite(o.first)||!Double.isFinite(o.second)))throw new IllegalArgumentException("scenario estimand requires --scenario-column, --first, and --second");if(o.estimand.equals("ame")&&o.ameColumn==null)throw new IllegalArgumentException("AME requires --ame-column");return o;}private static String value(String[] a,int i,String option){if(i>=a.length||a[i].isBlank())throw new IllegalArgumentException(option+" requires a value");return a[i];}}
}
