/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.regression.InstrumentStrengthDiagnostic;
import org.jlinalg.regression.InstrumentalVariableCovariance;
import org.jlinalg.regression.InstrumentalVariableOptions;
import org.jlinalg.regression.InstrumentalVariableRegression;
import org.jlinalg.regression.InstrumentalVariableResult;

/** Individual-level linear 2SLS command. */
final class InstrumentalVariableCli {
    private InstrumentalVariableCli() { }
    static int run(String[] args,PrintStream output,PrintStream error){try{Options o=Options.parse(args);if(o.help){output.println(help());return 0;}CliNumericTable table=CliNumericTable.read(o.input);double[] y=table.column(o.response);List<String> exogenousNames=new ArrayList<>();if(o.intercept)exogenousNames.add("(Intercept)");exogenousNames.addAll(o.exogenous);double[][] exogenous=table.design(o.exogenous,o.intercept),endogenous=table.matrix(o.endogenous),instruments=table.matrix(o.instruments);InstrumentalVariableOptions settings=new InstrumentalVariableOptions(o.covariance,o.confidence);InstrumentalVariableResult fit=o.covariance.requiresClusters()?InstrumentalVariableRegression.fit(y,exogenous,endogenous,instruments,table.integerGroups(o.cluster),settings,o.backend):InstrumentalVariableRegression.fit(y,exogenous,endogenous,instruments,settings,o.backend);List<Path> outputs=List.of(suffix(o.output,".coefficients.tsv"),suffix(o.output,".strength.tsv"),suffix(o.output,".manifest.tsv"));for(Path path:outputs)if(Files.exists(path)&&!o.overwrite)throw new IOException("output exists; use --overwrite: "+path);writeCoefficients(outputs.get(0),fit,concat(exogenousNames,o.endogenous));writeStrength(outputs.get(1),fit.strengthDiagnostics(),o.endogenous);writeManifest(outputs.get(2),fit,o);output.printf(Locale.ROOT,"observations=%d parameters=%d excluded_instruments=%d minimum_first_stage_f=%g backend=%s%n",fit.observations(),fit.parameters(),fit.excludedInstruments(),fit.minimumFirstStageFStatistic(),fit.backend().selectedBackend());return 0;}catch(IOException|IllegalArgumentException failure){error.println("jlinalg: "+failure.getMessage());return 2;}}
    private static void writeCoefficients(Path path,InstrumentalVariableResult fit,List<String> names)throws IOException{double[] beta=fit.coefficients(),se=fit.standardErrors(),lo=fit.confidenceLower(),hi=fit.confidenceUpper(),stat=fit.associationStatistics().statistics(),p=fit.associationStatistics().pValues();StringBuilder text=new StringBuilder("term\testimate\tstandard_error\tstatistic\tp_value\tlower\tupper\n");for(int i=0;i<beta.length;i++)text.append(names.get(i)).append('\t').append(beta[i]).append('\t').append(se[i]).append('\t').append(stat[i]).append('\t').append(p[i]).append('\t').append(lo[i]).append('\t').append(hi[i]).append('\n');write(path,text.toString());}
    private static void writeStrength(Path path,List<InstrumentStrengthDiagnostic> diagnostics,List<String> names)throws IOException{StringBuilder text=new StringBuilder("endogenous\tpartial_r2\tclassical_f\twald\teffective_f\tnumerator_df\tdenominator_df\tp_value\treference_distribution\n");for(var d:diagnostics)text.append(names.get(d.endogenousColumn())).append('\t').append(d.partialRSquared()).append('\t').append(d.classicalFStatistic()).append('\t').append(d.waldStatistic()).append('\t').append(d.effectiveFStatistic()).append('\t').append(d.numeratorDegreesOfFreedom()).append('\t').append(d.denominatorDegreesOfFreedom()).append('\t').append(d.pValue()).append('\t').append(d.referenceDistribution()).append('\n');write(path,text.toString());}
    private static void writeManifest(Path path,InstrumentalVariableResult fit,Options o)throws IOException{String text="key\tvalue\nmethod\t2sls\ncovariance\t"+o.covariance+"\nobservations\t"+fit.observations()+"\nstructural_rank\t"+fit.structuralRank()+"\ninstrument_rank\t"+fit.instrumentRank()+"\noveridentification_df\t"+fit.overidentificationDegreesOfFreedom()+"\nprojected_design_condition_number\t"+fit.projectedDesignConditionNumber()+"\nresidual_variance\t"+fit.residualVariance()+"\nbackend\t"+fit.backend().selectedBackend()+"\n";write(path,text);}
    private static void write(Path path,String text)throws IOException{Path parent=path.toAbsolutePath().getParent();if(parent!=null)Files.createDirectories(parent);Files.writeString(path,text,StandardCharsets.UTF_8,StandardOpenOption.CREATE,StandardOpenOption.TRUNCATE_EXISTING);}
    private static Path suffix(Path prefix,String suffix){return Path.of(prefix.toString()+suffix);}
    private static List<String> concat(List<String> first,List<String> second){List<String> result=new ArrayList<>(first);result.addAll(second);return result;}
    private static String help(){return """
        Usage: java -jar jlinalg.jar iv-regression --input TABLE --response Y
          --exogenous x1,x2 --endogenous exposure --instruments z1,z2 --out PREFIX
          [--covariance homoskedastic|hc0|hc1|cluster-cr0|cluster-cr1]
          [--cluster COLUMN] [--backend POLICY] [--confidence-level 0.95]
        The intercept is included unless --no-intercept is supplied.
        """;}
    private static final class Options{boolean help,overwrite,intercept=true;Path input,output;String response,cluster;List<String> exogenous=List.of(),endogenous=List.of(),instruments=List.of();InstrumentalVariableCovariance covariance=InstrumentalVariableCovariance.HC1;BackendPolicy backend=BackendPolicy.PREFERRED;double confidence=.95;static Options parse(String[] args){Options o=new Options();for(int i=0;i<args.length;i++){String a=args[i];switch(a){case "--help","-h"->o.help=true;case "--overwrite"->o.overwrite=true;case "--no-intercept"->o.intercept=false;case "--input"->o.input=Path.of(value(args,++i,a));case "--out"->o.output=Path.of(value(args,++i,a));case "--response"->o.response=value(args,++i,a);case "--exogenous"->o.exogenous=CliNumericTable.columns(value(args,++i,a));case "--endogenous"->o.endogenous=CliNumericTable.columns(value(args,++i,a));case "--instruments"->o.instruments=CliNumericTable.columns(value(args,++i,a));case "--cluster"->o.cluster=value(args,++i,a);case "--covariance"->o.covariance=InstrumentalVariableCovariance.valueOf(value(args,++i,a).toUpperCase(Locale.ROOT).replace('-','_'));case "--backend"->o.backend=BackendPolicy.valueOf(value(args,++i,a).toUpperCase(Locale.ROOT).replace('-','_'));case "--confidence-level"->o.confidence=Double.parseDouble(value(args,++i,a));default->throw new IllegalArgumentException("unknown option: "+a);}}if(o.help)return o;if(o.input==null||o.output==null||o.response==null||o.endogenous.isEmpty()||o.instruments.isEmpty())throw new IllegalArgumentException("--input, --response, --endogenous, --instruments, and --out are required");if(o.covariance.requiresClusters()&&o.cluster==null)throw new IllegalArgumentException("cluster covariance requires --cluster");if(!(o.confidence>0&&o.confidence<1))throw new IllegalArgumentException("confidence level must be in (0,1)");return o;}private static String value(String[] a,int i,String option){if(i>=a.length||a[i].isBlank())throw new IllegalArgumentException(option+" requires a value");return a[i];}}
}
