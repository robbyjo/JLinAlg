/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.confounding.AutoSva;
import org.jlinalg.confounding.ComBat;
import org.jlinalg.confounding.LatentFactorResult;
import org.jlinalg.confounding.Peer;
import org.jlinalg.confounding.PeerResult;
import org.jlinalg.confounding.PrincipalComponentConfounders;
import org.jlinalg.confounding.SurrogateVariableAnalysis;
import org.jlinalg.confounding.SvaResult;
import org.jlinalg.pipeline.DelimitedMatrixSource;
import org.jlinalg.pipeline.NumericBlock;
import org.jlinalg.pipeline.NumericBlockReader;

/** File-oriented latent-factor estimation and ComBat adjustment. */
final class ConfounderCli {
    private ConfounderCli() { }

    static int run(String command, String[] arguments, PrintStream output,
            PrintStream error) {
        try {
            Options options = Options.parse(command, arguments);
            if (options.help) { output.println(help(command)); return 0; }
            MatrixInput matrix = loadMatrix(options.omics, options.pheno,
                options.idColumn);
            if (command.equals("batch-adjust")) {
                runCombat(matrix, options, output);
            } else {
                runFactors(matrix, options, output);
            }
            return 0;
        } catch (IOException | IllegalArgumentException exception) {
            error.println("jlinalg: " + exception.getMessage());
            return 2;
        }
    }

    private static void runFactors(MatrixInput input, Options options,
            PrintStream output) throws IOException {
        double[][] full = input.design(options.fullColumns, true);
        double[][] reduced = input.design(options.nullColumns, true);
        BackendPolicy policy = BackendPolicy.valueOf(
            options.backend.toUpperCase(Locale.ROOT).replace('-', '_'));
        try (BackendContext context = BackendContext.select(policy)) {
            double[][] factors;
            double[][] loadings;
            double[][] adjusted;
            double[] variance;
            double[] weights = null;
            List<SvaResult.FactorSelection> selection = List.of();
            int iterations = 1;
            boolean converged = true;
            int selected;
            switch (options.method) {
                case "pca" -> {
                    selected = options.requireFactorCount();
                    double[][] source = options.fullColumns.isEmpty() ? input.values
                        : residualizeForPca(input.values, full);
                    LatentFactorResult result = PrincipalComponentConfounders.fit(source,
                        new PrincipalComponentConfounders.Options(selected,
                            options.center, options.scale), context.backend());
                    factors = result.factors(); loadings = result.loadings();
                    variance = result.varianceExplained();
                    adjusted = removeFromOriginal(input.values, factors);
                }
                case "sva" -> {
                    if (options.fullColumns.isEmpty())
                        throw new IllegalArgumentException("SVA requires --full-design");
                    selected = options.factorCount == null
                        ? SurrogateVariableAnalysis.estimateFactorCount(input.values,
                            full, options.permutations, options.seed, context.backend())
                        : options.factorCount;
                    if (selected == 0) throw new IllegalArgumentException(
                        "SVA estimated zero factors; no factor files were written");
                    SvaResult result = SurrogateVariableAnalysis.fit(input.values,
                        full, reduced, new SurrogateVariableAnalysis.Options(
                            selected, options.iterations), context.backend());
                    factors = result.surrogateVariables(); loadings = null;
                    variance = result.varianceExplained(); adjusted = result.adjusted();
                    weights = result.featureWeights(); iterations = result.iterations();
                }
                case "autosva" -> {
                    if (options.fullColumns.isEmpty())
                        throw new IllegalArgumentException("AutoSVA requires --full-design");
                    AutoSva.Options auto = options.factorCount == null
                        ? new AutoSva.Options(0, options.iterations, options.mrse,
                            options.maximumFactors, options.increment)
                        : new AutoSva.Options(options.factorCount, options.iterations,
                            options.mrse, 0, options.increment);
                    SvaResult result = AutoSva.fit(input.values, full, reduced,
                        full, auto, context.backend());
                    selected = result.factorCount(); factors = result.surrogateVariables();
                    loadings = null; variance = result.varianceExplained();
                    adjusted = result.adjusted(); weights = result.featureWeights();
                    selection = result.selectionPath(); iterations = result.iterations();
                    converged = result.converged();
                }
                case "peer" -> {
                    selected = options.requireFactorCount();
                    double[][] covariates = input.design(options.fullColumns, false);
                    if (options.fullColumns.isEmpty()) covariates = null;
                    Peer.Options peerOptions = new Peer.Options(selected,
                        options.maximumIterations, options.tolerance,
                        options.varianceTolerance, options.addMean, options.seed,
                        0.001, 0.1, 0.1, 10.0, 100.0);
                    PeerResult result = Peer.fit(input.values, covariates, peerOptions);
                    factors = result.factors(); loadings = result.loadings();
                    adjusted = result.residuals(); variance = new double[selected];
                    iterations = result.iterations(); converged = result.converged();
                }
                default -> throw new IllegalArgumentException(
                    "--method must be pca, sva, autosva, or peer");
            }
            List<Path> outputs = factorOutputs(options, loadings != null,
                weights != null, !selection.isEmpty());
            requireFresh(outputs, options.overwrite);
            writeFactors(outputs.get(0), input.sampleIds, factors, label(options.method));
            int index = 1;
            if (loadings != null) writeFeatureMatrix(outputs.get(index++), input.featureIds,
                loadings, label(options.method));
            if (weights != null) writeWeights(outputs.get(index++), input.featureIds, weights);
            if (!selection.isEmpty()) writeSelection(outputs.get(index++), selection);
            if (options.writeAdjusted) writeFeatureMatrix(outputs.get(index++),
                input.featureIds, adjusted, "sample", input.sampleIds);
            writeManifest(outputs.get(index), options, input, selected, iterations,
                converged, policy, variance);
            output.printf(Locale.ROOT,
                "method=%s aligned_samples=%d features=%d factors=%d backend=%s%n",
                options.method, input.sampleIds.size(), input.featureIds.size(),
                selected, context.provenance().selectedBackend());
        }
    }

    private static void runCombat(MatrixInput input, Options options,
            PrintStream output) throws IOException {
        if (!options.method.equals("combat"))
            throw new IllegalArgumentException("batch-adjust currently requires --method combat");
        if (options.batchColumn == null || input.phenotypeRows == null)
            throw new IllegalArgumentException("ComBat requires --pheno, --id, and --batch");
        String[] batch = input.strings(options.batchColumn);
        double[][] covariates = input.design(options.fullColumns, true);
        if (options.fullColumns.isEmpty()) covariates = null;
        ComBat.Options combatOptions = new ComBat.Options(options.parametric,
            options.meanOnly, options.referenceBatch, options.tolerance,
            options.maximumIterations);
        ComBat.Result result = ComBat.adjust(input.values, batch, covariates, combatOptions);
        Path adjusted = suffix(options.output, ".adjusted.tsv");
        Path manifest = suffix(options.output, ".manifest.tsv");
        requireFresh(List.of(adjusted, manifest), options.overwrite);
        writeFeatureMatrix(adjusted, input.featureIds, result.adjusted(),
            "sample", input.sampleIds);
        writeManifest(manifest, options, input, 0, 0, true,
            BackendPolicy.CPU, new double[0]);
        output.printf(Locale.ROOT,
            "method=combat aligned_samples=%d features=%d batches=%d mean_only=%s unadjusted_features=%d%n",
            input.sampleIds.size(), input.featureIds.size(), result.batches().size(),
            result.meanOnly(), result.unadjustedFeatures().length);
    }

    private static MatrixInput loadMatrix(Path path, Path phenotype,
            String idColumn) throws IOException {
        DelimitedMatrixSource source = DelimitedMatrixSource.open(path);
        List<String> sourceSamples = source.metadata().sampleIds();
        DelimitedData table = phenotype == null ? null : DelimitedData.read(phenotype);
        List<String[]> alignedRows = null;
        List<String> selectedSamples = sourceSamples;
        int[] order = null;
        if (table != null) {
            if (idColumn == null) throw new IllegalArgumentException("--id is required with --pheno");
            int id = table.column(idColumn);
            Map<String, String[]> byId = new LinkedHashMap<>();
            for (String[] row : table.rows()) {
                String key = row[id].trim();
                if (key.isEmpty() || byId.put(key, row) != null)
                    throw new IllegalArgumentException("phenotype IDs must be unique and nonblank");
            }
            List<Integer> selected = new ArrayList<>(); selectedSamples = new ArrayList<>(); alignedRows = new ArrayList<>();
            for (int i = 0; i < sourceSamples.size(); i++) {
                String[] row = byId.get(sourceSamples.get(i));
                if (row != null) { selected.add(i); selectedSamples.add(sourceSamples.get(i)); alignedRows.add(row); }
            }
            if (selectedSamples.isEmpty()) throw new IllegalArgumentException("omics and phenotype have no aligned sample IDs");
            order = selected.stream().mapToInt(Integer::intValue).toArray();
        }
        long rowCount = source.metadata().rowCount();
        if (rowCount > Integer.MAX_VALUE) throw new IllegalArgumentException("global factor methods require fewer than 2^31 features");
        List<String> featureIds = new ArrayList<>((int) rowCount);
        List<double[]> rows = new ArrayList<>((int) rowCount);
        try (NumericBlockReader reader = source.open(order)) {
            for (NumericBlock block; (block = reader.read(1024)) != null;)
                for (var row : block.rows()) { featureIds.add(row.id()); rows.add(row.values()); }
        }
        return new MatrixInput(featureIds, selectedSamples, rows.toArray(double[][]::new),
            table, alignedRows);
    }

    private record MatrixInput(List<String> featureIds, List<String> sampleIds,
            double[][] values, DelimitedData phenotype, List<String[]> phenotypeRows) {
        double[][] design(List<String> columns, boolean intercept) {
            if (columns.isEmpty() && !intercept) return new double[sampleIds.size()][0];
            if (columns.isEmpty()) {
                double[][] result = new double[sampleIds.size()][1];
                for (double[] row : result) row[0] = 1.0;
                return result;
            }
            if (!columns.isEmpty() && phenotype == null)
                throw new IllegalArgumentException("design columns require --pheno");
            double[][] result = new double[sampleIds.size()][columns.size() + (intercept ? 1 : 0)];
            int[] indices = columns.stream().mapToInt(phenotype::column).toArray();
            for (int row = 0; row < result.length; row++) {
                int out = 0; if (intercept) result[row][out++] = 1.0;
                for (int index : indices) {
                    try { result[row][out++] = Double.parseDouble(phenotypeRows.get(row)[index].trim()); }
                    catch (NumberFormatException failure) { throw new IllegalArgumentException("design value is not numeric for sample " + sampleIds.get(row), failure); }
                }
            }
            return result;
        }
        String[] strings(String column) { int index = phenotype.column(column); String[] result = new String[sampleIds.size()]; for (int i = 0; i < result.length; i++) result[i] = phenotypeRows.get(i)[index].trim(); return result; }
    }

    private static double[][] residualizeForPca(double[][] data, double[][] design) {
        double[][] q = qr(design), result = copy(data);
        for (int f = 0; f < data.length; f++) for (int k = 0; k < q[0].length; k++) { double c = 0; for (int j = 0; j < q.length; j++) c += data[f][j] * q[j][k]; for (int j = 0; j < q.length; j++) result[f][j] -= c * q[j][k]; }
        return result;
    }
    private static double[][] removeFromOriginal(double[][] data, double[][] factors) { return residualizeForPca(data, factors); }
    private static double[][] qr(double[][] x) { int n=x.length,p=x[0].length,r=0; double[][] q=new double[n][p]; for(int c=0;c<p;c++){double[] v=new double[n];for(int i=0;i<n;i++)v[i]=x[i][c];for(int k=0;k<r;k++){double z=0;for(int i=0;i<n;i++)z+=q[i][k]*v[i];for(int i=0;i<n;i++)v[i]-=z*q[i][k];}double norm=0;for(double z:v)norm+=z*z;norm=Math.sqrt(norm);if(norm>1e-12){for(int i=0;i<n;i++)q[i][r]=v[i]/norm;r++;}}double[][] out=new double[n][r];for(int i=0;i<n;i++)System.arraycopy(q[i],0,out[i],0,r);return out; }
    private static double[][] copy(double[][] x) { double[][] result=new double[x.length][];for(int i=0;i<x.length;i++)result[i]=x[i].clone();return result; }

    private static List<Path> factorOutputs(Options o, boolean loadings,
            boolean weights, boolean selection) {
        List<Path> result = new ArrayList<>(); result.add(suffix(o.output, ".factors.tsv"));
        if (loadings) result.add(suffix(o.output, ".loadings.tsv"));
        if (weights) result.add(suffix(o.output, ".feature-weights.tsv"));
        if (selection) result.add(suffix(o.output, ".factor-selection.tsv"));
        if (o.writeAdjusted) result.add(suffix(o.output, ".adjusted.tsv"));
        result.add(suffix(o.output, ".manifest.tsv")); return result;
    }
    private static Path suffix(Path prefix, String suffix) { return Path.of(prefix.toString() + suffix); }
    private static void requireFresh(List<Path> outputs, boolean overwrite) throws IOException { for (Path path : outputs) { Path parent=path.toAbsolutePath().getParent(); if(parent!=null)Files.createDirectories(parent); if(Files.exists(path)&&!overwrite)throw new IOException("output exists; use --overwrite: "+path); } }
    private static StandardOpenOption[] writeOptions(boolean overwrite) { return overwrite ? new StandardOpenOption[]{StandardOpenOption.CREATE,StandardOpenOption.TRUNCATE_EXISTING} : new StandardOpenOption[]{StandardOpenOption.CREATE_NEW}; }
    private static void writeFactors(Path path, List<String> ids, double[][] matrix, String label) throws IOException { writeMatrix(path,"sample_id",ids,matrix,label,null,true); }
    private static void writeFeatureMatrix(Path path, List<String> ids, double[][] matrix, String label) throws IOException { writeMatrix(path,"feature_id",ids,matrix,label,null,true); }
    private static void writeFeatureMatrix(Path path,List<String> ids,double[][] matrix,String label,List<String> columns) throws IOException { writeMatrix(path,"feature_id",ids,matrix,label,columns,true); }
    private static void writeMatrix(Path path,String idHeader,List<String> ids,double[][] matrix,String label,List<String> columns,boolean overwrite) throws IOException { StringBuilder text=new StringBuilder(idHeader); for(int j=0;j<matrix[0].length;j++)text.append('\t').append(columns==null?label+(j+1):columns.get(j));text.append('\n');for(int i=0;i<matrix.length;i++){text.append(ids.get(i));for(double value:matrix[i])text.append('\t').append(Double.isNaN(value)?"NA":Double.toString(value));text.append('\n');}Files.writeString(path,text,StandardCharsets.UTF_8,writeOptions(overwrite)); }
    private static void writeWeights(Path path,List<String> ids,double[] weights) throws IOException { StringBuilder text=new StringBuilder("feature_id\tweight\n");for(int i=0;i<ids.size();i++)text.append(ids.get(i)).append('\t').append(weights[i]).append('\n');Files.writeString(path,text,StandardCharsets.UTF_8,StandardOpenOption.CREATE,StandardOpenOption.TRUNCATE_EXISTING); }
    private static void writeSelection(Path path,List<SvaResult.FactorSelection> rows) throws IOException {StringBuilder text=new StringBuilder("factors\tf_ratio\tsucceeded\n");for(var row:rows)text.append(row.factors()).append('\t').append(row.fRatio()).append('\t').append(row.succeeded()).append('\n');Files.writeString(path,text,StandardCharsets.UTF_8,StandardOpenOption.CREATE,StandardOpenOption.TRUNCATE_EXISTING);}
    private static void writeManifest(Path path,Options o,MatrixInput input,int factors,int iterations,boolean converged,BackendPolicy backend,double[] variance) throws IOException {StringBuilder text=new StringBuilder("key\tvalue\n");text.append("method\t").append(o.method).append('\n').append("samples\t").append(input.sampleIds.size()).append('\n').append("features\t").append(input.featureIds.size()).append('\n').append("factors\t").append(factors).append('\n').append("iterations\t").append(iterations).append('\n').append("converged\t").append(converged).append('\n').append("backend_policy\t").append(backend.name().toLowerCase(Locale.ROOT)).append('\n').append("seed\t").append(o.seed).append('\n').append("full_design\t").append(String.join(",",o.fullColumns)).append('\n').append("null_design\t").append(String.join(",",o.nullColumns)).append('\n').append("variance_explained\t").append(Arrays.toString(variance)).append('\n');Files.writeString(path,text,StandardCharsets.UTF_8,StandardOpenOption.CREATE,StandardOpenOption.TRUNCATE_EXISTING);}
    private static String label(String method) { return switch(method){case "pca"->"PC";case "peer"->"PEER";default->"SV";}; }

    private static String help(String command) {
        return command.equals("batch-adjust") ? """
            Usage: java -jar jlinalg.jar batch-adjust --method combat --omics MATRIX
              --pheno TABLE --id COLUMN --batch COLUMN [--preserve c1,c2] --out PREFIX
              [--nonparametric] [--mean-only] [--reference-batch VALUE] [--overwrite]
            Input MATRIX is feature-by-sample CSV/TSV. ComBat expects normalized continuous data.
            """ : """
            Usage: java -jar jlinalg.jar confounders --method pca|sva|autosva|peer
              --omics MATRIX [--pheno TABLE --id COLUMN]
              [--full-design c1,c2] [--null-design c1,c2]
              --factors N|auto --out PREFIX [--write-adjusted] [--overwrite]
            PCA: --center/--no-center, --scale. SVA auto uses BE permutations.
            AutoSVA: auto searches the preserved-signal F ratio. PEER supports
              --max-iterations, --tolerance, --variance-tolerance, --add-mean, --seed.
            """;
    }

    private static final class Options {
        boolean help, overwrite, writeAdjusted, center=true, scale, addMean, parametric=true, meanOnly;
        String method, idColumn, batchColumn, referenceBatch, backend="preferred";
        Path omics, pheno, output; Integer factorCount; int iterations=20, permutations=20, maximumFactors, increment=5, maximumIterations=1000;
        long seed=1L; double mrse=0.001, tolerance=1e-3, varianceTolerance=1e-8;
        List<String> fullColumns=List.of(), nullColumns=List.of();
        int requireFactorCount(){if(factorCount==null)throw new IllegalArgumentException("this method requires numeric --factors");return factorCount;}
        static Options parse(String command,String[] args){Options o=new Options();for(int i=0;i<args.length;i++){String a=args[i];switch(a){case "--help","-h"->o.help=true;case "--overwrite"->o.overwrite=true;case "--write-adjusted"->o.writeAdjusted=true;case "--center"->o.center=true;case "--no-center"->o.center=false;case "--scale"->o.scale=true;case "--add-mean"->o.addMean=true;case "--nonparametric"->o.parametric=false;case "--mean-only"->o.meanOnly=true;case "--method"->o.method=value(args,++i,a).toLowerCase(Locale.ROOT);case "--omics"->o.omics=Path.of(value(args,++i,a));case "--pheno"->o.pheno=Path.of(value(args,++i,a));case "--id"->o.idColumn=value(args,++i,a);case "--batch"->o.batchColumn=value(args,++i,a);case "--reference-batch"->o.referenceBatch=value(args,++i,a);case "--out"->o.output=Path.of(value(args,++i,a));case "--backend"->o.backend=value(args,++i,a);case "--full-design","--preserve","--covariates"->o.fullColumns=columns(value(args,++i,a));case "--null-design"->o.nullColumns=columns(value(args,++i,a));case "--factors"->{String v=value(args,++i,a);o.factorCount=v.equalsIgnoreCase("auto")?null:positive(v,a);}case "--iterations"->o.iterations=positive(value(args,++i,a),a);case "--permutations"->o.permutations=positive(value(args,++i,a),a);case "--max-factors"->o.maximumFactors=positive(value(args,++i,a),a);case "--increment"->o.increment=positive(value(args,++i,a),a);case "--max-iterations"->o.maximumIterations=positive(value(args,++i,a),a);case "--seed"->o.seed=Long.parseLong(value(args,++i,a));case "--mrse"->o.mrse=nonnegative(value(args,++i,a),a);case "--tolerance"->o.tolerance=positiveDouble(value(args,++i,a),a);case "--variance-tolerance"->o.varianceTolerance=nonnegative(value(args,++i,a),a);default->throw new IllegalArgumentException("unknown option: "+a);}}if(o.help)return o;if(o.method==null)throw new IllegalArgumentException("--method is required");if(o.omics==null)throw new IllegalArgumentException("--omics is required");if(o.output==null)throw new IllegalArgumentException("--out is required");if(command.equals("batch-adjust")&&!o.method.equals("combat"))throw new IllegalArgumentException("batch-adjust requires --method combat");return o;}
        private static String value(String[] a,int i,String option){if(i>=a.length||a[i].isBlank())throw new IllegalArgumentException(option+" requires a value");return a[i];}
        private static List<String> columns(String x){List<String> r=Arrays.stream(x.split(",",-1)).map(String::trim).filter(v->!v.isEmpty()).toList();if(r.isEmpty())throw new IllegalArgumentException("design column list is empty");return r;}
        private static int positive(String x,String option){int r=Integer.parseInt(x);if(r<1)throw new IllegalArgumentException(option+" must be positive");return r;}
        private static double nonnegative(String x,String option){double r=Double.parseDouble(x);if(!(r>=0)||!Double.isFinite(r))throw new IllegalArgumentException(option+" must be finite and nonnegative");return r;}
        private static double positiveDouble(String x,String option){double r=Double.parseDouble(x);if(!(r>0)||!Double.isFinite(r))throw new IllegalArgumentException(option+" must be positive and finite");return r;}
    }
}
