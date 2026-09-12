/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.xwas.*;

/** LDSC, predicted molecular association and common-factor summary GWAS commands. */
final class SummaryXwasCli {
    private SummaryXwasCli() { }
    static int run(String command,String[] args,PrintStream out,PrintStream error) {
        try {
            if(Arrays.asList(args).contains("--help")){out.println(help(command));return 0;}
            switch(command) {
                case "ldsc" -> ldsc(args);
                case "twas", "pwas" -> molecular(args);
                case "genomic-factor" -> factor(args);
                default -> throw new IllegalArgumentException("unknown summary command");
            }
            out.println(command+" complete");return 0;
        } catch(IOException|RuntimeException ex){error.println("jlinalg: "+ex.getMessage());return 2;}
    }
    private static void ldsc(String[] args) throws IOException {
        Map<String,String> o=XwasFiles.options(args,"--input","--traits","--reference-variants","--blocks","--out");
        List<String> traits=XwasFiles.names(XwasFiles.required(o,"--traits"));
        for(String trait:traits)if(trait.contains(":"))throw new IllegalArgumentException("trait names cannot contain ':'");
        DelimitedData table=DelimitedData.read(XwasFiles.path(o,"--input"));
        XwasFiles.indexed(table,"variant");
        int n=table.rows().size(),t=traits.size(),b=Integer.parseInt(o.getOrDefault("--blocks","200"));
        if(b<3||b>n)throw new IllegalArgumentException("--blocks must be between 3 and the variant count");
        double[] ld=new double[n],wld=new double[n];double[][] z=new double[n][t],ns=new double[n][t];int[] blocks=new int[n];
        for(int i=0;i<n;i++) {
            String[] row=table.rows().get(i);
            String ea=row[table.column("ea")],oa=row[table.column("oa")];PredictedOmics.alleleSign(ea,oa,ea,oa);
            ld[i]=XwasFiles.number(table,row,"ld_score");wld[i]=XwasFiles.number(table,row,"weight_ld");
            blocks[i]=(int)((long)i*b/n);
            for(int j=0;j<t;j++){z[i][j]=XwasFiles.number(table,row,"z_"+traits.get(j));ns[i][j]=XwasFiles.number(table,row,"n_"+traits.get(j));}
        }
        var fit=LdScoreRegression.fit(ld,wld,z,ns,blocks,XwasFiles.number(XwasFiles.required(o,"--reference-variants")));
        double[] s=fit.geneticCovariance(),v=fit.samplingCovariance(),intercept=fit.intercepts(),ise=fit.interceptStandardErrors();
        double[] rg=fit.geneticCorrelations(),rgse=fit.correlationStandardErrors();
        List<String> pairs=XwasFiles.pairs(traits);int p=pairs.size(),k=0;
        StringBuilder text=new StringBuilder("trait1\ttrait2\testimate\tse\tp\tintercept\tintercept_se\trg\trg_se\trg_p\tstatus\n");
        for(int i=0;i<t;i++)for(int j=i;j<t;j++,k++) {
            double se=Math.sqrt(v[k*p+k]),rse=rgse[i*t+j];
            text.append(traits.get(i)).append('\t').append(traits.get(j)).append('\t').append(s[i*t+j]).append('\t').append(se)
                .append('\t').append(se>0?SummaryMath.p(s[i*t+j]/se):Double.NaN).append('\t').append(intercept[i*t+j]).append('\t').append(ise[i*t+j])
                .append('\t').append(rg[i*t+j]).append('\t').append(rse).append('\t').append(i!=j&&rse>0?SummaryMath.p(rg[i*t+j]/rse):Double.NaN)
                .append('\t').append(!Double.isFinite(rg[i*t+j])?"nonpositive_h2":!Double.isFinite(rse)?"invalid_rg_delete_h2":Math.abs(rg[i*t+j])>1?"rg_outside_unit_interval":"ok").append('\n');
        }
        Path output=XwasFiles.path(o,"--out");Map<Path,String> files=new LinkedHashMap<>();
        files.put(output,text.toString());files.put(Path.of(output+".S.tsv"),XwasFiles.matrix(traits,s));files.put(Path.of(output+".V.tsv"),XwasFiles.matrix(pairs,v));
        files.put(Path.of(output+".intercepts.tsv"),XwasFiles.matrix(traits,intercept));
        StringBuilder deleted=new StringBuilder("block\t"+String.join("\t",pairs)+"\n");
        double[][] d=fit.deleteValues();for(int r=0;r<b;r++){deleted.append(r);for(double value:d[r])deleted.append('\t').append(value);deleted.append('\n');}
        files.put(Path.of(output+".deletes.tsv"),deleted.toString());
        files.put(Path.of(output+".metadata.tsv"),"key\tvalue\nmethod\tone-step unpartitioned LDSC; two weight updates; final-weight block jackknife\nvariants\t"+n+"\nblocks\t"+b+"\nreference_variants\t"+o.get("--reference-variants")+"\nscale\tobserved standardized phenotype\nalignment\tall trait z columns must use row ea on one forward strand and genome build\n");
        XwasFiles.publish(files);
    }
    private static void molecular(String[] args) throws IOException {
        Map<String,String> o=XwasFiles.options(args,"--gwas","--weights","--reference","--ld","--joint","--out");
        if(o.containsKey("--joint")&&!Set.of("true","false").contains(o.get("--joint")))throw new IllegalArgumentException("--joint must be true or false");
        var reference=DelimitedData.read(XwasFiles.path(o,"--reference"));
        Map<String,String[]> refs=XwasFiles.indexed(reference,"variant");List<String> ids=new ArrayList<>(refs.keySet());
        var gwas=DelimitedData.read(XwasFiles.path(o,"--gwas"));Map<String,String[]> stats=XwasFiles.indexed(gwas,"variant");
        var weights=DelimitedData.read(XwasFiles.path(o,"--weights"));
        double[] ld=XwasFiles.matrix(XwasFiles.path(o,"--ld"),ids),z=new double[ids.size()],sd=new double[ids.size()];
        Map<String,Integer> index=new HashMap<>();
        for(int i=0;i<ids.size();i++) {
            String id=ids.get(i);index.put(id,i);String[] ref=refs.get(id),row=stats.get(id);
            if(row==null)throw new IllegalArgumentException("missing GWAS variant: "+id);
            z[i]=XwasFiles.number(gwas,row,"z")*PredictedOmics.alleleSign(row[gwas.column("ea")],row[gwas.column("oa")],ref[reference.column("ea")],ref[reference.column("oa")]);
            sd[i]=XwasFiles.number(reference,ref,"sd");
        }
        Map<String,double[]> models=new LinkedHashMap<>();Map<String,Set<String>> seen=new HashMap<>();
        for(String[] row:weights.rows()) {
            String model=XwasFiles.id(row[weights.column("model")]),id=XwasFiles.id(row[weights.column("variant")]);
            Integer i=index.get(id);if(i==null)throw new IllegalArgumentException("model variant absent from reference: "+id);
            if(!seen.computeIfAbsent(model,k->new HashSet<>()).add(id))throw new IllegalArgumentException("duplicate model variant: "+model+"/"+id);
            String[] ref=refs.get(id);
            double sign=PredictedOmics.alleleSign(row[weights.column("ea")],row[weights.column("oa")],ref[reference.column("ea")],ref[reference.column("oa")]);
            models.computeIfAbsent(model,k->new double[ids.size()])[i]=sign*XwasFiles.number(weights,row,"weight");
        }
        StringBuilder text=new StringBuilder("model\tvariants\tz\tp\tlog_p\tpredicted_variance\n");
        for(var entry:models.entrySet()) {
            var fit=PredictedOmics.test(z,entry.getValue(),sd,ld);
            text.append(entry.getKey()).append('\t').append(seen.get(entry.getKey()).size()).append('\t').append(fit.z()).append('\t').append(fit.pValue()).append('\t').append(fit.logPValue()).append('\t').append(fit.predictedVariance()).append('\n');
        }
        Path output=XwasFiles.path(o,"--out");Map<Path,String> files=new LinkedHashMap<>();files.put(output,text.toString());
        if(Boolean.parseBoolean(o.getOrDefault("--joint","false"))) {
            var joint=PredictedOmics.joint(z,models.values().toArray(double[][]::new),sd,ld);
            files.put(Path.of(output+".joint.tsv"),"chi_square\tdf\tp\n"+joint.chiSquare()+"\t"+joint.degreesOfFreedom()+"\t"+joint.pValue()+"\n");
        }
        files.put(Path.of(output+".metadata.tsv"),"key\tvalue\nweight_scale\traw effect-allele dosage\nld_scale\tcorrelation in reference ea orientation\ncoverage\tcomplete reference and model coverage required\nstrand\texact forward-strand match or swap; no complement inference\n");
        XwasFiles.publish(files);
    }
    private static void factor(String[] args) throws IOException {
        var o=XwasFiles.options(args,"--s","--v","--traits","--gwas","--sampling-correlation","--out");
        var traits=XwasFiles.names(XwasFiles.required(o,"--traits"));int t=traits.size();
        double[] s=XwasFiles.matrix(XwasFiles.path(o,"--s"),traits),v=XwasFiles.matrix(XwasFiles.path(o,"--v"),XwasFiles.pairs(traits));
        var fit=GenomicFactor.fit(s,t,v);double[] load=fit.loadings(),res=fit.residualVariances(),pc=fit.parameterCovariance();
        Path output=XwasFiles.path(o,"--out");Map<Path,String> files=new LinkedHashMap<>();
        StringBuilder text=new StringBuilder("trait\tloading\tloading_se\tresidual_variance\tresidual_se\n");
        for(int i=0;i<t;i++)text.append(traits.get(i)).append('\t').append(load[i]).append('\t').append(Math.sqrt(pc[i*2*t+i])).append('\t').append(res[i]).append('\t').append(Math.sqrt(pc[(t+i)*2*t+t+i])).append('\n');
        files.put(output,text.toString());
        files.put(Path.of(output+".fit.tsv"),"chi_square\tdf\tp\titerations\n"+fit.chiSquare()+"\t"+fit.degreesOfFreedom()+"\t"+fit.pValue()+"\t"+fit.iterations()+"\n");
        List<String> params=new ArrayList<>();for(String trait:traits)params.add("loading:"+trait);for(String trait:traits)params.add("residual:"+trait);
        files.put(Path.of(output+".parameter-covariance.tsv"),XwasFiles.matrix(params,pc));
        if(o.containsKey("--gwas")!=o.containsKey("--sampling-correlation"))throw new IllegalArgumentException("--gwas and --sampling-correlation must be supplied together");
        if(o.containsKey("--gwas")) {
            double[] correlation=XwasFiles.matrix(XwasFiles.path(o,"--sampling-correlation"),traits);PredictedOmics.correlation(correlation,t);
            var gwas=DelimitedData.read(XwasFiles.path(o,"--gwas"));XwasFiles.indexed(gwas,"variant");
            StringBuilder results=new StringBuilder("variant\tea\toa\tbeta\tse\tz\tp\tq_snp\tq_df\tq_p\n");
            for(String[] row:gwas.rows()) {
                double[] beta=new double[t],se=new double[t],c=new double[t*t];
                String ea=row[gwas.column("ea")],oa=row[gwas.column("oa")];PredictedOmics.alleleSign(ea,oa,ea,oa);
                for(int i=0;i<t;i++){beta[i]=XwasFiles.number(gwas,row,"beta_"+traits.get(i));se[i]=XwasFiles.number(gwas,row,"se_"+traits.get(i));if(!(se[i]>0))throw new IllegalArgumentException("GWAS SE must be positive");}
                for(int i=0;i<t;i++)for(int j=0;j<t;j++)c[i*t+j]=correlation[i*t+j]*se[i]*se[j];
                var a=GenomicFactor.associate(beta,c,load);
                results.append(row[gwas.column("variant")]).append('\t').append(ea).append('\t').append(oa).append('\t').append(a.beta()).append('\t').append(a.standardError()).append('\t').append(a.z()).append('\t').append(a.pValue()).append('\t').append(a.heterogeneity()).append('\t').append(a.heterogeneityDf()).append('\t').append(a.heterogeneityP()).append('\n');
            }
            files.put(Path.of(output+".gwas.tsv"),results.toString());
        }
        files.put(Path.of(output+".metadata.tsv"),"key\tvalue\nmethod\tfull WLS single factor; factor variance fixed at one\nsnp_inference\tconditional on fitted loadings; measurement uncertainty not propagated into SNP SE\nvariance_policy\tpositive residuals and positive definite sampling covariance; no matrix repair\n");
        XwasFiles.publish(files);
    }
    static String help(String command) {
        return switch(command) {
            case "ldsc" -> """
                Usage: jlinalg ldsc --input FILE --traits A,B --reference-variants M --out FILE [--blocks 200]
                Input columns: variant ea oa ld_score weight_ld z_A n_A z_B n_B ...
                All z columns must use the row effect allele; supply sorted, harmonized, QC-filtered summary rows.
                M counts reference variants used to calculate LD scores, not input rows. No LD calculation or munging.
                Produces estimates and OUT.S.tsv, OUT.V.tsv, OUT.intercepts.tsv, OUT.deletes.tsv, OUT.metadata.tsv.
                Observed scale, free intercept, two weight updates, equal contiguous final-weight jackknife blocks.
                """;
            case "twas", "pwas" -> """
                Usage: jlinalg twas|pwas --gwas FILE --weights FILE --reference FILE --ld FILE --out FILE [--joint true]
                GWAS: variant ea oa z. Weights: model variant ea oa weight (raw dosage coefficients).
                Reference: variant ea oa sd (genotype standard deviation). LD: labeled square correlation matrix.
                Complete coverage required; forward-strand exact match/swap only. Joint test requires linearly independent models.
                Models can represent genes, tissues, or proteins. Train weights with score-train or import external weights.
                """;
            default -> """
                Usage: jlinalg genomic-factor --s FILE --v FILE --traits A,B,C --out FILE
                  [--gwas FILE --sampling-correlation FILE]
                S: labeled trait genetic covariance. V: sampling covariance of upper triangle by rows (A:A,A:B,...).
                S/V can come directly from ldsc. Fits a full-WLS single factor with positive residual variances.
                GWAS: variant ea oa beta_A se_A beta_B se_B ... already aligned to the common row ea.
                Sampling correlation: labeled trait matrix accounting for overlap; identity requires independent errors.
                SNP effects must share the phenotype scale of S. SNP tests condition on fitted loadings.
                """;
        };
    }
}
