/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashSet;
import java.util.Set;
import java.util.UUID;
import java.util.function.Function;
import org.jlinalg.pipeline.*;
import org.jlinalg.settest.SetTestScoreState;

/** Cohort-side score exports; only aggregate values are written. Covariance
 * blocks are explicit coverage units, never a claim of zero cross-block LD. */
final class ConditionalGwasExport implements AutoCloseable {
    static final List<String> COLUMNS=List.of("score_variant_key","score_u","score_variance",
        "n_analyzed","n_cases","n_controls","effect_ac_cases","effect_ac_controls",
        "n_called_cases","n_called_controls","n_events","beta_score","se_score",
        "p_score_normal","p_score_calibrated","calibration_method","calibration_status",
        "null_model_id","conditioning_set_id","score_covariance_block");
    private final CliOptions options;
    private final Function<double[][],SetTestScoreState> scorer;
    private final double[] response;
    private final boolean binary;
    private final int events;
    private final String modelId=UUID.randomUUID().toString();
    private final String conditioningId;
    private final ManifestWriter metadata;
    private final BufferedWriter writer;
    private final Path temporary;
    private final List<VariantRecord> retained=new ArrayList<>();
    private final List<double[]> imputed=new ArrayList<>();
    private final List<Double> variances=new ArrayList<>();
    private long block=1,variants;
    private boolean closed;
    private String previousChromosome;
    private long previousPosition;
    private final Set<String> keysAtPosition=new HashSet<>();

    ConditionalGwasExport(CliOptions options,Function<double[][],SetTestScoreState> scorer,
            double[] response,boolean binary,int events,String model,String nullFamily,
            double dispersion,List<String> covariates,List<String> conditions,
            PhenotypeData.BinaryMapping coding) throws IOException {
        this.options=options;this.scorer=scorer;this.response=response.clone();this.binary=binary;this.events=events;
        if(binary && Arrays.stream(response).anyMatch(y->y!=0 && y!=1))
            throw new IllegalArgumentException("binary score export requires individual 0/1 responses; grouped-binomial summaries need a different count schema");
        conditioningId=conditions.isEmpty()?"none":modelId+":conditioning";
        metadata=new ManifestWriter().put("export_schema","jlinalg-conditional-score-v1")
            .put("status","complete").put("version",JLinAlgCli.version()).put("null_model_id",modelId)
            .put("conditioning_set_id",conditioningId).putStrings("conditioning_variants",conditions)
            .put("genome_build",options.scoreGenomeBuild).put("model",model).put("family",nullFamily)
            .put("formula",options.formula).putStrings("null_covariate_columns",covariates)
            .put("case_value",binary?(coding==null?"1":coding.caseValue()):null)
            .put("control_value",binary?(coding==null?"0":coding.controlValue()):null)
            .put("phenotype_transform","as-specified-in-formula").put("effect_allele","ALT")
            .put("genotype_coding","additive dosage").put("missing_genotypes","analysis-sample mean imputation")
            .put("case_control_allele_counts","ALT dosage sums before imputation; called counts accompany each sum")
            .put("analysis_sample","fixed complete phenotype/covariate cases aligned to genotype source")
            .put("n_analyzed",response.length).put("dispersion",Double.isFinite(dispersion)?dispersion:null)
            .put("null_converged",true).put("score_covariance_scale","unstandardized; includes inverse fitted-null dispersion for GLMs")
            .put("variance_method","model-based; nuisance covariates projected out")
            .put("relatedness_adjustment","none").put("ties",model.equals("cox")?options.ties:null)
            .put("tail_calibration","normal approximation only; no SPA or finite-sample guarantee")
            .put("covariance_coverage","complete within each exported block; cross-block covariance unavailable, not zero")
            .put("maximum_block_variants",options.scoreBlockSize)
            .put("conditional_scope",conditions.isEmpty()?"one-step approximation for new conditioning sets":"null refitted with listed variants; new conditioning sets still require refitting")
            .put("covariance_file",options.scoreCovariancePath().toAbsolutePath())
            .put("summary_file",options.output.toAbsolutePath());
        temporary=Path.of(options.scoreCovariancePath()+".partial");
        // A stale success manifest must not survive a failed overwrite run.
        if(options.overwrite)Files.deleteIfExists(options.scoreManifestPath());
        Path parent=temporary.toAbsolutePath().getParent();Files.createDirectories(parent);
        writer=Files.newBufferedWriter(temporary,options.overwrite
            ?new StandardOpenOption[]{StandardOpenOption.CREATE,StandardOpenOption.TRUNCATE_EXISTING}
            :new StandardOpenOption[]{StandardOpenOption.CREATE_NEW});
        writer.write("variant_i\tvariant_j\tscore_covariance\tnull_model_id\tscore_covariance_block\n");
    }

    List<String> summarize(VariantRecord variant) throws IOException {
        String key=key(variant);
        if(previousChromosome!=null) {
            int chromosomeOrder=org.jlinalg.raremetal.RareMetalStudy.compareChromosome(previousChromosome,variant.chromosome());
            if(chromosomeOrder>0 || chromosomeOrder==0 && previousPosition>variant.position())
                throw new IllegalArgumentException("score export requires chromosome/position sorted variants: "+key);
        }
        if(!variant.chromosome().equals(previousChromosome) || previousPosition!=variant.position())keysAtPosition.clear();
        if(!keysAtPosition.add(key))throw new IllegalArgumentException("duplicate score variant key: "+key);
        previousChromosome=variant.chromosome();previousPosition=variant.position();
        if(!retained.isEmpty() && (!retained.get(0).chromosome().equals(variant.chromosome()) || retained.size()>=options.scoreBlockSize))flush();
        double[] raw=variant.dosages(),g=impute(raw);int nc=0,nn=0,calledCases=0,calledControls=0;double ac=0,an=0;
        if(raw.length!=response.length)throw new IllegalArgumentException("summary genotype/phenotype alignment differs");
        if(binary)for(int i=0;i<raw.length;i++) {
            if(response[i]==1){nc++;if(Double.isFinite(raw[i])){calledCases++;ac+=raw[i];}}
            else {nn++;if(Double.isFinite(raw[i])){calledControls++;an+=raw[i];}}
        }
        SetTestScoreState state=scorer.apply(new double[][]{g});double u=state.scores()[0],v=state.information()[0];
        if(!Double.isFinite(u) || !(v>0) || !Double.isFinite(v))throw new IllegalArgumentException("invalid efficient score for "+key);
        if(!Double.isFinite(u/v) || !Double.isFinite(1/Math.sqrt(v)))
            throw new IllegalArgumentException("nonfinite one-step score estimate for "+key);
        retained.add(variant);imputed.add(g);variances.add(v);variants++;
        double z=u/Math.sqrt(v),p=2*jdistlib.Normal.cumulative(-Math.abs(z),0,1,true,false);
        return List.of(key,number(u),number(v),Integer.toString(response.length),binary?Integer.toString(nc):"",
            binary?Integer.toString(nn):"",binary?number(ac):"",binary?number(an):"",
            binary?Integer.toString(calledCases):"",binary?Integer.toString(calledControls):"",events<0?"":Integer.toString(events),
            number(u/v),number(1/Math.sqrt(v)),number(p),"","normal-score","normal_only",
            modelId,conditioningId,Long.toString(block));
    }

    private void flush() throws IOException {
        if(retained.isEmpty())return;
        double[] v=scorer.apply(imputed.toArray(double[][]::new)).information();int m=retained.size();
        for(int i=0;i<m;i++) {
            if(Math.abs(v[i*m+i]-variances.get(i))>1e-8*variances.get(i))throw new IllegalArgumentException("singleton and block score information disagree");
            for(int j=i;j<m;j++) {
                if(!Double.isFinite(v[i*m+j]))throw new IllegalArgumentException("nonfinite score covariance");
                writer.write(key(retained.get(i))+"\t"+key(retained.get(j))+"\t"+number(v[i*m+j])+"\t"+modelId+"\t"+block+"\n");
            }
        }
        retained.clear();imputed.clear();variances.clear();block++;
    }
    void finish() throws IOException {
        flush();close();
        if(options.overwrite)Files.move(temporary,options.scoreCovariancePath(),StandardCopyOption.REPLACE_EXISTING);
        else Files.move(temporary,options.scoreCovariancePath());
        Path manifestTemporary=Path.of(options.scoreManifestPath()+".partial");
        metadata.put("exported_variants",variants).put("covariance_blocks",block-1).write(manifestTemporary);
        if(options.overwrite)Files.move(manifestTemporary,options.scoreManifestPath(),StandardCopyOption.REPLACE_EXISTING);
        else Files.move(manifestTemporary,options.scoreManifestPath());
    }
    @Override public void close() throws IOException {if(!closed){closed=true;writer.close();}}
    static String key(VariantRecord variant) {
        if(!variant.genomic() || !variant.referenceAllele().matches("[ACGT]+") || !variant.alternateAllele().matches("[ACGT]+")
                || !variant.chromosome().matches("[A-Za-z0-9_.-]+"))
            throw new IllegalArgumentException("score export requires genomic coordinates and normalized REF/ALT alleles");
        return variant.chromosome()+":"+variant.position()+":"+variant.referenceAllele()+":"+variant.alternateAllele();
    }
    static double[] impute(double[] raw) {
        double[] g=raw.clone();double sum=0;int n=0;
        for(double value:g)if(Double.isFinite(value)) {
            if(value<0 || value>2)throw new IllegalArgumentException("additive dosage must be in [0,2]");sum+=value;n++;
        } else if(!Double.isNaN(value))throw new IllegalArgumentException("infinite dosage");
        if(n==0)throw new IllegalArgumentException("variant has no called samples");
        for(int i=0;i<g.length;i++)if(Double.isNaN(g[i]))g[i]=sum/n;return g;
    }
    /** A local prepass loads only explicitly requested conditioning variants. */
    static Conditioning conditioning(VariantSource source,List<String> samples,List<String> requested,double[][] x) throws IOException {
        if(requested.isEmpty())return new Conditioning(x,List.of());
        double[][] g=new double[requested.size()][];String[] keys=new String[requested.size()];
        try(VariantBlockReader reader=source.open(SampleAlignment.requireOrder(source.metadata().sampleIds(),samples))) {
            for(VariantBlock block;(block=reader.read(128))!=null;)for(VariantRecord variant:block.variants()) {
                for(int j=0;j<requested.size();j++)if(requested.get(j).equals(variant.id()) || (variant.genomic() && requested.get(j).equals(key(variant)))) {
                    if(g[j]!=null)throw new IllegalArgumentException("ambiguous conditioning variant: "+requested.get(j));
                    keys[j]=key(variant);g[j]=impute(variant.dosages());
                }
            }
        }
        for(int j=0;j<g.length;j++)if(g[j]==null)throw new IllegalArgumentException("conditioning variant not found: "+requested.get(j));
        if(Arrays.stream(keys).distinct().count()!=keys.length)throw new IllegalArgumentException("conditioning aliases resolve to the same variant");
        double[][] expanded=new double[x.length][x[0].length+g.length];
        for(int i=0;i<x.length;i++){System.arraycopy(x[i],0,expanded[i],0,x[i].length);for(int j=0;j<g.length;j++)expanded[i][x[i].length+j]=g[j][i];}
        return new Conditioning(expanded,List.of(keys));
    }
    record Conditioning(double[][] covariates,List<String> keys) { }
    static VariantSource omitConditioning(VariantSource source,List<String> keys) {
        if(keys.isEmpty())return source;
        return new VariantSource() {
            public VariantSourceMetadata metadata(){return source.metadata();}
            public VariantBlockReader open(int[] order) throws IOException {
                VariantBlockReader underlying=source.open(order);
                return new VariantBlockReader() {
                    private long next;
                    public VariantBlock read(int maximum) throws IOException {
                        for(VariantBlock block;(block=underlying.read(maximum))!=null;) {
                            List<VariantRecord> rows=block.variants().stream().filter(v->!keys.contains(key(v))).toList();
                            if(!rows.isEmpty()){long first=next;next+=rows.size();return new VariantBlock(first,rows);}
                        }
                        return null;
                    }
                    public void close() throws IOException {underlying.close();}
                };
            }
        };
    }
    private static String number(double value){return Double.isFinite(value)?Double.toString(value):"";}
}
