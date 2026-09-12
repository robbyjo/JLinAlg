/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.raremetal;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.tabix.*;
import org.jlinalg.pipeline.*;
import org.jlinalg.settest.LinearSetTestNullModel;

/** Quantitative-trait unrelated-sample Gaussian score exporter. Uses the fitted
 * null model's residual variance and original phenotype units. Retains only
 * the dosage residuals within the requested forward covariance window.
 * Missing dosages are mean imputed on a fixed complete-case phenotype sample.
 */
public final class GaussianScoreWriter {
    private GaussianScoreWriter() { }
    public static void write(VariantSource source,int[] sampleOrder,LinearSetTestNullModel model,
            long window,int maximumVariants,String genomeBuild,Path scoreFile,Path covarianceFile)throws IOException {
        write(source,sampleOrder,model,window,maximumVariants,genomeBuild,null,"original",scoreFile,covarianceFile);
    }
    /** Export a declared trait identity and units for strict model compatibility checks.
     * Units describe the supplied phenotype; no transformation or conversion occurs. */
    public static void write(VariantSource source,int[] sampleOrder,LinearSetTestNullModel model,
            long window,int maximumVariants,String genomeBuild,String traitId,String traitUnits,Path scoreFile,Path covarianceFile)throws IOException {
        double variance=model.residualVariance();
        if(!(variance>0)||!Double.isFinite(variance))throw new IllegalArgumentException("positive null residual variance required");
        double[] py=model.responseResiduals();for(int i=0;i<py.length;i++)py[i]/=variance;
        write(source,sampleOrder,model.observations(),py,g->{
            double[] pg=model.residualize(new double[][]{g})[0];
            for(int i=0;i<pg.length;i++)pg[i]/=variance;return pg;
        },variance,"unrelated-Gaussian",window,maximumVariants,genomeBuild,traitId,traitUnits,scoreFile,covarianceFile);
    }
    /** Export related-sample scores using the fitted REML projection P.
     * U=G'Py and V=G'PG; trait units remain original, disk covariance is V/N. */
    public static void write(VariantSource source,int[] sampleOrder,org.jlinalg.gwas.RemlAssociationScanner model,
            long window,int maximumVariants,String genomeBuild,Path scoreFile,Path covarianceFile)throws IOException {
        write(source,sampleOrder,model,window,maximumVariants,genomeBuild,null,"original",scoreFile,covarianceFile);
    }
    /** Export Gaussian REML scores with a declared trait identity and units. */
    public static void write(VariantSource source,int[] sampleOrder,org.jlinalg.gwas.RemlAssociationScanner model,
            long window,int maximumVariants,String genomeBuild,String traitId,String traitUnits,Path scoreFile,Path covarianceFile)throws IOException {
        write(source,sampleOrder,model.observations(),model.projectedResponse(),g->model.project(g,1),
            Double.NaN,"related-Gaussian-REML",window,maximumVariants,genomeBuild,traitId,traitUnits,scoreFile,covarianceFile);
    }
    private static void write(VariantSource source,int[] sampleOrder,int n,double[] residual,
            java.util.function.UnaryOperator<double[]> projection,double variance,String nullModel,
            long window,int maximumVariants,String genomeBuild,String traitId,String traitUnits,Path scoreFile,Path covarianceFile)throws IOException {
        if(window<1||maximumVariants<1||genomeBuild==null||genomeBuild.isBlank())throw new IllegalArgumentException("invalid score export options");
        if(!genomeBuild.equals(genomeBuild.trim())||genomeBuild.chars().anyMatch(Character::isISOControl))throw new IllegalArgumentException("genome build must not contain control characters or outer whitespace");
        String header="##ProgramName=JLinAlg (RareMetalWorker-compatible)\n##Version=4.14-format\n##AnalyzedSamples="+n+"\n##GenomeBuild="+genomeBuild
            +"\n##CovarianceWindow="+window+"\n##ResidualVariance="+variance+"\n"
            +ScoreModelMetadata.gaussianHeaders(traitId,traitUnits,nullModel)+"##MissingDosages=mean-imputed\n##HWE=exact-hard-calls-only\n";
        try(IndexedWriter scores=new IndexedWriter(scoreFile);IndexedWriter cov=new IndexedWriter(covarianceFile);
                VariantBlockReader reader=source.open(sampleOrder)) {
            scores.raw(header+"#CHROM\tPOS\tREF\tALT\tN_INFORMATIVE\tFOUNDER_AF\tALL_AF\tINFORMATIVE_ALT_AC\tCALL_RATE\tHWE_PVALUE\tN_REF\tN_HET\tN_ALT\tU_STAT\tSQRT_V_STAT\tALT_EFFSIZE\tPVALUE\n");
            cov.raw(header+"#CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES\n");
            ArrayDeque<Row> queue=new ArrayDeque<>();VariantRecord previous=null;
            for(VariantBlock block;(block=reader.read(1))!=null;) {
                VariantRecord variant=block.variants().get(0);
                String chr=RareMetalStudy.chromosome(variant.chromosome());long pos=variant.position();
                if(pos<1||pos>Integer.MAX_VALUE||!variant.referenceAllele().matches("[ACGT]+")||!variant.alternateAllele().matches("[ACGT]+"))throw new IOException("only normalized biallelic genomic variants are supported");
                if(previous!=null) {
                    int cmp=RareMetalStudy.compareChromosome(RareMetalStudy.chromosome(previous.chromosome()),chr);
                    if(cmp>0||cmp==0&&previous.position()>=pos)throw new IOException("variants must have unique sorted biallelic positions");
                }
                previous=variant;
                while(!queue.isEmpty()&&(!queue.peek().chromosome.equals(chr)||pos-queue.peek().position>window))flush(queue,cov,n);
                if(queue.size()>=maximumVariants)throw new IOException("covariance window exceeds maximum retained variants");
                double[] dosage=variant.dosages();if(dosage.length!=n)throw new IOException("sample alignment differs from null model");
                int called=0,homRef=0,heterozygous=0,homAlt=0;double sum=0;
                for(double d:dosage)if(Double.isFinite(d)){if(d<0||d>2)throw new IOException("dosages must lie in [0,2]");called++;sum+=d;if(d==0)homRef++;else if(d==1)heterozygous++;else if(d==2)homAlt++;}
                double hwe=org.jlinalg.genetics.HardyWeinberg.calculate(dosage,null,-1).pValue();
                double mean=called==0?0:sum/called;
                for(int i=0;i<n;i++)if(!Double.isFinite(dosage[i]))dosage[i]=mean;
                double[] projected=projection.apply(dosage);
                double u=dot(dosage,residual),v=dot(dosage,projected);
                // Exact monomorphic/imputed-constant variants have zero information.
                if(called==0||mean==0||mean==2||v<1e-24*n){u=0;v=0;Arrays.fill(projected,0);}
                String af=called==0?"NA":Double.toString(mean/2);
                double p=v>0?org.jlinalg.settest.SummarySetTests.singleVariant("variant",u,v).pValue():Double.NaN;
                String line=chr+"\t"+pos+"\t"+variant.referenceAllele()+"\t"+variant.alternateAllele()+"\t"+n+"\t"
                    +af+"\t"+af+"\t"+sum+"\t"+(double)called/n+"\t"+(Double.isNaN(hwe)?"NA":hwe)+"\t"+homRef+"\t"+heterozygous+"\t"+homAlt+"\t"+u+"\t"+Math.sqrt(v)
                    +"\t"+(v>0?Double.toString(u/v):"NA")+"\t"+(Double.isNaN(p)?"NA":Double.toString(p))+"\n";
                scores.row(chr,pos,line);queue.add(new Row(chr,pos,dosage,projected));
            }
            while(!queue.isEmpty())flush(queue,cov,n);
        }
    }
    private static void flush(ArrayDeque<Row> queue,IndexedWriter writer,int n)throws IOException {
        Row first=queue.peek();StringBuilder positions=new StringBuilder(),cov=new StringBuilder();
        for(Row row:queue){if(!positions.isEmpty()){positions.append(',');cov.append(',');}positions.append(row.position);cov.append(dot(first.dosage,row.residual)/n);}
        writer.row(first.chromosome,first.position,first.chromosome+"\t"+first.position+"\t"+positions+"\t"+cov+"\n");queue.remove();
    }
    private static double dot(double[] a,double[] b){double s=0;for(int i=0;i<a.length;i++)s+=a[i]*b[i];return s;}
    private record Row(String chromosome,long position,double[] dosage,double[] residual) { }
    private static final class IndexedWriter implements AutoCloseable {
        private final BlockCompressedOutputStream stream;private final TabixIndexCreator index;
        private final Path path;
        IndexedWriter(Path path)throws IOException {
            if(Files.exists(path)||Files.exists(Path.of(path+".tbi")))throw new IOException("output exists: "+path);
            this.path=path;stream=new BlockCompressedOutputStream(path.toFile());index=new TabixIndexCreator(new TabixFormat(0,1,2,2,'#',0));
        }
        void raw(String text)throws IOException{stream.write(text.getBytes(StandardCharsets.UTF_8));}
        void row(String chr,long pos,String text)throws IOException {
            index.addFeature(new Feature(){public String getContig(){return chr;}public int getStart(){return (int)pos;}public int getEnd(){return (int)pos;}},stream.getFilePointer());raw(text);
        }
        public void close()throws IOException {long end=stream.getFilePointer();stream.close();index.finalizeIndex(end).write(Path.of(path+".tbi"));}
    }
}
