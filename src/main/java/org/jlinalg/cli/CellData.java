/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import static org.jlinalg.cli.FollowupSupport.*;

/** Explicit sparse quantified-count contract shared by single-cell and tissue workflows. */
final class CellData {
    final DelimitedData cells,samples,features;
    final List<String> ids,genes;
    final Map<String,Integer> sampleIndex,featureIndex;
    final List<Map<Integer,Double>> counts;
    final double[] libraries,mitochondrial;
    final int[] detected;
    private CellData(DelimitedData cells,DelimitedData samples,DelimitedData features,
            List<String> ids,List<String> genes,Map<String,Integer> sampleIndex,
            List<Map<Integer,Double>> counts,double[] libraries,double[] mitochondrial,int[] detected) {
        this.cells=cells;this.samples=samples;this.features=features;this.ids=ids;this.genes=genes;
        this.sampleIndex=sampleIndex;this.featureIndex=new HashMap<>();for(int j=0;j<genes.size();j++)featureIndex.put(genes.get(j),j);
        this.counts=counts;this.libraries=libraries;this.mitochondrial=mitochondrial;this.detected=detected;
    }
    static Map<String,Integer> unique(DelimitedData data,String column) {
        int c=data.column(column);Map<String,Integer> result=new LinkedHashMap<>();
        for(int i=0;i<data.rows().size();i++) {
            String value=data.rows().get(i)[c];
            if(value.isBlank() || value.chars().anyMatch(Character::isISOControl) || result.put(value,i)!=null)throw new IllegalArgumentException("Duplicate/blank/invalid "+column);
        }
        return result;
    }
    static CellData read(Map<String,String> o) throws IOException {
        DelimitedData cells=bounded(Path.of(required(o,"cells"))),samples=bounded(Path.of(required(o,"samples"))),features=bounded(Path.of(required(o,"features")));
        for(DelimitedData table:List.of(cells,samples,features))for(String[] row:table.rows())for(String field:row)
            if(field.chars().anyMatch(Character::isISOControl))throw new IllegalArgumentException("Metadata fields may not contain control characters");
        Map<String,Integer> index=unique(cells,"obs_id"),si=unique(samples,"sample_id"),fi=unique(features,"feature_id");
        if(index.size()>100000 || fi.size()>50000 || si.size()>10000)throw new IllegalArgumentException("Limit: 100000 observations, 50000 features, 10000 samples");
        for(String col:List.of("donor_id","condition","batch","source","assay","visit")) {
            int c=samples.column(col);for(String[] row:samples.rows())if(row[c].isBlank())throw new IllegalArgumentException("Blank sample "+col);
        }
        for(String[] row:samples.rows())if(!Set.of("scRNA","spatialRNA").contains(row[samples.column("assay")]))
            throw new IllegalArgumentException("This count workflow supports assay scRNA or spatialRNA; protein/ATAC require distinct models");
        Set<String> used=new HashSet<>();
        for(String[] row:cells.rows()) {
            String sample=row[cells.column("sample_id")];
            if(!si.containsKey(sample))throw new IllegalArgumentException("Unknown cell sample_id: "+sample);
            used.add(sample);
            if(row[cells.column("cell_type")].isBlank())throw new IllegalArgumentException("Use explicit unknown cell_type instead of a blank annotation");
        }
        if(used.size()!=si.size())throw new IllegalArgumentException("Sample sheet contains samples without observations");
        boolean[] mito=new boolean[fi.size()];int mc=features.column("mitochondrial");
        for(int j=0;j<mito.length;j++) {
            String v=features.rows().get(j)[mc];
            if(!Set.of("true","false").contains(v))throw new IllegalArgumentException("mitochondrial must be true/false");
            mito[j]=Boolean.parseBoolean(v);
        }
        List<Map<Integer,Double>> values=new ArrayList<>();for(int i=0;i<index.size();i++)values.add(new TreeMap<>());
        double[] lib=new double[index.size()],mt=new double[index.size()];int[] detected=new int[index.size()];
        Path path=Path.of(required(o,"counts"));long entries=0;
        try(BufferedReader in=Files.newBufferedReader(path)) {
            String first=in.readLine();if(!"obs_id\tfeature_id\tcount".equals(first))throw new IllegalArgumentException("Sparse counts require TSV header obs_id, feature_id, count in that order");
            for(String line;(line=in.readLine())!=null;) {
                if(++entries>2_000_000)throw new IllegalArgumentException("Sparse input exceeds 2 million entries");
                List<String> fields=DelimitedData.parse(line,'\t',entries+1,path);
                if(fields.size()!=3)throw new IllegalArgumentException("Sparse count row must have three fields");
                Integer i=index.get(fields.get(0)),j=fi.get(fields.get(1));
                if(i==null || j==null)throw new IllegalArgumentException("Unknown observation/feature in counts");
                double value=finite(fields.get(2),"count");
                if(value<0 || value!=Math.rint(value) || value>1e12)throw new IllegalArgumentException("Raw nonnegative integer RNA counts <=1e12 required");
                if(values.get(i).put(j,value)!=null)throw new IllegalArgumentException("Duplicate observation-feature count");
                lib[i]+=value;if(lib[i]>9e15)throw new IllegalArgumentException("Library count exceeds exact integer bound");
                if(value>0)detected[i]++;if(mito[j])mt[i]+=value;
            }
        }
        return new CellData(cells,samples,features,List.copyOf(index.keySet()),List.copyOf(fi.keySet()),si,values,lib,mt,detected);
    }
    static DelimitedData bounded(Path path) throws IOException {
        if(Files.size(path)>32_000_000)throw new IllegalArgumentException("Metadata table exceeds 32 MB");
        return DelimitedData.read(path);
    }
    String cell(int i,String col){return cells.rows().get(i)[cells.column(col)];}
    String sample(String id,String col){return samples.rows().get(sampleIndex.get(id))[samples.column(col)];}
    double value(int i,int j){return counts.get(i).getOrDefault(j,0.0);}
    double normalized(int i,int j){return libraries[i]>0?Math.log1p(10000*value(i,j)/libraries[i]):0;}
    static void checkOutputSize(long rows) {
        if(rows>250_000)throw new IllegalArgumentException("Result/summary table exceeds 250000 rows; select a smaller hypothesis family");
    }
    void validateGroup(String column) {
        int c=cells.column(column);for(String[] row:cells.rows())if(row[c].isBlank())throw new IllegalArgumentException("Blank population/domain label in "+column);
    }

    boolean[] qc(Map<String,String> o,Path dir,Map<String,Object> manifest) throws IOException {
        int minFeatures=integer(o,"min-features",1,0,50000);double minCounts=number(o,"min-counts",1),maxMito=number(o,"max-mito",1);
        if(minCounts<0 || maxMito<0 || maxMito>1)throw new IllegalArgumentException("Invalid QC thresholds");
        boolean[] keep=new boolean[ids.size()];List<String[]> rows=new ArrayList<>();
        Map<String,int[]> tally=new LinkedHashMap<>();
        for(int i=0;i<ids.size();i++) {
            List<String> reason=new ArrayList<>();
            if(libraries[i]==0)reason.add("zero_library");
            if(libraries[i]<minCounts)reason.add("low_library");
            if(detected[i]<minFeatures)reason.add("low_features");
            if(libraries[i]>0 && mitochondrial[i]/libraries[i]>maxMito)reason.add("high_mitochondrial");
            if(cells.header().contains("exclude")) {
                String flag=cell(i,"exclude");if(!Set.of("true","false").contains(flag))throw new IllegalArgumentException("exclude must be true/false");
                if(flag.equals("true"))reason.add("imported_exclusion");
            }
            keep[i]=reason.isEmpty();int[] totals=tally.computeIfAbsent(cell(i,"sample_id"),s->new int[2]);totals[0]++;if(keep[i])totals[1]++;
            rows.add(new String[]{ids.get(i),cell(i,"sample_id"),cell(i,"cell_type"),Double.toString(libraries[i]),Integer.toString(detected[i]),libraries[i]>0?Double.toString(mitochondrial[i]/libraries[i]):"NA",Boolean.toString(keep[i]),String.join(";",reason)});
        }
        table(dir.resolve("cell-qc.tsv"),List.of("obs_id","sample_id","cell_type","library_count","detected_features","mitochondrial_fraction","retained","exclusion_reason"),rows);
        List<String[]> summary=new ArrayList<>();for(var e:tally.entrySet())summary.add(new String[]{e.getKey(),sample(e.getKey(),"donor_id"),Integer.toString(e.getValue()[0]),Integer.toString(e.getValue()[1])});
        table(dir.resolve("sample-qc.tsv"),List.of("sample_id","donor_id","input_observations","retained_observations"),summary);
        table(dir.resolve("input-cells.tsv"),cells.header(),cells.rows());table(dir.resolve("input-samples.tsv"),samples.header(),samples.rows());table(dir.resolve("input-features.tsv"),features.header(),features.rows());
        manifest.put("qc",Map.of("min_features",minFeatures,"min_counts",minCounts,"max_mito",maxMito));
        manifest.put("missing_evidence",List.of("Empty droplets, ambient RNA and doublets are not estimated; supply audited exclusions upstream.","Manual cell_type labels are imported, not validated reference mapping."));
        return keep;
    }
    record Aggregate(String id,String sample,String type,int cells,double library,double[] values) { }
    List<Aggregate> aggregate(boolean[] keep,String groupColumn) {
        Map<List<String>,double[]> sums=new LinkedHashMap<>();Map<List<String>,Integer> n=new LinkedHashMap<>();
        for(int i=0;i<ids.size();i++)if(keep[i]) {
            List<String> key=List.of(cell(i,"sample_id"),cell(i,groupColumn));
            if(key.get(1).isBlank())throw new IllegalArgumentException("Blank aggregate group");
            if(!sums.containsKey(key) && (long)(sums.size()+1)*genes.size()>10_000_000)throw new IllegalArgumentException("Aggregates exceed 10 million dense entries");
            double[] x=sums.computeIfAbsent(key,k->new double[genes.size()]);n.merge(key,1,Integer::sum);
            for(var v:counts.get(i).entrySet()){x[v.getKey()]+=v.getValue();if(x[v.getKey()]>9e15)throw new IllegalArgumentException("Aggregate exceeds exact integer bound");}
        }
        List<Aggregate> out=new ArrayList<>();for(var e:sums.entrySet()) {
            double total=Arrays.stream(e.getValue()).sum();if(total>9e15)throw new IllegalArgumentException("Aggregate library exceeds exact integer bound");
            out.add(new Aggregate("pb"+(out.size()+1),e.getKey().get(0),e.getKey().get(1),n.get(e.getKey()),total,e.getValue()));
        }
        return out;
    }
    void writeAggregates(List<Aggregate> aggregates,Path dir,int minCells) throws IOException {
        List<String[]> info=new ArrayList<>();
        for(Aggregate a:aggregates)info.add(new String[]{a.id,a.sample,sample(a.sample,"donor_id"),sample(a.sample,"condition"),a.type,Integer.toString(a.cells),Double.toString(a.library),a.cells>=minCells?"eligible":"insufficient_cells"});
        table(dir.resolve("pseudobulk-samples.tsv"),List.of("aggregate_id","sample_id","donor_id","condition","cell_type","cells","library_count","status"),info);
        try(var out=Files.newBufferedWriter(dir.resolve("pseudobulk-counts.tsv"))) {
            List<String> header=new ArrayList<>(List.of("feature_id"));aggregates.forEach(a->header.add(a.id));row(out,header.toArray(String[]::new));
            for(int j=0;j<genes.size();j++){String[] r=new String[aggregates.size()+1];r[0]=genes.get(j);for(int i=0;i<aggregates.size();i++)r[i+1]=Double.toString(aggregates.get(i).values[j]);row(out,r);}
        }
    }
}
