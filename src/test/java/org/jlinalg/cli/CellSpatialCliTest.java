/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class CellSpatialCliTest {
    @TempDir Path temp;
    @Test void exactlyFittedSampleResponseHasNoSpuriousPValue() {
        var design=org.jlinalg.singlecell.SampleInference.design(new String[]{"1","2","3","4","5","6","7","8"},
            new String[]{"a","a","a","a","b","b","b","b"},"a","b",new double[8][0],List.of(),false);
        String[] row=SingleCellCli.fitRow("type","gene",new double[]{.2,.2,.2,.2,.4,.4,.4,.4},design);
        assertEquals("NA",row[8]);
    }
    private List<String> args(String command,String method,Path out) {
        List<String> a=new ArrayList<>(List.of(command,"--no-config","--method",method,"--out",out.toString()));
        for(String f:List.of("counts","cells","samples","features"))a.addAll(List.of("--"+f,"examples/cell-spatial/"+f+".tsv"));
        if(command.equals("spatial"))a.addAll(List.of("--units","micrometer","--coordinate-system","synthetic-v1","--resolution","cell","--radius","10.1"));
        return a;
    }
    private int run(List<String> a){return JLinAlgCli.run(a.toArray(String[]::new),new PrintStream(OutputStream.nullOutputStream()),new PrintStream(OutputStream.nullOutputStream()));}
    @Test void pseudobulkPreservesRawTotalsAndReportsDonors() throws IOException {
        Path dir=temp.resolve("bulk");assertEquals(0,run(args("single-cell","pseudobulk",dir)));
        DelimitedData raw=DelimitedData.read(Path.of("examples/cell-spatial/counts.tsv")),bulk=DelimitedData.read(dir.resolve("pseudobulk-counts.tsv"));
        double expected=raw.rows().stream().mapToDouble(r->Double.parseDouble(r[2])).sum(),actual=0;
        for(String[] r:bulk.rows())for(int j=1;j<r.length;j++)actual+=Double.parseDouble(r[j]);assertEquals(expected,actual);
        assertEquals(16,DelimitedData.read(dir.resolve("pseudobulk-samples.tsv")).rows().size());
        assertTrue(Files.readString(dir.resolve("manifest.yaml")).contains("counts:"));
        assertTrue(Files.exists(Path.of(dir+".config.yaml")));assertEquals(2,run(args("single-cell","pseudobulk",dir)));
    }
    @Test void abundanceUsesSampleDesignAndCanFitCompletePairs() throws IOException {
        for(boolean paired:List.of(false,true)) {
            Path dir=temp.resolve("abundance"+paired);var a=args("single-cell","abundance",dir);
            a.addAll(List.of("--reference","control","--tested","case"));
            if(paired){a.set(a.indexOf("--samples")+1,"examples/cell-spatial/paired-samples.tsv");a.addAll(List.of("--paired","true"));}
            assertEquals(0,run(a));assertEquals(8,DelimitedData.read(dir.resolve("design.tsv")).rows().size());
            assertEquals(2,DelimitedData.read(dir.resolve("results.tsv")).rows().size());
        }
    }
    @Test void graphNeverConnectsOverlappingSectionsAndTestsAreReproducible() throws IOException {
        Path first=temp.resolve("spatial-a"),second=temp.resolve("spatial-b");
        for(Path dir:List.of(first,second)) {
            var a=args("spatial","autocorrelation",dir);a.addAll(List.of("--permutations","19","--feature-list","examples/cell-spatial/selected-features.tsv"));assertEquals(0,run(a));
            for(String[] r:DelimitedData.read(dir.resolve("edges.tsv")).rows())assertEquals(r[0].substring(0,4),r[1].substring(0,4));
        }
        assertEquals(Files.readString(first.resolve("autocorrelation.tsv")),Files.readString(second.resolve("autocorrelation.tsv")));
        assertTrue(Files.exists(first.resolve("section-1.svg")));
    }
    @Test void badSparseCountsAndUnknownOptionsLeaveNoPublishedOutput() throws IOException {
        Path bad=temp.resolve("bad.tsv");Files.writeString(bad,"obs_id\tfeature_id\tcount\nS1_0_0\tGENE01\t1\nS1_0_0\tGENE01\t2\n");
        var a=args("single-cell","qc",temp.resolve("bad-run"));a.set(a.indexOf("--counts")+1,bad.toString());assertEquals(2,run(a));assertFalse(Files.exists(temp.resolve("bad-run")));
        var option=args("single-cell","qc",temp.resolve("option"));option.addAll(List.of("--paired","true"));assertEquals(2,run(option));
    }
    @Test void rejectMixedSpotNeighborhoodsAndConfoundedDesign() throws IOException {
        var a=args("spatial","neighborhoods",temp.resolve("mixed"));a.set(a.indexOf("--resolution")+1,"spot");assertEquals(2,run(a));
        Path sheet=temp.resolve("samples.tsv");String content=Files.readString(Path.of("examples/cell-spatial/samples.tsv"));
        content=content.replace("control\tB1","control\tB0").replace("case\tB0","case\tB1");Files.writeString(sheet,content);
        a=args("single-cell","abundance",temp.resolve("confounded"));a.set(a.indexOf("--samples")+1,sheet.toString());a.addAll(List.of("--reference","control","--tested","case"));assertEquals(2,run(a));
    }
    @Test void sparseQcRetainsZeroLibraryExclusionAndFeatureUniverse() throws IOException {
        Path file=temp.resolve("sparse.tsv");Files.writeString(file,"obs_id\tfeature_id\tcount\nS1_0_0\tMT-GENE\t10\n");
        var a=args("single-cell","qc",temp.resolve("qc"));a.set(a.indexOf("--counts")+1,file.toString());a.addAll(List.of("--max-mito","0.5"));assertEquals(0,run(a));
        var table=DelimitedData.read(temp.resolve("qc/cell-qc.tsv"));assertEquals(320,table.rows().size());
        assertTrue(table.rows().get(0)[7].contains("high_mitochondrial"));assertTrue(table.rows().get(1)[7].contains("zero_library"));
    }
    @Test void partiallyNestedPairedBatchesKeepOnlyIdentifiableNuisanceBasis() throws IOException {
        Path input=temp.resolve("paired.tsv");var sheet=DelimitedData.read(Path.of("examples/cell-spatial/paired-samples.tsv"));
        String[] batch={"A","A","C","C","B","A","C","C"};List<String[]> rows=new ArrayList<>();
        for(int i=0;i<8;i++){String[] r=sheet.rows().get(i).clone();r[sheet.column("batch")]=batch[i];rows.add(r);}
        FollowupSupport.table(input,sheet.header(),rows);
        var a=args("single-cell","abundance",temp.resolve("nested"));a.set(a.indexOf("--samples")+1,input.toString());a.addAll(List.of("--reference","control","--tested","case","--paired","true"));
        assertEquals(0,run(a));var design=DelimitedData.read(temp.resolve("nested/design.tsv"));
        assertTrue(design.header().contains("batch:B"));assertFalse(design.header().contains("batch:C"));
        for(int i=0;i<8;i++)rows.get(i)[sheet.column("batch")]=i<4?"A":"B";
        FollowupSupport.table(input,sheet.header(),rows);
        a.set(a.indexOf("--out")+1,temp.resolve("confounded-pairs").toString());assertEquals(2,run(a));
    }
    @Test void invalidOptionsAndPathwaysCannotHideBehindUntestablePopulations() throws IOException {
        for(String[] invalid:List.of(new String[]{"--paired","maybe"},new String[]{"--covariates","absent"},new String[]{"--covariates","sample_id"})) {
            Path out=temp.resolve("invalid-"+invalid[1]);var a=args("single-cell","state",out);
            a.addAll(List.of("--reference","control","--tested","case","--min-cells","10000","--rscript","unused"));a.addAll(List.of(invalid));
            assertEquals(2,run(a));assertFalse(Files.exists(out));
        }
        Path bad=temp.resolve("bad-sets.tsv");Files.writeString(bad,"gene_set\tfeature_id\ns\tGENE01\ns\tGENE01\n");
        var a=args("single-cell","pathways",temp.resolve("invalid-pathway"));a.addAll(List.of("--reference","control","--tested","case","--min-cells","10000","--gene-sets",bad.toString()));assertEquals(2,run(a));
    }
    @Test void missingDomainLabelsAreRejectedForAbundance() throws IOException {
        var cells=DelimitedData.read(Path.of("examples/cell-spatial/cells.tsv"));List<String> header=new ArrayList<>(cells.header());header.add("domain");List<String[]> rows=new ArrayList<>();
        for(String[] cell:cells.rows()){String[] row=Arrays.copyOf(cell,cell.length+1);row[cell.length]="";rows.add(row);}
        Path input=temp.resolve("blank-domain.tsv");FollowupSupport.table(input,header,rows);
        var a=args("single-cell","abundance",temp.resolve("domain-run"));a.set(a.indexOf("--cells")+1,input.toString());a.addAll(List.of("--reference","control","--tested","case","--group-column","domain"));assertEquals(2,run(a));
    }
    @Test void resultCardinalityBoundsRejectBeforeAllocation() {
        assertThrows(IllegalArgumentException.class,()->CellData.checkOutputSize(250001));
        CellData.checkOutputSize(250000);
    }
}
