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
}
