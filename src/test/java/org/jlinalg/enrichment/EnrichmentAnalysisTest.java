/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.enrichment;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import org.junit.jupiter.api.Test;
import org.jlinalg.stats.StatisticalTests;
import org.jlinalg.stats.Alternative;

class EnrichmentAnalysisTest {
    @Test void exactOraAndZeroHitFamily() {
        Set<String> universe=new TreeSet<>(); for(int i=0;i<20;i++) universe.add("g"+i);
        var results=EnrichmentAnalysis.ora(universe,Set.of("g0","g1","g2","g3"),List.of(
            new GeneSet("hit","Hit",Set.of("g0","g1","g2","g4","g5")),
            new GeneSet("zero","Zero",Set.of("g10","g11"))),1,100,EnrichmentAnalysis.Fdr.BH);
        var r=results.get(0);
        // Exact enumeration: [C(5,3)C(15,1)+C(5,4)]/C(20,4).
        assertEquals(155.0/4845,r.pValue(),1e-14);
        assertEquals(2*r.pValue(),r.adjustedPValue(),1e-14);
        assertEquals(21,r.oddsRatio(),1e-14);
        assertEquals(1,results.get(1).pValue());
        assertEquals(StatisticalTests.fisherExact(new long[][]{{3,1},{2,14}},1,Alternative.GREATER,.95).pValue(),r.pValue(),1e-13);
    }
    @Test void correctionPreservesTinyTailsAndOrder() {
        double[] adjusted=EnrichmentAnalysis.adjustLog(new double[]{Math.log(.2),-1000,0,Math.log(.01)},EnrichmentAnalysis.Fdr.BY);
        assertEquals(-1000+Math.log(4*(1+.5+1.0/3+.25)),adjusted[1],1e-12);
        assertEquals(Math.log(.01*4/2*(1+.5+1.0/3+.25)),adjusted[3],1e-14);
    }
    @Test void biasedUrnAnalyticAndCentralLimits() {
        // Two red, one white, draw twice. P(two red)=4/5 * 2/3.
        assertEquals(8.0/15,Wallenius.upperTail(2,1,2,2,2),1e-15);
        assertEquals(1.0/3,Wallenius.upperTail(2,1,2,2,1),1e-15);
        assertEquals(0,Wallenius.upperTail(2,1,2,2,0));
        assertEquals(1,Wallenius.upperTail(2,1,2,2,Double.POSITIVE_INFINITY));
        assertThrows(IllegalArgumentException.class,()->Wallenius.upperTail(100,100,50,25,2,100));
    }
    @Test void referenceWalleniusAndSmootherFixtures() throws Exception {
        for(String line:Files.readAllLines(Path.of("src/test/resources/enrichment/wallenius.tsv")).subList(1,
                Files.readAllLines(Path.of("src/test/resources/enrichment/wallenius.tsv")).size())) {
            String[] x=line.split("\t");
            double expected=Double.parseDouble(x[5]);
            double actual=Wallenius.upperTail(Integer.parseInt(x[0]),Integer.parseInt(x[1]),Integer.parseInt(x[2]),Integer.parseInt(x[3]),Double.parseDouble(x[4]));
            assertEquals(expected,actual,Math.max(1e-13,expected*2e-7),line);
        }
        double[] input={0,0,1,1,0,1,0,0,0,1,1,1,0,0,1,0};
        double[] fitted=Gsameth.tricube(input);
        List<String> rows=Files.readAllLines(Path.of("src/test/resources/enrichment/smoother.tsv"));
        for(int i=0;i<input.length;i++) assertEquals(Double.parseDouble(rows.get(i+1).split("\t")[1]),fitted[i],2e-14);
    }
    @Test void methylationWeightsDeduplicateAndCap() {
        var w=Gsameth.weights(Map.of("p1",Set.of("A","B"),"p2",Set.of("A"),"p3",Set.of("C")),Set.of("p1","p2"));
        assertEquals(1.5,w.get("A").equivalentProbes());
        assertEquals(1,w.get("A").selectedWeight());
        assertEquals(.5,w.get("B").selectedWeight());
        assertEquals(0,w.get("C").selectedWeight());
    }
    @Test void highPrecisionTailsIncludingBiasedUrnCancellationCase() throws Exception {
        var lines=Files.readAllLines(Path.of("src/test/resources/enrichment/wallenius-precision.tsv"));
        for(String line:lines.subList(1,lines.size())) {
            String[] x=line.split("\t"); double expected=Double.parseDouble(x[5]);
            double actual=Wallenius.upperTail(Integer.parseInt(x[0]),Integer.parseInt(x[1]),Integer.parseInt(x[2]),Integer.parseInt(x[3]),Double.parseDouble(x[4]));
            assertEquals(expected,actual,expected*2e-12,line);
        }
    }
    @Test void repeatedFractionalMappingsDoNotCrossIntegerThresholds() {
        Map<String,Set<String>> mapping=new java.util.TreeMap<>();
        Set<String> genes=new TreeSet<>(); for(int i=0;i<10;i++) genes.add("g"+i);
        for(int i=0;i<10;i++) mapping.put("p"+i,genes);
        var weights=Gsameth.weights(mapping,mapping.keySet());
        weights.values().forEach(w -> { assertEquals(1,w.equivalentProbes()); assertEquals(1,w.selectedWeight()); });
    }
    @Test void completeGsamethReference() throws Exception {
        Map<String,Set<String>> mapping=new java.util.TreeMap<>();
        Set<String> selected=new TreeSet<>();
        for(String line:Files.readAllLines(Path.of("src/test/resources/enrichment/probes.tsv")).subList(1,
                Files.readAllLines(Path.of("src/test/resources/enrichment/probes.tsv")).size())) {
            String[] x=line.split("\t"); mapping.put(x[0],Set.of(x[1].split(";"))); if(x[2].equals("1")) selected.add(x[0]);
        }
        List<GeneSet> sets=new ArrayList<>();
        for(String line:Files.readAllLines(Path.of("src/test/resources/enrichment/sets.gmt"))) {
            String[] x=line.split("\t"); sets.add(new GeneSet(x[0],x[1],Set.of(java.util.Arrays.copyOfRange(x,2,x.length))));
        }
        var results=EnrichmentAnalysis.gsameth(mapping,selected,sets,1,100,EnrichmentAnalysis.Fdr.BH,1000000);
        List<String> fixture=Files.readAllLines(Path.of("src/test/resources/enrichment/gsameth.tsv"));
        for(String line:fixture.subList(1,fixture.size())) {
            String[] x=line.split("\t"); var r=results.stream().filter(a->a.id().equals(x[0])).findFirst().orElseThrow();
            assertEquals(Double.parseDouble(x[1]),r.weightedOverlap(),1e-13);
            assertEquals(Double.parseDouble(x[2]),r.selectionOdds(),1e-13);
            assertEquals(Double.parseDouble(x[3]),r.pValue(),2e-7);
            assertEquals(Double.parseDouble(x[4]),r.adjustedPValue(),2e-7);
        }
    }
}
