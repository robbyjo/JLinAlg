/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PlinkBedLdClumperTest {
    @TempDir Path temporaryDirectory;

    @Test
    void slowPhaseCannotAbortValidFortyOneFounderClump() throws Exception {
        checkPhaseTable(new int[][]{{4,0,2},{3,21,4},{3,3,1}}, .01, "slow");
    }

    @Test
    void midpointStationaryPointCannotHideClumpAboveDefaultThreshold() throws Exception {
        checkPhaseTable(new int[][]{{2,2,1},{1,13,3},{1,3,0}}, .001, "twin");
    }

    private void checkPhaseTable(int[][] counts, double threshold, String key) throws Exception {
        int n=java.util.Arrays.stream(counts).flatMapToInt(java.util.Arrays::stream).sum();
        int[] a=new int[n],b=new int[n];int index=0;
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<counts[i][j];k++) {a[index]=i;b[index++]=j;}
        Path prefix=temporaryDirectory.resolve(key);
        writeReference(prefix,new Variant("1","a",100,a),new Variant("1","b",200,b));
        var result=PlinkBedLdClumper.clump(prefix,List.of(new LdClumpCandidate("a",1e-9,"g"),
            new LdClumpCandidate("b",2e-9,"g")),new LdClumpOptions(10000,threshold,1));
        assertEquals(1,result.retained().size());assertEquals(1,result.exclusions().size());
        assertEquals("a",result.retained().get(0).variantId());
        var reference=new java.util.Properties();
        try(var in=getClass().getResourceAsStream("/r-reference/genetic-final-fixes.properties")) {reference.load(in);}
        double expected=Double.parseDouble(reference.getProperty("ld."+key+".rSquared"));
        assertEquals(expected,result.exclusions().get(0).rSquared(),expected*2e-10);
    }

    @Test
    void escapesIndependenceStationaryPointInUnphasedLikelihood() throws Exception {
        int[] a = new int[104], b = new int[104];
        java.util.Arrays.fill(a, 1); java.util.Arrays.fill(b, 1);
        a[0]=0; b[0]=0; a[1]=0; b[1]=2; a[2]=2; b[2]=0; a[3]=2; b[3]=2;
        Path prefix = temporaryDirectory.resolve("symmetric");
        writeReference(prefix,new Variant("1","a",100,a),new Variant("1","b",200,b));
        var result=PlinkBedLdClumper.clump(prefix,List.of(new LdClumpCandidate("a",1e-9,"g"),
            new LdClumpCandidate("b",2e-9,"g")),new LdClumpOptions(10000,.5,1));
        assertEquals(1,result.retained().size());
        // With four corner homozygotes and m double heterozygotes, writing
        // f00=f11=t, f01=f10=.5-t gives r^2=(m-4)/(m+4) at either maximum.
        assertEquals(96.0/104,result.exclusions().get(0).rSquared(),2e-11);
    }

    @Test
    void retainsLowestPVariantAndDropsMissingReferenceVariants() throws Exception {
        Path prefix = temporaryDirectory.resolve("reference");
        writeReference(prefix,
            new Variant("1", "rs1", 10_000,
                new int[] {0, 0, 1, 1, 2, 2, 0, 2}),
            new Variant("1", "rs2", 12_000,
                new int[] {0, 0, 1, 1, 2, 2, 0, 2}),
            new Variant("1", "rs3", 20_000_000,
                new int[] {0, 1, 2, 0, 1, 2, 1, 0}));
        List<LdClumpCandidate> candidates = List.of(
            new LdClumpCandidate("rs1", 2e-9, "BMI"),
            new LdClumpCandidate("rs2", 1e-9, "BMI"),
            new LdClumpCandidate("rs3", 3e-9, "BMI"),
            new LdClumpCandidate("rsMissing", 4e-9, "BMI"));

        LdClumpResult result = PlinkBedLdClumper.clump(
            prefix, candidates, LdClumpOptions.defaults());

        assertEquals(List.of("rs2", "rs3"), result.retained().stream()
            .map(LdClumpCandidate::variantId).toList());
        LdClumpExclusion linked = result.exclusions().stream()
            .filter(value -> value.candidate().variantId().equals("rs1"))
            .findFirst().orElseThrow();
        assertEquals(LdClumpExclusionReason.IN_LINKAGE_DISEQUILIBRIUM,
            linked.reason());
        assertEquals("rs2", linked.indexVariantId());
        assertEquals(1.0, linked.rSquared(), 1e-12);
        assertTrue(result.exclusions().stream().anyMatch(value ->
            value.candidate().variantId().equals("rsMissing")
                && value.reason()
                    == LdClumpExclusionReason.ABSENT_FROM_REFERENCE));
    }

    @Test
    void matchesIeugwasrSingleVariantGroupPassThrough() throws Exception {
        Path prefix = temporaryDirectory.resolve("reference");
        writeReference(prefix, new Variant("1", "rs1", 10_000,
            new int[] {0, 0, 1, 1, 2, 2, 0, 2}));
        LdClumpCandidate absent = new LdClumpCandidate(
            "not-in-reference", 1.0, "single");

        LdClumpResult result = PlinkBedLdClumper.clump(prefix,
            List.of(absent), new LdClumpOptions(10_000, 0.001, 0.01));

        assertEquals(List.of(absent), result.retained());
        assertTrue(result.exclusions().isEmpty());
    }

    private static void writeReference(Path prefix, Variant... variants)
            throws Exception {
        StringBuilder fam = new StringBuilder();
        for (int sample = 0; sample < variants[0].genotypes().length; sample++)
            fam.append("F").append(sample).append(" I").append(sample)
                .append(" 0 0 0 -9\n");
        Files.writeString(Path.of(prefix + ".fam"), fam);
        StringBuilder bim = new StringBuilder();
        ByteArrayOutputStream bed = new ByteArrayOutputStream();
        bed.write(new byte[] {(byte) 0x6c, (byte) 0x1b, 0x01});
        for (Variant variant : variants) {
            bim.append(variant.chromosome()).append('\t')
                .append(variant.id()).append("\t0\t")
                .append(variant.position()).append("\tA\tG\n");
            writePacked(bed, variant.genotypes());
        }
        Files.writeString(Path.of(prefix + ".bim"), bim);
        Files.write(Path.of(prefix + ".bed"), bed.toByteArray());
    }

    private static void writePacked(ByteArrayOutputStream output,
            int[] genotypes) {
        for (int start = 0; start < genotypes.length; start += 4) {
            int packed = 0;
            for (int within = 0; within < 4
                    && start + within < genotypes.length; within++) {
                int genotype = genotypes[start + within];
                int code = switch (genotype) {
                    case 0 -> 3;
                    case 1 -> 2;
                    case 2 -> 0;
                    default -> 1;
                };
                packed |= code << (2 * within);
            }
            output.write(packed);
        }
    }

    private record Variant(String chromosome, String id, long position,
            int[] genotypes) { }
}
