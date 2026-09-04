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
