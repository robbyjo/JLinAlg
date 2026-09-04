/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.ByteArrayOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;

/** Small PLINK 1 reference writer shared by CLI tests. */
final class PlinkTestReference {
    private PlinkTestReference() { }

    static void write(Path prefix, Variant... variants) throws Exception {
        StringBuilder fam = new StringBuilder();
        for (int sample = 0; sample < variants[0].genotypes().length; sample++)
            fam.append('F').append(sample).append(" I").append(sample)
                .append(" 0 0 0 -9\n");
        Files.writeString(Path.of(prefix + ".fam"), fam);
        StringBuilder bim = new StringBuilder();
        ByteArrayOutputStream bed = new ByteArrayOutputStream();
        bed.write(new byte[] {(byte) 0x6c, (byte) 0x1b, 0x01});
        for (Variant variant : variants) {
            bim.append(variant.chromosome()).append('\t')
                .append(variant.id()).append("\t0\t")
                .append(variant.position()).append("\tA\tG\n");
            for (int start = 0; start < variant.genotypes().length; start += 4) {
                int packed = 0;
                for (int within = 0; within < 4
                        && start + within < variant.genotypes().length; within++) {
                    int genotype = variant.genotypes()[start + within];
                    int code = switch (genotype) {
                        case 0 -> 3;
                        case 1 -> 2;
                        case 2 -> 0;
                        default -> 1;
                    };
                    packed |= code << (2 * within);
                }
                bed.write(packed);
            }
        }
        Files.writeString(Path.of(prefix + ".bim"), bim);
        Files.write(Path.of(prefix + ".bed"), bed.toByteArray());
    }

    record Variant(String chromosome, String id, long position,
            int[] genotypes) { }
}
