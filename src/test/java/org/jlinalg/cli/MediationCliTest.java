/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MediationCliTest {
    @TempDir Path temporaryDirectory;

    @Test
    void olsMediationOmitsOneCommonMissingRowAndWritesFiveEffects()
            throws Exception {
        Path input = temporaryDirectory.resolve("mediation.csv");
        Path output = temporaryDirectory.resolve("effects.csv");
        StringBuilder data = new StringBuilder("Y,X,M,age\n");
        for (int row = 0; row < 40; row++) {
            double x = ((row * 7) % 17) - 8.0;
            double age = 30 + row;
            double mediator = 1 + 0.6 * x + 0.02 * age
                + ((row * 7) % 5 - 2) * 0.08;
            double outcome = 2 + 0.25 * x + 0.9 * mediator
                + 0.01 * age + ((row * 11) % 7 - 3) * 0.06;
            data.append(outcome).append(',').append(x).append(',')
                .append(row == 9 ? "NA" : mediator).append(',')
                .append(age).append('\n');
        }
        Files.writeString(input, data);
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream console = new PrintStream(bytes);

        int status = JLinAlgCli.run(new String[] {
            "mediation", "--input", input.toString(),
            "--outcome", "Y", "--treatment", "X", "--mediator", "M",
            "--covariates", "age", "--out", output.toString()
        }, console, console);

        assertEquals(0, status, bytes.toString());
        DelimitedData effects = DelimitedData.read(output);
        assertEquals(5, effects.rows().size());
        assertTrue(effects.rows().stream().anyMatch(row ->
            row[effects.column("effect")].equals("indirect")));
        assertTrue(bytes.toString().contains(
            "Analysis samples: 39 (input=40, missing omitted=1)"));
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("model=ols"));
    }

    @Test
    void groupedMediationUsesReml() throws Exception {
        Path input = temporaryDirectory.resolve("grouped.tsv");
        Path output = temporaryDirectory.resolve("grouped-effects.tsv");
        StringBuilder data = new StringBuilder("Y\tX\tM\tage\tgroup\n");
        for (int row = 0; row < 60; row++) {
            int group = row % 12;
            double x = (row % 10) - 4.5;
            double age = 40 + (row % 17);
            double groupEffect = (group - 5.5) * 0.07;
            double mediator = 1 + 0.5 * x + 0.01 * age + groupEffect
                + ((row * 3) % 7 - 3) * 0.04;
            double outcome = 2 + 0.2 * x + 0.8 * mediator
                + 0.01 * age + groupEffect
                + ((row * 5) % 11 - 5) * 0.03;
            data.append(outcome).append('\t').append(x).append('\t')
                .append(mediator).append('\t').append(age).append('\t')
                .append("G").append(group).append('\n');
        }
        Files.writeString(input, data);

        int status = JLinAlgCli.run(new String[] {
            "mediation", "--input", input.toString(),
            "--outcome", "Y", "--treatment", "X", "--mediator", "M",
            "--covariates", "age", "--group", "group",
            "--backend", "cpu", "--out", output.toString()
        });

        assertEquals(0, status);
        assertEquals(6, Files.readAllLines(output).size());
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("model=grouped-reml"));
    }

    @Test
    void pedigreeMediationRetainsAndReportsRepeatedMissingSubject()
            throws Exception {
        Path input = temporaryDirectory.resolve("pedigree-mediation.tsv");
        Path pedigree = temporaryDirectory.resolve("pedigree.tsv");
        Path output = temporaryDirectory.resolve("pedigree-effects.tsv");
        StringBuilder data =
            new StringBuilder("Y\tX\tM\tage\tsubject\n");
        for (int row = 0; row < 80; row++) {
            int subject = row % 20;
            double x = ((row * 7) % 17) - 8.0;
            double age = 30 + (row % 23);
            double subjectEffect = (subject - 9.5) * 0.05;
            double mediator = 1 + 0.5 * x + 0.01 * age
                + subjectEffect + ((row * 3) % 7 - 3) * 0.04;
            double outcome = 2 + 0.2 * x + 0.8 * mediator
                + 0.01 * age + subjectEffect
                + ((row * 5) % 11 - 5) * 0.03;
            data.append(outcome).append('\t').append(x).append('\t')
                .append(mediator).append('\t').append(age).append('\t')
                .append("S").append(subject).append('\n');
        }
        Files.writeString(input, data);
        StringBuilder pedigreeData =
            new StringBuilder("family\tmember\tsire\tdam\n");
        for (int subject = 0; subject < 19; subject++) {
            pedigreeData.append("F").append(subject).append("\tS")
                .append(subject).append("\t0\t0\n");
        }
        Files.writeString(pedigree, pedigreeData);
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream console = new PrintStream(bytes);

        int status = JLinAlgCli.run(new String[] {
            "mediation", "--input", input.toString(),
            "--outcome", "Y", "--treatment", "X", "--mediator", "M",
            "--covariates", "age", "--individual-id", "subject",
            "--pedigree", pedigree.toString(), "--pedigree-id", "member",
            "--sire-id", "sire", "--dam-id", "dam",
            "--pedigree-family-id", "family",
            "--backend", "cpu", "--out", output.toString()
        }, console, console);

        assertEquals(0, status, bytes.toString());
        assertTrue(bytes.toString().contains(
            "Pedigree singletons: 1 families across 4 observations"));
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("model=pedigree-reml"));
        assertTrue(log.contains("pedigree_file_members=19"));
        assertTrue(log.contains(
            "pedigree_file_observations_matched=76"));
        assertTrue(log.contains(
            "pedigree_unqualified_aliases_resolved=76"));
        assertTrue(log.contains("pedigree_singletons_added=1"));
        assertTrue(log.contains("pedigree_singleton_observations=4"));
    }

    @Test
    void pedigreeMediationRejectsAnAllSingletonIdMismatch()
            throws Exception {
        Path input = temporaryDirectory.resolve("mismatched-mediation.tsv");
        Path pedigree = temporaryDirectory.resolve("mismatched-pedigree.tsv");
        Path output = temporaryDirectory.resolve("mismatched-effects.tsv");
        Files.writeString(input,
            "Y\tX\tM\tsubject\n"
            + "1\t0\t0.5\tU1\n2\t1\t1.2\tU2\n3\t2\t2.1\tU3\n");
        Files.writeString(pedigree,
            "family\tmember\tsire\tdam\nF1\tP1\t0\t0\n");
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream console = new PrintStream(bytes);

        int status = JLinAlgCli.run(new String[] {
            "mediation", "--input", input.toString(),
            "--outcome", "Y", "--treatment", "X", "--mediator", "M",
            "--individual-id", "subject", "--pedigree", pedigree.toString(),
            "--pedigree-id", "member", "--sire-id", "sire",
            "--dam-id", "dam", "--pedigree-family-id", "family",
            "--out", output.toString()
        }, console, console);

        assertEquals(2, status);
        assertTrue(bytes.toString().contains(
            "no input observations match pedigree members; check "
                + "--individual-id and pedigree ID qualification"));
    }
}
