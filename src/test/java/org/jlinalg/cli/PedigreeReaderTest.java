/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PedigreeReaderTest {
    @TempDir Path temporaryDirectory;

    @Test
    void sparseReaderInbreedingMatchesTabularRelationship() throws Exception {
        Path input = temporaryDirectory.resolve("pedigree.tsv");
        Files.writeString(input,
            "id\tsire\tdam\n"
            + "E\tC\tD\n"
            + "A\t0\t0\n"
            + "D\tA\tB\n"
            + "B\t0\t0\n"
            + "C\tA\tB\n");

        PedigreeReader.Loaded loaded =
            PedigreeReader.read(input, "id", "sire", "dam", null);
        double[] expected = Pedigree.of(loaded.individuals())
            .inbreedingCoefficients();

        assertArrayEquals(expected, loaded.inbreedingCoefficients(), 1e-12);
    }

    @Test
    void uniqueRawIdsResolveAndConnectParentsAcrossFamilyLabels()
            throws Exception {
        Path input = temporaryDirectory.resolve("family-pedigree.tsv");
        Files.writeString(input,
            "family\tid\tsire\tdam\n"
            + "F0\tP\t0\t0\n"
            + "F1\tC1\tP\t0\n"
            + "F2\tC2\tP\t0\n");

        PedigreeReader.Loaded loaded = PedigreeReader.read(
            input, "id", "sire", "dam", "family");
        Pedigree pedigree = Pedigree.of(loaded.individuals());

        assertEquals("F0:P", loaded.resolveObservationId("P"));
        assertEquals("F1:C1", loaded.resolveObservationId("C1"));
        assertEquals("F2:C2", loaded.resolveObservationId("F2:C2"));
        assertEquals(0.25, pedigree.relationship("F1:C1", "F2:C2"),
            1e-12);
    }

    @Test
    void repeatedRawIdRequiresQualifiedPhenotypeKey() throws Exception {
        Path input = temporaryDirectory.resolve("repeated-family-id.tsv");
        Files.writeString(input,
            "family\tid\tsire\tdam\n"
            + "F1\tP\t0\t0\n"
            + "F2\tP\t0\t0\n");

        PedigreeReader.Loaded loaded = PedigreeReader.read(
            input, "id", "sire", "dam", "family");

        assertEquals("F1:P", loaded.resolveObservationId("F1:P"));
        IllegalArgumentException failure = assertThrows(
            IllegalArgumentException.class,
            () -> loaded.resolveObservationId("P"));
        assertEquals("ambiguous unqualified pedigree individual P; "
            + "use family:individual in the phenotype matching column",
            failure.getMessage());
    }

    @Test
    void qualifiedParentDisambiguatesRepeatedRawIdAcrossFamilies()
            throws Exception {
        Path input = temporaryDirectory.resolve("qualified-parent.tsv");
        Files.writeString(input,
            "family\tid\tsire\tdam\n"
            + "F1\tP\t0\t0\n"
            + "F2\tP\t0\t0\n"
            + "F3\tC\tF1:P\t0\n");

        PedigreeReader.Loaded loaded = PedigreeReader.read(
            input, "id", "sire", "dam", "family");
        Pedigree pedigree = Pedigree.of(loaded.individuals());

        assertEquals(0.5, pedigree.relationship("F1:P", "F3:C"), 1e-12);
        assertEquals(0.0, pedigree.relationship("F2:P", "F3:C"), 1e-12);
    }

    @Test
    void ambiguousParentAcrossFamiliesFailsClearly() throws Exception {
        Path input = temporaryDirectory.resolve("ambiguous-parent.tsv");
        Files.writeString(input,
            "family\tid\tsire\tdam\n"
            + "F1\tP\t0\t0\n"
            + "F2\tP\t0\t0\n"
            + "F3\tC\tP\t0\n");

        IllegalArgumentException failure = assertThrows(
            IllegalArgumentException.class,
            () -> PedigreeReader.read(
                input, "id", "sire", "dam", "family"));
        assertEquals("parent P of F3:C is ambiguous across pedigree families",
            failure.getMessage());
    }

    @Test
    void rawAndQualifiedNamespaceCollisionFailsClearly() throws Exception {
        Path input = temporaryDirectory.resolve("colliding-ids.tsv");
        Files.writeString(input,
            "family\tid\tsire\tdam\n"
            + "F1\tP\t0\t0\n"
            + "F2\tF1:P\t0\t0\n");

        IllegalArgumentException failure = assertThrows(
            IllegalArgumentException.class,
            () -> PedigreeReader.read(
                input, "id", "sire", "dam", "family"));
        assertEquals("pedigree ID namespace collision: F1:P is both a "
            + "qualified key and a different raw individual ID",
            failure.getMessage());
    }

    @Test
    void multigenerationHierarchyMatchesPedigreemmWithCousinMating()
            throws Exception {
        Properties reference = new Properties();
        try (var reader = Files.newBufferedReader(Path.of(
                "src/test/resources/r-reference/"
                    + "pedigreemm-hierarchy.properties"))) {
            reference.load(reader);
        }
        String[] ids = reference.getProperty("ids").split(",");
        String[] sires = reference.getProperty("sire").split(",");
        String[] dams = reference.getProperty("dam").split(",");
        Path input = temporaryDirectory.resolve("pedigreemm-hierarchy.tsv");
        StringBuilder table = new StringBuilder("family\tid\tsire\tdam\n");
        List<String> canonicalIds = new ArrayList<>();
        int[] order = {10, 5, 11, 0, 8, 2, 9, 6, 3, 7, 1, 4};
        for (int source : order) {
            String family = "ROW" + source;
            canonicalIds.add(family + ":" + ids[source]);
            table.append(family).append('\t').append(ids[source])
                .append('\t').append(sires[source])
                .append('\t').append(dams[source]).append('\n');
        }
        Files.writeString(input, table);

        PedigreeReader.Loaded loaded = PedigreeReader.read(
            input, "id", "sire", "dam", "family");
        Pedigree pedigree = Pedigree.of(loaded.individuals());
        PedigreeRandomEffectTerm sparse =
            PedigreeRandomEffectTerm.ofSparse("pedigree", canonicalIds,
                loaded.individuals(), loaded.inbreedingCoefficients());

        for (int index = 0; index < ids.length; index++)
            assertEquals("ROW" + index + ":" + ids[index],
                loaded.resolveObservationId(ids[index]));
        assertArrayEquals(reorder(numbers(reference, "F"), order),
            loaded.inbreedingCoefficients(), 1e-15);
        assertArrayEquals(reorderSquare(numbers(reference, "A"), order),
            pedigree.relationshipMatrix(), 1e-15);
        assertArrayEquals(reorderSquare(numbers(reference, "Ainv"), order),
            dense(sparse.precision()), 2e-15);
        assertEquals(0.5, pedigree.relationship("ROW5:S1", "ROW6:S2"));
        assertEquals(0.25, pedigree.relationship("ROW5:S1", "ROW7:H1"));
        assertEquals(0.125, pedigree.relationship("ROW8:C1", "ROW9:C2"));
        assertEquals(0.0625, loaded.inbreedingCoefficients()[0], 0.0);
    }

    private static double[] numbers(Properties properties, String key) {
        return Arrays.stream(properties.getProperty(key).split(","))
            .mapToDouble(Double::parseDouble).toArray();
    }

    private static double[] reorder(double[] values, int[] order) {
        double[] result = new double[order.length];
        for (int index = 0; index < order.length; index++)
            result[index] = values[order[index]];
        return result;
    }

    private static double[] reorderSquare(double[] values, int[] order) {
        double[] result = new double[values.length];
        int size = order.length;
        for (int row = 0; row < size; row++) {
            for (int column = 0; column < size; column++) {
                result[row * size + column] =
                    values[order[row] * size + order[column]];
            }
        }
        return result;
    }

    private static double[] dense(SparsePrecisionMatrix matrix) {
        double[] result = new double[
            matrix.dimension() * matrix.dimension()];
        int[] starts = matrix.rowStarts();
        int[] columns = matrix.columnIndices();
        double[] values = matrix.values();
        for (int row = 0; row < matrix.dimension(); row++) {
            for (int index = starts[row]; index < starts[row + 1]; index++) {
                result[row * matrix.dimension() + columns[index]] =
                    values[index];
            }
        }
        return result;
    }
}
