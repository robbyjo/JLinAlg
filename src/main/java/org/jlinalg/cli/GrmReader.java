/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.genetics.GenomicRelationshipMatrix;

/** Reads labeled dense and GCTA binary genomic relationship matrices. */
final class GrmReader {
    private GrmReader() { }

    static Loaded read(Path path) throws IOException {
        if (gcta(path)) return readGcta(path);
        return readDense(path);
    }

    private static Loaded readDense(Path path) throws IOException {
        DelimitedData data = DelimitedData.read(path);
        if (data.header().size() < 2)
            throw new IOException(
                "a dense GRM needs a row-ID column and matrix columns: " + path);
        List<String> ids = data.header().subList(1, data.header().size());
        int size = ids.size();
        if (data.rows().size() != size)
            throw new IOException("dense GRM must be square: " + path);
        Map<String, String[]> byId = new HashMap<>();
        for (String[] row : data.rows()) {
            String id = row[0].trim();
            if (id.isEmpty() || byId.put(id, row) != null)
                throw new IOException(
                    "dense GRM row IDs must be unique and nonblank: " + path);
        }
        if (!byId.keySet().equals(new HashSet<>(ids)))
            throw new IOException(
                "dense GRM row and column ID sets differ: " + path);
        double[] matrix = new double[size * size];
        for (int row = 0; row < size; row++) {
            String[] values = byId.get(ids.get(row));
            for (int column = 0; column < size; column++) {
                String token = values[column + 1].trim();
                try {
                    matrix[row * size + column] = Double.parseDouble(token);
                } catch (NumberFormatException exception) {
                    throw new IOException(
                        "invalid dense GRM value for " + ids.get(row) + ","
                            + ids.get(column) + ": " + token,
                        exception);
                }
            }
        }
        try {
            return new Loaded(
                new GenomicRelationshipMatrix(ids, matrix), "dense");
        } catch (IllegalArgumentException exception) {
            throw new IOException("invalid dense GRM: " + exception.getMessage(),
                exception);
        }
    }

    private static Loaded readGcta(Path requested) throws IOException {
        Paths paths = gctaPaths(requested);
        List<String[]> pairs = new ArrayList<>();
        for (String line : Files.readAllLines(
                paths.ids(), StandardCharsets.UTF_8)) {
            if (line.isBlank()) continue;
            String[] fields = line.trim().split("\\s+");
            if (fields.length < 2)
                throw new IOException(
                    "GCTA GRM ID rows require FID and IID: " + paths.ids());
            pairs.add(new String[] {fields[0], fields[1]});
        }
        if (pairs.isEmpty())
            throw new IOException("GCTA GRM ID file is empty: " + paths.ids());
        Set<String> individualIds = new HashSet<>();
        boolean uniqueIndividuals = true;
        for (String[] pair : pairs)
            uniqueIndividuals &= individualIds.add(pair[1]);
        List<String> ids = new ArrayList<>(pairs.size());
        for (String[] pair : pairs)
            ids.add(uniqueIndividuals ? pair[1] : pair[0] + ":" + pair[1]);

        int size = ids.size();
        long values = (long) size * (size + 1) / 2;
        long expectedBytes = Math.multiplyExact(values, Float.BYTES);
        long actualBytes = Files.size(paths.binary());
        if (actualBytes != expectedBytes)
            throw new IOException("GCTA GRM binary size is " + actualBytes
                + " bytes; expected " + expectedBytes + " for " + size
                + " samples: " + paths.binary());
        double[] matrix = new double[Math.multiplyExact(size, size)];
        try (DataInputStream input = new DataInputStream(
                new BufferedInputStream(Files.newInputStream(paths.binary())))) {
            for (int row = 0; row < size; row++)
                for (int column = 0; column <= row; column++) {
                    int bits;
                    try {
                        bits = Integer.reverseBytes(input.readInt());
                    } catch (EOFException exception) {
                        throw new IOException(
                            "truncated GCTA GRM binary: " + paths.binary(),
                            exception);
                    }
                    double value = Float.intBitsToFloat(bits);
                    matrix[row * size + column] = value;
                    matrix[column * size + row] = value;
                }
        }
        try {
            return new Loaded(
                new GenomicRelationshipMatrix(ids, matrix), "gcta-binary");
        } catch (IllegalArgumentException exception) {
            throw new IOException("invalid GCTA GRM: " + exception.getMessage(),
                exception);
        }
    }

    private static boolean gcta(Path path) {
        String name = path.getFileName().toString()
            .toLowerCase(Locale.ROOT);
        return name.endsWith(".grm.bin") || name.endsWith(".grm.id")
            || Files.exists(Path.of(path.toString() + ".grm.bin"));
    }

    private static Paths gctaPaths(Path path) {
        String value = path.toString();
        String lower = value.toLowerCase(Locale.ROOT);
        String prefix;
        if (lower.endsWith(".grm.bin"))
            prefix = value.substring(0, value.length() - ".grm.bin".length());
        else if (lower.endsWith(".grm.id"))
            prefix = value.substring(0, value.length() - ".grm.id".length());
        else prefix = value;
        return new Paths(Path.of(prefix + ".grm.bin"),
            Path.of(prefix + ".grm.id"));
    }

    record Loaded(GenomicRelationshipMatrix matrix, String format) { }
    private record Paths(Path binary, Path ids) { }
}
