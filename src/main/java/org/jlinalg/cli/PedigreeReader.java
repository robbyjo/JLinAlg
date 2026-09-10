/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.pedigree.PedigreeIndividual;

/** Reads pedigree ancestry and derives exact inbreeding without dense A. */
final class PedigreeReader {
    private PedigreeReader() { }

    static Loaded read(
            Path path, String individualColumn, String sireColumn,
            String damColumn, String familyColumn) throws IOException {
        DelimitedData table = DelimitedData.read(path);
        int individual = table.column(individualColumn);
        int sire = table.column(sireColumn);
        int dam = table.column(damColumn);
        int family = familyColumn == null ? -1 : table.column(familyColumn);
        List<PedigreeIndividual> entries =
            new ArrayList<>(table.rows().size());
        Map<String, Integer> ids = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            String familyId = family < 0 ? null
                : required(row[family], familyColumn);
            String id = qualify(familyId,
                required(row[individual], individualColumn));
            if (ids.putIfAbsent(id, ids.size()) != null)
                throw new IllegalArgumentException(
                    "duplicate pedigree individual: " + id);
            entries.add(new PedigreeIndividual(id,
                qualify(familyId, parent(row[sire])),
                qualify(familyId, parent(row[dam]))));
        }
        return new Loaded(entries, inbreeding(entries, ids));
    }

    private static double[] inbreeding(
            List<PedigreeIndividual> entries,
            Map<String, Integer> index) {
        int size = entries.size();
        int[] pending = new int[size];
        List<List<Integer>> children = new ArrayList<>(size);
        for (int member = 0; member < size; member++)
            children.add(new ArrayList<>());
        for (int child = 0; child < size; child++) {
            PedigreeIndividual value = entries.get(child);
            for (String parent : List.of(
                    value.sireId() == null ? "" : value.sireId(),
                    value.damId() == null ? "" : value.damId())) {
                if (parent.isEmpty()) continue;
                Integer parentIndex = index.get(parent);
                if (parentIndex == null)
                    throw new IllegalArgumentException(
                        "parent " + parent + " of " + value.id()
                            + " is absent from the pedigree");
                pending[child]++;
                children.get(parentIndex).add(child);
            }
        }
        ArrayDeque<Integer> ready = new ArrayDeque<>();
        for (int member = 0; member < size; member++)
            if (pending[member] == 0) ready.add(member);
        int[] order = new int[size];
        int[] rank = new int[size];
        int count = 0;
        while (!ready.isEmpty()) {
            int parent = ready.remove();
            rank[parent] = count;
            order[count++] = parent;
            for (int child : children.get(parent))
                if (--pending[child] == 0) ready.add(child);
        }
        if (count != size)
            throw new IllegalArgumentException(
                "pedigree contains an ancestry cycle");

        double[] result = new double[size];
        Map<Long, Double> relationships = new HashMap<>();
        for (int position = 0; position < size; position++) {
            int member = order[position];
            PedigreeIndividual value = entries.get(member);
            int sire = parentIndex(value.sireId(), index);
            int dam = parentIndex(value.damId(), index);
            result[member] = sire < 0 || dam < 0 ? 0.0
                : 0.5 * relationship(sire, dam, entries, index,
                    rank, result, relationships);
        }
        return result;
    }

    private static double relationship(
            int first, int second, List<PedigreeIndividual> entries,
            Map<String, Integer> index, int[] rank, double[] inbreeding,
            Map<Long, Double> cache) {
        if (first == second) return 1.0 + inbreeding[first];
        int later = rank[first] > rank[second] ? first : second;
        int earlier = later == first ? second : first;
        long key = ((long) Math.min(first, second) << 32)
            | Integer.toUnsignedLong(Math.max(first, second));
        Double retained = cache.get(key);
        if (retained != null) return retained;
        PedigreeIndividual value = entries.get(later);
        int sire = parentIndex(value.sireId(), index);
        int dam = parentIndex(value.damId(), index);
        double result = 0.5 * (
            (sire < 0 ? 0.0 : relationship(sire, earlier, entries,
                index, rank, inbreeding, cache))
            + (dam < 0 ? 0.0 : relationship(dam, earlier, entries,
                index, rank, inbreeding, cache)));
        cache.put(key, result);
        return result;
    }

    private static int parentIndex(
            String id, Map<String, Integer> index) {
        if (id == null) return -1;
        Integer result = index.get(id);
        if (result == null)
            throw new IllegalArgumentException(
                "parent is absent from the pedigree: " + id);
        return result;
    }

    private static String required(String value, String column) {
        String result = value == null ? "" : value.trim();
        if (result.isEmpty())
            throw new IllegalArgumentException(
                "missing pedigree value in column " + column);
        return result;
    }

    private static String parent(String value) {
        if (value == null) return null;
        String result = value.trim();
        String lower = result.toLowerCase(Locale.ROOT);
        return result.isEmpty() || lower.equals("na") || lower.equals("nan")
                || lower.equals("null") || result.equals("0")
                || result.equals(".") || result.equals("-9")
            ? null : result;
    }

    private static String qualify(String family, String individual) {
        if (individual == null) return null;
        return family == null ? individual : family + ":" + individual;
    }

    record Loaded(
            List<PedigreeIndividual> individuals,
            double[] inbreedingCoefficients) {
        Loaded {
            individuals = List.copyOf(individuals);
            inbreedingCoefficients = inbreedingCoefficients.clone();
        }
        @Override public double[] inbreedingCoefficients() {
            return inbreedingCoefficients.clone();
        }
    }
}
