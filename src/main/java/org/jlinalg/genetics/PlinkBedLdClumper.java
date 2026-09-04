/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Greedy PLINK-style clumping over a variant-major PLINK 1 BED reference.
 *
 * <p>Pairwise LD is the unphased maximum-likelihood haplotype r-squared used
 * by PLINK's {@code --clump}, evaluated over founders with nonmissing calls.
 * Index variants are considered by ascending p-value independently per
 * group.</p>
 */
public final class PlinkBedLdClumper {
    private static final byte[] BED_HEADER = {
        (byte) 0x6c, (byte) 0x1b, (byte) 0x01};
    private static final int EM_ITERATIONS = 1_000;
    private static final double EM_TOLERANCE = 1e-12;

    private PlinkBedLdClumper() { }

    /**
     * Clumps unique group/variant candidates against a PLINK fileset prefix.
     *
     * @param plinkPrefix path without {@code .bed/.bim/.fam}
     * @param candidates candidate variants in their desired output order
     * @param options clumping thresholds
     * @return retained variants and exclusions
     * @throws IOException when the reference is malformed or unreadable
     */
    public static LdClumpResult clump(Path plinkPrefix,
            List<LdClumpCandidate> candidates, LdClumpOptions options)
            throws IOException {
        if (plinkPrefix == null) throw new NullPointerException("plinkPrefix");
        if (candidates == null) throw new NullPointerException("candidates");
        if (options == null) throw new NullPointerException("options");
        validateUnique(candidates);
        Reference reference = Reference.open(plinkPrefix,
            candidates.stream().map(LdClumpCandidate::variantId)
                .collect(java.util.stream.Collectors.toSet()));
        Map<String, List<IndexedCandidate>> groups = new LinkedHashMap<>();
        for (int index = 0; index < candidates.size(); index++) {
            LdClumpCandidate candidate = candidates.get(index);
            groups.computeIfAbsent(candidate.group(), ignored -> new ArrayList<>())
                .add(new IndexedCandidate(candidate, index));
        }

        Set<Integer> retainedIndexes = new HashSet<>();
        List<IndexedExclusion> exclusions = new ArrayList<>();
        try (FileChannel bed = reference.openBed()) {
            Map<Long, byte[]> genotypes = new HashMap<>();
            for (List<IndexedCandidate> group : groups.values())
                clumpGroup(group, reference, bed, genotypes, options,
                    retainedIndexes, exclusions);
        }

        List<LdClumpCandidate> retained = new ArrayList<>();
        for (int index = 0; index < candidates.size(); index++)
            if (retainedIndexes.contains(index)) retained.add(candidates.get(index));
        exclusions.sort(Comparator.comparingInt(IndexedExclusion::index));
        return new LdClumpResult(retained, exclusions.stream()
            .map(IndexedExclusion::exclusion).toList());
    }

    private static void clumpGroup(List<IndexedCandidate> group,
            Reference reference, FileChannel bed, Map<Long, byte[]> genotypes,
            LdClumpOptions options, Set<Integer> retained,
            List<IndexedExclusion> exclusions) throws IOException {
        // ieugwasr::ld_clump() returns a one-row group without consulting PLINK.
        if (group.size() == 1) {
            retained.add(group.get(0).index());
            return;
        }
        List<IndexedCandidate> ordered = new ArrayList<>();
        for (IndexedCandidate indexed : group) {
            Variant variant = reference.variants().get(
                indexed.candidate().variantId());
            if (variant == null) {
                exclusions.add(new IndexedExclusion(indexed.index(),
                    new LdClumpExclusion(indexed.candidate(),
                        LdClumpExclusionReason.ABSENT_FROM_REFERENCE,
                        null, Double.NaN)));
            } else if (indexed.candidate().pValue()
                    > options.indexPValueThreshold()) {
                exclusions.add(new IndexedExclusion(indexed.index(),
                    new LdClumpExclusion(indexed.candidate(),
                        LdClumpExclusionReason.ABOVE_INDEX_P_VALUE_THRESHOLD,
                        null, Double.NaN)));
            } else ordered.add(indexed);
        }
        ordered.sort(Comparator
            .comparingDouble((IndexedCandidate value) ->
                value.candidate().pValue())
            .thenComparingInt(IndexedCandidate::index));
        Set<String> assigned = new LinkedHashSet<>();
        long windowBases = Math.multiplyExact(
            (long) options.windowKilobases(), 1_000L);
        for (int index = 0; index < ordered.size(); index++) {
            IndexedCandidate lead = ordered.get(index);
            if (!assigned.add(lead.candidate().variantId())) continue;
            retained.add(lead.index());
            Variant leadVariant = reference.variants().get(
                lead.candidate().variantId());
            byte[] leadGenotypes = genotypes.computeIfAbsent(
                leadVariant.bedIndex(), ignored -> uncheckedRead(
                    reference, bed, leadVariant.bedIndex()));
            for (int otherIndex = index + 1;
                    otherIndex < ordered.size(); otherIndex++) {
                IndexedCandidate other = ordered.get(otherIndex);
                if (assigned.contains(other.candidate().variantId())) continue;
                Variant otherVariant = reference.variants().get(
                    other.candidate().variantId());
                if (!leadVariant.chromosome().equals(otherVariant.chromosome())
                        || distance(leadVariant.position(), otherVariant.position())
                            > windowBases) continue;
                byte[] otherGenotypes = genotypes.computeIfAbsent(
                    otherVariant.bedIndex(), ignored -> uncheckedRead(
                        reference, bed, otherVariant.bedIndex()));
                double rSquared = haplotypeRSquared(leadGenotypes,
                    otherGenotypes, reference.founders());
                if (Double.isFinite(rSquared)
                        && rSquared > options.rSquaredThreshold()) {
                    assigned.add(other.candidate().variantId());
                    exclusions.add(new IndexedExclusion(other.index(),
                        new LdClumpExclusion(other.candidate(),
                            LdClumpExclusionReason.IN_LINKAGE_DISEQUILIBRIUM,
                            lead.candidate().variantId(), rSquared)));
                }
            }
        }
    }

    private static byte[] uncheckedRead(Reference reference, FileChannel bed,
            long bedIndex) {
        try {
            return reference.readGenotypes(bed, bedIndex);
        } catch (IOException exception) {
            throw new GenotypeReadFailure(exception);
        }
    }

    private static long distance(long first, long second) {
        long difference = first - second;
        if (difference == Long.MIN_VALUE) return Long.MAX_VALUE;
        return Math.abs(difference);
    }

    private static double haplotypeRSquared(byte[] first, byte[] second,
            boolean[] founders) {
        long[][] counts = new long[3][3];
        long observations = 0;
        for (int sample = 0; sample < founders.length; sample++) {
            if (!founders[sample]) continue;
            int a = first[sample];
            int b = second[sample];
            if (a < 0 || b < 0) continue;
            counts[a][b]++;
            observations++;
        }
        if (observations == 0) return Double.NaN;

        double alleleA = 0.0;
        double alleleB = 0.0;
        for (int a = 0; a <= 2; a++) for (int b = 0; b <= 2; b++) {
            alleleA += a * counts[a][b];
            alleleB += b * counts[a][b];
        }
        alleleA /= 2.0 * observations;
        alleleB /= 2.0 * observations;
        if (!(alleleA > 0.0 && alleleA < 1.0
                && alleleB > 0.0 && alleleB < 1.0)) return Double.NaN;

        double[] frequency = {
            (1.0 - alleleA) * (1.0 - alleleB),
            (1.0 - alleleA) * alleleB,
            alleleA * (1.0 - alleleB),
            alleleA * alleleB};
        double[] fixed = fixedHaplotypeCounts(counts);
        double doubleHeterozygotes = counts[1][1];
        for (int iteration = 0; iteration < EM_ITERATIONS; iteration++) {
            double coupling = frequency[0] * frequency[3];
            double repulsion = frequency[1] * frequency[2];
            double probability = coupling + repulsion == 0.0
                ? 0.5 : coupling / (coupling + repulsion);
            double[] next = fixed.clone();
            next[0] += doubleHeterozygotes * probability;
            next[3] += doubleHeterozygotes * probability;
            next[1] += doubleHeterozygotes * (1.0 - probability);
            next[2] += doubleHeterozygotes * (1.0 - probability);
            for (int haplotype = 0; haplotype < next.length; haplotype++)
                next[haplotype] /= 2.0 * observations;
            double difference = 0.0;
            for (int haplotype = 0; haplotype < next.length; haplotype++)
                difference = Math.max(difference,
                    Math.abs(next[haplotype] - frequency[haplotype]));
            frequency = next;
            if (difference < EM_TOLERANCE) break;
        }
        double pA = frequency[2] + frequency[3];
        double pB = frequency[1] + frequency[3];
        double denominator = pA * (1.0 - pA) * pB * (1.0 - pB);
        if (!(denominator > 0.0)) return Double.NaN;
        double disequilibrium = frequency[3] - pA * pB;
        return Math.min(1.0, disequilibrium * disequilibrium / denominator);
    }

    private static double[] fixedHaplotypeCounts(long[][] counts) {
        double[] result = new double[4];
        result[0] += 2.0 * counts[0][0];
        result[0] += counts[0][1] + counts[1][0];
        result[1] += counts[0][1] + 2.0 * counts[0][2] + counts[1][2];
        result[2] += counts[1][0] + 2.0 * counts[2][0] + counts[2][1];
        result[3] += counts[1][2] + counts[2][1] + 2.0 * counts[2][2];
        return result;
    }

    private static void validateUnique(List<LdClumpCandidate> candidates) {
        Set<String> seen = new HashSet<>();
        for (LdClumpCandidate candidate : candidates) {
            String key = candidate.group() + '\u0000' + candidate.variantId();
            if (!seen.add(key)) throw new IllegalArgumentException(
                "candidate variants must be unique within each group: "
                    + candidate.variantId());
        }
    }

    private record IndexedCandidate(LdClumpCandidate candidate, int index) { }
    private record IndexedExclusion(int index, LdClumpExclusion exclusion) { }
    private record Variant(String chromosome, String id, long position,
            long bedIndex) { }

    private static final class GenotypeReadFailure extends RuntimeException {
        private static final long serialVersionUID = 1L;
        GenotypeReadFailure(IOException cause) { super(cause); }
    }

    private record Reference(Path bed, boolean[] founders,
            Map<String, Variant> variants, int sampleCount,
            int bytesPerVariant) {
        static Reference open(Path prefix, Set<String> requested)
                throws IOException {
            Path bed = Path.of(prefix + ".bed");
            Path bim = Path.of(prefix + ".bim");
            Path fam = Path.of(prefix + ".fam");
            boolean[] founders = readFounders(fam);
            Map<String, Variant> variants = readVariants(bim, requested);
            int bytesPerVariant = (founders.length + 3) / 4;
            long variantCount;
            try (java.util.stream.Stream<String> lines = Files.lines(bim,
                    StandardCharsets.UTF_8)) {
                variantCount = lines.filter(line -> !line.isBlank()).count();
            }
            long expected = Math.addExact(3L,
                Math.multiplyExact(variantCount, bytesPerVariant));
            if (Files.size(bed) != expected) throw new IOException(
                "PLINK BED size does not match BIM/FAM dimensions: " + bed);
            try (FileChannel channel = FileChannel.open(bed,
                    StandardOpenOption.READ)) {
                ByteBuffer header = ByteBuffer.allocate(3);
                readFully(channel, header, 0L);
                if (!Arrays.equals(header.array(), BED_HEADER))
                    throw new IOException(
                        "reference is not a variant-major PLINK 1 BED: " + bed);
            }
            return new Reference(bed, founders, variants, founders.length,
                bytesPerVariant);
        }

        FileChannel openBed() throws IOException {
            return FileChannel.open(bed, StandardOpenOption.READ);
        }

        byte[] readGenotypes(FileChannel channel, long variantIndex)
                throws IOException {
            ByteBuffer packed = ByteBuffer.allocate(bytesPerVariant);
            long offset = Math.addExact(3L,
                Math.multiplyExact(variantIndex, bytesPerVariant));
            readFully(channel, packed, offset);
            byte[] result = new byte[sampleCount];
            for (int sample = 0; sample < sampleCount; sample++) {
                int code = (packed.array()[sample / 4] >>>
                    (2 * (sample % 4))) & 3;
                result[sample] = switch (code) {
                    case 0 -> 2;
                    case 2 -> 1;
                    case 3 -> 0;
                    default -> -1;
                };
            }
            return result;
        }

        private static void readFully(FileChannel channel, ByteBuffer target,
                long position) throws IOException {
            while (target.hasRemaining()) {
                int count = channel.read(target, position + target.position());
                if (count < 0) throw new EOFException("truncated PLINK BED");
                if (count == 0) continue;
            }
        }

        private static boolean[] readFounders(Path fam) throws IOException {
            List<Boolean> values = new ArrayList<>();
            try (BufferedReader input = Files.newBufferedReader(fam,
                    StandardCharsets.UTF_8)) {
                for (String line; (line = input.readLine()) != null;) {
                    if (line.isBlank()) continue;
                    String[] fields = line.trim().split("\\s+");
                    if (fields.length < 6) throw new IOException(
                        "PLINK FAM row has fewer than six fields: " + fam);
                    values.add(fields[2].equals("0") && fields[3].equals("0"));
                }
            }
            if (values.isEmpty()) throw new IOException("PLINK FAM is empty: " + fam);
            boolean[] result = new boolean[values.size()];
            for (int index = 0; index < result.length; index++)
                result[index] = values.get(index);
            if (values.stream().noneMatch(Boolean::booleanValue))
                throw new IOException("PLINK reference contains no founders: " + fam);
            return result;
        }

        private static Map<String, Variant> readVariants(Path bim,
                Set<String> requested) throws IOException {
            Map<String, Variant> result = new HashMap<>();
            Set<String> allIds = new HashSet<>();
            long index = 0;
            try (BufferedReader input = Files.newBufferedReader(bim,
                    StandardCharsets.UTF_8)) {
                for (String line; (line = input.readLine()) != null;) {
                    if (line.isBlank()) continue;
                    String[] fields = line.trim().split("\\s+");
                    if (fields.length < 6) throw new IOException(
                        "PLINK BIM row has fewer than six fields: " + bim);
                    if (!allIds.add(fields[1])) throw new IOException(
                        "duplicate variant ID in PLINK BIM: " + fields[1]);
                    if (requested.contains(fields[1])) {
                        long position;
                        try {
                            position = Long.parseLong(fields[3]);
                        } catch (NumberFormatException exception) {
                            throw new IOException("invalid PLINK BIM position for "
                                + fields[1], exception);
                        }
                        result.put(fields[1], new Variant(fields[0], fields[1],
                            position, index));
                    }
                    index++;
                }
            }
            if (index == 0) throw new IOException("PLINK BIM is empty: " + bim);
            return Map.copyOf(result);
        }
    }
}
