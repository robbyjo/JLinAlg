/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import com.github.luben.zstd.ZstdInputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.InflaterInputStream;

/**
 * BGEN layout-2 additive dosage reader.
 *
 * <p>BGEN 1.2 zlib, BGEN 1.3 zstd, and uncompressed layout-2 blocks are
 * supported at 1--32 probability bits. Biallelic phased/unphased and variable
 * ploidy records are decoded. A file without embedded sample IDs must supply
 * them explicitly.</p>
 */
public final class BgenVariantSource implements VariantSource {
    private final Path path;
    private final Header header;
    private final VariantSourceMetadata metadata;

    public BgenVariantSource(Path path) throws IOException {
        this(path, null);
    }

    public BgenVariantSource(Path path, List<String> externalSampleIds)
            throws IOException {
        this.path = path.toAbsolutePath().normalize();
        try (RandomAccessFile input = new RandomAccessFile(this.path.toFile(), "r")) {
            header = readHeader(input, externalSampleIds);
        }
        metadata = new VariantSourceMetadata(this.path, header.variantCount(),
            header.sampleIds(), DataFormat.BGEN);
    }

    @Override public VariantSourceMetadata metadata() { return metadata; }

    @Override
    public VariantBlockReader open(int[] requestedOrder) throws IOException {
        int[] order = DelimitedVariantSource.normalizeOrder(requestedOrder,
            metadata.sampleIds().size());
        RandomAccessFile input = new RandomAccessFile(path.toFile(), "r");
        try {
            input.seek(header.variantStart());
            return new Reader(input, order);
        } catch (IOException | RuntimeException exception) {
            input.close();
            throw exception;
        }
    }

    private final class Reader implements VariantBlockReader {
        private final RandomAccessFile input;
        private final int[] order;
        private long variantIndex;
        private boolean closed;

        private Reader(RandomAccessFile input, int[] order) {
            this.input = input;
            this.order = order;
        }

        @Override
        public VariantBlock read(int maximumVariants) throws IOException {
            if (maximumVariants < 1)
                throw new IllegalArgumentException("maximum variants must be positive");
            if (closed) throw new IOException("BGEN reader is closed");
            if (variantIndex >= header.variantCount()) return null;
            long first = variantIndex;
            List<VariantRecord> variants = new ArrayList<>(maximumVariants);
            while (variants.size() < maximumVariants
                    && variantIndex < header.variantCount()) {
                variants.add(readVariant());
                variantIndex++;
            }
            return new VariantBlock(first, variants);
        }

        private VariantRecord readVariant() throws IOException {
            String variantId = readString16(input);
            String rsid = readString16(input);
            String chromosome = readString16(input);
            long position = unsignedInt(input);
            int alleleCount = unsignedShort(input);
            if (alleleCount < 2)
                throw new IOException("BGEN variant has fewer than two alleles");
            List<String> alleles = new ArrayList<>(alleleCount);
            for (int index = 0; index < alleleCount; index++)
                alleles.add(readString32(input));
            int storedLength = arrayLength(unsignedInt(input), "BGEN block length");
            byte[] probabilityData;
            if (header.compression() == 0) {
                probabilityData = bytes(input, storedLength);
            } else {
                int uncompressedLength = arrayLength(
                    unsignedInt(input), "BGEN uncompressed block length");
                if (storedLength < 4)
                    throw new IOException("compressed BGEN block is shorter than 4 bytes");
                byte[] compressed = bytes(input, storedLength - 4);
                probabilityData = decompress(compressed, uncompressedLength,
                    header.compression());
            }
            if (alleleCount != 2)
                throw new IOException(
                    "BGEN multiallelic records are not yet supported: "
                    + displayId(variantId, rsid, chromosome, position, alleles));
            double[] sourceDosages = decode(probabilityData,
                header.sampleCount(), alleleCount);
            double[] selected = new double[order.length];
            for (int index = 0; index < order.length; index++)
                selected[index] = sourceDosages[order[index]];
            String id = displayId(variantId, rsid, chromosome, position, alleles);
            return new VariantRecord(id, chromosome, position,
                alleles.get(0), alleles.get(1), selected, Double.NaN);
        }

        @Override public void close() throws IOException {
            if (!closed) {
                closed = true;
                input.close();
            }
        }
    }

    private static Header readHeader(
            RandomAccessFile input, List<String> externalSampleIds)
            throws IOException {
        long offset = unsignedInt(input);
        long headerLength = unsignedInt(input);
        long variantCount = unsignedInt(input);
        int sampleCount = arrayLength(unsignedInt(input), "BGEN sample count");
        byte[] magic = bytes(input, 4);
        boolean validMagic = magic[0] == 'b' && magic[1] == 'g'
            && magic[2] == 'e' && magic[3] == 'n';
        boolean zeroMagic = magic[0] == 0 && magic[1] == 0
            && magic[2] == 0 && magic[3] == 0;
        if (!validMagic && !zeroMagic)
            throw new IOException("file does not contain a BGEN header");
        if (headerLength < 20 || headerLength > offset)
            throw new IOException("invalid BGEN header length");
        input.seek(4 + headerLength - 4);
        long flags = unsignedInt(input);
        int compression = (int) (flags & 3);
        int layout = (int) ((flags >>> 2) & 15);
        boolean embeddedSamples = (flags & 0x80000000L) != 0;
        if (compression > 2)
            throw new IOException("unsupported BGEN compression code: " + compression);
        if (layout != 2)
            throw new IOException(
                "only BGEN layout 2 (v1.2/v1.3) is supported; found layout "
                + layout);

        List<String> sampleIds;
        if (embeddedSamples) {
            input.seek(4 + headerLength);
            long sampleBlockLength = unsignedInt(input);
            int blockSamples = arrayLength(
                unsignedInt(input), "BGEN sample block count");
            if (blockSamples != sampleCount)
                throw new IOException(
                    "BGEN sample block count does not match header");
            List<String> parsed = new ArrayList<>(sampleCount);
            for (int index = 0; index < sampleCount; index++)
                parsed.add(readString16(input));
            long expectedEnd = 4 + headerLength + sampleBlockLength;
            if (input.getFilePointer() > expectedEnd)
                throw new IOException("BGEN sample block length is invalid");
            sampleIds = List.copyOf(parsed);
        } else {
            if (externalSampleIds == null)
                throw new IOException(
                    "BGEN file has no embedded sample IDs; supply them explicitly");
            if (externalSampleIds.size() != sampleCount)
                throw new IOException(
                    "external BGEN sample count does not match header");
            sampleIds = List.copyOf(externalSampleIds);
        }
        return new Header(variantCount, sampleCount, sampleIds,
            offset + 4, compression);
    }

    private static double[] decode(
            byte[] block, int expectedSamples, int expectedAlleles)
            throws IOException {
        Bytes input = new Bytes(block);
        int samples = arrayLength(input.unsignedInt(), "BGEN probability sample count");
        int alleles = input.unsignedShort();
        if (samples != expectedSamples || alleles != expectedAlleles)
            throw new IOException(
                "BGEN probability dimensions do not match header");
        input.unsignedByte(); // minimum ploidy, retained for validation by rows
        input.unsignedByte(); // maximum ploidy
        int[] ploidy = new int[samples];
        boolean[] missing = new boolean[samples];
        for (int sample = 0; sample < samples; sample++) {
            int value = input.unsignedByte();
            ploidy[sample] = value & 63;
            missing[sample] = (value & 128) != 0;
        }
        int phased = input.unsignedByte();
        int bits = input.unsignedByte();
        if ((phased != 0 && phased != 1) || bits < 1 || bits > 32)
            throw new IOException("invalid BGEN phased flag or probability precision");
        BitReader probabilities = new BitReader(block, input.position(), bits);
        double denominator = bits == 32
            ? 4294967295.0 : (1L << bits) - 1.0;
        double[] dosages = new double[samples];
        for (int sample = 0; sample < samples; sample++) {
            int copies = ploidy[sample];
            if (copies < 1)
                throw new IOException("BGEN sample ploidy must be positive");
            double dosage;
            if (phased == 1) {
                dosage = 0;
                for (int copy = 0; copy < copies; copy++) {
                    double referenceProbability =
                        probabilities.next() / denominator;
                    dosage += 1.0 - referenceProbability;
                }
            } else {
                double storedSum = 0;
                dosage = 0;
                for (int alternateCopies = 0;
                        alternateCopies < copies; alternateCopies++) {
                    double probability = probabilities.next() / denominator;
                    storedSum += probability;
                    dosage += alternateCopies * probability;
                }
                double last = 1.0 - storedSum;
                if (last < -1e-6)
                    throw new IOException(
                        "BGEN genotype probabilities sum above one");
                dosage += copies * Math.max(0, last);
            }
            dosages[sample] = missing[sample] ? Double.NaN : dosage;
        }
        return dosages;
    }

    private static byte[] decompress(
            byte[] compressed, int expectedLength, int compression)
            throws IOException {
        InputStream base = new ByteArrayInputStream(compressed);
        try (InputStream decoded = compression == 1
                ? new InflaterInputStream(base)
                : new ZstdInputStream(base);
             ByteArrayOutputStream output =
                 new ByteArrayOutputStream(expectedLength)) {
            decoded.transferTo(output);
            byte[] result = output.toByteArray();
            if (result.length != expectedLength)
                throw new IOException("BGEN decompressed block length mismatch");
            return result;
        }
    }

    private static String displayId(
            String variantId, String rsid, String chromosome, long position,
            List<String> alleles) {
        if (!variantId.isBlank()) return variantId;
        if (!rsid.isBlank()) return rsid;
        return chromosome + ":" + position + ":" + alleles.get(0)
            + ":" + alleles.get(1);
    }

    private static String readString16(RandomAccessFile input) throws IOException {
        return new String(bytes(input, unsignedShort(input)),
            StandardCharsets.UTF_8);
    }

    private static String readString32(RandomAccessFile input) throws IOException {
        return new String(bytes(input,
            arrayLength(unsignedInt(input), "BGEN string length")),
            StandardCharsets.UTF_8);
    }

    private static int unsignedShort(RandomAccessFile input) throws IOException {
        int low = input.read();
        int high = input.read();
        if ((low | high) < 0) throw new EOFException("unexpected end of BGEN file");
        return low | (high << 8);
    }

    private static long unsignedInt(RandomAccessFile input) throws IOException {
        long first = input.read();
        long second = input.read();
        long third = input.read();
        long fourth = input.read();
        if ((first | second | third | fourth) < 0)
            throw new EOFException("unexpected end of BGEN file");
        return first | (second << 8) | (third << 16) | (fourth << 24);
    }

    private static byte[] bytes(RandomAccessFile input, int length)
            throws IOException {
        byte[] result = new byte[length];
        input.readFully(result);
        return result;
    }

    private static int arrayLength(long value, String name) throws IOException {
        if (value < 0 || value > Integer.MAX_VALUE)
            throw new IOException(name + " exceeds the Java array limit");
        return (int) value;
    }

    private record Header(long variantCount, int sampleCount,
        List<String> sampleIds, long variantStart, int compression) { }

    private static final class Bytes {
        private final byte[] values;
        private int position;
        private Bytes(byte[] values) { this.values = values; }
        private int position() { return position; }
        private int unsignedByte() throws IOException {
            if (position >= values.length)
                throw new EOFException("truncated BGEN probability block");
            return values[position++] & 255;
        }
        private int unsignedShort() throws IOException {
            return unsignedByte() | (unsignedByte() << 8);
        }
        private long unsignedInt() throws IOException {
            return unsignedByte() | ((long) unsignedByte() << 8)
                | ((long) unsignedByte() << 16)
                | ((long) unsignedByte() << 24);
        }
    }

    /** BGEN probabilities are packed least-significant bit first. */
    private static final class BitReader {
        private final byte[] values;
        private final int bits;
        private long bitOffset;
        private BitReader(byte[] values, int byteOffset, int bits) {
            this.values = values;
            this.bits = bits;
            bitOffset = byteOffset * 8L;
        }
        private long next() throws IOException {
            long result = 0;
            for (int bit = 0; bit < bits; bit++, bitOffset++) {
                int byteIndex = (int) (bitOffset >>> 3);
                if (byteIndex >= values.length)
                    throw new EOFException("truncated BGEN packed probabilities");
                int shift = (int) (bitOffset & 7);
                result |= (long) ((values[byteIndex] >>> shift) & 1) << bit;
            }
            return result;
        }
    }
}
