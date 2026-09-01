/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.zip.DeflaterOutputStream;
import java.util.zip.GZIPOutputStream;
import com.github.luben.zstd.Zstd;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gwas.AssociationScanOptions;
import org.jlinalg.gwas.GenotypeMissingPolicy;
import org.jlinalg.gwas.RemlAssociationScanner;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PipelineInputTest {
    @TempDir Path temporaryDirectory;

    @Test
    void delimitedSourceReadsGenomicRowsAndReordersSamples() throws Exception {
        Path input = temporaryDirectory.resolve("variants.csv");
        Files.writeString(input,
            "id,chromosome,position,ref,alt,S1,S2,S3\n"
            + "\"v,1\",1,100,A,G,0,1,2\n"
            + "v2,1,200,C,T,.,0.5,1.5\n");
        VariantSource source = VariantSources.open(input);

        assertEquals(DataFormat.CSV, source.metadata().format());
        assertEquals(List.of("S1", "S2", "S3"),
            source.metadata().sampleIds());
        assertEquals(2, source.metadata().variantCount());
        try (VariantBlockReader reader = source.open(new int[] {2, 0})) {
            VariantBlock block = reader.read(1);
            assertEquals("v,1", block.variants().get(0).id());
            assertArrayEquals(new double[] {2, 0},
                block.variants().get(0).dosages(), 0);
            VariantRecord second = reader.read(2).variants().get(0);
            assertArrayEquals(new double[] {1.5, Double.NaN},
                second.dosages());
            assertNull(reader.read(1));
        }
    }

    @Test
    void vcfAndVcfGzPreferDosagesAndExpandAlternateAlleles() throws Exception {
        String text = "##fileformat=VCFv4.2\n"
            + "##contig=<ID=1,length=1000>\n"
            + "##INFO=<ID=R2,Number=1,Type=Float,Description=\"quality\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n"
            + "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"dosage\">\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            + "1\t10\trs1\tA\tG\t.\tPASS\tR2=0.9\tGT:DS\t0/0:0\t0/1:1\t1/1:2\n"
            + "1\t20\trs2\tC\tT,G\t.\tPASS\t.\tGT:DS\t0/1:1,0\t1/2:1,1\t2/2:0,2\n";
        Path vcf = temporaryDirectory.resolve("input.vcf");
        Files.writeString(vcf, text);
        Path gzip = temporaryDirectory.resolve("input.vcf.gz");
        try (OutputStream output = new GZIPOutputStream(
                Files.newOutputStream(gzip))) {
            output.write(text.getBytes(StandardCharsets.UTF_8));
        }
        for (Path path : List.of(vcf, gzip)) {
            VariantSource source = VariantSources.open(path);
            try (VariantBlockReader reader = source.open(new int[] {2, 0, 1})) {
                VariantBlock block = reader.read(10);
                assertEquals(3, block.variants().size());
                assertArrayEquals(new double[] {2, 0, 1},
                    block.variants().get(0).dosages(), 0);
                assertEquals(0.9,
                    block.variants().get(0).imputationQuality(), 1e-6);
                assertArrayEquals(new double[] {0, 1, 1},
                    block.variants().get(1).dosages(), 0);
                assertArrayEquals(new double[] {2, 0, 1},
                    block.variants().get(2).dosages(), 0);
            }
        }
        Path bcf = temporaryDirectory.resolve("input.bcf");
        try (VCFIterator input = new VCFIteratorBuilder().open(vcf);
             VariantContextWriter writer = new VariantContextWriterBuilder()
                 .setOutputPath(bcf)
                 .unsetOption(Options.INDEX_ON_THE_FLY)
                 .setOption(Options.FORCE_BCF)
                 .build()) {
            writer.writeHeader(input.getHeader());
            while (input.hasNext()) writer.add(input.next());
        }
        VariantSource bcfSource = VariantSources.open(bcf);
        try (VariantBlockReader reader = bcfSource.open()) {
            VariantBlock block = reader.read(10);
            assertEquals(3, block.variants().size());
            assertArrayEquals(new double[] {0, 1, 2},
                block.variants().get(0).dosages(), 0);
        }
    }

    @Test
    void bgenLayoutTwoDecodesUncompressedZlibAndZstdProbabilities()
            throws Exception {
        for (int compression = 0; compression <= 2; compression++) {
            Path bgen = temporaryDirectory.resolve(
                "input-" + compression + ".bgen");
            writeBgen(bgen, compression);
            VariantSource source = VariantSources.open(bgen);

            assertEquals(List.of("S1", "S2", "S3"),
                source.metadata().sampleIds());
            try (VariantBlockReader reader = source.open(
                    new int[] {2, 1, 0})) {
                VariantRecord variant = reader.read(10).variants().get(0);
                assertEquals("rs10", variant.id());
                assertEquals("1", variant.chromosome());
                assertEquals(10, variant.position());
                assertArrayEquals(new double[] {2, 1, 0},
                    variant.dosages(), 1e-12);
                assertNull(reader.read(1));
            }
        }
    }

    @Test
    void alignedFrequencyFiltersRetainStructuredReasons() {
        VariantRecord variant = new VariantRecord("rare", "1", 1, "A", "G",
            new double[] {0, 0, 1, Double.NaN}, 0.8);
        VariantFilterOptions options = VariantFilterOptions.builder()
            .minimumMac(2).maximumMissingRate(0.2)
            .minimumImputationQuality(0.9).build();
        VariantFilterResult result = VariantFilters.evaluate(variant, options);

        assertFalse(result.included());
        assertEquals(1, result.statistics().minorAlleleCount(), 0);
        assertEquals(1.0 / 6, result.statistics().minorAlleleFrequency(), 1e-12);
        assertTrue(result.reasons().contains(
            VariantFilterReason.BELOW_MINIMUM_MAC));
        assertTrue(result.reasons().contains(
            VariantFilterReason.TOO_MANY_MISSING));
        assertTrue(result.reasons().contains(
            VariantFilterReason.BELOW_IMPUTATION_QUALITY));
    }

    @Test
    void omicsTransformsAreComposableTieAwareAndNonMutating() {
        double[] source = {0, 1, 1, 100, Double.NaN};
        double[] transformed = OmicsTransforms.compose(
            OmicsTransforms.winsorize(0, 0.75),
            OmicsTransforms.rankInverseNormal()).apply(source);

        assertArrayEquals(new double[] {0, 1, 1, 100, Double.NaN}, source);
        assertEquals(transformed[1], transformed[2], 0);
        assertTrue(transformed[0] < transformed[1]);
        assertTrue(transformed[3] > transformed[1]);
        assertTrue(Double.isNaN(transformed[4]));
        double[] z = OmicsTransforms.zScore().apply(new double[] {1, 2, 3});
        assertEquals(0, z[1], 1e-12);
        assertEquals(0, (z[0] + z[1] + z[2]) / 3, 1e-12);
    }

    @Test
    void streamedFastOlsFiltersAndReturnsOrderedEffects() throws Exception {
        Path input = temporaryDirectory.resolve("scan.tsv");
        Files.writeString(input,
            "id\tS1\tS2\tS3\tS4\tS5\tS6\n"
            + "signal\t0\t0\t1\t1\t2\t2\n"
            + "constant\t1\t1\t1\t1\t1\t1\n"
            + "noise\t0\t1\t0\t1\t0\t1\n");
        VariantSource source = VariantSources.open(input);
        double[] response = {1, 1.2, 2, 2.1, 3, 3.2};
        double[][] covariates = {{1}, {1}, {1}, {1}, {1}, {1}};
        AssociationPipelineResult result = StreamingAssociationPipeline.fastOls(
            source, source.metadata().sampleIds(), response, covariates,
            null, null, OlsOptions.defaults(),
            AssociationEngineOptions.acceleratedSerial()
                .withBackendPolicy(BackendPolicy.CPU),
            new AssociationPipelineOptions(1,
                VariantFilterOptions.defaults()));

        assertEquals(3, result.sourceVariants());
        assertEquals(2, result.testedVariants());
        assertEquals(List.of("signal", "noise"), result.estimates().stream()
            .map(value -> value.variant().id()).toList());
        assertEquals(1, result.excludedVariantCount());
        assertTrue(result.estimates().get(0).beta() > 0.8);

        Path output = temporaryDirectory.resolve("scan-results.tsv");
        AssociationPipelineSummary summary;
        try (DelimitedAssociationWriter writer =
                new DelimitedAssociationWriter(output, '\t')) {
            summary = StreamingAssociationPipeline.fastOlsTo(
                source, source.metadata().sampleIds(), response, covariates,
                null, null, OlsOptions.defaults(),
                AssociationEngineOptions.acceleratedSerial()
                    .withBackendPolicy(BackendPolicy.CPU),
                new AssociationPipelineOptions(2,
                    VariantFilterOptions.defaults()), writer);
        }
        List<String> lines = Files.readAllLines(output);
        assertEquals(3, summary.sourceVariants());
        assertEquals(2, summary.testedVariants());
        assertEquals(1, summary.excludedVariants());
        assertEquals(4, lines.size());
        assertTrue(lines.get(0).contains("negative_log10_p"));
        assertTrue(lines.stream().anyMatch(line ->
            line.startsWith("tested\tsignal\t")));
        assertTrue(lines.stream().anyMatch(line ->
            line.startsWith("excluded\tconstant\t")));
    }

    @Test
    void omicsMatrixStreamsTransformedChangingResponses() throws Exception {
        Path input = temporaryDirectory.resolve("expression.csv");
        Files.writeString(input,
            "feature,S1,S2,S3,S4,S5,S6\n"
            + "geneA,1,1.2,2,2.1,3,3.2\n"
            + "geneB,4,3.9,4.1,4,4.2,4.1\n"
            + "geneMissing,1,.,2,2,3,3\n");
        NumericMatrixSource source = DelimitedMatrixSource.open(input);
        double[][] design = {
            {1, 0}, {1, 0}, {1, 1},
            {1, 1}, {1, 2}, {1, 2}
        };

        OmicsAssociationResult result =
            StreamingOmicsAssociationPipeline.scanResponses(
                source, source.metadata().sampleIds(), design, 1,
                OmicsTransforms.winsorize(0, 1),
                OmicsMissingPolicy.MEAN_IMPUTE, 2,
                null, null, OlsOptions.defaults(),
                AssociationEngineOptions.acceleratedSerial()
                    .withBackendPolicy(BackendPolicy.CPU));

        assertEquals(3, result.sourceFeatures());
        assertEquals(List.of("geneA", "geneB", "geneMissing"),
            result.estimates().stream()
                .map(OmicsAssociationEstimate::featureId).toList());
        assertTrue(result.estimates().get(0).beta() > 0.8);
        assertTrue(Math.abs(result.estimates().get(1).beta()) < 0.2);
        assertTrue(result.failures().isEmpty());
    }

    @Test
    void streamedP3dFitsOneMixedNullAndScansBlocks() throws Exception {
        Path input = temporaryDirectory.resolve("mixed.csv");
        Files.writeString(input,
            "id,S1,S2,S3,S4,S5,S6,S7,S8\n"
            + "g1,0,0,1,1,2,2,1,2\n"
            + "g2,0,1,0,1,0,1,0,1\n");
        VariantSource source = VariantSources.open(input);
        double[] response = {1, 2, 2, 4, 5, 7, 7, 9};
        double[][] covariates = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        RemlAssociationScanner scanner = RemlAssociationScanner.prepare(
            response, covariates,
            List.of(VarianceComponent.identity(
                "residual", response.length)),
            RemlOptions.builder().initialVariances(2).build(),
            BackendPolicy.CPU);

        AssociationPipelineResult result =
            StreamingAssociationPipeline.remlP3d(
                source, source.metadata().sampleIds(), scanner,
                new AssociationScanOptions(
                    1, GenotypeMissingPolicy.MEAN_IMPUTE, 1),
                new AssociationPipelineOptions(
                    1, VariantFilterOptions.defaults()));

        assertEquals(2, result.estimates().size());
        assertEquals(List.of("g1", "g2"), result.estimates().stream()
            .map(value -> value.variant().id()).toList());
        assertTrue(result.estimates().stream()
            .allMatch(value -> Double.isFinite(value.pValue())));
    }

    private static void writeBgen(Path path, int compression)
            throws IOException {
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        int headerLength = 20;
        int sampleBlockLength = 8 + 3 * 4;
        int variantStart = 4 + headerLength + sampleBlockLength;
        littleInt(output, variantStart - 4);
        littleInt(output, headerLength);
        littleInt(output, 1);
        littleInt(output, 3);
        output.writeBytes("bgen".getBytes(StandardCharsets.US_ASCII));
        littleInt(output, (2 << 2) | compression | 0x80000000);
        littleInt(output, sampleBlockLength);
        littleInt(output, 3);
        for (String sample : List.of("S1", "S2", "S3"))
            string16(output, sample);

        string16(output, "rs10");
        string16(output, "rs10");
        string16(output, "1");
        littleInt(output, 10);
        littleShort(output, 2);
        string32(output, "A");
        string32(output, "G");
        ByteArrayOutputStream probabilities = new ByteArrayOutputStream();
        littleInt(probabilities, 3);
        littleShort(probabilities, 2);
        probabilities.write(2);
        probabilities.write(2);
        probabilities.write(new byte[] {2, 2, 2});
        probabilities.write(0);
        probabilities.write(8);
        probabilities.write(new byte[] {
            (byte) 255, 0,
            0, (byte) 255,
            0, 0
        });
        byte[] raw = probabilities.toByteArray();
        if (compression == 0) {
            littleInt(output, raw.length);
            output.writeBytes(raw);
        } else {
            byte[] compressed;
            if (compression == 1) {
                ByteArrayOutputStream buffer = new ByteArrayOutputStream();
                try (DeflaterOutputStream deflater =
                        new DeflaterOutputStream(buffer)) {
                    deflater.write(raw);
                }
                compressed = buffer.toByteArray();
            } else {
                compressed = Zstd.compress(raw);
            }
            littleInt(output, compressed.length + 4);
            littleInt(output, raw.length);
            output.writeBytes(compressed);
        }
        Files.write(path, output.toByteArray());
    }

    private static void string16(ByteArrayOutputStream output, String value) {
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        littleShort(output, bytes.length);
        output.writeBytes(bytes);
    }

    private static void string32(ByteArrayOutputStream output, String value) {
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        littleInt(output, bytes.length);
        output.writeBytes(bytes);
    }

    private static void littleShort(ByteArrayOutputStream output, int value) {
        output.write(value & 255);
        output.write((value >>> 8) & 255);
    }

    private static void littleInt(ByteArrayOutputStream output, int value) {
        output.write(value & 255);
        output.write((value >>> 8) & 255);
        output.write((value >>> 16) & 255);
        output.write((value >>> 24) & 255);
    }
}
