/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/** One-row-per-event CSV/TSV writer for results, QC exclusions, and failures. */
public final class DelimitedAssociationWriter
        implements AssociationPipelineSink, Closeable {
    private static final String[] HEADER = {
        "status", "id", "chromosome", "position", "ref", "alt",
        "called", "missing", "missing_rate", "mac", "maf", "info",
        "beta", "se", "statistic", "df", "p", "log10_p",
        "negative_log10_p", "filter_reason", "error_type", "message"
    };
    private final BufferedWriter writer;
    private final char delimiter;

    public DelimitedAssociationWriter(Path path, char delimiter)
            throws IOException {
        if (path == null || (delimiter != ',' && delimiter != '\t'))
            throw new IllegalArgumentException(
                "output path and comma/tab delimiter are required");
        writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8);
        this.delimiter = delimiter;
        write(HEADER);
    }

    @Override
    public void acceptEstimate(AssociationPipelineEstimate estimate)
            throws IOException {
        VariantRecord variant = estimate.variant();
        VariantStatistics qc = estimate.variantStatistics();
        write(new String[] {"tested", variant.id(), variant.chromosome(),
            Long.toString(variant.position()), variant.referenceAllele(),
            variant.alternateAllele(), Integer.toString(qc.calledSamples()),
            Integer.toString(qc.missingSamples()), value(qc.missingRate()),
            value(qc.minorAlleleCount()), value(qc.minorAlleleFrequency()),
            value(variant.imputationQuality()), value(estimate.beta()),
            value(estimate.standardError()), value(estimate.statistic()),
            value(estimate.degreesOfFreedom()), value(estimate.pValue()),
            value(estimate.log10PValue()),
            value(estimate.negativeLog10PValue()), "", "", ""});
    }

    @Override
    public void acceptExcluded(VariantFilterResult excluded)
            throws IOException {
        VariantRecord variant = excluded.variant();
        VariantStatistics qc = excluded.statistics();
        write(new String[] {"excluded", variant.id(), variant.chromosome(),
            Long.toString(variant.position()), variant.referenceAllele(),
            variant.alternateAllele(), Integer.toString(qc.calledSamples()),
            Integer.toString(qc.missingSamples()), value(qc.missingRate()),
            value(qc.minorAlleleCount()), value(qc.minorAlleleFrequency()),
            value(variant.imputationQuality()), "", "", "", "", "", "",
            "", excluded.reasons().stream().map(Enum::name)
                .collect(java.util.stream.Collectors.joining(";")), "", ""});
    }

    @Override
    public void acceptFailure(AssociationPipelineFailure failure)
            throws IOException {
        write(new String[] {"failure", failure.variantId(), "", "", "", "",
            "", "", "", "", "", "", "", "", "", "", "", "", "",
            "", failure.exceptionType(), failure.message()});
    }

    @Override public void close() throws IOException { writer.close(); }

    private void write(String[] fields) throws IOException {
        for (int index = 0; index < fields.length; index++) {
            if (index > 0) writer.write(delimiter);
            writer.write(escape(fields[index]));
        }
        writer.newLine();
    }

    private String escape(String value) {
        String text = value == null ? "" : value;
        if (text.indexOf(delimiter) < 0 && text.indexOf('"') < 0
                && text.indexOf('\n') < 0 && text.indexOf('\r') < 0)
            return text;
        return '"' + text.replace("\"", "\"\"") + '"';
    }

    private static String value(double number) {
        return Double.toString(number);
    }
}
