/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.jlinalg.pipeline.AssociationPipelineEstimate;
import org.jlinalg.pipeline.AssociationPipelineFailure;
import org.jlinalg.pipeline.AssociationPipelineSink;
import org.jlinalg.pipeline.OmicsAssociationEstimate;
import org.jlinalg.pipeline.OmicsAssociationSink;
import org.jlinalg.pipeline.VariantFilterResult;
import org.jlinalg.pipeline.VariantRecord;
import org.jlinalg.pipeline.VariantStatistics;

/** Common streaming result table for genotype and generic omics scans. */
final class CliResultSink
        implements AssociationPipelineSink, OmicsAssociationSink, AutoCloseable {
    private static final List<String> BASE_HEADER = List.of(
        "status", "omics_type", "id", "chromosome", "position",
        "reference_allele", "alternate_allele", "effect_allele",
        "other_allele", "effect_allele_frequency",
        "minor_allele_frequency", "minor_allele_count",
        "called", "missing", "missing_rate", "imputation_info",
        "hwe_p_all", "hwe_n_all", "hwe_p_cases", "hwe_n_cases",
        "hwe_p_controls", "hwe_n_controls", "hwe_method",
        "beta", "standard_error", "statistic", "statistic_type",
        "df_numerator", "df_denominator", "df_method",
        "partial_r2", "partial_r2_method", "p_value",
        "filter_reason", "error_type", "message");
    private final String omicsType;
    private final String statisticType;
    private final String dfMethod;
    private final AnnotationLookup annotation;
    private final int[] caseControlGroups;
    private final ExternalBh output;

    CliResultSink(
            Path path, boolean overwrite, String omicsType,
            String statisticType, String dfMethod,
            AnnotationLookup annotation, int[] caseControlGroups)
            throws IOException {
        this.omicsType = omicsType;
        this.statisticType = statisticType;
        this.dfMethod = dfMethod;
        this.annotation = annotation;
        this.caseControlGroups = caseControlGroups == null ? null
            : caseControlGroups.clone();
        output = new ExternalBh(path, overwrite);
        List<String> header = new ArrayList<>(BASE_HEADER);
        for (String column : annotation.columns())
            header.add("annot_" + column);
        output.writeHeader(header);
    }

    @Override
    public void acceptEstimate(AssociationPipelineEstimate estimate)
            throws IOException {
        VariantRecord variant = estimate.variant();
        VariantStatistics qc = estimate.variantStatistics();
        Hwe hwe = hwe(variant);
        write("ok", variant.id(), variant.chromosome(),
            variant.position() == 0 ? "" : Long.toString(variant.position()),
            variant.referenceAllele(), variant.alternateAllele(),
            variant.alternateAllele(), variant.referenceAllele(),
            number(qc.alternateAlleleFrequency()),
            number(qc.minorAlleleFrequency()),
            number(qc.minorAlleleCount()),
            Integer.toString(qc.calledSamples()),
            Integer.toString(qc.missingSamples()),
            number(qc.missingRate()), number(variant.imputationQuality()), hwe,
            estimate.beta(), estimate.standardError(), estimate.statistic(),
            estimate.degreesOfFreedom(), estimate.pValue(), "", "", "");
    }

    @Override
    public void acceptExcluded(VariantFilterResult excluded)
            throws IOException {
        VariantRecord variant = excluded.variant();
        VariantStatistics qc = excluded.statistics();
        Hwe hwe = hwe(variant);
        write("filtered", variant.id(), variant.chromosome(),
            variant.position() == 0 ? "" : Long.toString(variant.position()),
            variant.referenceAllele(), variant.alternateAllele(),
            variant.alternateAllele(), variant.referenceAllele(),
            number(qc.alternateAlleleFrequency()),
            number(qc.minorAlleleFrequency()),
            number(qc.minorAlleleCount()),
            Integer.toString(qc.calledSamples()),
            Integer.toString(qc.missingSamples()),
            number(qc.missingRate()), number(variant.imputationQuality()), hwe,
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            excluded.reasons().stream().map(Enum::name)
                .collect(java.util.stream.Collectors.joining(";")), "", "");
    }

    @Override
    public void acceptEstimate(OmicsAssociationEstimate estimate)
            throws IOException {
        write("ok", estimate.featureId(), "", "", "", "", "", "",
            "", "", "", "", "", "", "", Hwe.empty(),
            estimate.beta(), estimate.standardError(), estimate.statistic(),
            estimate.degreesOfFreedom(), estimate.pValue(), "", "", "");
    }

    @Override
    public void acceptFailure(AssociationPipelineFailure failure)
            throws IOException {
        write("failed", failure.variantId(), "", "", "", "", "", "",
            "", "", "", "", "", "", "", Hwe.empty(),
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            "", failure.exceptionType(), failure.message());
    }

    void finish() throws IOException { output.finish(); }
    long adjustedTests() { return output.tests(); }
    @Override public void close() throws IOException { output.close(); }

    private void write(
            String status, String id, String chromosome, String position,
            String reference, String alternate, String effect, String other,
            String eaf, String maf, String mac, String called, String missing,
            String missingRate, String info, Hwe hwe,
            double beta, double standardError, double statistic,
            double degreesOfFreedom, double pValue,
            String filter, String error, String message) throws IOException {
        double partial = statisticType.equals("t")
            ? partialR2(statistic, degreesOfFreedom) : Double.NaN;
        List<String> fields = new ArrayList<>(BASE_HEADER.size()
            + annotation.columns().size());
        fields.addAll(List.of(status, omicsType, id, chromosome, position,
            reference, alternate, effect, other, eaf, maf, mac, called,
            missing, missingRate, info, number(hwe.all().pValue()),
            integer(hwe.all().samples()), number(hwe.cases().pValue()),
            integer(hwe.cases().samples()), number(hwe.controls().pValue()),
            integer(hwe.controls().samples()), hwe.method(),
            number(beta), number(standardError), number(statistic),
            statisticType, "1", number(degreesOfFreedom), dfMethod,
            number(partial), Double.isFinite(partial)
                ? "test-statistic" : "", number(pValue),
            filter, error, message));
        fields.addAll(List.of(annotation.values(id)));
        output.write(fields, pValue);
    }

    private Hwe hwe(VariantRecord variant) {
        double[] dosage = variant.dosages();
        HardyWeinberg.Result all = HardyWeinberg.calculate(
            dosage, caseControlGroups, -1);
        HardyWeinberg.Result cases = caseControlGroups == null
            ? new HardyWeinberg.Result(Double.NaN, 0, "")
            : HardyWeinberg.calculate(dosage, caseControlGroups, 1);
        HardyWeinberg.Result controls = caseControlGroups == null
            ? new HardyWeinberg.Result(Double.NaN, 0, "")
            : HardyWeinberg.calculate(dosage, caseControlGroups, 0);
        return new Hwe(all, cases, controls, all.method());
    }

    private static double partialR2(double statistic, double df) {
        if (!Double.isFinite(statistic) || !Double.isFinite(df) || df <= 0)
            return Double.NaN;
        double squared = statistic * statistic;
        return squared / (squared + df);
    }
    private static String integer(int value) {
        return value == 0 ? "" : Integer.toString(value);
    }
    private static String number(double value) {
        return Double.isFinite(value) ? Double.toString(value) : "";
    }

    private record Hwe(
            HardyWeinberg.Result all,
            HardyWeinberg.Result cases,
            HardyWeinberg.Result controls,
            String method) {
        private static Hwe empty() {
            HardyWeinberg.Result empty =
                new HardyWeinberg.Result(Double.NaN, 0, "");
            return new Hwe(empty, empty, empty, "");
        }
    }
}
