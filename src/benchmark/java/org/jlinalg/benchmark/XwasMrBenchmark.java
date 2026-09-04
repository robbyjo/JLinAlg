/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.benchmark;

import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.jlinalg.mr.AlleleHarmonizer;
import org.jlinalg.mr.HarmonizedInstrument;
import org.jlinalg.mr.MendelianRandomization;
import org.jlinalg.mr.MrEstimate;
import org.jlinalg.mr.MrOptions;
import org.jlinalg.mr.SummaryAssociation;
import org.jlinalg.mr.XwasMrBatchResult;
import org.jlinalg.mr.XwasMrExposure;
import org.jlinalg.mr.XwasMrOptions;
import org.jlinalg.mr.XwasMrOutcome;
import org.jlinalg.mr.XwasMrPipeline;
import org.jlinalg.mr.XwasMrScreeningMethod;
import org.jlinalg.mr.XwasMrSignificanceFilter;

/** Reproducible synthetic exposure-by-phenotype IVW benchmark shared with R. */
public final class XwasMrBenchmark {
    private static volatile long checksum;

    private XwasMrBenchmark() { }

    public static void main(String[] arguments) throws Exception {
        int exposureCount = integer("jlinalg.benchmark.exposures", 300);
        int outcomeCount = integer("jlinalg.benchmark.outcomes", 150);
        int instrumentCount = integer("jlinalg.benchmark.instruments", 10);
        int threads = integer("jlinalg.benchmark.parallelism",
            Math.min(8, Runtime.getRuntime().availableProcessors()));
        int warmups = integer("jlinalg.benchmark.warmups", 2);
        int measurements = integer("jlinalg.benchmark.measurements", 5);
        Path output = Path.of(System.getProperty("jlinalg.benchmark.output",
            "build/benchmarks/xwas-mr")).toAbsolutePath().normalize();
        if (exposureCount < 1 || outcomeCount < 1 || instrumentCount < 3
                || threads < 1 || warmups < 0 || measurements < 1)
            throw new IllegalArgumentException("invalid benchmark dimensions");
        Files.createDirectories(output);

        Data data = data(exposureCount, outcomeCount, instrumentCount);
        XwasMrPipeline pipeline = XwasMrPipeline.prepare(
            data.exposures(), data.outcomes());
        XwasMrOptions options = new XwasMrOptions(threads,
            Math.max(32, threads * 16), XwasMrScreeningMethod.IVW_MULTIPLICATIVE_RANDOM,
            XwasMrSignificanceFilter.pValueAtMost(0.0),
            new MrOptions(0.95, 2, 20260904L));

        for (int warmup = 0; warmup < warmups; warmup++) consume(
            pipeline.scan(options));
        List<Double> seconds = new ArrayList<>();
        for (int measurement = 1; measurement <= measurements; measurement++) {
            System.gc();
            long started = System.nanoTime();
            XwasMrBatchResult result = pipeline.scan(options);
            double elapsed = (System.nanoTime() - started) / 1e9;
            consume(result);
            seconds.add(elapsed);
            System.out.printf(Locale.ROOT,
                "JLinAlg measurement=%d exposures=%d outcomes=%d pairs=%d "
                    + "instruments=%d threads=%d seconds=%.9f pairs_per_second=%.3f%n",
                measurement, exposureCount, outcomeCount, result.totalPairs(),
                instrumentCount, threads, elapsed, result.totalPairs() / elapsed);
        }
        writeTimings(output.resolve("jlinalg_timings.csv"), seconds,
            exposureCount, outcomeCount, instrumentCount, threads);
        writeValidation(output.resolve("jlinalg_validation.csv"),
            Math.min(5, exposureCount), Math.min(5, outcomeCount),
            instrumentCount);
        if (checksum == 0L)
            throw new IllegalStateException("benchmark result was not consumed");
    }

    private static void consume(XwasMrBatchResult result) {
        checksum += result.totalPairs() + result.analyzablePairs()
            + result.belowThresholdPairs() + result.hits().size();
    }

    private static Data data(int exposures, int outcomes, int instruments) {
        List<XwasMrExposure> exposureValues = new ArrayList<>(exposures);
        for (int exposure = 0; exposure < exposures; exposure++)
            exposureValues.add(new XwasMrExposure("E" + exposure,
                "Exposure " + exposure, exposureAssociations(exposure, instruments)));
        List<XwasMrOutcome> outcomeValues = new ArrayList<>(outcomes);
        for (int outcome = 0; outcome < outcomes; outcome++) {
            List<SummaryAssociation> associations = new ArrayList<>(
                exposures * instruments);
            for (int exposure = 0; exposure < exposures; exposure++)
                for (int instrument = 0; instrument < instruments; instrument++)
                    associations.add(outcomeAssociation(exposure, outcome,
                        instrument));
            outcomeValues.add(new XwasMrOutcome("O" + outcome,
                "Outcome " + outcome, category(outcome), associations));
        }
        return new Data(exposureValues, outcomeValues);
    }

    private static List<SummaryAssociation> exposureAssociations(
            int exposure, int instruments) {
        List<SummaryAssociation> result = new ArrayList<>(instruments);
        for (int instrument = 0; instrument < instruments; instrument++)
            result.add(new SummaryAssociation(variant(exposure, instrument),
                "A", "C", exposureEffect(exposure, instrument),
                exposureStandardError(instrument), 0.2));
        return result;
    }

    private static SummaryAssociation outcomeAssociation(int exposure,
            int outcome, int instrument) {
        return new SummaryAssociation(variant(exposure, instrument), "A", "C",
            outcomeEffect(exposure, outcome, instrument),
            outcomeStandardError(instrument), 0.2);
    }

    private static double exposureEffect(int exposure, int instrument) {
        return 0.05 + 0.004 * instrument + 0.0001 * (exposure % 17);
    }

    private static double exposureStandardError(int instrument) {
        return 0.01 + 0.0002 * (instrument % 3);
    }

    private static double outcomeEffect(int exposure, int outcome,
            int instrument) {
        double causal = ((exposure * 13 + outcome * 7) % 11 - 5) * 0.02;
        double noise = 0.01 * Math.sin((exposure + 1) * 0.17
            + (outcome + 1) * 0.11 + (instrument + 1) * 0.37);
        return causal * exposureEffect(exposure, instrument) + noise;
    }

    private static double outcomeStandardError(int instrument) {
        return 0.04 + 0.001 * (instrument % 4);
    }

    private static String variant(int exposure, int instrument) {
        return "v" + exposure + '_' + instrument;
    }

    private static String category(int outcome) {
        return switch (outcome % 3) {
            case 0 -> "cardiovascular";
            case 1 -> "kidney";
            default -> "lung";
        };
    }

    private static void writeValidation(Path path, int exposures, int outcomes,
            int instruments) throws Exception {
        try (BufferedWriter output = Files.newBufferedWriter(path,
                StandardCharsets.UTF_8)) {
            output.write("exposure,outcome,beta,se,statistic,p_value,cochran_q,");
            output.write("heterogeneity_p_value,i_squared\n");
            for (int exposure = 0; exposure < exposures; exposure++)
                for (int outcome = 0; outcome < outcomes; outcome++) {
                    List<SummaryAssociation> x = exposureAssociations(
                        exposure, instruments);
                    List<SummaryAssociation> y = new ArrayList<>(instruments);
                    for (int instrument = 0; instrument < instruments; instrument++)
                        y.add(outcomeAssociation(exposure, outcome, instrument));
                    List<HarmonizedInstrument> values = AlleleHarmonizer.harmonize(
                        x, y).instruments();
                    MrEstimate estimate = MendelianRandomization.ivw(
                        values, true, 0.95);
                    output.write(String.format(Locale.ROOT,
                        "%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g%n",
                        exposure, outcome, estimate.estimate(),
                        estimate.standardError(), estimate.statistic(),
                        estimate.pValue(), estimate.cochranQ(),
                        estimate.heterogeneityPValue(), estimate.iSquared()));
                }
        }
    }

    private static void writeTimings(Path path, List<Double> seconds,
            int exposures, int outcomes, int instruments, int threads)
            throws Exception {
        try (BufferedWriter output = Files.newBufferedWriter(path,
                StandardCharsets.UTF_8)) {
            output.write("runtime,measurement,exposures,outcomes,pairs,instruments,");
            output.write("threads,seconds,pairs_per_second\n");
            long pairs = (long) exposures * outcomes;
            for (int index = 0; index < seconds.size(); index++) {
                double elapsed = seconds.get(index);
                output.write(String.format(Locale.ROOT,
                    "JLinAlg,%d,%d,%d,%d,%d,%d,%.9f,%.6f%n",
                    index + 1, exposures, outcomes, pairs, instruments, threads,
                    elapsed, pairs / elapsed));
            }
        }
    }

    private static int integer(String name, int fallback) {
        int value = Integer.getInteger(name, fallback);
        if (value < 0) throw new IllegalArgumentException(name + " is invalid");
        return value;
    }

    private record Data(List<XwasMrExposure> exposures,
            List<XwasMrOutcome> outcomes) { }
}
