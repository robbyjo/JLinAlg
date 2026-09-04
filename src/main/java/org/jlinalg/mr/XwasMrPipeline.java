/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

/**
 * Bounded, parallel exposure-by-outcome MR screening for xWAS applications.
 *
 * <p>Exposure instruments are assumed to have already been selected and
 * clumped. Outcome associations are indexed by variant once. Pair evaluations
 * are processed in bounded blocks, and full MR diagnostics are calculated only
 * for pairs passing the requested primary IVW threshold. Results are emitted in
 * exposure-major then outcome-major order independently of thread scheduling.</p>
 */
public final class XwasMrPipeline {
    private static final int MINIMUM_DIAGNOSTIC_INSTRUMENTS = 3;

    private final List<XwasMrExposure> exposures;
    private final List<XwasMrOutcome> outcomes;
    private final List<Map<String, List<SummaryAssociation>>> outcomeIndexes;

    private XwasMrPipeline(List<XwasMrExposure> exposures,
            List<XwasMrOutcome> outcomes) {
        this.exposures = List.copyOf(exposures);
        this.outcomes = List.copyOf(outcomes);
        List<Map<String, List<SummaryAssociation>>> indexes =
            new ArrayList<>(outcomes.size());
        for (XwasMrOutcome outcome : outcomes)
            indexes.add(index(outcome.associations()));
        outcomeIndexes = List.copyOf(indexes);
    }

    /** Validates and prepares reusable exposure and phenotype inputs. */
    public static XwasMrPipeline prepare(
            List<XwasMrExposure> exposures,
            List<XwasMrOutcome> outcomes) {
        if (exposures == null || exposures.isEmpty()
                || exposures.stream().anyMatch(value -> value == null))
            throw new IllegalArgumentException(
                "at least one nonnull exposure is required");
        if (outcomes == null || outcomes.isEmpty()
                || outcomes.stream().anyMatch(value -> value == null))
            throw new IllegalArgumentException(
                "at least one nonnull outcome is required");
        uniqueIds(exposures.stream().map(XwasMrExposure::id).toList(),
            "exposure");
        uniqueIds(outcomes.stream().map(XwasMrOutcome::id).toList(),
            "outcome");
        return new XwasMrPipeline(exposures, outcomes);
    }

    /** Returns the number of prepared molecular or conventional exposures. */
    public int exposureCount() { return exposures.size(); }

    /** Returns the number of prepared outcome phenotypes. */
    public int outcomeCount() { return outcomes.size(); }

    /**
     * Screens the complete exposure-by-outcome grid with bounded memory.
     *
     * <p>Pairs with fewer than three harmonized instruments are counted and
     * skipped. Fixed or multiplicative-random IVW is used for the inexpensive
     * first stage. Only significant pairs receive the full
     * {@link MendelianRandomization#analyze(List, MrOptions)} analysis.</p>
     */
    public XwasMrBatchResult scan(XwasMrOptions options) {
        return scan(options, ignored -> { });
    }

    /**
     * Screens the grid and streams every successful primary IVW result.
     *
     * <p>The sink is called serially in exposure-major, outcome-major order
     * after each bounded parallel block. It therefore does not need to be
     * thread safe. Pairs that cannot be screened are represented by the batch
     * accounting and failure collections rather than sink calls.</p>
     */
    public XwasMrBatchResult scan(XwasMrOptions options,
            XwasMrScreeningSink sink) {
        if (options == null)
            throw new IllegalArgumentException("xWAS MR options are required");
        Objects.requireNonNull(sink, "sink");
        long totalPairs = Math.multiplyExact((long) exposures.size(),
            outcomes.size());
        int workers = (int) Math.min((long) options.parallelism(), totalPairs);
        long started = System.nanoTime();
        List<XwasMrHit> hits = new ArrayList<>();
        List<XwasMrFailure> failures = new ArrayList<>();
        long analyzable = 0L;
        long below = 0L;
        long insufficient = 0L;

        ForkJoinPool pool = new ForkJoinPool(workers);
        try {
            for (long from = 0L; from < totalPairs;
                    from += options.pairBlockSize()) {
                int size = (int) Math.min(options.pairBlockSize(),
                    totalPairs - from);
                PairEvaluation[] block = new PairEvaluation[size];
                evaluateBlock(pool, from, block, options);
                for (PairEvaluation evaluation : block) {
                    if (evaluation.screening() != null)
                        sink.accept(evaluation.screening());
                    switch (evaluation.status()) {
                        case HIT -> {
                            analyzable++;
                            hits.add(evaluation.hit());
                        }
                        case BELOW_THRESHOLD -> {
                            analyzable++;
                            below++;
                        }
                        case INSUFFICIENT_INSTRUMENTS -> insufficient++;
                        case FAILED -> {
                            analyzable++;
                            failures.add(evaluation.failure());
                        }
                    }
                }
            }
        } finally {
            pool.shutdown();
        }
        return new XwasMrBatchResult(hits, failures, totalPairs, analyzable,
            below, insufficient, workers, System.nanoTime() - started);
    }

    private void evaluateBlock(ForkJoinPool pool, long from,
            PairEvaluation[] block, XwasMrOptions options) {
        try {
            pool.submit(() -> IntStream.range(0, block.length).parallel()
                .forEach(offset -> block[offset] = evaluate(from + offset,
                    options))).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("xWAS MR scan was interrupted",
                exception);
        } catch (ExecutionException exception) {
            Throwable cause = exception.getCause();
            if (cause instanceof RuntimeException runtime) throw runtime;
            throw new IllegalStateException("xWAS MR scan failed", cause);
        }
    }

    private PairEvaluation evaluate(long pairIndex, XwasMrOptions options) {
        int exposureIndex = (int) (pairIndex / outcomes.size());
        int outcomeIndex = (int) (pairIndex % outcomes.size());
        XwasMrExposure exposure = exposures.get(exposureIndex);
        XwasMrOutcome outcome = outcomes.get(outcomeIndex);
        List<SummaryAssociation> relevantOutcome = relevant(
            exposure.clumpedInstruments(), outcomeIndexes.get(outcomeIndex));
        HarmonizationResult harmonized = AlleleHarmonizer.harmonize(
            exposure.clumpedInstruments(), relevantOutcome);
        if (harmonized.instruments().size()
                < MINIMUM_DIAGNOSTIC_INSTRUMENTS)
            return PairEvaluation.insufficient();
        MrEstimate screening;
        try {
            screening = screen(harmonized.instruments(), options);
        } catch (RuntimeException exception) {
            return PairEvaluation.failed(failure(exposureIndex, outcomeIndex,
                exposure, outcome, exception), null);
        }
        boolean thresholdPassed = options.significanceFilter().includes(
            screening.pValue());
        double negativeLog10 = screening.pValue() == 0.0
            ? Double.POSITIVE_INFINITY : -Math.log10(screening.pValue());
        XwasMrScreeningResult screened = new XwasMrScreeningResult(
            exposureIndex, outcomeIndex, exposure.id(), exposure.label(),
            outcome.id(), outcome.label(), outcome.category(), screening,
            negativeLog10, thresholdPassed);
        if (!thresholdPassed) return PairEvaluation.below(screened);
        try {
            MrAnalysisResult analysis = MendelianRandomization.analyze(
                harmonized.instruments(), options.diagnosticOptions());
            return PairEvaluation.hit(new XwasMrHit(exposureIndex,
                outcomeIndex, exposure.id(), exposure.label(), outcome.id(),
                outcome.label(), outcome.category(), screening,
                negativeLog10, analysis, harmonized.exclusions()), screened);
        } catch (RuntimeException exception) {
            return PairEvaluation.failed(failure(exposureIndex, outcomeIndex,
                exposure, outcome, exception), screened);
        }
    }

    private static XwasMrFailure failure(int exposureIndex, int outcomeIndex,
            XwasMrExposure exposure, XwasMrOutcome outcome,
            RuntimeException exception) {
        String message = exception.getMessage() == null
            ? "" : exception.getMessage();
        return new XwasMrFailure(exposureIndex, outcomeIndex, exposure.id(),
            outcome.id(), exception.getClass().getSimpleName(), message);
    }

    private static MrEstimate screen(List<HarmonizedInstrument> instruments,
            XwasMrOptions options) {
        return switch (options.screeningMethod()) {
            case IVW_FIXED -> MendelianRandomization.ivw(
                instruments, false,
                options.diagnosticOptions().confidenceLevel());
            case IVW_MULTIPLICATIVE_RANDOM -> MendelianRandomization.ivw(
                instruments, true,
                options.diagnosticOptions().confidenceLevel());
        };
    }

    private static List<SummaryAssociation> relevant(
            List<SummaryAssociation> instruments,
            Map<String, List<SummaryAssociation>> outcomeIndex) {
        List<SummaryAssociation> result = new ArrayList<>();
        Set<String> requested = new HashSet<>();
        for (SummaryAssociation instrument : instruments)
            if (requested.add(instrument.variantId())) {
                List<SummaryAssociation> found = outcomeIndex.get(
                    instrument.variantId());
                if (found != null) result.addAll(found);
            }
        return result;
    }

    private static Map<String, List<SummaryAssociation>> index(
            List<SummaryAssociation> associations) {
        Map<String, List<SummaryAssociation>> mutable = new LinkedHashMap<>();
        for (SummaryAssociation association : associations)
            mutable.computeIfAbsent(association.variantId(),
                ignored -> new ArrayList<>()).add(association);
        Map<String, List<SummaryAssociation>> result = new HashMap<>();
        for (Map.Entry<String, List<SummaryAssociation>> entry
                : mutable.entrySet())
            result.put(entry.getKey(), List.copyOf(entry.getValue()));
        return Map.copyOf(result);
    }

    private static void uniqueIds(List<String> ids, String kind) {
        Set<String> unique = new HashSet<>();
        for (String id : ids) if (!unique.add(id))
            throw new IllegalArgumentException(
                "duplicate " + kind + " id: " + id);
    }

    private enum PairStatus {
        HIT,
        BELOW_THRESHOLD,
        INSUFFICIENT_INSTRUMENTS,
        FAILED
    }

    private record PairEvaluation(PairStatus status, XwasMrHit hit,
            XwasMrFailure failure, XwasMrScreeningResult screening) {
        static PairEvaluation hit(XwasMrHit hit,
                XwasMrScreeningResult screening) {
            return new PairEvaluation(PairStatus.HIT, hit, null, screening);
        }
        static PairEvaluation below(XwasMrScreeningResult screening) {
            return new PairEvaluation(PairStatus.BELOW_THRESHOLD, null, null,
                screening);
        }
        static PairEvaluation insufficient() {
            return new PairEvaluation(
                PairStatus.INSUFFICIENT_INSTRUMENTS, null, null, null);
        }
        static PairEvaluation failed(XwasMrFailure failure,
                XwasMrScreeningResult screening) {
            return new PairEvaluation(PairStatus.FAILED, null, failure,
                screening);
        }
    }
}
