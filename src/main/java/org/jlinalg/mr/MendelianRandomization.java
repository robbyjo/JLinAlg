/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import jdistlib.ChiSquare;
import jdistlib.Normal;

/** Core two-sample summary-data MR estimators for independent instruments. */
public final class MendelianRandomization {
    private MendelianRandomization() {
    }

    /** Runs the core MR analysis with reproducible defaults. */
    public static MrAnalysisResult analyze(List<HarmonizedInstrument> instruments) {
        return analyze(instruments, MrOptions.defaults());
    }

    /**
     * Runs Wald ratios, fixed and multiplicative-random IVW, MR-Egger,
     * bootstrapped weighted median, strength, heterogeneity, and leave-one-out.
     */
    public static MrAnalysisResult analyze(
            List<HarmonizedInstrument> instruments, MrOptions options) {
        List<HarmonizedInstrument> values = validated(instruments, 3);
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        List<WaldRatio> ratios = waldRatios(values);
        MrEstimate fixed = ivw(values, false, options.confidenceLevel());
        MrEstimate random = ivw(values, true, options.confidenceLevel());
        MrEggerResult egger = egger(values, options.confidenceLevel());
        MrEstimate median = weightedMedian(values, options);
        List<LeaveOneOutEstimate> leaveOneOut = leaveOneOut(
            values, options.confidenceLevel());
        double meanF = ratios.stream().mapToDouble(WaldRatio::fStatistic)
            .average().orElse(Double.NaN);
        return new MrAnalysisResult(ratios, fixed, random, egger,
            median, leaveOneOut, meanF,
            warnings(ratios, fixed, egger));
    }

    /** Returns first-order single-variant Wald ratio estimates. */
    public static List<WaldRatio> waldRatios(
            List<HarmonizedInstrument> instruments) {
        List<HarmonizedInstrument> values = validated(instruments, 1);
        List<WaldRatio> result = new ArrayList<>(values.size());
        for (HarmonizedInstrument instrument : values) {
            double estimate = instrument.outcomeEffect()
                / instrument.exposureEffect();
            double standardError = instrument.outcomeStandardError()
                / Math.abs(instrument.exposureEffect());
            double statistic = estimate / standardError;
            result.add(new WaldRatio(instrument.variantId(), estimate,
                standardError, statistic, normalPValue(statistic),
                square(instrument.exposureEffect()
                    / instrument.exposureStandardError())));
        }
        return List.copyOf(result);
    }

    /** Returns an IVW estimate for independent instruments. */
    public static MrEstimate ivw(
            List<HarmonizedInstrument> instruments,
            boolean multiplicativeRandomEffects,
            double confidenceLevel) {
        List<HarmonizedInstrument> values = validated(instruments, 2);
        validateConfidence(confidenceLevel);
        double numerator = 0.0;
        double denominator = 0.0;
        for (HarmonizedInstrument instrument : values) {
            double weight = 1.0 / square(instrument.outcomeStandardError());
            numerator += weight * instrument.exposureEffect()
                * instrument.outcomeEffect();
            denominator += weight * square(instrument.exposureEffect());
        }
        if (!(denominator > 0.0) || !Double.isFinite(denominator)) {
            throw new IllegalArgumentException("IVW information is not finite and positive");
        }
        double estimate = numerator / denominator;
        double q = 0.0;
        for (HarmonizedInstrument instrument : values) {
            double residual = instrument.outcomeEffect()
                - estimate * instrument.exposureEffect();
            q += square(residual / instrument.outcomeStandardError());
        }
        int degreesOfFreedom = values.size() - 1;
        double dispersion = multiplicativeRandomEffects
            ? Math.max(1.0, q / degreesOfFreedom) : 1.0;
        double standardError = Math.sqrt(dispersion / denominator);
        MrMethod method = multiplicativeRandomEffects
            ? MrMethod.IVW_MULTIPLICATIVE_RANDOM : MrMethod.IVW_FIXED;
        return estimate(method, estimate, standardError, confidenceLevel,
            q, degreesOfFreedom, dispersion, values.size());
    }

    /** Returns MR-Egger after orienting every exposure association positively. */
    public static MrEggerResult egger(
            List<HarmonizedInstrument> instruments, double confidenceLevel) {
        List<HarmonizedInstrument> values = validated(instruments, 3);
        validateConfidence(confidenceLevel);
        double sumWeight = 0.0;
        double sumX = 0.0;
        double sumY = 0.0;
        double sumXX = 0.0;
        double sumXY = 0.0;
        for (HarmonizedInstrument instrument : values) {
            double sign = Math.copySign(1.0, instrument.exposureEffect());
            double x = Math.abs(instrument.exposureEffect());
            double y = sign * instrument.outcomeEffect();
            double weight = 1.0 / square(instrument.outcomeStandardError());
            sumWeight += weight;
            sumX += weight * x;
            sumY += weight * y;
            sumXX += weight * x * x;
            sumXY += weight * x * y;
        }
        double determinant = sumWeight * sumXX - sumX * sumX;
        if (!(determinant > 0.0) || !Double.isFinite(determinant)) {
            throw new IllegalArgumentException(
                "MR-Egger requires varying exposure associations");
        }
        double slope = (sumWeight * sumXY - sumX * sumY) / determinant;
        double intercept = (sumY - slope * sumX) / sumWeight;
        double q = 0.0;
        for (HarmonizedInstrument instrument : values) {
            double sign = Math.copySign(1.0, instrument.exposureEffect());
            double residual = sign * instrument.outcomeEffect()
                - intercept - slope * Math.abs(instrument.exposureEffect());
            q += square(residual / instrument.outcomeStandardError());
        }
        int degreesOfFreedom = values.size() - 2;
        double dispersion = Math.max(1.0, q / degreesOfFreedom);
        double slopeStandardError = Math.sqrt(dispersion * sumWeight / determinant);
        double interceptStandardError = Math.sqrt(dispersion * sumXX / determinant);
        MrEstimate slopeEstimate = estimate(MrMethod.MR_EGGER,
            slope, slopeStandardError, confidenceLevel,
            q, degreesOfFreedom, dispersion, values.size());
        double interceptStatistic = intercept / interceptStandardError;
        double critical = normalCritical(confidenceLevel);
        return new MrEggerResult(slopeEstimate, intercept,
            interceptStandardError, interceptStatistic,
            normalPValue(interceptStatistic),
            intercept - critical * interceptStandardError,
            intercept + critical * interceptStandardError,
            iSquaredGx(values));
    }

    /** Returns the weighted-median estimate with parametric-bootstrap inference. */
    public static MrEstimate weightedMedian(
            List<HarmonizedInstrument> instruments, MrOptions options) {
        List<HarmonizedInstrument> values = validated(instruments, 3);
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        double point = weightedMedianPoint(values, null, null);
        Random random = new Random(options.randomSeed());
        double[] sampledX = new double[values.size()];
        double[] sampledY = new double[values.size()];
        double[] replicates = new double[options.weightedMedianBootstrapReplicates()];
        for (int replicate = 0; replicate < replicates.length; replicate++) {
            for (int index = 0; index < values.size(); index++) {
                HarmonizedInstrument instrument = values.get(index);
                sampledX[index] = instrument.exposureEffect()
                    + instrument.exposureStandardError() * random.nextGaussian();
                sampledY[index] = instrument.outcomeEffect()
                    + instrument.outcomeStandardError() * random.nextGaussian();
                if (sampledX[index] == 0.0) {
                    sampledX[index] = Math.copySign(Double.MIN_NORMAL,
                        instrument.exposureEffect());
                }
            }
            replicates[replicate] = weightedMedianPoint(
                values, sampledX, sampledY);
        }
        double mean = Arrays.stream(replicates).average().orElseThrow();
        double sumSquares = 0.0;
        for (double value : replicates) {
            sumSquares += square(value - mean);
        }
        double standardError = Math.sqrt(sumSquares / (replicates.length - 1));
        return estimate(MrMethod.WEIGHTED_MEDIAN, point, standardError,
            options.confidenceLevel(), Double.NaN, 0, Double.NaN, values.size());
    }

    private static double weightedMedianPoint(
            List<HarmonizedInstrument> instruments,
            double[] sampledX,
            double[] sampledY) {
        WeightedRatio[] ordered = new WeightedRatio[instruments.size()];
        double totalWeight = 0.0;
        for (int index = 0; index < instruments.size(); index++) {
            HarmonizedInstrument instrument = instruments.get(index);
            double x = sampledX == null ? instrument.exposureEffect() : sampledX[index];
            double y = sampledY == null ? instrument.outcomeEffect() : sampledY[index];
            double weight = square(x / instrument.outcomeStandardError());
            ordered[index] = new WeightedRatio(y / x, weight);
            totalWeight += weight;
        }
        if (!(totalWeight > 0.0) || !Double.isFinite(totalWeight)) {
            throw new IllegalArgumentException(
                "weighted-median weights are not finite and positive");
        }
        Arrays.sort(ordered, Comparator.comparingDouble(WeightedRatio::ratio));
        double[] probability = new double[ordered.length];
        double cumulative = 0.0;
        for (int index = 0; index < ordered.length; index++) {
            double normalized = ordered[index].weight() / totalWeight;
            probability[index] = cumulative + 0.5 * normalized;
            cumulative += normalized;
        }
        if (0.5 <= probability[0]) {
            return ordered[0].ratio();
        }
        for (int index = 1; index < ordered.length; index++) {
            if (0.5 <= probability[index]) {
                double fraction = (0.5 - probability[index - 1])
                    / (probability[index] - probability[index - 1]);
                return ordered[index - 1].ratio() + fraction
                    * (ordered[index].ratio() - ordered[index - 1].ratio());
            }
        }
        return ordered[ordered.length - 1].ratio();
    }

    private static List<LeaveOneOutEstimate> leaveOneOut(
            List<HarmonizedInstrument> instruments, double confidenceLevel) {
        List<LeaveOneOutEstimate> result = new ArrayList<>(instruments.size());
        for (int omitted = 0; omitted < instruments.size(); omitted++) {
            List<HarmonizedInstrument> retained = new ArrayList<>(instruments);
            HarmonizedInstrument removed = retained.remove(omitted);
            result.add(new LeaveOneOutEstimate(removed.variantId(),
                ivw(retained, false, confidenceLevel)));
        }
        return List.copyOf(result);
    }

    private static double iSquaredGx(List<HarmonizedInstrument> instruments) {
        double sumWeight = 0.0;
        double weightedSum = 0.0;
        for (HarmonizedInstrument instrument : instruments) {
            double weight = 1.0 / square(instrument.exposureStandardError());
            sumWeight += weight;
            weightedSum += weight * Math.abs(instrument.exposureEffect());
        }
        double mean = weightedSum / sumWeight;
        double q = 0.0;
        for (HarmonizedInstrument instrument : instruments) {
            double difference = Math.abs(instrument.exposureEffect()) - mean;
            q += square(difference / instrument.exposureStandardError());
        }
        return q == 0.0 ? 0.0
            : Math.max(0.0, (q - (instruments.size() - 1.0)) / q);
    }

    private static List<String> warnings(
            List<WaldRatio> ratios,
            MrEstimate fixed,
            MrEggerResult egger) {
        List<String> result = new ArrayList<>();
        long weak = ratios.stream().filter(value -> value.fStatistic() < 10.0).count();
        if (weak > 0) {
            result.add(weak + " instrument(s) have an approximate F statistic below 10");
        }
        if (egger.iSquaredGx() < 0.90) {
            result.add("MR-Egger I-squared GX is below 0.90; dilution from exposure measurement error may be material");
        }
        if (fixed.heterogeneityPValue() < 0.05) {
            result.add("IVW Cochran Q rejects homogeneity at the 0.05 level");
        }
        if (egger.interceptPValue() < 0.05) {
            result.add("MR-Egger intercept differs from zero at the 0.05 level");
        }
        return List.copyOf(result);
    }

    static MrEstimate estimate(
            MrMethod method,
            double estimate,
            double standardError,
            double confidenceLevel,
            double q,
            int degreesOfFreedom,
            double dispersion,
            int instrumentCount) {
        double statistic = standardError == 0.0
            ? Math.copySign(Double.POSITIVE_INFINITY, estimate)
            : estimate / standardError;
        double critical = normalCritical(confidenceLevel);
        double heterogeneityP = degreesOfFreedom > 0 && Double.isFinite(q)
            ? ChiSquare.cumulative(q, degreesOfFreedom, false, false)
            : Double.NaN;
        double iSquared = degreesOfFreedom > 0 && q > 0.0
            ? Math.max(0.0, (q - degreesOfFreedom) / q) : 0.0;
        return new MrEstimate(method, estimate, standardError,
            statistic, normalPValue(statistic),
            estimate - critical * standardError,
            estimate + critical * standardError,
            q, degreesOfFreedom, heterogeneityP, iSquared,
            dispersion, instrumentCount);
    }

    static List<HarmonizedInstrument> validated(
            List<HarmonizedInstrument> instruments, int minimumCount) {
        if (instruments == null || instruments.size() < minimumCount) {
            throw new IllegalArgumentException(
                "at least " + minimumCount + " instruments are required");
        }
        SetBuilder variants = new SetBuilder();
        for (HarmonizedInstrument instrument : instruments) {
            if (instrument == null) {
                throw new IllegalArgumentException("instruments must not contain null");
            }
            if (!variants.add(instrument.variantId())) {
                throw new IllegalArgumentException(
                    "instrument variant identifiers must be unique");
            }
            if (!Double.isFinite(instrument.exposureEffect())
                    || instrument.exposureEffect() == 0.0
                    || !(instrument.exposureStandardError() > 0.0)
                    || !Double.isFinite(instrument.exposureStandardError())
                    || !Double.isFinite(instrument.outcomeEffect())
                    || !(instrument.outcomeStandardError() > 0.0)
                    || !Double.isFinite(instrument.outcomeStandardError())) {
                throw new IllegalArgumentException(
                    "instrument effects and standard errors are invalid");
            }
        }
        return List.copyOf(instruments);
    }

    static void validateConfidence(double confidenceLevel) {
        if (!(confidenceLevel > 0.0) || !(confidenceLevel < 1.0)
                || !Double.isFinite(confidenceLevel)) {
            throw new IllegalArgumentException(
                "confidenceLevel must be finite and lie in (0, 1)");
        }
    }

    static double normalCritical(double confidenceLevel) {
        return Normal.quantile(0.5 + confidenceLevel / 2.0,
            0.0, 1.0, true, false);
    }

    static double normalPValue(double statistic) {
        return Math.min(1.0, 2.0 * Normal.cumulative(
            Math.abs(statistic), 0.0, 1.0, false, false));
    }

    private static double square(double value) {
        return value * value;
    }

    private record WeightedRatio(double ratio, double weight) {
    }

    private static final class SetBuilder {
        private final java.util.Set<String> values = new java.util.HashSet<>();

        private boolean add(String value) {
            return value != null && values.add(value);
        }
    }
}
