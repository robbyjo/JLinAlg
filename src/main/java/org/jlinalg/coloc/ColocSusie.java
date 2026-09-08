/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.coloc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.susie.SusieResult;

/**
 * Multi-signal colocalization compatible with {@code coloc::coloc.susie}.
 */
public final class ColocSusie {
    private ColocSusie() { }

    /** Colocalizes the credible signals in two JLinAlg SuSiE fits. */
    public static ColocSusieResult analyze(
            SusieResult trait1, SusieResult trait2) {
        return analyze(trait1, trait2, ColocOptions.defaults());
    }

    /** Colocalizes the credible signals in two JLinAlg SuSiE fits. */
    public static ColocSusieResult analyze(
            SusieResult trait1, SusieResult trait2, ColocOptions options) {
        return analyze(ColocSusieInput.from(trait1),
            ColocSusieInput.from(trait2), options);
    }

    /** Colocalizes arbitrary SuSiE log-Bayes-factor matrices. */
    public static ColocSusieResult analyze(
            ColocSusieInput trait1, ColocSusieInput trait2) {
        return analyze(trait1, trait2, ColocOptions.defaults());
    }

    /**
     * Colocalizes every retained pair of credible signals. Variant identifiers
     * are intersected in trait-1 order, matching R's {@code intersect} rule.
     */
    public static ColocSusieResult analyze(
            ColocSusieInput trait1, ColocSusieInput trait2,
            ColocOptions options) {
        if (trait1 == null || trait2 == null || options == null) {
            throw new IllegalArgumentException(
                "two SuSiE inputs and options are required");
        }
        validateWeights(options.rawTrait1PriorWeights(), trait1.variants(), 1);
        validateWeights(options.rawTrait2PriorWeights(), trait2.variants(), 2);

        Match match = match(trait1.variantNames(), trait2.variantNames());
        if (match.names().isEmpty() || trait1.signals() == 0
                || trait2.signals() == 0) {
            return new ColocSusieResult(match.names(), List.of(),
                new double[0], options, 0);
        }

        int common = match.names().size();
        double[] first = restrict(trait1, match.first(), common);
        double[] second = restrict(trait2, match.second(), common);
        double[] overlap1 = posteriorOverlap(trait1, match.first());
        double[] overlap2 = posteriorOverlap(trait2, match.second());
        int[] lead1 = rowMaxima(first, trait1.signals(), common);
        int[] lead2 = rowMaxima(second, trait2.signals(), common);

        Weights weights = weights(options, match);
        int totalPairs;
        int posteriorLength;
        try {
            totalPairs = Math.multiplyExact(trait1.signals(), trait2.signals());
            posteriorLength = Math.multiplyExact(totalPairs, common);
        } catch (ArithmeticException exception) {
            throw new IllegalArgumentException(
                "colocalization result is too large", exception);
        }
        List<ColocSignalPair> pairs = new ArrayList<>(totalPairs);
        double[] posterior = new double[posteriorLength];
        double[] hypotheses = new double[5];
        int retained = 0;

        // expand.grid(i, j) in coloc varies i fastest, so trait 2 is outermost.
        for (int signal2 = 0; signal2 < trait2.signals(); signal2++) {
            for (int signal1 = 0; signal1 < trait1.signals(); signal1++) {
                if (options.trimByPosterior()
                        && (overlap1[signal1] < options.minimumPosteriorOverlap()
                            || overlap2[signal2]
                                < options.minimumPosteriorOverlap())) {
                    continue;
                }
                int firstOffset = signal1 * common;
                int secondOffset = signal2 * common;
                int outputOffset = retained * common;
                hypotheses(first, firstOffset, second, secondOffset,
                    weights, options, common, hypotheses, posterior, outputOffset);
                pairs.add(new ColocSignalPair(
                    trait1.rawEffectIndices()[signal1],
                    trait2.rawEffectIndices()[signal2], common,
                    match.names().get(lead1[signal1]),
                    match.names().get(lead2[signal2]),
                    hypotheses[0], hypotheses[1], hypotheses[2],
                    hypotheses[3], hypotheses[4]));
                retained++;
            }
        }
        double[] compact = retained == totalPairs ? posterior
            : java.util.Arrays.copyOf(posterior, retained * common);
        return new ColocSusieResult(match.names(), pairs, compact,
            options, totalPairs - retained);
    }

    private static void hypotheses(
            double[] first, int firstOffset,
            double[] second, int secondOffset,
            Weights weights, ColocOptions options, int variants,
            double[] output, double[] posterior, int posteriorOffset) {
        double firstSum = Double.NEGATIVE_INFINITY;
        double secondSum = Double.NEGATIVE_INFINITY;
        double distinctSum = Double.NEGATIVE_INFINITY;
        double sharedSum = Double.NEGATIVE_INFINITY;
        double logPrior1 = Math.log(options.trait1Prior());
        double logPrior2 = Math.log(options.trait2Prior());
        double multiplier = Math.log(options.sharedPrior())
            - Math.log(options.trait1Prior())
            - Math.log(options.trait2Prior());
        for (int variant = 0; variant < variants; variant++) {
            double weighted1 = (weights == null ? logPrior1 : weights.logTrait1()[variant])
                + first[firstOffset + variant];
            double weighted2 = (weights == null ? logPrior2 : weights.logTrait2()[variant])
                + second[secondOffset + variant];
            // Enumerate i != j incrementally, without subtracting nearly equal
            // product and diagonal sums. Each ordered pair enters exactly once.
            distinctSum = logAdd(distinctSum,
                logAdd(weighted1 + secondSum, weighted2 + firstSum));
            firstSum = logAdd(firstSum, weighted1);
            secondSum = logAdd(secondSum, weighted2);
            double both = weighted1 + weighted2;
            sharedSum = logAdd(sharedSum, multiplier + both);
            posterior[posteriorOffset + variant] = multiplier + both;
        }
        output[0] = 0.0;
        output[1] = firstSum;
        output[2] = secondSum;
        output[3] = distinctSum;
        output[4] = sharedSum;
        for (double value : output) {
            if (Double.isNaN(value) || value == Double.POSITIVE_INFINITY)
                throw new IllegalArgumentException("log Bayes factors exceed numerical range");
        }
        double posteriorSum = 0;
        for (int variant = 0; variant < variants; variant++) {
            // If H4 has zero support its conditional distribution is undefined;
            // use the documented zero vector, but retain the other hypotheses.
            posterior[posteriorOffset + variant] = sharedSum == Double.NEGATIVE_INFINITY
                ? 0.0 : Math.exp(posterior[posteriorOffset + variant] - sharedSum);
            posteriorSum += posterior[posteriorOffset + variant];
        }
        if (posteriorSum > 0)
            for (int variant = 0; variant < variants; variant++)
                posterior[posteriorOffset + variant] /= posteriorSum;
        normalizeLogs(output);
    }

    private static double[] posteriorOverlap(
            ColocSusieInput input, int[] commonIndices) {
        int variables = input.variants();
        double[] values = input.rawLogBayesFactors();
        double[] result = new double[input.signals()];
        for (int signal = 0; signal < input.signals(); signal++) {
            int offset = signal * variables;
            double all = logSum(values, offset, variables);
            double common = Double.NEGATIVE_INFINITY;
            for (int index : commonIndices) {
                common = logAdd(common, values[offset + index]);
            }
            result[signal] = Math.exp(common - all);
        }
        return result;
    }

    private static double[] restrict(
            ColocSusieInput input, int[] indices, int common) {
        double[] source = input.rawLogBayesFactors();
        double[] result = new double[input.signals() * common];
        for (int signal = 0; signal < input.signals(); signal++) {
            int sourceOffset = signal * input.variants();
            int targetOffset = signal * common;
            for (int variant = 0; variant < common; variant++) {
                result[targetOffset + variant] =
                    source[sourceOffset + indices[variant]];
            }
        }
        return result;
    }

    private static int[] rowMaxima(
            double[] values, int rows, int columns) {
        int[] result = new int[rows];
        for (int row = 0; row < rows; row++) {
            int offset = row * columns;
            int maximum = 0;
            for (int column = 1; column < columns; column++) {
                if (values[offset + column] > values[offset + maximum]) {
                    maximum = column;
                }
            }
            result[row] = maximum;
        }
        return result;
    }

    private static Match match(List<String> first, List<String> second) {
        Map<String, Integer> secondIndex = new HashMap<>(second.size() * 2);
        for (int index = 0; index < second.size(); index++) {
            secondIndex.put(second.get(index), index);
        }
        List<String> names = new ArrayList<>();
        int[] firstBuffer = new int[Math.min(first.size(), second.size())];
        int[] secondBuffer = new int[firstBuffer.length];
        int found = 0;
        for (int index = 0; index < first.size(); index++) {
            Integer other = secondIndex.get(first.get(index));
            if (other != null) {
                names.add(first.get(index));
                firstBuffer[found] = index;
                secondBuffer[found] = other;
                found++;
            }
        }
        return new Match(List.copyOf(names),
            java.util.Arrays.copyOf(firstBuffer, found),
            java.util.Arrays.copyOf(secondBuffer, found));
    }

    private static Weights weights(ColocOptions options, Match match) {
        double[] supplied1 = options.rawTrait1PriorWeights();
        double[] supplied2 = options.rawTrait2PriorWeights();
        if (supplied1 == null && supplied2 == null) return null;
        int common = match.names().size();
        double[] log1 = new double[common];
        double[] log2 = new double[common];
        double sum1 = Double.NEGATIVE_INFINITY;
        double sum2 = Double.NEGATIVE_INFINITY;
        for (int variant = 0; variant < common; variant++) {
            log1[variant] = supplied1 == null ? 0.0 : Math.log(supplied1[match.first()[variant]]);
            log2[variant] = supplied2 == null ? 0.0 : Math.log(supplied2[match.second()[variant]]);
            sum1 = logAdd(sum1, log1[variant]);
            sum2 = logAdd(sum2, log2[variant]);
        }
        for (int variant = 0; variant < common; variant++) {
            log1[variant] = Math.log(common) + Math.log(options.trait1Prior())
                + (log1[variant] - sum1);
            log2[variant] = Math.log(common) + Math.log(options.trait2Prior())
                + (log2[variant] - sum2);
        }
        return new Weights(log1, log2);
    }

    private static void validateWeights(
            double[] weights, int variants, int trait) {
        if (weights != null && weights.length != variants) {
            throw new IllegalArgumentException("trait " + trait
                + " prior weight count must equal its variant count");
        }
    }

    private static double logSum(double[] values, int offset, int length) {
        double maximum = Double.NEGATIVE_INFINITY;
        for (int index = 0; index < length; index++) {
            maximum = Math.max(maximum, values[offset + index]);
        }
        if (maximum == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }
        double sum = 0.0;
        for (int index = 0; index < length; index++) {
            sum += Math.exp(values[offset + index] - maximum);
        }
        return maximum + Math.log(sum);
    }

    private static double logAdd(double first, double second) {
        if (first == Double.NEGATIVE_INFINITY) return second;
        if (second == Double.NEGATIVE_INFINITY) return first;
        double maximum = Math.max(first, second);
        return maximum + Math.log1p(Math.exp(Math.min(first, second) - maximum));
    }

    private static void normalizeLogs(double[] values) {
        double denominator = logSum(values, 0, values.length);
        for (int index = 0; index < values.length; index++) {
            values[index] = Math.exp(values[index] - denominator);
        }
    }

    private record Match(List<String> names, int[] first, int[] second) { }
    private record Weights(double[] logTrait1, double[] logTrait2) { }
}
