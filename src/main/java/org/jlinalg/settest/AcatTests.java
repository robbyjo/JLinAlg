/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import jdistlib.Normal;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.pipeline.VariantFilterResult;

/** ACAT-V and canonical six-component ACAT-O score tests. */
public final class AcatTests {
    public static final double DEFAULT_MAC_THRESHOLD = 10;

    private AcatTests() { }

    /** Published ACAT-V with Beta(MAF; a1,a2) variant coefficients. */
    public static AcatVResult acatV(
            String id, SetTestScoreState state, double[] minorAlleleFrequencies,
            double[] minorAlleleCounts, double shape1, double shape2,
            double minorAlleleCountThreshold) {
        validateInputs(state, minorAlleleFrequencies, minorAlleleCounts,
            minorAlleleCountThreshold);
        return acatVBeta(id, state.variants(), state,
            minorAlleleFrequencies, minorAlleleCounts, shape1, shape2,
            minorAlleleCountThreshold, List.of(), AcatTests::normalPValue);
    }

    /** Published default ACAT-V with Beta(MAF; 1,25) coefficients. */
    public static AcatVResult acatV(
            String id, SetTestScoreState state, double[] minorAlleleFrequencies,
            double[] minorAlleleCounts) {
        return acatV(id, state, minorAlleleFrequencies, minorAlleleCounts,
            1, 25, DEFAULT_MAC_THRESHOLD);
    }

    /**
     * Custom ACAT-V. The same supplied weights are used as burden coefficients
     * and as the Cauchy weights, matching the reference ACAT R interface.
     */
    public static AcatVResult acatV(
            String id, SetTestScoreState state, double[] minorAlleleFrequencies,
            double[] minorAlleleCounts, double[] weights,
            double minorAlleleCountThreshold) {
        validateInputs(state, minorAlleleFrequencies, minorAlleleCounts,
            minorAlleleCountThreshold);
        validateWeights(weights, state.variants());
        return acatVCustom(id, state.variants(), state,
            minorAlleleFrequencies, minorAlleleCounts, weights, weights,
            minorAlleleCountThreshold, List.of(), AcatTests::normalPValue);
    }

    /** Canonical ACAT-O combining Burden, SKAT, and ACAT-V at two weights. */
    public static AcatOResult acatO(
            String id, SetTestScoreState state, double[] minorAlleleFrequencies,
            double[] minorAlleleCounts, double minorAlleleCountThreshold) {
        validateInputs(state, minorAlleleFrequencies, minorAlleleCounts,
            minorAlleleCountThreshold);
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            return acatO(id, state.variants(), state,
                minorAlleleFrequencies, minorAlleleCounts,
                minorAlleleCountThreshold, List.of(), AcatTests::normalPValue,
                context.backend());
        }
    }

    /** Canonical ACAT-O with the published MAC threshold of ten. */
    public static AcatOResult acatO(
            String id, SetTestScoreState state, double[] minorAlleleFrequencies,
            double[] minorAlleleCounts) {
        return acatO(id, state, minorAlleleFrequencies, minorAlleleCounts,
            DEFAULT_MAC_THRESHOLD);
    }

    static AcatVResult acatVBeta(
            String id, int requested, SetTestScoreState state,
            double[] mafs, double[] macs, double shape1, double shape2,
            double threshold, List<VariantFilterResult> excluded,
            DoubleUnaryOperator tail) {
        validateInputs(state, mafs, macs, threshold);
        double[] coefficients = new double[state.variants()];
        double[] pWeights = new double[state.variants()];
        for (int index = 0; index < coefficients.length; index++) {
            coefficients[index] = VariantWeights.betaBurden(
                mafs[index], shape1, shape2);
            pWeights[index] = coefficientPWeight(
                coefficients[index], mafs[index]);
        }
        int collapsed = collapsed(macs, threshold);
        double collapsedWeight = Double.NaN;
        if (collapsed > 0) {
            double meanMaf = 0;
            for (int index = 0; index < mafs.length; index++)
                if (macs[index] <= threshold) meanMaf += mafs[index];
            meanMaf /= collapsed;
            double coefficient = VariantWeights.betaBurden(
                meanMaf, shape1, shape2);
            collapsedWeight = coefficientPWeight(coefficient, meanMaf);
        }
        return acatV(id, requested, state, mafs, macs, coefficients,
            pWeights, collapsedWeight, threshold, excluded, tail,
            "beta(" + shape1 + "," + shape2 + ")", null);
    }

    static AcatVResult acatVCustom(
            String id, int requested, SetTestScoreState state,
            double[] mafs, double[] macs, double[] burdenWeights,
            double[] pWeights, double threshold,
            List<VariantFilterResult> excluded, DoubleUnaryOperator tail) {
        validateInputs(state, mafs, macs, threshold);
        validateWeights(burdenWeights, state.variants());
        validateWeights(pWeights, state.variants());
        int collapsed = collapsed(macs, threshold);
        double collapsedWeight = 0;
        if (collapsed > 0) {
            for (int index = 0; index < pWeights.length; index++)
                if (macs[index] <= threshold)
                    collapsedWeight += pWeights[index];
            collapsedWeight /= collapsed;
        }
        return acatV(id, requested, state, mafs, macs, burdenWeights,
            pWeights, collapsedWeight, threshold, excluded, tail, "custom", null);
    }

    static AcatOResult acatO(
            String id, int requested, SetTestScoreState state,
            double[] mafs, double[] macs, double threshold,
            List<VariantFilterResult> excluded, DoubleUnaryOperator tail,
            ComputeBackend backend) {
        validateInputs(state, mafs, macs, threshold);
        if (backend == null) throw new IllegalArgumentException("compute backend is required");
        List<AcatOResult.Component> components = new ArrayList<>(6);
        AcatVResult[] acatV = new AcatVResult[2];
        double[][] shapes = {{1, 25}, {1, 1}};
        double[] skatP = new double[2], burdenP = new double[2];
        for (int scheme = 0; scheme < shapes.length; scheme++) {
            double[] coefficients = new double[state.variants()];
            double[] pWeights = new double[state.variants()];
            for (int index = 0; index < coefficients.length; index++)
            {
                coefficients[index] = VariantWeights.betaBurden(
                    mafs[index], shapes[scheme][0], shapes[scheme][1]);
                pWeights[index] = coefficientPWeight(
                    coefficients[index], mafs[index]);
            }
            SetTestScoreState weighted = weighted(state, coefficients);
            burdenP[scheme] = burdenPValue(weighted, null, tail);
            skatP[scheme] = SetTests.kernelResult(id,
                "skat-beta-" + (int) shapes[scheme][0] + "-"
                    + (int) shapes[scheme][1],
                requested, state.variants(), weighted, excluded, backend).pValue();
            int collapsed = collapsed(macs, threshold);
            double collapsedWeight = Double.NaN;
            if (collapsed > 0) {
                double meanMaf = 0;
                for (int index = 0; index < mafs.length; index++)
                    if (macs[index] <= threshold) meanMaf += mafs[index];
                meanMaf /= collapsed;
                double coefficient = VariantWeights.betaBurden(meanMaf,
                    shapes[scheme][0], shapes[scheme][1]);
                collapsedWeight = coefficientPWeight(coefficient, meanMaf);
            }
            acatV[scheme] = acatV(id, requested, state, mafs, macs,
                coefficients, pWeights, collapsedWeight, threshold, excluded,
                tail, "beta(" + shapes[scheme][0] + ","
                    + shapes[scheme][1] + ")", weighted);
        }
        components.add(new AcatOResult.Component("skat-beta-1-25", skatP[0]));
        components.add(new AcatOResult.Component("skat-beta-1-1", skatP[1]));
        components.add(new AcatOResult.Component("burden-beta-1-25", burdenP[0]));
        components.add(new AcatOResult.Component("burden-beta-1-1", burdenP[1]));
        components.add(new AcatOResult.Component("acat-v-beta-1-25", acatV[0].pValue()));
        components.add(new AcatOResult.Component("acat-v-beta-1-1", acatV[1].pValue()));
        double[] pValues = components.stream()
            .mapToDouble(AcatOResult.Component::pValue).toArray();
        Acat.Result result = Acat.combine(pValues);
        return new AcatOResult(id, requested, state.variants(), threshold,
            components, result.statistic(), result.pValue(),
            result.log10PValue(), excluded);
    }

    private static AcatVResult acatV(
            String id, int requested, SetTestScoreState state,
            double[] mafs, double[] macs, double[] burdenWeights,
            double[] pWeights, double collapsedWeight, double threshold,
            List<VariantFilterResult> excluded, DoubleUnaryOperator tail,
            String weighting, SetTestScoreState preparedWeighted) {
        SetTestScoreState weighted = preparedWeighted == null
            ? weighted(state, burdenWeights) : preparedWeighted;
        List<String> names = new ArrayList<>();
        List<Integer> counts = new ArrayList<>();
        List<Double> pValues = new ArrayList<>(), combinationWeights = new ArrayList<>();
        int collapsed = collapsed(macs, threshold);
        if (collapsed > 0) {
            boolean[] selected = new boolean[state.variants()];
            for (int index = 0; index < selected.length; index++)
                selected[index] = macs[index] <= threshold;
            double p = burdenPValue(weighted, selected, tail);
            names.add("ultra-rare-burden");
            counts.add(collapsed);
            pValues.add(p);
            combinationWeights.add(collapsedWeight);
        }
        double[] scores = state.scoresView(), information = state.informationView();
        int dimension = state.variants();
        for (int index = 0; index < dimension; index++) {
            if (macs[index] <= threshold) continue;
            double variance = information[index * dimension + index];
            if (!(variance > 0)) {
                if (scores[index] == 0) continue;
                throw new IllegalArgumentException(
                    "ACAT-V score has nonpositive marginal information");
            }
            names.add("variant-" + (index + 1));
            counts.add(1);
            pValues.add(tail.applyAsDouble(scores[index] / Math.sqrt(variance)));
            combinationWeights.add(pWeights[index]);
        }
        // The reference ACAT-V implementation removes exact-one components
        // in the mixed burden-plus-marginal branch before combining them.
        if (collapsed > 0 && collapsed < dimension) {
            for (int index = pValues.size() - 1; index >= 0; index--) {
                if (pValues.get(index) != 1) continue;
                names.remove(index);
                counts.remove(index);
                pValues.remove(index);
                combinationWeights.remove(index);
            }
        }
        if (pValues.isEmpty())
            throw new IllegalArgumentException("ACAT-V has no informative components");
        double[] numericP = pValues.stream().mapToDouble(Double::doubleValue).toArray();
        double[] numericWeights = combinationWeights.stream()
            .mapToDouble(Double::doubleValue).toArray();
        double[] normalized = Acat.normalizedWeights(numericP, numericWeights);
        List<AcatVResult.Component> components = new ArrayList<>(numericP.length);
        for (int index = 0; index < numericP.length; index++)
            components.add(new AcatVResult.Component(names.get(index),
                counts.get(index), numericP[index], normalized[index]));
        Acat.Result result = Acat.combine(numericP, numericWeights);
        return new AcatVResult(id, requested, state.variants(), weighting,
            threshold, collapsed, components, result.statistic(),
            result.pValue(), result.log10PValue(), excluded);
    }

    private static SetTestScoreState weighted(
            SetTestScoreState state, double[] weights) {
        int dimension = state.variants();
        double[] scores = state.scores(), information = state.information();
        for (int row = 0; row < dimension; row++) {
            scores[row] *= weights[row];
            for (int column = 0; column < dimension; column++)
                information[row * dimension + column] *=
                    weights[row] * weights[column];
        }
        return new SetTestScoreState(scores, information, dimension);
    }

    private static double burdenPValue(
            SetTestScoreState weighted, boolean[] selected,
            DoubleUnaryOperator tail) {
        int dimension = weighted.variants();
        double score = 0, variance = 0, scale = 0;
        for (int row = 0; row < dimension; row++) {
            if (selected != null && !selected[row]) continue;
            score += weighted.scoresView()[row];
            for (int column = 0; column < dimension; column++) {
                if (selected != null && !selected[column]) continue;
                double value = weighted.informationView()[row * dimension + column];
                variance += value;
                scale += Math.abs(value);
            }
        }
        if (!(variance > 1e-12 * scale)) {
            if (Math.abs(score) <= 1e-12 * Math.sqrt(Math.max(1, scale))) return 1;
            throw new IllegalArgumentException("ACAT burden has no positive information");
        }
        double p = tail.applyAsDouble(score / Math.sqrt(variance));
        if (!Double.isFinite(p) || p < 0 || p > 1)
            throw new IllegalArgumentException("ACAT tail probability is invalid");
        return p;
    }

    private static double normalPValue(double statistic) {
        return 2 * Normal.cumulative(-Math.abs(statistic), 0, 1, true, false);
    }

    private static double coefficientPWeight(double coefficient, double maf) {
        return coefficient * coefficient * maf * (1 - maf);
    }

    private static int collapsed(double[] macs, double threshold) {
        int result = 0;
        for (double mac : macs) if (mac <= threshold) result++;
        return result;
    }

    private static void validateInputs(
            SetTestScoreState state, double[] mafs, double[] macs,
            double threshold) {
        if (state == null || mafs == null || macs == null
                || mafs.length != state.variants()
                || macs.length != state.variants())
            throw new IllegalArgumentException(
                "ACAT requires one MAF and MAC per score statistic");
        if (!Double.isFinite(threshold) || threshold < 0)
            throw new IllegalArgumentException("ACAT MAC threshold must be nonnegative");
        for (int index = 0; index < mafs.length; index++) {
            if (!Double.isFinite(mafs[index]) || mafs[index] <= 0
                    || mafs[index] > 0.5 || !Double.isFinite(macs[index])
                    || macs[index] <= 0)
                throw new IllegalArgumentException(
                    "ACAT MAFs must be in (0,0.5] and MACs must be positive");
        }
        validateScoreMoments(state);
    }

    // O(m^2) checks preserve ACAT-V's fast path; SKAT additionally performs
    // a full eigensystem/PSD check when it is part of ACAT-O.
    private static void validateScoreMoments(SetTestScoreState state) {
        int dimension = state.variants();
        double[] scores = state.scoresView(), information = state.informationView();
        double scale = 0;
        for (double value : information) {
            if (!Double.isFinite(value))
                throw new IllegalArgumentException("ACAT score covariance is nonfinite");
            scale = Math.max(scale, Math.abs(value));
        }
        for (double score : scores) if (!Double.isFinite(score))
            throw new IllegalArgumentException("ACAT score is nonfinite");
        if (!(scale > 0)) throw new IllegalArgumentException("ACAT score has no information");
        double tolerance = 1e-12 * scale;
        for (int row = 0; row < dimension; row++) {
            double diagonal = information[row * dimension + row];
            if (diagonal < -tolerance)
                throw new IllegalArgumentException("ACAT score variance is negative");
            if (diagonal <= tolerance && Math.abs(scores[row]) > 1e-7 * Math.sqrt(scale))
                throw new IllegalArgumentException(
                    "ACAT score is inconsistent with zero marginal information");
            for (int column = 0; column < row; column++) {
                double left = information[row * dimension + column];
                double right = information[column * dimension + row];
                if (Math.abs(left - right) > tolerance)
                    throw new IllegalArgumentException(
                        "ACAT score covariance is not symmetric");
                double bound = Math.sqrt(Math.max(0, diagonal)
                    * Math.max(0, information[column * dimension + column]));
                if (Math.abs(0.5 * (left + right)) > bound + tolerance)
                    throw new IllegalArgumentException(
                        "ACAT score covariance violates a principal minor");
            }
        }
    }

    private static void validateWeights(double[] weights, int variants) {
        if (weights == null || weights.length != variants)
            throw new IllegalArgumentException("ACAT requires one weight per variant");
        boolean positive = false;
        for (double weight : weights) {
            if (!Double.isFinite(weight) || weight < 0)
                throw new IllegalArgumentException(
                    "ACAT weights must be finite and nonnegative");
            positive |= weight > 0;
        }
        if (!positive) throw new IllegalArgumentException(
            "at least one ACAT weight must be positive");
    }
}
