/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import jdistlib.T;

/** Multiplicative-dispersion PET/PEESE and Duval-Tweedie trim-and-fill. */
public final class MetaBiasCorrections {
    private MetaBiasCorrections() { }
    public enum Side { AUTO, LEFT, RIGHT }
    public enum TrimEstimator { L0, R0, Q0 }

    public static Result pet(List<MetaStudy> studies) { return regression(studies, false, 0.95); }
    public static Result peese(List<MetaStudy> studies) { return regression(studies, true, 0.95); }
    public static Result pet(List<MetaStudy> studies, double confidenceLevel) { return regression(studies, false, confidenceLevel); }
    public static Result peese(List<MetaStudy> studies, double confidenceLevel) { return regression(studies, true, confidenceLevel); }

    private static Result regression(List<MetaStudy> studies, boolean variancePredictor, double level) {
        MetaMath.Data data = MetaMath.data(studies);
        if (studies.size() < 3 || !(level > 0 && level < 1))
            throw new IllegalArgumentException("bias regression needs at least three studies and confidence in (0,1)");
        int n = studies.size();
        double[] x = new double[n * 2];
        for (int i = 0; i < n; i++) {
            x[2 * i] = 1;
            x[2 * i + 1] = variancePredictor ? data.variances()[i] : Math.sqrt(data.variances()[i]);
        }
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            MetaMath.Fit fit = MetaMath.fit(data, x, 2, 0, context.backend());
            double scale = fit.qe() / (n - 2), critical = T.quantile((1 + level) / 2, n - 2, true, false);
            double[] cov = fit.covariance().clone(), se = new double[2], lo = new double[2], hi = new double[2];
            for (int j = 0; j < cov.length; j++) cov[j] *= scale;
            for (int j = 0; j < 2; j++) {
                se[j] = Math.sqrt(cov[j * 2 + j]);
                lo[j] = fit.beta()[j] - critical * se[j]; hi[j] = fit.beta()[j] + critical * se[j];
            }
            return new Result(variancePredictor ? "PEESE" : "PET", fit.beta()[0], fit.beta()[1], n,
                AssociationStatistics.studentT(fit.beta(), se, n - 2, DegreesOfFreedomMethod.RESIDUAL), cov, lo, hi);
        }
    }

    /** Defaults match trimfill(rma(..., method="REML")): L0 and regression-selected side. */
    public static TrimFillResult trimAndFill(List<MetaStudy> studies) {
        return trimAndFill(studies, MetaAnalysisOptions.randomEffects(), Side.AUTO, TrimEstimator.L0, BackendPolicy.CPU);
    }

    public static TrimFillResult trimAndFill(List<MetaStudy> studies, MetaAnalysisOptions options,
            Side side, TrimEstimator estimator, BackendPolicy policy) {
        MetaMath.data(studies);
        if (studies.size() < 3 || options == null || side == null || estimator == null || policy == null)
            throw new IllegalArgumentException("trim-and-fill requires three studies and explicit controls");
        int n = studies.size();
        if (side == Side.AUTO) {
            double[][] se = new double[n][1];
            for (int i = 0; i < n; i++) se[i][0] = studies.get(i).standardError();
            double slope = MetaRegression.fit(studies, se, List.of("SE"), true, options, policy).beta()[1];
            side = slope < 0 ? Side.RIGHT : Side.LEFT;
        }
        double sign = side == Side.RIGHT ? -1 : 1;
        List<MetaStudy> sorted = new ArrayList<>();
        for (MetaStudy s : studies) sorted.add(new MetaStudy(s.name(), sign * s.effectSize(), s.standardError()));
        sorted.sort(Comparator.comparingDouble(MetaStudy::effectSize));
        int missing = 0, previous = -1, iteration = 0;
        double center = 0;
        while (missing != previous) {
            if (++iteration > options.maximumIterations()) throw new ArithmeticException("trim-and-fill did not converge");
            previous = missing;
            if (n - missing < 2) throw new ArithmeticException("trim-and-fill retained fewer than two studies");
            center = MetaAnalysis.fit(sorted.subList(0, n - missing), options, policy).pooledEffectSize();
            final double location = center;
            Integer[] order = new Integer[n];
            for (int i = 0; i < n; i++) order[i] = i;
            // Stable ordering reproduces ties.method="first" after sorting by effect.
            Arrays.sort(order, Comparator.comparingDouble(i -> Math.abs(sorted.get(i).effectSize() - location)));
            double positiveRankSum = 0, maxNegativeRank = 0;
            for (int rank = 0; rank < n; rank++) {
                double deviation = sorted.get(order[rank]).effectSize() - center;
                if (deviation > 0) positiveRankSum += rank + 1;
                if (deviation < 0) maxNegativeRank = rank + 1;
            }
            double count = switch (estimator) {
                case L0 -> (4 * positiveRankSum - (double) n * (n + 1)) / (2 * n - 1);
                case R0 -> n - maxNegativeRank - 1;
                case Q0 -> n - 0.5 - Math.sqrt(2.0 * n * n - 4 * positiveRankSum + 0.25);
            };
            if (!Double.isFinite(count)) throw new ArithmeticException("trim-and-fill count is undefined for these data");
            missing = (int) Math.max(0, Math.rint(count));
        }
        List<MetaStudy> augmented = new ArrayList<>(studies);
        for (int i = n - missing; i < n; i++) augmented.add(new MetaStudy("Filled " + (i - n + missing + 1),
            sign * (2 * center - sorted.get(i).effectSize()), sorted.get(i).standardError()));
        MetaAnalysisResult adjusted = MetaAnalysis.fit(augmented, options, policy);
        return new TrimFillResult(adjusted.pooledEffectSize(), missing, estimator.name(), side,
            iteration, augmented, adjusted);
    }

    public record Result(String method, double intercept, double slope, int studyCount,
            AssociationStatistics associationStatistics, double[] covariance,
            double[] confidenceLower, double[] confidenceUpper) {
        public Result { covariance = covariance.clone(); confidenceLower = confidenceLower.clone(); confidenceUpper = confidenceUpper.clone(); }
        @Override public double[] covariance() { return covariance.clone(); }
        @Override public double[] confidenceLower() { return confidenceLower.clone(); }
        @Override public double[] confidenceUpper() { return confidenceUpper.clone(); }
        public double[] standardErrors() { return associationStatistics.standardErrors(); }
        public double[] pValues() { return associationStatistics.pValues(); }
    }
    public record TrimFillResult(double adjustedEffect, int imputedStudyCount, String method,
            Side side, int iterations, List<MetaStudy> augmentedStudies, MetaAnalysisResult adjustedFit) {
        public TrimFillResult { augmentedStudies = List.copyOf(augmentedStudies); }
    }
}
