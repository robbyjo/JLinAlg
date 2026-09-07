/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.accelerator.CholeskyFactor;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;

/** Sandwich inference for clustered study effects and meta-regression. */
public final class MetaClusterRobust {
    private MetaClusterRobust() { }

    public static MetaClusterRobustResult fit(List<MetaStudy> studies,
            String[] clusters, double[][] moderators, List<String> names,
            MetaAnalysisOptions options, BackendPolicy backendPolicy) {
        MetaMath.Data data = MetaMath.data(studies);
        if (clusters == null || clusters.length != studies.size()
                || options == null || backendPolicy == null)
            throw new IllegalArgumentException("clusters, options, and backend policy are required");
        int moderatorsCount = moderators == null || moderators.length == 0 ? 0 : moderators[0].length;
        if (moderators == null || moderators.length != studies.size() || names == null
                || names.size() != moderatorsCount)
            throw new IllegalArgumentException("moderators and names must match studies");
        int columns = moderatorsCount + 1;
        double[] design = new double[studies.size() * columns];
        List<String> coefficientNames = new ArrayList<>(); coefficientNames.add("(Intercept)");
        coefficientNames.addAll(names);
        for (int i = 0; i < studies.size(); i++) {
            design[i * columns] = 1.0;
            if (clusters[i] == null || clusters[i].isBlank()) throw new IllegalArgumentException("cluster labels must not be blank");
            if (moderators[i] == null || moderators[i].length != moderatorsCount) throw new IllegalArgumentException("moderator rows must have equal width");
            System.arraycopy(moderators[i], 0, design, i * columns + 1, moderatorsCount);
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            double tau = MetaMath.estimateTauSquared(data, design, columns, options, context.backend());
            MetaMath.Fit fit = MetaMath.fit(data, design, columns, tau, context.backend());
            double[] meat = new double[columns * columns];
            Map<String, double[]> byCluster = new LinkedHashMap<>();
            for (int i = 0; i < studies.size(); i++) {
                double fitted = 0.0;
                for (int j = 0; j < columns; j++) fitted += design[i * columns + j] * fit.beta()[j];
                double score = fit.weights()[i] * (data.effects()[i] - fitted);
                double[] sum = byCluster.computeIfAbsent(clusters[i], key -> new double[columns]);
                for (int j = 0; j < columns; j++) sum[j] += design[i * columns + j] * score;
            }
            for (double[] score : byCluster.values()) for (int r = 0; r < columns; r++)
                for (int c = 0; c < columns; c++) meat[r * columns + c] += score[r] * score[c];
            double[] covariance = sandwich(fit.covariance(), meat, columns);
            double df = Math.max(1.0, byCluster.size() - 1.0);
            double[] errors = new double[columns];
            for (int i = 0; i < columns; i++) errors[i] = Math.sqrt(Math.max(0.0, covariance[i * columns + i]));
            AssociationStatistics stats = options.inferenceMethod() == MetaInferenceMethod.NORMAL
                ? AssociationStatistics.normal(fit.beta(), errors)
                : AssociationStatistics.studentT(fit.beta(), errors, df, DegreesOfFreedomMethod.RESIDUAL);
            return new MetaClusterRobustResult(coefficientNames, stats, covariance, tau,
                byCluster.size(), context.provenance());
        }
    }

    private static double[] sandwich(double[] bread, double[] meat, int n) {
        CholeskyFactor factor = null;
        double[] result = new double[n * n];
        // bread is symmetric; direct multiplication avoids an extra matrix object.
        for (int r = 0; r < n; r++) for (int c = 0; c < n; c++) {
            double value = 0.0;
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
                value += bread[r * n + i] * meat[i * n + j] * bread[j * n + c];
            result[r * n + c] = value;
        }
        return result;
    }
}
