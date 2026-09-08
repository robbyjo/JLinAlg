/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import jdistlib.Normal;

/**
 * Summary regression under z ~ N(R theta, R), theta_j = jointBeta_j / marginalSE_j.
 * This is exact for a Gaussian score model with known common residual variance
 * and SEs proportional to inverse genotype norms. GWAS Wald SEs give a summary
 * approximation; heterogeneous sample sizes and estimated LD add uncertainty
 * not represented here. No residual variance re-estimation or ridge is applied.
 */
public final class ConditionalAssociationModel {
    private final List<ConditionalAssociation> associations;
    private final double[][] correlation;
    private final double[] z;

    public ConditionalAssociationModel(List<ConditionalAssociation> associations,
            double[][] alleleAlignedCorrelation) {
        if (associations == null || associations.isEmpty())
            throw new IllegalArgumentException("associations are required");
        this.associations = List.copyOf(associations);
        int n = associations.size();
        if (alleleAlignedCorrelation == null || alleleAlignedCorrelation.length != n)
            throw new IllegalArgumentException("LD dimensions must match associations");
        correlation = new double[n][n];
        z = new double[n];
        var ids = new HashSet<String>();
        for (int i = 0; i < n; i++) {
            if (!ids.add(associations.get(i).variantId()))
                throw new IllegalArgumentException("variant IDs must be unique");
            if (alleleAlignedCorrelation[i] == null || alleleAlignedCorrelation[i].length != n)
                throw new IllegalArgumentException("LD must be square");
            correlation[i] = alleleAlignedCorrelation[i].clone();
            z[i] = associations.get(i).marginalEffect() / associations.get(i).marginalStandardError();
        }
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            double r = correlation[i][j];
            if (!Double.isFinite(r) || Math.abs(r) > 1.0
                    || Math.abs(r - correlation[j][i]) > 1e-12
                    || (i == j && Math.abs(r - 1.0) > 1e-12))
                throw new IllegalArgumentException("LD must be finite, symmetric, bounded and have unit diagonal");
        }
        factor(correlation); // Reject singular/indefinite LD explicitly, never silently regularize.
    }

    public List<ConditionalAssociation> associations() { return associations; }

    /** Tests each SNP given C; a member of C is tested given C minus itself. */
    public List<ConditionalAssociationResult> condition(int[] conditioningIndices,
            double confidenceLevel) {
        if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0))
            throw new IllegalArgumentException("confidence must be in (0,1)");
        if (conditioningIndices == null) throw new IllegalArgumentException("conditioning indices are required");
        var unique = new HashSet<Integer>();
        for (int index : conditioningIndices) if (index < 0 || index >= z.length || !unique.add(index))
            throw new IllegalArgumentException("conditioning indices must be unique and in range");
        int k = conditioningIndices.length;
        double[][] sub = new double[k][k];
        double[] zs = new double[k];
        for (int i = 0; i < k; i++) {
            zs[i] = z[conditioningIndices[i]];
            for (int j = 0; j < k; j++) sub[i][j] = correlation[conditioningIndices[i]][conditioningIndices[j]];
        }
        double[][] lower = factor(sub);
        double[] joint = solve(lower, zs);
        double critical = Normal.quantile(0.5 + confidenceLevel / 2.0, 0, 1, true, false);
        List<ConditionalAssociationResult> result = new ArrayList<>();
        for (int target = 0; target < z.length; target++) {
            int position = -1;
            List<String> ids = new ArrayList<>();
            for (int j = 0; j < k; j++) {
                if (conditioningIndices[j] == target) position = j;
                else ids.add(associations.get(conditioningIndices[j]).variantId());
            }
            double residualVariance, theta;
            if (position >= 0) {
                double[] unit = new double[k]; unit[position] = 1;
                residualVariance = 1.0 / solve(lower, unit)[position];
                theta = joint[position];
            } else {
                double[] r = new double[k];
                for (int j = 0; j < k; j++) r[j] = correlation[target][conditioningIndices[j]];
                double[] projected = solve(lower, r);
                residualVariance = 1.0;
                double residualScore = z[target];
                for (int j = 0; j < k; j++) {
                    residualVariance -= r[j] * projected[j];
                    residualScore -= r[j] * joint[j];
                }
                theta = residualScore / residualVariance;
            }
            if (!(residualVariance > 1e-12))
                throw new IllegalArgumentException("conditional genotype information is numerically singular");
            double scale = associations.get(target).marginalStandardError();
            double beta = scale * theta, se = scale / Math.sqrt(residualVariance);
            double statistic = theta * Math.sqrt(residualVariance);
            if (!Double.isFinite(beta) || !Double.isFinite(se)
                    || !Double.isFinite(beta - critical * se) || !Double.isFinite(beta + critical * se))
                throw new IllegalArgumentException("conditional estimates or intervals exceed finite numeric range");
            result.add(new ConditionalAssociationResult(associations.get(target).variantId(),
                beta, se, statistic, 2 * Normal.cumulative(-Math.abs(statistic), 0, 1, true, false),
                beta - critical * se, beta + critical * se, residualVariance, ids));
        }
        return List.copyOf(result);
    }

    private static double[][] factor(double[][] matrix) {
        int n = matrix.length;
        double[][] lower = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
            double value = matrix[i][j];
            for (int k = 0; k < j; k++) value -= lower[i][k] * lower[j][k];
            if (i == j) {
                if (!(value > 1e-12)) throw new IllegalArgumentException("LD must be positive definite and numerically nonsingular");
                lower[i][j] = Math.sqrt(value);
            } else lower[i][j] = value / lower[j][j];
        }
        return lower;
    }

    private static double[] solve(double[][] lower, double[] rhs) {
        double[] x = rhs.clone();
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < i; j++) x[i] -= lower[i][j] * x[j];
            x[i] /= lower[i][i];
        }
        for (int i = x.length - 1; i >= 0; i--) {
            for (int j = i + 1; j < x.length; j++) x[i] -= lower[j][i] * x[j];
            x[i] /= lower[i][i];
        }
        return x;
    }
}
