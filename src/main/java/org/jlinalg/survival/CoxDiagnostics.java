/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import jdistlib.ChiSquare;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Post-fit diagnostics for fixed-effect Cox proportional-hazards models. */
public final class CoxDiagnostics {
    private CoxDiagnostics() { }

    /**
     * Computes cluster-sandwich inference, martingale/deviance/score/dfbeta
     * residuals, Schoenfeld residuals, and log-time PH score tests.
     * A null cluster list treats each observation as an independent cluster.
     */
    public static CoxDiagnosticsResult analyze(
            CoxSurvivalData survival, double[][] covariates,
            CoxResult fit, double[] offset, List<String> clusterIds) {
        if (survival == null || covariates == null || fit == null
                || covariates.length != survival.observations())
            throw new IllegalArgumentException(
                "survival, covariates, and fitted Cox result are required");
        int rows = survival.observations();
        int columns = fit.beta().length;
        double[] design = MatrixOps.rowMajor(covariates, rows);
        if (design.length != rows * columns)
            throw new IllegalArgumentException(
                "diagnostic covariates must match fitted coefficients");
        double[] offsets = offset == null ? new double[rows] : offset.clone();
        if (offsets.length != rows)
            throw new IllegalArgumentException(
                "one diagnostic offset is required per observation");
        for (double value : offsets) if (!Double.isFinite(value))
            throw new IllegalArgumentException("offset must be finite");
        List<String> clusters = clusters(clusterIds, rows);
        Accumulation values = accumulate(survival, design, columns,
            fit.beta(), offsets, fit.ties());
        double[] bread = fit.covariance();
        double[] dfbeta = multiplyRows(values.score(), rows, columns, bread);
        double[] robust = sandwich(values.score(), rows, columns,
            clusters, bread);
        double[] robustSe = new double[columns];
        for (int column = 0; column < columns; column++)
            robustSe[column] = Math.sqrt(Math.max(0.0,
                robust[column * columns + column]));
        double[] deviance = new double[rows];
        for (int row = 0; row < rows; row++) {
            double martingale = values.martingale()[row];
            double event = survival.eventView()[row] ? 1.0 : 0.0;
            double inside = -2.0 * (martingale
                + (event == 0.0 ? 0.0
                    : Math.log(Math.max(1e-300, event - martingale))));
            deviance[row] = Math.copySign(Math.sqrt(Math.max(0.0, inside)),
                martingale);
        }
        CoxResiduals residuals = new CoxResiduals(values.martingale(),
            deviance, values.score(), dfbeta, rows, columns,
            values.schoenfeldRows(), values.schoenfeldTimes(),
            values.schoenfeld());
        return new CoxDiagnosticsResult(robust, robustSe,
            new java.util.HashSet<>(clusters).size(), residuals,
            proportionalHazards(values.groups(), columns));
    }

    static double[] scoreResiduals(
            CoxSurvivalData survival, double[] design, int columns,
            double[] beta, double[] offset, CoxTies ties) {
        return accumulate(survival, design, columns, beta, offset, ties)
            .score().clone();
    }

    private static Accumulation accumulate(
            CoxSurvivalData survival, double[] design, int columns,
            double[] beta, double[] offset, CoxTies ties) {
        int rows = survival.observations();
        double[] eta = offset.clone();
        for (int row = 0; row < rows; row++)
            for (int column = 0; column < columns; column++)
                eta[row] += design[row * columns + column] * beta[column];
        double[] martingale = new double[rows];
        double[] score = new double[rows * columns];
        List<Integer> schoenfeldRows = new ArrayList<>();
        List<Double> schoenfeldTimes = new ArrayList<>();
        List<Double> schoenfeld = new ArrayList<>();
        List<EventGroup> groups = new ArrayList<>();
        Set<Integer> strata = new TreeSet<>();
        for (int value : survival.strataView()) strata.add(value);
        for (int stratum : strata) {
            Set<Double> times = new TreeSet<>();
            for (int row = 0; row < rows; row++)
                if (survival.strataView()[row] == stratum
                        && survival.eventView()[row])
                    times.add(survival.stopView()[row]);
            for (double time : times) {
                List<Integer> riskRows = new ArrayList<>();
                List<Integer> deaths = new ArrayList<>();
                double maximumEta = Double.NEGATIVE_INFINITY;
                for (int row = 0; row < rows; row++) {
                    if (survival.strataView()[row] != stratum) continue;
                    if (survival.startView()[row] < time
                            && survival.stopView()[row] >= time) {
                        riskRows.add(row);
                        maximumEta = Math.max(maximumEta, eta[row]);
                    }
                    if (survival.eventView()[row]
                            && survival.stopView()[row] == time) deaths.add(row);
                }
                double risk0 = 0.0;
                double death0 = 0.0;
                double[] risk1 = new double[columns];
                double[] death1 = new double[columns];
                double[] risk2 = new double[columns * columns];
                double[] death2 = new double[columns * columns];
                boolean[] isDeath = new boolean[rows];
                for (int row : deaths) isDeath[row] = true;
                for (int row : riskRows) {
                    double weight = Math.exp(eta[row] - maximumEta);
                    risk0 += weight;
                    moments(design, row, columns, weight, risk1, risk2);
                    if (isDeath[row]) {
                        death0 += weight;
                        moments(design, row, columns, weight,
                            death1, death2);
                    }
                }
                int steps = ties == CoxTies.EFRON ? deaths.size() : 1;
                double[] meanSum = new double[columns];
                double[] groupInformation = new double[columns * columns];
                for (int step = 0; step < steps; step++) {
                    double fraction = ties == CoxTies.EFRON
                        ? step / (double) deaths.size() : 0.0;
                    double multiplier = ties == CoxTies.EFRON ? 1.0
                        : deaths.size();
                    double denominator = risk0 - fraction * death0;
                    double[] mean = new double[columns];
                    for (int column = 0; column < columns; column++) {
                        mean[column] = (risk1[column]
                            - fraction * death1[column]) / denominator;
                        meanSum[column] += multiplier * mean[column];
                    }
                    for (int left = 0; left < columns; left++)
                        for (int right = 0; right < columns; right++) {
                            double second = (risk2[left * columns + right]
                                - fraction * death2[left * columns + right])
                                / denominator;
                            groupInformation[left * columns + right] +=
                                multiplier * (second
                                    - mean[left] * mean[right]);
                        }
                    for (int row : riskRows) {
                        double adjustment = isDeath[row] ? 1.0 - fraction : 1.0;
                        double expected = multiplier * adjustment
                            * Math.exp(eta[row] - maximumEta) / denominator;
                        martingale[row] -= expected;
                        for (int column = 0; column < columns; column++)
                            score[row * columns + column] -= expected
                                * (design[row * columns + column]
                                    - mean[column]);
                    }
                }
                double[] groupResidual = new double[columns];
                for (int row : deaths) {
                    martingale[row] += 1.0;
                    schoenfeldRows.add(row);
                    schoenfeldTimes.add(time);
                    for (int column = 0; column < columns; column++) {
                        double residual = design[row * columns + column]
                            - meanSum[column] / deaths.size();
                        score[row * columns + column] += residual;
                        groupResidual[column] += residual;
                        schoenfeld.add(residual);
                    }
                }
                groups.add(new EventGroup(time, deaths.size(),
                    groupResidual, groupInformation));
            }
        }
        return new Accumulation(martingale, score,
            schoenfeldRows.stream().mapToInt(Integer::intValue).toArray(),
            schoenfeldTimes.stream().mapToDouble(Double::doubleValue).toArray(),
            schoenfeld.stream().mapToDouble(Double::doubleValue).toArray(),
            List.copyOf(groups));
    }

    private static CoxProportionalHazardsTest proportionalHazards(
            List<EventGroup> groups, int columns) {
        int events = groups.stream().mapToInt(EventGroup::deaths).sum();
        double mean = groups.stream().mapToDouble(group ->
            group.deaths() * Math.log(group.time())).sum() / events;
        double[] score = new double[columns];
        double[] information = new double[columns * columns];
        for (EventGroup group : groups) {
            double transformed = Math.log(group.time()) - mean;
            for (int column = 0; column < columns; column++)
                score[column] += transformed * group.residualSum()[column];
            double squared = transformed * transformed;
            for (int index = 0; index < information.length; index++)
                information[index] += squared * group.information()[index];
        }
        double[] statistics = new double[columns];
        double[] pValues = new double[columns];
        for (int column = 0; column < columns; column++) {
            double variance = information[column * columns + column];
            statistics[column] = variance > 0.0
                ? score[column] * score[column] / variance : Double.NaN;
            pValues[column] = Double.isFinite(statistics[column])
                ? ChiSquare.cumulative(statistics[column], 1.0, false, false)
                : Double.NaN;
        }
        double global;
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            double[] solved = CoxMath.solvePositive(context.backend(),
                information, columns, score, 1e-10);
            global = 0.0;
            for (int column = 0; column < columns; column++)
                global += score[column] * solved[column];
        }
        return new CoxProportionalHazardsTest(statistics, pValues, global,
            columns, ChiSquare.cumulative(global, columns, false, false),
            events, "log(time)");
    }

    private static double[] sandwich(double[] score, int rows, int columns,
            List<String> clusters, double[] bread) {
        Map<String, double[]> sums = new LinkedHashMap<>();
        for (int row = 0; row < rows; row++) {
            double[] value = sums.computeIfAbsent(clusters.get(row),
                ignored -> new double[columns]);
            for (int column = 0; column < columns; column++)
                value[column] += score[row * columns + column];
        }
        double[] meat = new double[columns * columns];
        for (double[] value : sums.values())
            for (int left = 0; left < columns; left++)
                for (int right = 0; right < columns; right++)
                    meat[left * columns + right] += value[left] * value[right];
        double[] left = matrixMultiply(bread, meat, columns);
        return matrixMultiply(left, bread, columns);
    }

    private static double[] multiplyRows(double[] rows, int rowCount,
            int columns, double[] matrix) {
        double[] result = new double[rows.length];
        for (int row = 0; row < rowCount; row++)
            for (int destination = 0; destination < columns; destination++)
                for (int source = 0; source < columns; source++)
                    result[row * columns + destination] +=
                        rows[row * columns + source]
                            * matrix[source * columns + destination];
        return result;
    }

    private static double[] matrixMultiply(
            double[] left, double[] right, int dimension) {
        double[] result = new double[dimension * dimension];
        for (int row = 0; row < dimension; row++)
            for (int column = 0; column < dimension; column++)
                for (int inner = 0; inner < dimension; inner++)
                    result[row * dimension + column] +=
                        left[row * dimension + inner]
                            * right[inner * dimension + column];
        return result;
    }

    private static void moments(double[] design, int row, int columns,
            double weight, double[] first, double[] second) {
        int start = row * columns;
        for (int left = 0; left < columns; left++) {
            first[left] += weight * design[start + left];
            for (int right = 0; right < columns; right++)
                second[left * columns + right] += weight
                    * design[start + left] * design[start + right];
        }
    }

    private static List<String> clusters(List<String> supplied, int rows) {
        if (supplied == null) {
            String[] result = new String[rows];
            for (int row = 0; row < rows; row++)
                result[row] = Integer.toString(row);
            return List.of(result);
        }
        if (supplied.size() != rows
                || supplied.stream().anyMatch(value -> value == null))
            throw new IllegalArgumentException(
                "one nonnull cluster id is required per observation");
        return List.copyOf(supplied);
    }

    private record EventGroup(double time, int deaths,
        double[] residualSum, double[] information) { }
    private record Accumulation(double[] martingale, double[] score,
        int[] schoenfeldRows, double[] schoenfeldTimes, double[] schoenfeld,
        List<EventGroup> groups) { }
}
