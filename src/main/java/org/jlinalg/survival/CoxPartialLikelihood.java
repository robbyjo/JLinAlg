/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/** Shared Cox partial-likelihood risk-set accumulation. */
final class CoxPartialLikelihood {
    private CoxPartialLikelihood() { }

    record Evaluation(double logLikelihood, double[] score, double[] information) { }

    static Evaluation evaluate(
            CoxSurvivalData survival,
            double[] design, int columns,
            double[] coefficients, double[] offset,
            CoxTies ties) {
        int rows = survival.observations();
        double[] eta = linearPredictor(
            design, rows, columns, coefficients, offset);
        if (rightCensored(survival))
            return evaluateRightCensored(
                survival, design, columns, eta, ties);
        double logLikelihood = 0;
        double[] score = new double[columns];
        double[] information = new double[columns * columns];
        int[] strata = survival.strataView();
        for (int stratum : distinctStrata(strata)) {
            for (double time : eventTimes(survival, stratum, false)) {
                int deaths = 0;
                double maximumEta = Double.NEGATIVE_INFINITY;
                for (int row = 0; row < rows; row++) {
                    if (strata[row] == stratum && atRisk(survival, row, time))
                        maximumEta = Math.max(maximumEta, eta[row]);
                    if (strata[row] == stratum
                            && survival.eventView()[row]
                            && survival.stopView()[row] == time)
                        deaths++;
                }
                if (deaths == 0 || !Double.isFinite(maximumEta)) continue;
                double risk0 = 0;
                double death0 = 0;
                double[] risk1 = new double[columns];
                double[] death1 = new double[columns];
                double[] risk2 = new double[columns * columns];
                double[] death2 = new double[columns * columns];
                for (int row = 0; row < rows; row++) {
                    if (strata[row] != stratum) continue;
                    boolean death = survival.eventView()[row]
                        && survival.stopView()[row] == time;
                    if (death) {
                        logLikelihood += eta[row];
                        for (int column = 0; column < columns; column++)
                            score[column] += design[row * columns + column];
                    }
                    if (!atRisk(survival, row, time)) continue;
                    double risk = Math.exp(eta[row] - maximumEta);
                    risk0 += risk;
                    addMoments(design, row, columns, risk, risk1, risk2);
                    if (death) {
                        death0 += risk;
                        addMoments(design, row, columns, risk,
                            death1, death2);
                    }
                }
                int steps = ties == CoxTies.EFRON ? deaths : 1;
                for (int step = 0; step < steps; step++) {
                    double fraction = ties == CoxTies.EFRON
                        ? step / (double) deaths : 0;
                    double multiplier = ties == CoxTies.EFRON ? 1 : deaths;
                    double denominator = risk0 - fraction * death0;
                    if (!(denominator > 0) || !Double.isFinite(denominator))
                        throw new IllegalArgumentException(
                            "Cox risk-set denominator is nonpositive");
                    logLikelihood -= multiplier
                        * (maximumEta + Math.log(denominator));
                    double[] mean = new double[columns];
                    for (int column = 0; column < columns; column++) {
                        mean[column] = (risk1[column]
                            - fraction * death1[column]) / denominator;
                        score[column] -= multiplier * mean[column];
                    }
                    for (int left = 0; left < columns; left++) {
                        for (int right = 0; right < columns; right++) {
                            double second = (risk2[left * columns + right]
                                - fraction * death2[left * columns + right])
                                / denominator;
                            information[left * columns + right] += multiplier
                                * (second - mean[left] * mean[right]);
                        }
                    }
                }
            }
        }
        return new Evaluation(logLikelihood, score, information);
    }

    static List<BaselineHazardPoint> baseline(
            CoxSurvivalData survival,
            double[] design, int columns,
            double[] coefficients, double[] offset,
            CoxTies ties) {
        int rows = survival.observations();
        double[] eta = linearPredictor(
            design, rows, columns, coefficients, offset);
        if (rightCensored(survival))
            return baselineRightCensored(survival, eta, ties);
        List<BaselineHazardPoint> result = new ArrayList<>();
        int[] strata = survival.strataView();
        for (int stratum : distinctStrata(strata)) {
            double cumulative = 0;
            for (double time : eventTimes(survival, stratum, true)) {
                int deaths = 0;
                double maximumEta = Double.NEGATIVE_INFINITY;
                for (int row = 0; row < rows; row++) {
                    if (strata[row] == stratum && atRisk(survival, row, time))
                        maximumEta = Math.max(maximumEta, eta[row]);
                    if (strata[row] == stratum
                            && survival.eventView()[row]
                            && survival.stopView()[row] == time)
                        deaths++;
                }
                double risk0 = 0;
                double death0 = 0;
                for (int row = 0; row < rows; row++) {
                    if (strata[row] != stratum
                            || !atRisk(survival, row, time)) continue;
                    double value = Math.exp(eta[row] - maximumEta);
                    risk0 += value;
                    if (survival.eventView()[row]
                            && survival.stopView()[row] == time)
                        death0 += value;
                }
                double increment = 0;
                if (ties == CoxTies.EFRON) {
                    for (int step = 0; step < deaths; step++)
                        increment += Math.exp(-maximumEta)
                            / (risk0 - step / (double) deaths * death0);
                } else {
                    increment = deaths * Math.exp(-maximumEta) / risk0;
                }
                cumulative += increment;
                result.add(new BaselineHazardPoint(stratum, time, deaths,
                    increment, cumulative, Math.exp(-cumulative)));
            }
        }
        return List.copyOf(result);
    }

    private static Evaluation evaluateRightCensored(
            CoxSurvivalData survival,
            double[] design, int columns,
            double[] eta, CoxTies ties) {
        double logLikelihood = 0;
        double[] score = new double[columns];
        double[] information = new double[columns * columns];
        int[] strata = survival.strataView();
        for (int stratum : distinctStrata(strata)) {
            List<Integer> ordered = rowsInStratum(survival, stratum);
            ordered.sort(Comparator.comparingDouble(
                (Integer row) -> survival.stopView()[row]).reversed());
            Map<Double, List<Integer>> deaths = deathsByTime(
                survival, stratum);
            List<Double> times = new ArrayList<>(deaths.keySet());
            times.sort(Comparator.reverseOrder());
            int riskPointer = 0;
            double maximumEta = Double.NEGATIVE_INFINITY;
            double risk0 = 0;
            double[] risk1 = new double[columns];
            double[] risk2 = new double[columns * columns];
            for (double time : times) {
                while (riskPointer < ordered.size()
                        && survival.stopView()[ordered.get(riskPointer)] >= time) {
                    int row = ordered.get(riskPointer++);
                    if (eta[row] > maximumEta) {
                        double factor = Double.isFinite(maximumEta)
                            ? Math.exp(maximumEta - eta[row]) : 0;
                        risk0 *= factor;
                        scale(risk1, factor);
                        scale(risk2, factor);
                        maximumEta = eta[row];
                    }
                    double risk = Math.exp(eta[row] - maximumEta);
                    risk0 += risk;
                    addMoments(design, row, columns, risk, risk1, risk2);
                }
                List<Integer> eventRows = deaths.get(time);
                int deathCount = eventRows.size();
                double death0 = 0;
                double[] death1 = new double[columns];
                double[] death2 = new double[columns * columns];
                for (int row : eventRows) {
                    logLikelihood += eta[row];
                    for (int column = 0; column < columns; column++)
                        score[column] += design[row * columns + column];
                    double risk = Math.exp(eta[row] - maximumEta);
                    death0 += risk;
                    addMoments(design, row, columns, risk, death1, death2);
                }
                int steps = ties == CoxTies.EFRON ? deathCount : 1;
                for (int step = 0; step < steps; step++) {
                    double fraction = ties == CoxTies.EFRON
                        ? step / (double) deathCount : 0;
                    double multiplier = ties == CoxTies.EFRON ? 1 : deathCount;
                    double denominator = risk0 - fraction * death0;
                    logLikelihood -= multiplier
                        * (maximumEta + Math.log(denominator));
                    double[] mean = new double[columns];
                    for (int column = 0; column < columns; column++) {
                        mean[column] = (risk1[column]
                            - fraction * death1[column]) / denominator;
                        score[column] -= multiplier * mean[column];
                    }
                    for (int left = 0; left < columns; left++)
                        for (int right = 0; right < columns; right++) {
                            double second = (risk2[left * columns + right]
                                - fraction * death2[left * columns + right])
                                / denominator;
                            information[left * columns + right] += multiplier
                                * (second - mean[left] * mean[right]);
                        }
                }
            }
        }
        return new Evaluation(logLikelihood, score, information);
    }

    private static List<BaselineHazardPoint> baselineRightCensored(
            CoxSurvivalData survival, double[] eta, CoxTies ties) {
        List<BaselineHazardPoint> result = new ArrayList<>();
        int[] strata = survival.strataView();
        for (int stratum : distinctStrata(strata)) {
            List<Integer> ordered = rowsInStratum(survival, stratum);
            ordered.sort(Comparator.comparingDouble(
                (Integer row) -> survival.stopView()[row]).reversed());
            Map<Double, List<Integer>> deaths = deathsByTime(
                survival, stratum);
            List<Double> times = new ArrayList<>(deaths.keySet());
            times.sort(Comparator.reverseOrder());
            int riskPointer = 0;
            double maximumEta = Double.NEGATIVE_INFINITY;
            double risk0 = 0;
            List<RawHazard> descending = new ArrayList<>();
            for (double time : times) {
                while (riskPointer < ordered.size()
                        && survival.stopView()[ordered.get(riskPointer)] >= time) {
                    int row = ordered.get(riskPointer++);
                    if (eta[row] > maximumEta) {
                        risk0 *= Double.isFinite(maximumEta)
                            ? Math.exp(maximumEta - eta[row]) : 0;
                        maximumEta = eta[row];
                    }
                    risk0 += Math.exp(eta[row] - maximumEta);
                }
                List<Integer> eventRows = deaths.get(time);
                int deathCount = eventRows.size();
                double death0 = 0;
                for (int row : eventRows)
                    death0 += Math.exp(eta[row] - maximumEta);
                double increment = 0;
                if (ties == CoxTies.EFRON)
                    for (int step = 0; step < deathCount; step++)
                        increment += Math.exp(-maximumEta)
                            / (risk0 - step / (double) deathCount * death0);
                else
                    increment = deathCount * Math.exp(-maximumEta) / risk0;
                descending.add(new RawHazard(time, deathCount, increment));
            }
            double cumulative = 0;
            for (int index = descending.size() - 1; index >= 0; index--) {
                RawHazard value = descending.get(index);
                cumulative += value.increment();
                result.add(new BaselineHazardPoint(stratum, value.time(),
                    value.events(), value.increment(), cumulative,
                    Math.exp(-cumulative)));
            }
        }
        return List.copyOf(result);
    }

    private static List<Integer> rowsInStratum(
            CoxSurvivalData survival, int stratum) {
        List<Integer> result = new ArrayList<>();
        for (int row = 0; row < survival.observations(); row++)
            if (survival.strataView()[row] == stratum) result.add(row);
        return result;
    }

    private static Map<Double, List<Integer>> deathsByTime(
            CoxSurvivalData survival, int stratum) {
        Map<Double, List<Integer>> result = new HashMap<>();
        for (int row = 0; row < survival.observations(); row++)
            if (survival.strataView()[row] == stratum
                    && survival.eventView()[row])
                result.computeIfAbsent(survival.stopView()[row],
                    ignored -> new ArrayList<>()).add(row);
        return result;
    }

    private static boolean rightCensored(CoxSurvivalData survival) {
        for (double value : survival.startView()) if (value != 0) return false;
        return true;
    }

    private static void scale(double[] values, double factor) {
        for (int index = 0; index < values.length; index++)
            values[index] *= factor;
    }

    private static double[] linearPredictor(
            double[] design, int rows, int columns,
            double[] coefficients, double[] offset) {
        double[] result = offset.clone();
        for (int row = 0; row < rows; row++)
            for (int column = 0; column < columns; column++)
                result[row] += design[row * columns + column]
                    * coefficients[column];
        return result;
    }

    private static void addMoments(
            double[] design, int row, int columns, double weight,
            double[] first, double[] second) {
        for (int left = 0; left < columns; left++) {
            double leftValue = design[row * columns + left];
            first[left] += weight * leftValue;
            for (int right = 0; right < columns; right++)
                second[left * columns + right] += weight * leftValue
                    * design[row * columns + right];
        }
    }

    private static boolean atRisk(
            CoxSurvivalData survival, int row, double time) {
        return survival.startView()[row] < time
            && survival.stopView()[row] >= time;
    }

    private static Set<Integer> distinctStrata(int[] strata) {
        Set<Integer> result = new TreeSet<>();
        for (int stratum : strata) result.add(stratum);
        return result;
    }

    private static Set<Double> eventTimes(
            CoxSurvivalData survival, int stratum, boolean ascending) {
        TreeSet<Double> result = new TreeSet<>();
        for (int row = 0; row < survival.observations(); row++)
            if (survival.strataView()[row] == stratum
                    && survival.eventView()[row])
                result.add(survival.stopView()[row]);
        return ascending ? result : result.descendingSet();
    }

    private record RawHazard(double time, int events, double increment) { }
}
