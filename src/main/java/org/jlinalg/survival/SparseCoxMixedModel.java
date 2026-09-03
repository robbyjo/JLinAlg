/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.TreeMap;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.matrix.CsrMatrix;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

/**
 * Block-plus-dense Laplace Cox model for a pedigree frailty and small ordinary
 * random-effect terms. Pedigree precision is factored directly in coefficient
 * space without materializing an observation-scale covariance.
 *
 * <p>The sparse kernel currently supports one-stratum right-censored data,
 * distinct event times, and one unit-valued pedigree incidence per row.</p>
 */
public final class SparseCoxMixedModel {
    private SparseCoxMixedModel() { }

    public static CoxMixedResult fit(
            CoxSurvivalData survival, double[][] fixedEffects,
            PedigreeRandomEffectTerm pedigreeEffect,
            List<RandomEffectTerm> ordinaryEffects,
            double[] offset, CoxMixedOptions options,
            BackendPolicy backendPolicy) {
        try (Prepared prepared = prepare(survival, pedigreeEffect,
                ordinaryEffects, options, backendPolicy)) {
            return prepared.fit(fixedEffects, offset);
        }
    }

    public static Prepared prepare(
            CoxSurvivalData survival,
            PedigreeRandomEffectTerm pedigreeEffect,
            List<RandomEffectTerm> ordinaryEffects,
            CoxMixedOptions options, BackendPolicy backendPolicy) {
        return new Prepared(survival, pedigreeEffect, ordinaryEffects,
            options, backendPolicy);
    }

    /** Reusable single-worker sparse pedigree Cox scan. */
    public static final class Prepared implements AutoCloseable {
        private final CoxSurvivalData survival;
        private final CoxRiskSetPlan riskSets;
        private final SurvivalPlan survivalPlan;
        private final RandomEffectTerm pedigreeDesign;
        private final SparsePrecisionMatrix precision;
        private final List<RandomEffectTerm> ordinaryEffects;
        private final CoxMixedOptions options;
        private final int rows;
        private final int sparseColumns;
        private final int ordinaryColumns;
        private final int[] observationCoefficient;
        private final int[] coefficientObservation;
        private final int[] ordinaryStarts;
        private final double[] ordinaryDesign;
        private final SparsePattern pattern;
        private final BackendContext context;
        private final ComputeBackend backend;
        private final PreparedSparseCholesky factor;
        private final double precisionLogDeterminant;
        private double[] sharedLogVariances;
        private boolean closed;

        private Prepared(
                CoxSurvivalData survival,
                PedigreeRandomEffectTerm pedigreeEffect,
                List<RandomEffectTerm> ordinaryEffects,
                CoxMixedOptions options, BackendPolicy backendPolicy) {
            if (survival == null || pedigreeEffect == null
                    || ordinaryEffects == null || options == null
                    || backendPolicy == null)
                throw new IllegalArgumentException(
                    "sparse pedigree Cox inputs are required");
            rows = survival.observations();
            this.survival = survival;
            riskSets = CoxRiskSetPlan.prepare(survival);
            if (riskSets == null)
                throw new IllegalArgumentException(
                    "sparse pedigree Cox currently requires right-censored data");
            survivalPlan = SurvivalPlan.prepare(survival);
            pedigreeDesign = pedigreeEffect.randomEffect();
            precision = pedigreeEffect.precision();
            this.ordinaryEffects = List.copyOf(ordinaryEffects);
            this.options = options;
            if (pedigreeDesign.observations() != rows
                    || !pedigreeDesign.sparse())
                throw new IllegalArgumentException(
                    "pedigree incidence must be sparse and row-aligned");
            for (RandomEffectTerm term : this.ordinaryEffects)
                if (term == null || term.observations() != rows)
                    throw new IllegalArgumentException(
                        "ordinary frailty terms must be row-aligned");
            if (options.initialVariances().length
                    != 1 + this.ordinaryEffects.size())
                throw new IllegalArgumentException(
                    "one variance is required per sparse and ordinary frailty term");
            sparseColumns = pedigreeDesign.coefficients();
            observationCoefficient = incidence(pedigreeDesign, rows);
            coefficientObservation = inverseIncidence(
                observationCoefficient, sparseColumns);
            ordinaryStarts = new int[this.ordinaryEffects.size()];
            int ordinaryCount = 0;
            for (int term = 0; term < this.ordinaryEffects.size(); term++) {
                ordinaryStarts[term] = ordinaryCount;
                ordinaryCount += this.ordinaryEffects.get(term).coefficients();
            }
            ordinaryColumns = ordinaryCount;
            ordinaryDesign = ordinaryDesign(this.ordinaryEffects,
                ordinaryStarts, rows, ordinaryColumns);
            pattern = SparsePattern.prepare(precision);
            BackendContext selected = BackendContext.select(backendPolicy);
            PreparedSparseCholesky selectedFactor = null;
            double selectedLogDeterminant;
            try {
                backend = selected.backend();
                selectedLogDeterminant = backend.dcsrpotrf(
                    lowerPrecision(precision), MatrixTriangle.LOWER,
                    SparseOrdering.MINIMUM_DEGREE).logDeterminant();
                selectedFactor = backend.prepareDcsrpotrf(
                    pattern.matrix(new double[pattern.nonzeros()],
                        options.initialVariances()[0]), MatrixTriangle.LOWER,
                    SparseOrdering.MINIMUM_DEGREE);
            } catch (RuntimeException | Error failure) {
                if (selectedFactor != null) selectedFactor.close();
                selected.close();
                throw failure;
            }
            context = selected;
            factor = selectedFactor;
            precisionLogDeterminant = selectedLogDeterminant;
            sharedLogVariances = Arrays.stream(options.initialVariances())
                .map(Math::log).toArray();
        }

        public CoxMixedResult fit(double[][] fixedEffects, double[] offset) {
            if (closed)
                throw new IllegalStateException("prepared sparse Cox fit is closed");
            if (fixedEffects == null || fixedEffects.length != rows
                    || fixedEffects[0] == null || fixedEffects[0].length == 0)
                throw new IllegalArgumentException(
                    "fixed effects must be nonempty and row-aligned");
            int fixedColumns = fixedEffects[0].length;
            double[] fixed = MatrixOps.rowMajor(fixedEffects, rows);
            validateFixed(fixed, rows, fixedColumns);
            double[] offsets = offset == null ? new double[rows]
                : MatrixOps.finiteCopy(offset, "offset");
            if (offsets.length != rows)
                throw new IllegalArgumentException(
                    "one offset is required per observation");
            Fit fit = optimize(fixed, fixedColumns, offsets,
                sharedLogVariances.clone());
            sharedLogVariances = fit.logVariances().clone();
            return result(fit, fixedColumns);
        }

        /** Fits conditional modes with caller-fixed frailty variances. */
        public CoxMixedResult fitAtVariances(
                double[][] fixedEffects, double[] offset,
                double... variances) {
            if (closed)
                throw new IllegalStateException("prepared sparse Cox fit is closed");
            if (fixedEffects == null || fixedEffects.length != rows
                    || fixedEffects[0] == null || fixedEffects[0].length == 0
                    || variances == null
                    || variances.length != sharedLogVariances.length)
                throw new IllegalArgumentException(
                    "fixed effects and one variance per frailty are required");
            double[] logs = new double[variances.length];
            for (int index = 0; index < variances.length; index++) {
                if (!(variances[index] >= options.minimumVariance()
                        && variances[index] <= options.maximumVariance()))
                    throw new IllegalArgumentException(
                        "fixed frailty variances must lie within bounds");
                logs[index] = Math.log(variances[index]);
            }
            int fixedColumns = fixedEffects[0].length;
            double[] fixed = MatrixOps.rowMajor(fixedEffects, rows);
            validateFixed(fixed, rows, fixedColumns);
            double[] offsets = offset == null ? new double[rows]
                : MatrixOps.finiteCopy(offset, "offset");
            if (offsets.length != rows)
                throw new IllegalArgumentException(
                    "one offset is required per observation");
            Mode mode = mode(fixed, fixedColumns, offsets, logs, null);
            return result(new Fit(mode, logs, 0, mode.converged()), fixedColumns);
        }

        private Fit optimize(
                double[] fixed, int fixedColumns, double[] offset,
                double[] logVariances) {
            Mode current = mode(fixed, fixedColumns, offset, logVariances, null);
            double step = 1.0;
            boolean converged = false;
            int iterations = 0;
            for (int iteration = 1;
                    iteration <= options.maximumVarianceIterations(); iteration++) {
                iterations = iteration;
                boolean improved = false;
                boolean stationary = true;
                for (int term = 0; term < logVariances.length; term++) {
                    double original = logVariances[term];
                    double leftLog = clamp(original - step);
                    double rightLog = clamp(original + step);
                    Mode left = leftLog == original ? null
                        : varianceMode(fixed, fixedColumns, offset,
                            logVariances, term, leftLog, current);
                    Mode right = rightLog == original ? null
                        : varianceMode(fixed, fixedColumns, offset,
                            logVariances, term, rightLog, current);
                    Mode best = current;
                    double bestLog = original;
                    if (left != null && left.laplace() > best.laplace()) {
                        best = left;
                        bestLog = leftLog;
                    }
                    if (right != null && right.laplace() > best.laplace()) {
                        best = right;
                        bestLog = rightLog;
                    }
                    if (best == current && left != null && right != null) {
                        VariancePoint refined = refineVariance(fixed,
                            fixedColumns, offset, logVariances, term,
                            original, current, leftLog, left, rightLog, right);
                        best = refined.mode();
                        bestLog = refined.logVariance();
                    }
                    if (best != current) {
                        current = best;
                        logVariances[term] = bestLog;
                        improved = true;
                        stationary = false;
                    } else if (left != null && right != null
                            && (left.laplace() > current.laplace()
                                || right.laplace() > current.laplace())) {
                        stationary = false;
                    }
                }
                if (stationary) {
                    converged = true;
                    break;
                }
                if (!improved) step *= 0.5;
                if (step <= options.logVarianceTolerance()) {
                    converged = true;
                    break;
                }
            }
            current = mode(fixed, fixedColumns, offset, logVariances, current);
            return new Fit(current, logVariances, iterations,
                converged && current.converged());
        }

        private VariancePoint refineVariance(
                double[] fixed, int fixedColumns, double[] offset,
                double[] logs, int term, double centerLog, Mode center,
                double leftLog, Mode left, double rightLog, Mode right) {
            double curvature = left.laplace() - 2 * center.laplace()
                + right.laplace();
            if (!(curvature < 0)
                    || Math.abs(centerLog - leftLog
                        - (rightLog - centerLog)) > 1e-10)
                return new VariancePoint(centerLog, center);
            double radius = rightLog - centerLog;
            double displacement = 0.5 * radius
                * (left.laplace() - right.laplace()) / curvature;
            double candidateLog = Math.max(leftLog,
                Math.min(rightLog, centerLog + displacement));
            Mode candidate = varianceMode(fixed, fixedColumns, offset,
                logs, term, candidateLog, center);
            double selectedLog = candidate.laplace() > center.laplace()
                ? candidateLog : centerLog;
            Mode selected = candidate.laplace() > center.laplace()
                ? candidate : center;
            for (double radius2 : new double[] {0.2, 0.002}) {
                double low = clamp(selectedLog - radius2);
                double high = clamp(selectedLog + radius2);
                if (low == selectedLog || high == selectedLog) break;
                Mode lowMode = varianceMode(fixed, fixedColumns, offset,
                    logs, term, low, selected);
                Mode highMode = varianceMode(fixed, fixedColumns, offset,
                    logs, term, high, selected);
                double curve = lowMode.laplace() - 2 * selected.laplace()
                    + highMode.laplace();
                if (curve < 0) {
                    double shift = 0.5 * radius2
                        * (lowMode.laplace() - highMode.laplace()) / curve;
                    double vertex = Math.max(low,
                        Math.min(high, selectedLog + shift));
                    Mode vertexMode = varianceMode(fixed, fixedColumns,
                        offset, logs, term, vertex, selected);
                    if (vertexMode.laplace() > selected.laplace()) {
                        selectedLog = vertex;
                        selected = vertexMode;
                    }
                }
            }
            return new VariancePoint(selectedLog, selected);
        }

        private Mode varianceMode(
                double[] fixed, int fixedColumns, double[] offset,
                double[] logs, int term, double value, Mode start) {
            double[] candidate = logs.clone();
            candidate[term] = value;
            return mode(fixed, fixedColumns, offset, candidate, start);
        }

        private Mode mode(
                double[] fixed, int fixedColumns, double[] offset,
                double[] logVariances, Mode start) {
            int denseColumns = fixedColumns + ordinaryColumns;
            double[] dense = denseDesign(fixed, fixedColumns);
            double[] denseCoefficients = start == null
                ? new double[denseColumns] : start.dense().clone();
            double[] sparseCoefficients = start == null
                ? new double[sparseColumns] : start.sparse().clone();
            Evaluation evaluation = evaluate(dense, denseColumns,
                fixedColumns, denseCoefficients, sparseCoefficients,
                offset, logVariances);
            boolean converged = false;
            int iterations = 0;
            for (int iteration = 1;
                    iteration <= options.coxOptions().maximumIterations(); iteration++) {
                iterations = iteration;
                if (maximum(evaluation.denseScore(), evaluation.sparseScore())
                        <= options.coxOptions().scoreTolerance()) {
                    converged = true;
                    break;
                }
                SystemSolution solution = solve(evaluation, denseColumns, false);
                double scale = 1.0;
                Evaluation candidate = null;
                double[] nextDense = null;
                double[] nextSparse = null;
                for (int halving = 0;
                        halving <= options.coxOptions().maximumStepHalvings();
                        halving++) {
                    nextDense = add(denseCoefficients, solution.dense(), scale);
                    nextSparse = add(sparseCoefficients, solution.sparse(), scale);
                    candidate = evaluate(dense, denseColumns, fixedColumns,
                        nextDense, nextSparse, offset, logVariances);
                    if (candidate.penalized() >= evaluation.penalized() - 1e-12)
                        break;
                    scale *= 0.5;
                }
                if (candidate.penalized() < evaluation.penalized() - 1e-10)
                    break;
                double change = Math.max(relative(denseCoefficients, nextDense),
                    relative(sparseCoefficients, nextSparse));
                denseCoefficients = nextDense;
                sparseCoefficients = nextSparse;
                evaluation = candidate;
                if (change <= options.coxOptions().relativeTolerance()
                        && maximum(evaluation.denseScore(),
                            evaluation.sparseScore())
                            <= Math.sqrt(options.coxOptions().scoreTolerance())) {
                    converged = true;
                    break;
                }
            }
            evaluation = evaluate(dense, denseColumns, fixedColumns,
                denseCoefficients, sparseCoefficients, offset, logVariances);
            SystemSolution finalSystem = solve(evaluation, denseColumns, true);
            double laplace = laplace(evaluation, fixedColumns, denseColumns,
                logVariances);
            return new Mode(denseCoefficients, sparseCoefficients,
                finalSystem.denseCovariance(), evaluation, laplace,
                iterations, converged);
        }

        private Evaluation evaluate(
                double[] dense, int denseColumns, int fixedColumns,
                double[] denseCoefficients, double[] sparseCoefficients,
                double[] offset, double[] logVariances) {
            double[] eta = offset.clone();
            for (int row = 0; row < rows; row++) {
                int start = row * denseColumns;
                for (int column = 0; column < denseColumns; column++)
                    eta[row] += dense[start + column]
                        * denseCoefficients[column];
                eta[row] += sparseCoefficients[observationCoefficient[row]];
            }
            double maximumEta = Arrays.stream(eta).max().orElse(0.0);
            double[] risk = new double[rows];
            double risk0 = 0.0;
            double[] risk1 = new double[denseColumns];
            double[] risk2 = new double[denseColumns * denseColumns];
            for (int row = 0; row < rows; row++) {
                double weight = Math.exp(eta[row] - maximumEta);
                risk[row] = weight;
                risk0 += weight;
                addDenseMoments(dense, row, denseColumns, weight,
                    risk1, risk2);
            }
            int events = survivalPlan.eventRows().length;
            double[] prefixA = new double[events];
            double[] prefixB = new double[events];
            double[] prefixC = new double[events * denseColumns];
            double[] denseScore = new double[denseColumns];
            double[] denseInformation =
                new double[denseColumns * denseColumns];
            double partial = 0.0;
            int removal = 0;
            double cumulativeA = 0.0;
            double cumulativeB = 0.0;
            double[] cumulativeC = new double[denseColumns];
            double[] mean = new double[denseColumns];
            for (int event = 0; event < events; event++) {
                int eventRow = survivalPlan.eventRows()[event];
                double time = survival.stopView()[eventRow];
                while (removal < survivalPlan.removalRows().length
                        && survival.stopView()[survivalPlan.removalRows()[removal]]
                            < time) {
                    int row = survivalPlan.removalRows()[removal++];
                    double weight = risk[row];
                    risk0 -= weight;
                    subtractDenseMoments(dense, row, denseColumns, weight,
                        risk1, risk2);
                }
                partial += eta[eventRow] - maximumEta - Math.log(risk0);
                int eventStart = eventRow * denseColumns;
                for (int column = 0; column < denseColumns; column++) {
                    mean[column] = risk1[column] / risk0;
                    denseScore[column] += dense[eventStart + column]
                        - mean[column];
                }
                for (int left = 0; left < denseColumns; left++)
                    for (int right = 0; right <= left; right++)
                        denseInformation[left * denseColumns + right] +=
                            risk2[left * denseColumns + right] / risk0
                            - mean[left] * mean[right];
                cumulativeA += 1.0 / risk0;
                cumulativeB += 1.0 / (risk0 * risk0);
                prefixA[event] = cumulativeA;
                prefixB[event] = cumulativeB;
                for (int column = 0; column < denseColumns; column++) {
                    cumulativeC[column] += risk1[column] / (risk0 * risk0);
                    prefixC[event * denseColumns + column] = cumulativeC[column];
                }
            }
            symmetrizeLower(denseInformation, denseColumns);
            double[] sparseScore = new double[sparseColumns];
            double[] cross = new double[sparseColumns * denseColumns];
            double[] sparseInformation = new double[pattern.nonzeros()];
            for (int coefficient = 0; coefficient < sparseColumns; coefficient++) {
                int row = coefficientObservation[coefficient];
                if (row < 0) continue;
                int last = survivalPlan.lastEvent()[row];
                if (last < 0) continue;
                double weight = risk[row];
                double a = prefixA[last];
                double b = prefixB[last];
                sparseScore[coefficient] =
                    (survival.eventView()[row] ? 1.0 : 0.0) - weight * a;
                sparseInformation[pattern.diagonalPositions()[coefficient]] +=
                    weight * a - weight * weight * b;
                int denseStart = row * denseColumns;
                for (int column = 0; column < denseColumns; column++)
                    cross[coefficient * denseColumns + column] = weight
                        * (dense[denseStart + column] * a
                            - prefixC[last * denseColumns + column]);
            }
            double sparseVariance = Math.exp(logVariances[0]);
            double[] precisionProduct = multiply(precision, sparseCoefficients);
            double penalty = 0.0;
            for (int coefficient = 0;
                    coefficient < sparseColumns; coefficient++) {
                sparseScore[coefficient] -=
                    precisionProduct[coefficient] / sparseVariance;
                penalty += 0.5 * sparseCoefficients[coefficient]
                    * precisionProduct[coefficient] / sparseVariance;
            }
            for (int term = 0; term < ordinaryEffects.size(); term++) {
                double variance = Math.exp(logVariances[term + 1]);
                int start = fixedColumns + ordinaryStarts[term];
                for (int column = 0;
                        column < ordinaryEffects.get(term).coefficients();
                        column++) {
                    int index = start + column;
                    denseScore[index] -= denseCoefficients[index] / variance;
                    denseInformation[index * denseColumns + index] +=
                        1.0 / variance;
                    penalty += 0.5 * denseCoefficients[index]
                        * denseCoefficients[index] / variance;
                }
            }
            return new Evaluation(partial, partial - penalty, eta,
                denseScore, sparseScore, denseInformation, cross,
                sparseInformation, sparseVariance);
        }

        private SystemSolution solve(
                Evaluation evaluation, int denseColumns, boolean covariance) {
            factor.refactor(pattern.matrix(evaluation.sparseInformation(),
                evaluation.sparseVariance()));
            double[] solved = new double[sparseColumns * (denseColumns + 1)];
            for (int row = 0; row < sparseColumns; row++) {
                solved[row * (denseColumns + 1)] =
                    evaluation.sparseScore()[row];
                for (int column = 0; column < denseColumns; column++)
                    solved[row * (denseColumns + 1) + column + 1] =
                        evaluation.cross()[row * denseColumns + column];
            }
            factor.solveInPlace(solved, denseColumns + 1);
            double[] schur = evaluation.denseInformation().clone();
            double[] right = evaluation.denseScore().clone();
            for (int left = 0; left < denseColumns; left++)
                for (int sparse = 0; sparse < sparseColumns; sparse++) {
                    double value = evaluation.cross()[sparse * denseColumns + left];
                    right[left] -= value * solved[sparse * (denseColumns + 1)];
                    for (int column = 0; column < denseColumns; column++)
                        schur[left * denseColumns + column] -= value
                            * solved[sparse * (denseColumns + 1) + column + 1];
                }
            symmetrize(schur, denseColumns);
            CholeskyFactor denseFactor = CoxMath.factor(backend, schur,
                denseColumns, options.coxOptions().informationRidge());
            double[] denseStep = denseFactor.solve(right);
            double[] sparseStep = new double[sparseColumns];
            for (int row = 0; row < sparseColumns; row++) {
                double value = solved[row * (denseColumns + 1)];
                for (int column = 0; column < denseColumns; column++)
                    value -= solved[row * (denseColumns + 1) + column + 1]
                        * denseStep[column];
                sparseStep[row] = value;
            }
            double[] denseCovariance = covariance
                ? denseFactor.solve(MatrixOps.identity(denseColumns), denseColumns)
                : null;
            return new SystemSolution(denseStep, sparseStep, denseCovariance);
        }

        private double laplace(
                Evaluation evaluation, int fixedColumns, int denseColumns,
                double[] logVariances) {
            factor.refactor(pattern.matrix(evaluation.sparseInformation(),
                evaluation.sparseVariance()));
            double randomLogDeterminant = factor.logDeterminant();
            if (ordinaryColumns > 0) {
                double[] cross = new double[sparseColumns * ordinaryColumns];
                for (int row = 0; row < sparseColumns; row++)
                    for (int column = 0; column < ordinaryColumns; column++)
                        cross[row * ordinaryColumns + column] =
                            evaluation.cross()[row * denseColumns
                                + fixedColumns + column];
                double[] solved = cross.clone();
                factor.solveInPlace(solved, ordinaryColumns);
                double[] schur = new double[ordinaryColumns * ordinaryColumns];
                for (int left = 0; left < ordinaryColumns; left++)
                    for (int right = 0; right < ordinaryColumns; right++) {
                        double value = evaluation.denseInformation()[
                            (fixedColumns + left) * denseColumns
                                + fixedColumns + right];
                        for (int row = 0; row < sparseColumns; row++)
                            value -= cross[row * ordinaryColumns + left]
                                * solved[row * ordinaryColumns + right];
                        schur[left * ordinaryColumns + right] = value;
                    }
                randomLogDeterminant += CoxMath.factor(backend, schur,
                    ordinaryColumns,
                    options.coxOptions().informationRidge()).logDeterminant();
            }
            double precisionLog = precisionLogDeterminant
                - sparseColumns * logVariances[0];
            for (int term = 0; term < ordinaryEffects.size(); term++)
                precisionLog -= ordinaryEffects.get(term).coefficients()
                    * logVariances[term + 1];
            return evaluation.penalized() + 0.5 * precisionLog
                - 0.5 * randomLogDeterminant;
        }

        private CoxMixedResult result(Fit fit, int fixedColumns) {
            Mode mode = fit.mode();
            double[] beta = Arrays.copyOf(mode.dense(), fixedColumns);
            double[] covariance = new double[fixedColumns * fixedColumns];
            int denseColumns = mode.dense().length;
            for (int row = 0; row < fixedColumns; row++)
                for (int column = 0; column < fixedColumns; column++)
                    covariance[row * fixedColumns + column] =
                        mode.denseCovariance()[row * denseColumns + column];
            List<CoxRandomEffectEstimates> effects = new ArrayList<>();
            effects.add(new CoxRandomEffectEstimates(pedigreeDesign.name(),
                pedigreeDesign.coefficientNames(),
                Math.exp(fit.logVariances()[0]), mode.sparse()));
            for (int term = 0; term < ordinaryEffects.size(); term++) {
                RandomEffectTerm value = ordinaryEffects.get(term);
                int start = fixedColumns + ordinaryStarts[term];
                effects.add(new CoxRandomEffectEstimates(value.name(),
                    value.coefficientNames(),
                    Math.exp(fit.logVariances()[term + 1]),
                    Arrays.copyOfRange(mode.dense(), start,
                        start + value.coefficients())));
            }
            return new CoxMixedResult(beta, covariance, effects,
                CoxPartialLikelihood.baselineFromLinearPredictor(survival,
                    mode.evaluation().eta(), options.coxOptions().ties(), riskSets),
                mode.evaluation().partial(), mode.evaluation().penalized(),
                mode.laplace(), options, fit.varianceIterations(),
                fit.converged(), fit.converged()
                    ? "sparse Laplace variance and mode tolerances reached"
                    : "sparse Laplace optimization did not converge",
                context.provenance());
        }

        private double[] denseDesign(double[] fixed, int fixedColumns) {
            int columns = fixedColumns + ordinaryColumns;
            double[] result = new double[rows * columns];
            for (int row = 0; row < rows; row++) {
                System.arraycopy(fixed, row * fixedColumns,
                    result, row * columns, fixedColumns);
                System.arraycopy(ordinaryDesign, row * ordinaryColumns,
                    result, row * columns + fixedColumns, ordinaryColumns);
            }
            return result;
        }

        private double clamp(double value) {
            return Math.max(Math.log(options.minimumVariance()),
                Math.min(Math.log(options.maximumVariance()), value));
        }

        @Override public void close() {
            if (closed) return;
            closed = true;
            factor.close();
            context.close();
        }
    }

    private static int[] incidence(RandomEffectTerm term, int rows) {
        int[] starts = term.rowPointers();
        int[] columns = term.columnIndices();
        double[] values = term.sparseValues();
        int[] result = new int[rows];
        for (int row = 0; row < rows; row++) {
            if (starts[row + 1] - starts[row] != 1
                    || values[starts[row]] != 1.0)
                throw new IllegalArgumentException(
                    "sparse pedigree Cox requires one unit incidence per row");
            result[row] = columns[starts[row]];
        }
        return result;
    }

    private static int[] inverseIncidence(int[] incidence, int columns) {
        int[] result = new int[columns];
        Arrays.fill(result, -1);
        for (int row = 0; row < incidence.length; row++) {
            int column = incidence[row];
            if (result[column] >= 0)
                throw new IllegalArgumentException(
                    "sparse pedigree Cox requires at most one row per individual");
            result[column] = row;
        }
        return result;
    }

    private static double[] ordinaryDesign(
            List<RandomEffectTerm> terms, int[] starts,
            int rows, int columns) {
        double[] result = new double[rows * columns];
        for (int term = 0; term < terms.size(); term++) {
            RandomEffectTerm value = terms.get(term);
            double[] design = value.design();
            for (int row = 0; row < rows; row++)
                System.arraycopy(design, row * value.coefficients(), result,
                    row * columns + starts[term], value.coefficients());
        }
        return result;
    }

    private static CsrMatrix lowerPrecision(SparsePrecisionMatrix precision) {
        int[] starts = precision.rowStarts();
        int[] columns = precision.columnIndices();
        double[] values = precision.values();
        int count = 0;
        for (int row = 0; row < precision.dimension(); row++)
            for (int index = starts[row]; index < starts[row + 1]; index++)
                if (columns[index] <= row) count++;
        int[] resultStarts = new int[precision.dimension() + 1];
        int[] resultColumns = new int[count];
        double[] resultValues = new double[count];
        int position = 0;
        for (int row = 0; row < precision.dimension(); row++) {
            resultStarts[row] = position + 1;
            for (int index = starts[row]; index < starts[row + 1]; index++)
                if (columns[index] <= row) {
                    resultColumns[position] = columns[index] + 1;
                    resultValues[position++] = values[index];
                }
        }
        resultStarts[precision.dimension()] = position + 1;
        return new CsrMatrix(precision.dimension(), precision.dimension(),
            resultValues, resultColumns, resultStarts);
    }

    private static double[] multiply(
            SparsePrecisionMatrix matrix, double[] vector) {
        int[] starts = matrix.rowStarts();
        int[] columns = matrix.columnIndices();
        double[] values = matrix.values();
        double[] result = new double[matrix.dimension()];
        for (int row = 0; row < matrix.dimension(); row++)
            for (int index = starts[row]; index < starts[row + 1]; index++)
                result[row] += values[index] * vector[columns[index]];
        return result;
    }

    private static void addDenseMoments(
            double[] design, int row, int columns, double weight,
            double[] first, double[] second) {
        int start = row * columns;
        for (int left = 0; left < columns; left++) {
            double x = design[start + left];
            first[left] += weight * x;
            for (int right = 0; right <= left; right++)
                second[left * columns + right] +=
                    weight * x * design[start + right];
        }
    }

    private static void subtractDenseMoments(
            double[] design, int row, int columns, double weight,
            double[] first, double[] second) {
        addDenseMoments(design, row, columns, -weight, first, second);
    }

    private static double[] add(double[] values, double[] step, double scale) {
        double[] result = values.clone();
        for (int index = 0; index < result.length; index++)
            result[index] += scale * step[index];
        return result;
    }

    private static double relative(double[] oldValues, double[] newValues) {
        double result = 0.0;
        for (int index = 0; index < oldValues.length; index++)
            result = Math.max(result, Math.abs(newValues[index] - oldValues[index])
                / (1.0 + Math.abs(oldValues[index])));
        return result;
    }

    private static double maximum(double[] first, double[] second) {
        return Math.max(CoxMath.maximumAbsolute(first),
            CoxMath.maximumAbsolute(second));
    }

    private static void validateFixed(
            double[] fixed, int rows, int columns) {
        for (int column = 0; column < columns; column++) {
            double minimum = Double.POSITIVE_INFINITY;
            double maximum = Double.NEGATIVE_INFINITY;
            for (int row = 0; row < rows; row++) {
                minimum = Math.min(minimum, fixed[row * columns + column]);
                maximum = Math.max(maximum, fixed[row * columns + column]);
            }
            if (!(maximum > minimum))
                throw new IllegalArgumentException(
                    "Cox fixed effects must not include a constant column");
        }
    }

    private static void symmetrizeLower(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++)
            for (int column = 0; column < row; column++)
                matrix[column * dimension + row] =
                    matrix[row * dimension + column];
    }

    private static void symmetrize(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++)
            for (int column = 0; column < row; column++) {
                double value = 0.5 * (matrix[row * dimension + column]
                    + matrix[column * dimension + row]);
                matrix[row * dimension + column] = value;
                matrix[column * dimension + row] = value;
            }
    }

    private record SurvivalPlan(
            int[] eventRows, int[] removalRows, int[] lastEvent) {
        static SurvivalPlan prepare(CoxSurvivalData survival) {
            int firstStratum = survival.strataView()[0];
            List<Integer> events = new ArrayList<>();
            List<Integer> rows = new ArrayList<>();
            for (int row = 0; row < survival.observations(); row++) {
                if (survival.strataView()[row] != firstStratum)
                    throw new IllegalArgumentException(
                        "sparse pedigree Cox currently supports one stratum");
                rows.add(row);
                if (survival.eventView()[row]) events.add(row);
            }
            Comparator<Integer> ascending = Comparator.comparingDouble(
                row -> survival.stopView()[row]);
            events.sort(ascending);
            rows.sort(ascending);
            double previous = Double.NEGATIVE_INFINITY;
            for (int row : events) {
                double time = survival.stopView()[row];
                if (time == previous)
                    throw new IllegalArgumentException(
                        "sparse pedigree Cox currently requires distinct event times");
                previous = time;
            }
            int[] eventRows = events.stream().mapToInt(Integer::intValue).toArray();
            double[] times = events.stream().mapToDouble(
                row -> survival.stopView()[row]).toArray();
            int[] last = new int[survival.observations()];
            for (int row = 0; row < last.length; row++) {
                int position = Arrays.binarySearch(times, survival.stopView()[row]);
                last[row] = position >= 0 ? position : -position - 2;
            }
            return new SurvivalPlan(eventRows,
                rows.stream().mapToInt(Integer::intValue).toArray(), last);
        }
    }

    private record SparsePattern(
            int dimension, int[] rowStarts, int[] columns,
            int[] diagonalPositions, int[] precisionPositions,
            double[] precisionValues) {
        static SparsePattern prepare(SparsePrecisionMatrix precision) {
            int dimension = precision.dimension();
            int[] starts = precision.rowStarts();
            int[] columns = precision.columnIndices();
            double[] values = precision.values();
            TreeMap<Long, Boolean> keys = new TreeMap<>();
            for (int row = 0; row < dimension; row++) {
                keys.put((long) row * dimension + row, Boolean.TRUE);
                for (int index = starts[row]; index < starts[row + 1]; index++) {
                    int column = columns[index];
                    if (column <= row)
                        keys.put((long) row * dimension + column, Boolean.TRUE);
                }
            }
            int[] resultStarts = new int[dimension + 1];
            int[] resultColumns = new int[keys.size()];
            int position = 0;
            int currentRow = 0;
            for (long key : keys.keySet()) {
                int row = (int) (key / dimension);
                while (currentRow < row)
                    resultStarts[++currentRow] = position;
                resultColumns[position++] = (int) (key % dimension);
            }
            while (currentRow < dimension)
                resultStarts[++currentRow] = position;
            int[] diagonal = new int[dimension];
            for (int row = 0; row < dimension; row++)
                diagonal[row] = Arrays.binarySearch(resultColumns,
                    resultStarts[row], resultStarts[row + 1], row);
            List<Integer> precisionPositions = new ArrayList<>();
            List<Double> precisionValues = new ArrayList<>();
            for (int row = 0; row < dimension; row++)
                for (int index = starts[row]; index < starts[row + 1]; index++) {
                    int column = columns[index];
                    if (column > row) continue;
                    precisionPositions.add(Arrays.binarySearch(resultColumns,
                        resultStarts[row], resultStarts[row + 1], column));
                    precisionValues.add(values[index]);
                }
            return new SparsePattern(dimension, resultStarts, resultColumns,
                diagonal, precisionPositions.stream().mapToInt(
                    Integer::intValue).toArray(),
                precisionValues.stream().mapToDouble(Double::doubleValue).toArray());
        }

        int nonzeros() { return columns.length; }

        CsrMatrix matrix(double[] data, double variance) {
            double[] numeric = data.clone();
            for (int index = 0; index < precisionPositions.length; index++)
                numeric[precisionPositions[index]] +=
                    precisionValues[index] / variance;
            int[] csrStarts = rowStarts.clone();
            int[] csrColumns = columns.clone();
            for (int index = 0; index < csrStarts.length; index++)
                csrStarts[index]++;
            for (int index = 0; index < csrColumns.length; index++)
                csrColumns[index]++;
            return new CsrMatrix(dimension, dimension, numeric,
                csrColumns, csrStarts);
        }

    }
    private record Evaluation(
        double partial, double penalized, double[] eta,
        double[] denseScore, double[] sparseScore,
        double[] denseInformation, double[] cross,
        double[] sparseInformation, double sparseVariance) { }
    private record SystemSolution(
        double[] dense, double[] sparse, double[] denseCovariance) { }
    private record Mode(
        double[] dense, double[] sparse, double[] denseCovariance,
        Evaluation evaluation, double laplace, int iterations,
        boolean converged) { }
    private record Fit(
        Mode mode, double[] logVariances, int varianceIterations,
        boolean converged) { }
    private record VariancePoint(double logVariance, Mode mode) { }
}
