/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;
import jdistlib.Normal;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.model.MissingDataPolicy;

/** Cluster-streaming generalized estimating equations. */
public final class Gee {
    private static final double MINIMUM_VARIANCE = 1e-12;
    private static final double MAXIMUM_CORRELATION = 0.98;

    private Gee() { }

    /** Fits a GEE with unit weights, zero offset, and default controls. */
    public static GeeResult fit(
            double[] response,
            double[][] design,
            int[] cluster,
            GlmFamily family) {
        return fit(response, design, cluster, null, family,
            null, null, GeeOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a GEE from a conventional rectangular Java matrix. */
    public static GeeResult fit(
            double[] response,
            double[][] design,
            int[] cluster,
            int[] repeated,
            GlmFamily family,
            double[] weights,
            double[] offset,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajorUnchecked(design, response.length);
        return fit(response, rowMajor, response.length, design[0].length,
            cluster, repeated, family, weights, offset, options, backendPolicy);
    }

    /** Fits a GEE from a contiguous row-major design matrix. */
    public static GeeResult fit(
            double[] response,
            double[] design,
            int rows,
            int columns,
            int[] cluster,
            int[] repeated,
            GlmFamily family,
            double[] weights,
            double[] offset,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (family == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "family, options, and backendPolicy are required");
        }
        Prepared data = prepare(response, design, rows, columns,
            cluster, repeated, weights, offset, options, family);
        validateOptions(data, options, family);

        double[] coefficients = startingCoefficients(
            data, family, options, backendPolicy);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return fitPrepared(data, family, options, coefficients,
                context.backend(), context.provenance());
        }
    }

    private static GeeResult fitPrepared(
            Prepared data,
            GlmFamily family,
            GeeOptions options,
            double[] coefficients,
            ComputeBackend backend,
            BackendProvenance provenance) {
        double[] associationParameters = initialAssociation(data, options);
        Scale scale = new Scale(options.dispersion(),
            new double[0], constant(data.rows(), options.dispersion()));
        State state = state(data, family, coefficients, scale.observationScale(), backend);
        boolean converged = false;
        String message = "maximum iterations reached";
        int iterations = 0;

        for (int iteration = 1;
                iteration <= options.maximumIterations(); iteration++) {
            iterations = iteration;
            scale = estimateScale(data, family, state, options, backend);
            state = state(data, family, coefficients,
                scale.observationScale(), backend);
            double[] updatedAssociation = estimateAssociation(
                data, state, scale.observationScale(), options, backend);
            if (iteration > 1 && updatedAssociation.length == associationParameters.length) {
                for (int index = 0; index < updatedAssociation.length; index++) {
                    updatedAssociation[index] = 0.5 * associationParameters[index]
                        + 0.5 * updatedAssociation[index];
                }
            }
            double associationChange = relativeMaximumChange(
                associationParameters, updatedAssociation);
            associationParameters = updatedAssociation;

            Accumulation accumulation = accumulate(data, state,
                scale.observationScale(), associationParameters,
                options, backend);
            double[] score = accumulation.score().clone();
            double[] inverseBread = inverseSymmetric(
                accumulation.bread(), data.columns(), backend);
            if (options.method() == GeeMethod.BIAS_REDUCED) {
                score = adjustedScore(data, state, scale.observationScale(),
                    associationParameters, options, inverseBread, backend);
            } else if (options.method() == GeeMethod.JEFFREYS) {
                addInPlace(score, jeffreysAdjustment(data, family, coefficients,
                    scale.observationScale(), associationParameters,
                    options, backend));
            }
            double[] step = multiply(inverseBread, data.columns(), score);
            double[] candidate = null;
            State candidateState = null;
            double baselineScoreNorm = norm(score);
            double multiplier = 1.0;
            for (int attempt = 0; attempt < 24; attempt++) {
                double[] trial = addScaled(coefficients, step, multiplier);
                try {
                    State trialState = state(data, family, trial,
                        scale.observationScale(), backend);
                    Accumulation trialAccumulation = accumulate(data, trialState,
                        scale.observationScale(), associationParameters,
                        options, backend);
                    if (norm(trialAccumulation.score()) <= baselineScoreNorm
                            + 1e-10 * (1.0 + baselineScoreNorm)) {
                        candidate = trial;
                        candidateState = trialState;
                        break;
                    }
                } catch (IllegalArgumentException exception) {
                    // A smaller scoring step may restore the mean domain.
                }
                multiplier *= 0.5;
            }
            if (candidate == null) {
                if (norm(step) <= options.relativeTolerance()
                        * (1.0 + norm(coefficients))) {
                    converged = true;
                    message = "estimating-equation step tolerance reached";
                } else {
                    message = "step-halving could not reduce the estimating equation";
                }
                break;
            }
            double coefficientChange = relativeMaximumChange(coefficients, candidate);
            coefficients = candidate;
            state = candidateState;
            if (coefficientChange <= options.relativeTolerance()
                    && associationChange <= Math.sqrt(options.relativeTolerance())) {
                converged = true;
                message = "coefficient and association tolerances reached";
                break;
            }
        }

        scale = estimateScale(data, family, state, options, backend);
        state = state(data, family, coefficients, scale.observationScale(), backend);
        associationParameters = estimateAssociation(
            data, state, scale.observationScale(), options, backend);

        if (options.method() == GeeMethod.BIAS_CORRECTED
                || options.method() == GeeMethod.ONE_STEP_JEFFREYS
                || options.method() == GeeMethod.HYBRID_JEFFREYS) {
            Accumulation preliminary = accumulate(data, state,
                scale.observationScale(), associationParameters, options, backend);
            double[] inverseBread = inverseSymmetric(
                preliminary.bread(), data.columns(), backend);
            double[] adjusted;
            if (options.method() == GeeMethod.BIAS_CORRECTED) {
                adjusted = adjustedScore(data, state, scale.observationScale(),
                    associationParameters, options, inverseBread, backend);
            } else {
                adjusted = preliminary.score().clone();
                addInPlace(adjusted, jeffreysAdjustment(data, family, coefficients,
                    scale.observationScale(), associationParameters,
                    options, backend));
            }
            double multiplier = options.method() == GeeMethod.HYBRID_JEFFREYS
                ? 0.5 : 1.0;
            coefficients = addScaled(coefficients,
                multiply(inverseBread, data.columns(), adjusted), multiplier);
            state = state(data, family, coefficients,
                scale.observationScale(), backend);
            message += "; one adjusted step applied";
        }

        scale = estimateScale(data, family, state, options, backend);
        state = state(data, family, coefficients, scale.observationScale(), backend);
        associationParameters = estimateAssociation(
            data, state, scale.observationScale(), options, backend);
        Accumulation finalAccumulation = accumulate(data, state,
            scale.observationScale(), associationParameters, options, backend);
        int p = data.columns();
        double[] naive = inverseSymmetric(finalAccumulation.bread(), p, backend);
        double[] robust = sandwich(naive, finalAccumulation.meat(), p);
        double correction = data.clusters() > p
            ? (double) data.clusters() / (data.clusters() - p) : Double.NaN;
        double[] dfAdjusted = robust.clone();
        if (Double.isFinite(correction)) {
            scaleInPlace(dfAdjusted, correction);
        } else {
            Arrays.fill(dfAdjusted, Double.NaN);
        }
        double[] correctedMeat = biasCorrectedMeat(data, state,
            scale.observationScale(), associationParameters,
            options, naive, backend);
        double[] biasCorrected = sandwich(naive, correctedMeat, p);
        double[] selected = switch (options.covariance()) {
            case NAIVE -> naive;
            case ROBUST -> robust;
            case DF_ADJUSTED -> dfAdjusted;
            case BIAS_CORRECTED -> biasCorrected;
        };
        Inference inference = inference(coefficients, selected,
            options.confidenceLevel());
        GeeCriteria criteria = criteria(data, state, scale.observationScale(),
            robust, family, p, backend);

        return new GeeResult(
            family.name(), options.correlation(), options.association(),
            options.covariance(), options.method(), coefficients, selected,
            naive, robust, dfAdjusted, biasCorrected,
            inference.standardErrors(), inference.statistics(), inference.pValues(),
            inference.lower(), inference.upper(),
            data.output(state.linearPredictor()), data.output(state.means()),
            data.output(state.pearsonResiduals()), associationParameters,
            scale.average(), scale.coefficients(), criteria,
            data.rows(), data.clusters(), data.minimumClusterSize(),
            data.maximumClusterSize(), p, iterations, converged, message,
            data.retainedRows(), data.originalRows(), provenance);
    }

    private static Prepared prepare(
            double[] response,
            double[] design,
            int rows,
            int columns,
            int[] cluster,
            int[] repeated,
            double[] suppliedWeights,
            double[] suppliedOffset,
            GeeOptions options,
            GlmFamily family) {
        if (rows < 1 || columns < 1 || response == null || response.length != rows
                || design == null || design.length != rows * columns
                || cluster == null || cluster.length != rows
                || repeated != null && repeated.length != rows) {
            throw new IllegalArgumentException("GEE input dimensions are invalid");
        }
        if (suppliedWeights != null && suppliedWeights.length != rows
                || suppliedOffset != null && suppliedOffset.length != rows) {
            throw new IllegalArgumentException("weight or offset length is invalid");
        }
        double[] optionScale = options.scaleDesign();
        if (optionScale != null && options.scaleDesignRows() != rows) {
            throw new IllegalArgumentException(
                "scaleDesign rows must equal response length");
        }
        double[] weights = suppliedWeights == null ? constant(rows, 1.0)
            : suppliedWeights.clone();
        double[] offset = suppliedOffset == null ? new double[rows]
            : suppliedOffset.clone();
        List<Integer> retained = new ArrayList<>();
        for (int row = 0; row < rows; row++) {
            boolean finite = Double.isFinite(response[row])
                && Double.isFinite(weights[row]) && Double.isFinite(offset[row]);
            for (int column = 0; finite && column < columns; column++) {
                finite = Double.isFinite(design[row * columns + column]);
            }
            for (int column = 0; finite && optionScale != null
                    && column < options.scaleDesignColumns(); column++) {
                finite = Double.isFinite(
                    optionScale[row * options.scaleDesignColumns() + column]);
            }
            if (!finite && options.missingDataPolicy() == MissingDataPolicy.ERROR) {
                throw new IllegalArgumentException(
                    "model inputs must be finite; non-finite row " + row);
            }
            if (finite) {
                if (!(weights[row] > 0.0)) {
                    throw new IllegalArgumentException(
                        "weights must be strictly positive");
                }
                family.validateResponse(response[row], weights[row]);
                retained.add(row);
            }
        }
        if (retained.size() <= columns) {
            throw new IllegalArgumentException(
                "GEE requires more complete observations than parameters");
        }
        Integer[] order = retained.toArray(Integer[]::new);
        Comparator<Integer> comparator = Comparator.comparingInt(
            (Integer row) -> cluster[row]);
        if (repeated != null) {
            comparator = comparator.thenComparingInt(row -> repeated[row]);
        } else {
            comparator = comparator.thenComparingInt(row -> row);
        }
        Arrays.sort(order, comparator);

        int n = order.length;
        int[] globalWaves = repeated == null ? null : recodeWaves(repeated, retained);
        double[] keptResponse = new double[n];
        double[] keptDesign = new double[n * columns];
        double[] keptWeights = new double[n];
        double[] keptOffset = new double[n];
        double[] keptScale = optionScale == null ? null
            : new double[n * options.scaleDesignColumns()];
        int[] keptClusters = new int[n];
        int[] waves = new int[n];
        int[] sortedOriginal = new int[n];
        int currentCluster = 0;
        int withinCluster = 0;
        for (int destination = 0; destination < n; destination++) {
            int source = order[destination];
            if (destination == 0 || cluster[source] != currentCluster) {
                currentCluster = cluster[source];
                withinCluster = 0;
            } else {
                withinCluster++;
            }
            keptResponse[destination] = response[source];
            System.arraycopy(design, source * columns,
                keptDesign, destination * columns, columns);
            keptWeights[destination] = weights[source];
            keptOffset[destination] = offset[source];
            if (keptScale != null) {
                System.arraycopy(optionScale,
                    source * options.scaleDesignColumns(), keptScale,
                    destination * options.scaleDesignColumns(),
                    options.scaleDesignColumns());
            }
            keptClusters[destination] = cluster[source];
            waves[destination] = repeated == null
                ? withinCluster : globalWaves[source];
            sortedOriginal[destination] = source;
            if (destination > 0 && keptClusters[destination] == keptClusters[destination - 1]
                    && waves[destination] == waves[destination - 1]) {
                throw new IllegalArgumentException(
                    "repeated values must be unique within each cluster");
            }
        }
        int[] starts = clusterStarts(keptClusters);
        int[] retainedRows = retained.stream().mapToInt(Integer::intValue).toArray();
        int[] outputPosition = new int[n];
        int[] rank = new int[rows];
        Arrays.fill(rank, -1);
        for (int index = 0; index < retainedRows.length; index++) {
            rank[retainedRows[index]] = index;
        }
        for (int index = 0; index < n; index++) {
            outputPosition[index] = rank[sortedOriginal[index]];
        }
        int maximumWave = Arrays.stream(waves).max().orElse(0) + 1;
        return new Prepared(keptResponse, keptDesign, keptWeights, keptOffset,
            keptScale, options.scaleDesignColumns(), n, columns,
            keptClusters, waves, starts, maximumWave,
            retainedRows, outputPosition, rows);
    }

    private static int[] recodeWaves(int[] repeated, List<Integer> retained) {
        TreeSet<Integer> levels = new TreeSet<>();
        for (int row : retained) levels.add(repeated[row]);
        int[] result = new int[repeated.length];
        int level = 0;
        for (int value : levels) {
            for (int row : retained) {
                if (repeated[row] == value) result[row] = level;
            }
            level++;
        }
        return result;
    }

    private static int[] clusterStarts(int[] cluster) {
        int count = 1;
        for (int row = 1; row < cluster.length; row++) {
            if (cluster[row] != cluster[row - 1]) count++;
        }
        int[] starts = new int[count + 1];
        int next = 1;
        for (int row = 1; row < cluster.length; row++) {
            if (cluster[row] != cluster[row - 1]) starts[next++] = row;
        }
        starts[count] = cluster.length;
        return starts;
    }

    private static void validateOptions(
            Prepared data, GeeOptions options, GlmFamily family) {
        if (data.clusters() < 2) {
            throw new IllegalArgumentException("GEE requires at least two clusters");
        }
        if (options.initialCoefficients() != null
                && options.initialCoefficients().length != data.columns()) {
            throw new IllegalArgumentException(
                "initial coefficient count must equal design columns");
        }
        if (options.correlation() == GeeCorrelation.FIXED
                && options.fixedAssociationDimension() < data.maximumWave()) {
            throw new IllegalArgumentException(
                "fixed association matrix is smaller than the wave count");
        }
        int pairs = data.maximumWave() * (data.maximumWave() - 1) / 2;
        if (options.correlation() == GeeCorrelation.USER_DEFINED
                && options.correlationDesignRows() != pairs) {
            throw new IllegalArgumentException(
                "correlationDesign must have choose(maximum waves, 2) rows");
        }
        if (options.association() == GeeAssociation.ODDS_RATIO) {
            if (!family.name().toLowerCase(java.util.Locale.ROOT).contains("binomial")) {
                throw new IllegalArgumentException(
                    "odds-ratio association requires a binomial family");
            }
            for (double value : data.response()) {
                if (value != 0.0 && value != 1.0) {
                    throw new IllegalArgumentException(
                        "odds-ratio association requires binary 0/1 responses");
                }
            }
            if (options.correlation() != GeeCorrelation.INDEPENDENCE
                    && options.correlation() != GeeCorrelation.EXCHANGEABLE
                    && options.correlation() != GeeCorrelation.UNSTRUCTURED
                    && options.correlation() != GeeCorrelation.FIXED
                    && options.correlation() != GeeCorrelation.USER_DEFINED) {
                throw new IllegalArgumentException(
                    "odds-ratio association supports independence, exchangeable, "
                    + "unstructured, fixed, and user-defined structures");
            }
        }
    }

    private static double[] startingCoefficients(
            Prepared data,
            GlmFamily family,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        double[] supplied = options.initialCoefficients();
        if (supplied != null) return supplied;
        GlmOptions glmOptions = GlmOptions.builder()
            .maximumIterations(options.maximumIterations())
            .relativeTolerance(options.relativeTolerance())
            .confidenceLevel(options.confidenceLevel())
            .missingDataPolicy(MissingDataPolicy.ERROR)
            .build();
        GlmResult result = Glm.fit(data.response(), data.design(),
            data.rows(), data.columns(), family, data.weights(), data.offset(),
            glmOptions, backendPolicy);
        return result.coefficients();
    }

    private static State state(
            Prepared data,
            GlmFamily family,
            double[] coefficients,
            double[] observationScale,
            ComputeBackend backend) {
        double[] eta = MatrixOps.multiply(backend, data.design(),
            data.rows(), data.columns(), coefficients);
        double[] means = new double[data.rows()];
        double[] derivatives = new double[data.rows()];
        double[] variances = new double[data.rows()];
        double[] residuals = new double[data.rows()];
        double[] pearson = new double[data.rows()];
        double deviance = 0.0;
        for (int row = 0; row < data.rows(); row++) {
            eta[row] += data.offset()[row];
            means[row] = family.inverseLink(eta[row]);
            derivatives[row] = family.meanDerivative(eta[row]);
            variances[row] = Math.max(MINIMUM_VARIANCE,
                family.variance(means[row]));
            if (!Double.isFinite(means[row]) || !Double.isFinite(derivatives[row])
                    || derivatives[row] == 0.0
                    || !Double.isFinite(variances[row])) {
                throw new IllegalArgumentException(
                    "family produced invalid GEE mean, derivative, or variance");
            }
            residuals[row] = data.response()[row] - means[row];
            pearson[row] = residuals[row] * Math.sqrt(data.weights()[row]
                / (variances[row] * Math.max(MINIMUM_VARIANCE,
                    observationScale[row])));
            deviance += data.weights()[row]
                * family.unitDeviance(data.response()[row], means[row]);
        }
        return new State(eta, means, derivatives, variances,
            residuals, pearson, deviance);
    }

    private static Scale estimateScale(
            Prepared data,
            GlmFamily family,
            State state,
            GeeOptions options,
            ComputeBackend backend) {
        if (options.fixedDispersion() || family.fixedDispersion()) {
            double value = options.fixedDispersion() ? options.dispersion() : 1.0;
            return new Scale(value, new double[0], constant(data.rows(), value));
        }
        double sum = 0.0;
        for (int row = 0; row < data.rows(); row++) {
            sum += data.weights()[row] * state.residuals()[row]
                * state.residuals()[row] / state.variances()[row];
        }
        double global = Math.max(MINIMUM_VARIANCE,
            sum / Math.max(1, data.rows() - data.columns()));
        if (data.scaleDesign() == null) {
            return new Scale(global, new double[0], constant(data.rows(), global));
        }
        int q = data.scaleColumns();
        double[] gamma = new double[q];
        for (int iteration = 0; iteration < 30; iteration++) {
            double[] eta = MatrixOps.multiply(backend, data.scaleDesign(),
                data.rows(), q, gamma);
            double[] target = new double[data.rows()];
            double[] weightedDesign = new double[data.scaleDesign().length];
            double[] weightedTarget = new double[data.rows()];
            double rootWeight = Math.sqrt(0.5);
            for (int row = 0; row < data.rows(); row++) {
                double phi = Math.exp(Math.max(-30.0, Math.min(30.0, eta[row])));
                double squared = data.weights()[row] * state.residuals()[row]
                    * state.residuals()[row] / state.variances()[row];
                target[row] = eta[row] + (squared - phi) / phi;
                weightedTarget[row] = rootWeight * target[row];
                for (int column = 0; column < q; column++) {
                    weightedDesign[row * q + column] = rootWeight
                        * data.scaleDesign()[row * q + column];
                }
            }
            double[] next = LeastSquaresSolver.solve(weightedDesign,
                weightedTarget, data.rows(), q, false, backend).coefficients();
            if (relativeMaximumChange(gamma, next) < 1e-8) {
                gamma = next;
                break;
            }
            gamma = next;
        }
        double[] eta = MatrixOps.multiply(backend, data.scaleDesign(),
            data.rows(), q, gamma);
        double[] phi = new double[data.rows()];
        double average = 0.0;
        for (int row = 0; row < data.rows(); row++) {
            phi[row] = Math.exp(Math.max(-30.0, Math.min(30.0, eta[row])));
            average += phi[row];
        }
        return new Scale(average / data.rows(), gamma, phi);
    }

    private static double[] initialAssociation(Prepared data, GeeOptions options) {
        int count = associationParameterCount(data, options);
        double[] result = new double[count];
        if (options.association() == GeeAssociation.ODDS_RATIO) {
            Arrays.fill(result, 1.0);
        }
        if (options.correlation() == GeeCorrelation.FIXED) {
            return lowerTriangle(options.fixedAssociation(),
                options.fixedAssociationDimension(), data.maximumWave());
        }
        return result;
    }

    private static int associationParameterCount(Prepared data, GeeOptions options) {
        return switch (options.correlation()) {
            case INDEPENDENCE -> 0;
            case EXCHANGEABLE, AR1 -> 1;
            case M_DEPENDENT -> options.dependenceOrder();
            case TOEPLITZ -> Math.max(0, data.maximumWave() - 1);
            case UNSTRUCTURED, FIXED -> data.maximumWave()
                * (data.maximumWave() - 1) / 2;
            case USER_DEFINED -> options.correlationDesignColumns();
        };
    }

    private static double[] estimateAssociation(
            Prepared data,
            State state,
            double[] phi,
            GeeOptions options,
            ComputeBackend backend) {
        if (options.correlation() == GeeCorrelation.INDEPENDENCE) {
            return new double[0];
        }
        if (options.correlation() == GeeCorrelation.FIXED) {
            return lowerTriangle(options.fixedAssociation(),
                options.fixedAssociationDimension(), data.maximumWave());
        }
        if (options.association() == GeeAssociation.ODDS_RATIO) {
            return estimateOddsRatios(data, options, backend);
        }
        int count = associationParameterCount(data, options);
        if (options.correlation() == GeeCorrelation.USER_DEFINED) {
            return estimateUserCorrelation(data, state, options, backend);
        }
        double[] sums = new double[count];
        double[] counts = new double[count];
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            int start = data.starts()[cluster];
            int end = data.starts()[cluster + 1];
            for (int first = start; first < end; first++) {
                for (int second = start; second < first; second++) {
                    if (options.correlation() == GeeCorrelation.AR1
                            && Math.abs(data.waves()[first]
                                - data.waves()[second]) != 1) {
                        continue;
                    }
                    int index = associationIndex(data.waves()[first],
                        data.waves()[second], data.maximumWave(), options);
                    if (index >= 0 && index < count) {
                        double firstPearson = state.residuals()[first]
                            * Math.sqrt(data.weights()[first]
                                / (state.variances()[first] * phi[first]));
                        double secondPearson = state.residuals()[second]
                            * Math.sqrt(data.weights()[second]
                                / (state.variances()[second] * phi[second]));
                        sums[index] += firstPearson * secondPearson;
                        counts[index]++;
                    }
                }
            }
        }
        double lowerExchange = -1.0
            / Math.max(1.0, data.maximumClusterSize() - 1.0) + 1e-6;
        for (int index = 0; index < count; index++) {
            double denominator = counts[index] - data.columns();
            double value = denominator <= 0.0 ? 0.0 : sums[index] / denominator;
            double lower = options.correlation() == GeeCorrelation.EXCHANGEABLE
                ? lowerExchange : -MAXIMUM_CORRELATION;
            sums[index] = clamp(value, lower, MAXIMUM_CORRELATION);
        }
        return sums;
    }

    private static double[] estimateUserCorrelation(
            Prepared data,
            State state,
            GeeOptions options,
            ComputeBackend backend) {
        int q = options.correlationDesignColumns();
        int pairCount = 0;
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            int size = data.starts()[cluster + 1] - data.starts()[cluster];
            pairCount += size * (size - 1) / 2;
        }
        if (pairCount < q) return new double[q];
        double[] design = new double[pairCount * q];
        double[] target = new double[pairCount];
        double[] z = options.correlationDesign();
        int row = 0;
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            int start = data.starts()[cluster];
            int end = data.starts()[cluster + 1];
            for (int first = start; first < end; first++) {
                for (int second = start; second < first; second++) {
                    int pair = pairIndex(data.waves()[first], data.waves()[second]);
                    System.arraycopy(z, pair * q, design, row * q, q);
                    target[row] = state.pearsonResiduals()[first]
                        * state.pearsonResiduals()[second];
                    row++;
                }
            }
        }
        try {
            return LeastSquaresSolver.solve(
                design, target, pairCount, q, false, backend).coefficients();
        } catch (IllegalArgumentException exception) {
            return new double[q];
        }
    }

    private static double[] estimateOddsRatios(
            Prepared data, GeeOptions options, ComputeBackend backend) {
        if (options.correlation() == GeeCorrelation.USER_DEFINED) {
            return estimateUserOddsRatios(data, options, backend);
        }
        int count = associationParameterCount(data, options);
        double[][] cells = new double[count][4];
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            int start = data.starts()[cluster];
            int end = data.starts()[cluster + 1];
            for (int first = start; first < end; first++) {
                for (int second = start; second < first; second++) {
                    int index = options.correlation() == GeeCorrelation.EXCHANGEABLE
                        ? 0 : pairIndex(data.waves()[first], data.waves()[second]);
                    int firstValue = (int) data.response()[first];
                    int secondValue = (int) data.response()[second];
                    int cell = firstValue == 1
                        ? (secondValue == 1 ? 0 : 1)
                        : (secondValue == 1 ? 2 : 3);
                    cells[index][cell]++;
                }
            }
        }
        double[] result = new double[count];
        double adding = options.oddsRatioContinuityCorrection();
        for (int index = 0; index < count; index++) {
            result[index] = (cells[index][0] + adding) * (cells[index][3] + adding)
                / ((cells[index][1] + adding) * (cells[index][2] + adding));
            result[index] = clamp(result[index], 1e-4, 1e4);
        }
        return result;
    }

    private static double[] estimateUserOddsRatios(
            Prepared data, GeeOptions options, ComputeBackend backend) {
        int pairs = data.maximumWave() * (data.maximumWave() - 1) / 2;
        double[][] cells = new double[pairs][4];
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            int start = data.starts()[cluster];
            int end = data.starts()[cluster + 1];
            for (int first = start; first < end; first++) {
                for (int second = start; second < first; second++) {
                    int pair = pairIndex(data.waves()[first], data.waves()[second]);
                    int firstValue = (int) data.response()[first];
                    int secondValue = (int) data.response()[second];
                    int cell = firstValue == 1
                        ? (secondValue == 1 ? 0 : 1)
                        : (secondValue == 1 ? 2 : 3);
                    cells[pair][cell]++;
                }
            }
        }
        double adding = options.oddsRatioContinuityCorrection();
        double[] target = new double[pairs];
        for (int pair = 0; pair < pairs; pair++) {
            target[pair] = Math.log((cells[pair][0] + adding)
                * (cells[pair][3] + adding)
                / ((cells[pair][1] + adding) * (cells[pair][2] + adding)));
        }
        try {
            return LeastSquaresSolver.solve(options.correlationDesign(), target,
                pairs, options.correlationDesignColumns(), false, backend).coefficients();
        } catch (IllegalArgumentException exception) {
            return new double[options.correlationDesignColumns()];
        }
    }

    private static int associationIndex(
            int firstWave, int secondWave, int maximumWave, GeeOptions options) {
        int difference = Math.abs(firstWave - secondWave);
        return switch (options.correlation()) {
            case EXCHANGEABLE, AR1 -> 0;
            case M_DEPENDENT -> difference >= 1
                && difference <= options.dependenceOrder() ? difference - 1 : -1;
            case TOEPLITZ -> difference - 1;
            case UNSTRUCTURED, FIXED -> pairIndex(firstWave, secondWave);
            case INDEPENDENCE, USER_DEFINED -> -1;
        };
    }

    private static Accumulation accumulate(
            Prepared data,
            State state,
            double[] phi,
            double[] associationParameters,
            GeeOptions options,
            ComputeBackend backend) {
        int p = data.columns();
        double[] bread = new double[p * p];
        double[] meat = new double[p * p];
        double[] score = new double[p];
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            ClusterComponents components = clusterComponents(data, state, phi,
                associationParameters, options, cluster, backend);
            double[] clusterScore = new double[p];
            for (int column = 0; column < p; column++) {
                for (int row = 0; row < components.size(); row++) {
                    clusterScore[column] += components.derivativeDesign()[row * p + column]
                        * components.solvedResidual()[row];
                }
                score[column] += clusterScore[column];
            }
            for (int first = 0; first < p; first++) {
                for (int second = 0; second < p; second++) {
                    double value = 0.0;
                    for (int row = 0; row < components.size(); row++) {
                        value += components.derivativeDesign()[row * p + first]
                            * components.solvedDerivativeDesign()[row * p + second];
                    }
                    bread[first * p + second] += value;
                    meat[first * p + second] += clusterScore[first]
                        * clusterScore[second];
                }
            }
        }
        symmetrize(bread, p);
        symmetrize(meat, p);
        return new Accumulation(bread, meat, score);
    }

    private static ClusterComponents clusterComponents(
            Prepared data,
            State state,
            double[] phi,
            double[] associationParameters,
            GeeOptions options,
            int cluster,
            ComputeBackend backend) {
        int start = data.starts()[cluster];
        int end = data.starts()[cluster + 1];
        int size = end - start;
        int p = data.columns();
        double[] derivativeDesign = new double[size * p];
        double[] residual = new double[size];
        for (int local = 0; local < size; local++) {
            int row = start + local;
            residual[local] = state.residuals()[row];
            for (int column = 0; column < p; column++) {
                derivativeDesign[local * p + column] = state.derivatives()[row]
                    * data.design()[row * p + column];
            }
        }
        double[] covariance = covarianceMatrix(data, state, phi,
            associationParameters, options, start, end);
        CholeskyFactor factor = covarianceFactor(covariance, size, backend);
        double[] solvedResidual = factor.solve(residual);
        double[] solvedDesign = new double[size * p];
        double[] rightSide = new double[size];
        for (int column = 0; column < p; column++) {
            for (int row = 0; row < size; row++) {
                rightSide[row] = derivativeDesign[row * p + column];
            }
            double[] solved = factor.solve(rightSide);
            for (int row = 0; row < size; row++) {
                solvedDesign[row * p + column] = solved[row];
            }
        }
        return new ClusterComponents(size, derivativeDesign, residual,
            solvedResidual, solvedDesign, covariance);
    }

    private static double[] covarianceMatrix(
            Prepared data,
            State state,
            double[] phi,
            double[] associationParameters,
            GeeOptions options,
            int start,
            int end) {
        int size = end - start;
        double[] covariance = new double[size * size];
        double[] standardDeviation = new double[size];
        for (int local = 0; local < size; local++) {
            int row = start + local;
            standardDeviation[local] = Math.sqrt(phi[row]
                * state.variances()[row] / data.weights()[row]);
            covariance[local * size + local] = standardDeviation[local]
                * standardDeviation[local];
        }
        for (int first = 1; first < size; first++) {
            for (int second = 0; second < first; second++) {
                int firstRow = start + first;
                int secondRow = start + second;
                double correlation = pairCorrelation(data, state,
                    associationParameters, options, firstRow, secondRow);
                double value = correlation * standardDeviation[first]
                    * standardDeviation[second];
                covariance[first * size + second] = value;
                covariance[second * size + first] = value;
            }
        }
        return covariance;
    }

    private static double pairCorrelation(
            Prepared data,
            State state,
            double[] parameters,
            GeeOptions options,
            int firstRow,
            int secondRow) {
        if (options.correlation() == GeeCorrelation.INDEPENDENCE) return 0.0;
        int firstWave = data.waves()[firstRow];
        int secondWave = data.waves()[secondRow];
        double association;
        if (options.correlation() == GeeCorrelation.FIXED) {
            int dimension = options.fixedAssociationDimension();
            association = options.fixedAssociation()[firstWave * dimension + secondWave];
        } else if (options.correlation() == GeeCorrelation.USER_DEFINED) {
            int pair = pairIndex(firstWave, secondWave);
            association = 0.0;
            double[] design = options.correlationDesign();
            for (int parameter = 0; parameter < parameters.length; parameter++) {
                association += design[pair * parameters.length + parameter]
                    * parameters[parameter];
            }
            if (options.association() == GeeAssociation.ODDS_RATIO) {
                association = Math.exp(clamp(association, -20.0, 20.0));
            }
        } else {
            int index = associationIndex(firstWave, secondWave,
                data.maximumWave(), options);
            if (index < 0) return 0.0;
            association = options.correlation() == GeeCorrelation.AR1
                ? Math.pow(parameters[0], Math.abs(firstWave - secondWave))
                : parameters[index];
        }
        if (options.association() == GeeAssociation.CORRELATION) {
            return clamp(association, -MAXIMUM_CORRELATION, MAXIMUM_CORRELATION);
        }
        double joint = binaryJointProbability(state.means()[firstRow],
            state.means()[secondRow], association);
        double denominator = Math.sqrt(state.variances()[firstRow]
            * state.variances()[secondRow]);
        return clamp((joint - state.means()[firstRow] * state.means()[secondRow])
            / denominator, -MAXIMUM_CORRELATION, MAXIMUM_CORRELATION);
    }

    private static double binaryJointProbability(
            double first, double second, double oddsRatio) {
        if (Math.abs(oddsRatio - 1.0) < 1e-8) return first * second;
        double a = oddsRatio - 1.0;
        double b = 1.0 + a * (first + second);
        double discriminant = Math.max(0.0,
            b * b - 4.0 * oddsRatio * a * first * second);
        double joint = (b - Math.sqrt(discriminant)) / (2.0 * a);
        double lower = Math.max(0.0, first + second - 1.0);
        double upper = Math.min(first, second);
        return clamp(joint, lower, upper);
    }

    private static CholeskyFactor covarianceFactor(
            double[] covariance, int size, ComputeBackend backend) {
        double[] candidate = covariance.clone();
        double maximumDiagonal = 0.0;
        for (int index = 0; index < size; index++) {
            maximumDiagonal = Math.max(maximumDiagonal,
                Math.abs(covariance[index * size + index]));
        }
        double jitter = Math.max(1e-12, maximumDiagonal * 1e-12);
        IllegalArgumentException failure = null;
        for (int attempt = 0; attempt < 10; attempt++) {
            try {
                return backend.dpotrf(candidate, size);
            } catch (IllegalArgumentException exception) {
                failure = exception;
                candidate = covariance.clone();
                for (int index = 0; index < size; index++) {
                    candidate[index * size + index] += jitter;
                }
                jitter *= 10.0;
            }
        }
        throw new IllegalArgumentException(
            "working covariance is not positive definite", failure);
    }

    private static double[] adjustedScore(
            Prepared data,
            State state,
            double[] phi,
            double[] associationParameters,
            GeeOptions options,
            double[] inverseBread,
            ComputeBackend backend) {
        int p = data.columns();
        double[] score = new double[p];
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            double[] clusterScore = biasAdjustedClusterScore(data, state, phi,
                associationParameters, options, inverseBread, cluster, backend);
            addInPlace(score, clusterScore);
        }
        return score;
    }

    private static double[] biasCorrectedMeat(
            Prepared data,
            State state,
            double[] phi,
            double[] associationParameters,
            GeeOptions options,
            double[] inverseBread,
            ComputeBackend backend) {
        int p = data.columns();
        double[] meat = new double[p * p];
        for (int cluster = 0; cluster < data.clusters(); cluster++) {
            double[] score = biasAdjustedClusterScore(data, state, phi,
                associationParameters, options, inverseBread, cluster, backend);
            for (int first = 0; first < p; first++) {
                for (int second = 0; second < p; second++) {
                    meat[first * p + second] += score[first] * score[second];
                }
            }
        }
        return meat;
    }

    private static double[] biasAdjustedClusterScore(
            Prepared data,
            State state,
            double[] phi,
            double[] associationParameters,
            GeeOptions options,
            double[] inverseBread,
            int cluster,
            ComputeBackend backend) {
        ClusterComponents components = clusterComponents(data, state, phi,
            associationParameters, options, cluster, backend);
        int size = components.size();
        int p = data.columns();
        double[] dA = multiplyRectangular(
            components.derivativeDesign(), size, p, inverseBread, p);
        double[] leverage = multiplyTransposeRight(
            dA, size, p, components.solvedDerivativeDesign(), size);
        double[] identityMinus = new double[size * size];
        for (int row = 0; row < size; row++) {
            for (int column = 0; column < size; column++) {
                identityMinus[row * size + column] = (row == column ? 1.0 : 0.0)
                    - leverage[row * size + column];
            }
        }
        double[] adjustedResidual;
        try {
            adjustedResidual = backend.dgetrf(identityMinus, size)
                .solve(components.residual());
        } catch (IllegalArgumentException exception) {
            adjustedResidual = components.residual();
        }
        CholeskyFactor factor = covarianceFactor(
            components.covariance(), size, backend);
        double[] solved = factor.solve(adjustedResidual);
        double[] score = new double[p];
        for (int column = 0; column < p; column++) {
            for (int row = 0; row < size; row++) {
                score[column] += components.derivativeDesign()[row * p + column]
                    * solved[row];
            }
        }
        return score;
    }

    private static double[] jeffreysAdjustment(
            Prepared data,
            GlmFamily family,
            double[] coefficients,
            double[] phi,
            double[] associationParameters,
            GeeOptions options,
            ComputeBackend backend) {
        int p = data.columns();
        double[] result = new double[p];
        for (int column = 0; column < p; column++) {
            double step = Math.cbrt(Math.ulp(1.0))
                * (1.0 + Math.abs(coefficients[column]));
            double[] lower = coefficients.clone();
            double[] upper = coefficients.clone();
            lower[column] -= step;
            upper[column] += step;
            try {
                State lowerState = state(data, family, lower, phi, backend);
                State upperState = state(data, family, upper, phi, backend);
                double lowerLogDet = logDeterminant(accumulate(data, lowerState,
                    phi, associationParameters, options, backend).bread(), p, backend);
                double upperLogDet = logDeterminant(accumulate(data, upperState,
                    phi, associationParameters, options, backend).bread(), p, backend);
                result[column] = options.jeffreysPower()
                    * (upperLogDet - lowerLogDet) / (2.0 * step);
            } catch (IllegalArgumentException exception) {
                result[column] = 0.0;
            }
        }
        return result;
    }

    private static GeeCriteria criteria(
            Prepared data,
            State state,
            double[] phi,
            double[] robustCovariance,
            GlmFamily family,
            int parameters,
            ComputeBackend backend) {
        GeeOptions independence = GeeOptions.builder()
            .correlation(GeeCorrelation.INDEPENDENCE)
            .fixedDispersion(1.0)
            .build();
        Accumulation independent = accumulate(data, state, phi,
            new double[0], independence, backend);
        double trace = traceProduct(independent.bread(), robustCovariance, parameters);
        double quasiLikelihood = -0.5 * state.deviance();
        double qic = -2.0 * quasiLikelihood + 2.0 * trace;
        double qicu = -2.0 * quasiLikelihood + 2.0 * parameters;
        double qicc = data.clusters() > parameters + 1
            ? qic + 2.0 * parameters * (parameters + 1.0)
                / (data.clusters() - parameters - 1.0)
            : Double.NaN;
        return new GeeCriteria(quasiLikelihood, qic, qicu,
            trace, qicc, parameters);
    }

    private static Inference inference(
            double[] coefficients, double[] covariance, double confidenceLevel) {
        int p = coefficients.length;
        double[] standardErrors = new double[p];
        double[] statistics = new double[p];
        double[] pValues = new double[p];
        double[] lower = new double[p];
        double[] upper = new double[p];
        double critical = Normal.quantile(0.5 + confidenceLevel / 2.0,
            0.0, 1.0, true, false);
        for (int column = 0; column < p; column++) {
            double variance = covariance[column * p + column];
            standardErrors[column] = variance >= 0.0
                ? Math.sqrt(variance) : Double.NaN;
            statistics[column] = coefficients[column] / standardErrors[column];
            pValues[column] = Double.isFinite(statistics[column])
                ? 2.0 * Normal.cumulative(Math.abs(statistics[column]),
                    0.0, 1.0, false, false)
                : Double.NaN;
            lower[column] = coefficients[column] - critical * standardErrors[column];
            upper[column] = coefficients[column] + critical * standardErrors[column];
        }
        return new Inference(standardErrors, statistics, pValues, lower, upper);
    }

    private static double[] inverseSymmetric(
            double[] matrix, int dimension, ComputeBackend backend) {
        double[] identity = new double[dimension * dimension];
        for (int index = 0; index < dimension; index++) {
            identity[index * dimension + index] = 1.0;
        }
        try {
            CholeskyFactor factor = backend.dpotrf(matrix, dimension);
            double[] inverse = new double[dimension * dimension];
            for (int column = 0; column < dimension; column++) {
                double[] right = new double[dimension];
                right[column] = 1.0;
                double[] solved = factor.solve(right);
                for (int row = 0; row < dimension; row++) {
                    inverse[row * dimension + column] = solved[row];
                }
            }
            symmetrize(inverse, dimension);
            return inverse;
        } catch (IllegalArgumentException exception) {
            try {
                return backend.dsytrf(matrix, dimension).solve(identity, dimension);
            } catch (IllegalArgumentException fallback) {
                throw new IllegalArgumentException(
                    "GEE sensitivity matrix is singular", fallback);
            }
        }
    }

    private static double logDeterminant(
            double[] matrix, int dimension, ComputeBackend backend) {
        return backend.dpotrf(matrix, dimension).logDeterminant();
    }

    private static double[] sandwich(
            double[] inverseBread, double[] meat, int dimension) {
        double[] temporary = multiplyRectangular(
            inverseBread, dimension, dimension, meat, dimension);
        double[] result = multiplyTransposeRight(
            temporary, dimension, dimension, inverseBread, dimension);
        symmetrize(result, dimension);
        return result;
    }

    private static double[] multiply(
            double[] matrix, int dimension, double[] vector) {
        double[] result = new double[dimension];
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                result[row] += matrix[row * dimension + column] * vector[column];
            }
        }
        return result;
    }

    private static double[] multiplyRectangular(
            double[] left, int rows, int shared,
            double[] right, int columns) {
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < columns; column++) {
                double value = 0.0;
                for (int inner = 0; inner < shared; inner++) {
                    value += left[row * shared + inner]
                        * right[inner * columns + column];
                }
                result[row * columns + column] = value;
            }
        }
        return result;
    }

    private static double[] multiplyTransposeRight(
            double[] left, int rows, int shared,
            double[] right, int rightRows) {
        double[] result = new double[rows * rightRows];
        for (int row = 0; row < rows; row++) {
            for (int other = 0; other < rightRows; other++) {
                double value = 0.0;
                for (int inner = 0; inner < shared; inner++) {
                    value += left[row * shared + inner]
                        * right[other * shared + inner];
                }
                result[row * rightRows + other] = value;
            }
        }
        return result;
    }

    private static double traceProduct(
            double[] left, double[] right, int dimension) {
        double result = 0.0;
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                result += left[row * dimension + column]
                    * right[column * dimension + row];
            }
        }
        return result;
    }

    private static double[] lowerTriangle(
            double[] matrix, int dimension, int usedDimension) {
        double[] result = new double[usedDimension * (usedDimension - 1) / 2];
        for (int first = 1; first < usedDimension; first++) {
            for (int second = 0; second < first; second++) {
                result[pairIndex(first, second)] =
                    matrix[first * dimension + second];
            }
        }
        return result;
    }

    private static int pairIndex(int first, int second) {
        int high = Math.max(first, second);
        int low = Math.min(first, second);
        return high * (high - 1) / 2 + low;
    }

    private static double[] addScaled(
            double[] values, double[] step, double multiplier) {
        double[] result = values.clone();
        for (int index = 0; index < result.length; index++) {
            result[index] += multiplier * step[index];
        }
        return result;
    }

    private static void addInPlace(double[] destination, double[] source) {
        for (int index = 0; index < destination.length; index++) {
            destination[index] += source[index];
        }
    }

    private static void scaleInPlace(double[] values, double multiplier) {
        for (int index = 0; index < values.length; index++) {
            values[index] *= multiplier;
        }
    }

    private static void symmetrize(double[] matrix, int dimension) {
        for (int first = 1; first < dimension; first++) {
            for (int second = 0; second < first; second++) {
                double average = 0.5 * (matrix[first * dimension + second]
                    + matrix[second * dimension + first]);
                matrix[first * dimension + second] = average;
                matrix[second * dimension + first] = average;
            }
        }
    }

    private static double relativeMaximumChange(
            double[] current, double[] candidate) {
        if (current.length != candidate.length) return Double.POSITIVE_INFINITY;
        double maximum = 0.0;
        for (int index = 0; index < current.length; index++) {
            maximum = Math.max(maximum,
                Math.abs(candidate[index] - current[index])
                    / (1.0 + Math.abs(current[index])));
        }
        return maximum;
    }

    private static double norm(double[] values) {
        double sum = 0.0;
        for (double value : values) sum += value * value;
        return Math.sqrt(sum);
    }

    private static double clamp(double value, double lower, double upper) {
        return Math.max(lower, Math.min(upper, value));
    }

    private static double[] constant(int length, double value) {
        double[] result = new double[length];
        Arrays.fill(result, value);
        return result;
    }

    private record Prepared(
            double[] response,
            double[] design,
            double[] weights,
            double[] offset,
            double[] scaleDesign,
            int scaleColumns,
            int rows,
            int columns,
            int[] cluster,
            int[] waves,
            int[] starts,
            int maximumWave,
            int[] retainedRows,
            int[] outputPosition,
            int originalRows) {
        int clusters() { return starts.length - 1; }
        int minimumClusterSize() {
            int result = Integer.MAX_VALUE;
            for (int index = 0; index < clusters(); index++) {
                result = Math.min(result, starts[index + 1] - starts[index]);
            }
            return result;
        }
        int maximumClusterSize() {
            int result = 0;
            for (int index = 0; index < clusters(); index++) {
                result = Math.max(result, starts[index + 1] - starts[index]);
            }
            return result;
        }
        double[] output(double[] sorted) {
            double[] result = new double[rows];
            for (int index = 0; index < rows; index++) {
                result[outputPosition[index]] = sorted[index];
            }
            return result;
        }
    }

    private record State(
            double[] linearPredictor,
            double[] means,
            double[] derivatives,
            double[] variances,
            double[] residuals,
            double[] pearsonResiduals,
            double deviance) {
    }

    private record Scale(
            double average,
            double[] coefficients,
            double[] observationScale) {
    }

    private record Accumulation(
            double[] bread,
            double[] meat,
            double[] score) {
    }

    private record ClusterComponents(
            int size,
            double[] derivativeDesign,
            double[] residual,
            double[] solvedResidual,
            double[] solvedDerivativeDesign,
            double[] covariance) {
    }

    private record Inference(
            double[] standardErrors,
            double[] statistics,
            double[] pValues,
            double[] lower,
            double[] upper) {
    }
}
