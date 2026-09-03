/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.LinearMixedModelResult;
import org.jlinalg.mixed.RandomEffectEstimates;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.reml.RemlOptions;

/**
 * Retains batch/pedigree sparse structure and symbolic factorization while a
 * P-spline predictor changes across a high-throughput Gaussian GAMM scan.
 */
public final class PreparedGammPredictorScan implements AutoCloseable {
    private static final String SMOOTH_TERM = "__scan_smooth__";

    private final double[] response;
    private final double[] parametricDesign;
    private final int observations;
    private final int parametricColumns;
    private final int basisDimension;
    private final List<RandomEffectTerm> staticTerms;
    private final BackendContext context;
    private final SparseLinearMixedModel.PreparedChangingTerm prepared;
    private volatile boolean closed;

    /**
     * Creates a general prepared scan. Precision bases correspond one-for-one
     * with the static random-effect terms.
     */
    public PreparedGammPredictorScan(
            double[] response,
            double[] parametricDesign,
            int observations,
            int parametricColumns,
            List<RandomEffectTerm> staticTerms,
            List<SparsePrecisionMatrix> staticPrecisions,
            int basisDimension,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(
            response, parametricDesign, observations, parametricColumns);
        if (staticTerms == null || staticTerms.isEmpty()
                || staticPrecisions == null
                || staticTerms.size() != staticPrecisions.size()
                || basisDimension <= 2 || options == null
                || backendPolicy == null) {
            throw new IllegalArgumentException(
                "static random effects, precision bases, and GAMM controls are required");
        }
        for (RandomEffectTerm term : staticTerms) {
            if (term == null || term.observations() != observations
                    || SMOOTH_TERM.equals(term.name())) {
                throw new IllegalArgumentException(
                    "static random-effect rows and names are invalid");
            }
        }
        this.response = response.clone();
        this.parametricDesign = parametricDesign.clone();
        this.observations = observations;
        this.parametricColumns = parametricColumns;
        this.basisDimension = basisDimension;
        this.staticTerms = List.copyOf(staticTerms);
        this.context = BackendContext.select(backendPolicy);
        try {
            int randomColumns = basisDimension - 2;
            List<String> names = new ArrayList<>(randomColumns);
            for (int column = 0; column < randomColumns; column++) {
                names.add("smooth.pen" + (column + 1));
            }
            this.prepared =
                SparseLinearMixedModel.prepareChangingTermWithBackend(
                    observations, staticTerms, staticPrecisions,
                    SMOOTH_TERM, randomColumns, names, options,
                    context.backend(), context.provenance());
        } catch (RuntimeException | Error failure) {
            context.close();
            throw failure;
        }
    }

    /** Creates a prepared batch random-intercept GAMM scan. */
    public static PreparedGammPredictorScan batch(
            double[] response,
            double[] parametricDesign,
            int observations,
            int parametricColumns,
            String batchName,
            List<String> batches,
            int basisDimension,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        RandomEffectTerm batch =
            RandomEffectTerm.randomIntercept(batchName, batches);
        return new PreparedGammPredictorScan(response, parametricDesign,
            observations, parametricColumns, List.of(batch),
            List.of(SparsePrecisionMatrix.identity(batch.coefficients())),
            basisDimension, options, backendPolicy);
    }

    /** Creates a prepared additive-pedigree plus batch GAMM scan. */
    public static PreparedGammPredictorScan pedigreeBatch(
            double[] response,
            double[] parametricDesign,
            int observations,
            int parametricColumns,
            String pedigreeName,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            String batchName,
            List<String> batches,
            int basisDimension,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        PedigreeRandomEffectTerm pedigreeTerm =
            PedigreeRandomEffectTerm.of(
                pedigreeName, observationIndividualIds, pedigree);
        RandomEffectTerm batch =
            RandomEffectTerm.randomIntercept(batchName, batches);
        return new PreparedGammPredictorScan(response, parametricDesign,
            observations, parametricColumns,
            List.of(pedigreeTerm.randomEffect(), batch),
            List.of(pedigreeTerm.precision(),
                SparsePrecisionMatrix.identity(batch.coefficients())),
            basisDimension, options, backendPolicy);
    }

    /** Fits one changing cubic, second-difference P-spline predictor. */
    public GammScanResult fit(String termName, double[] covariate) {
        requireOpen();
        Objects.requireNonNull(termName, "termName");
        PSplineMixedModelCompiler.Compiled compiled =
            PSplineMixedModelCompiler.compile(
                parametricDesign, observations, parametricColumns,
                List.of(PSplineTerm.of(
                    termName, covariate, basisDimension)),
                context.backend());
        PSplineMixedModelCompiler.Term smooth = compiled.terms().get(0);
        if (smooth.randomColumns() != prepared.changingCoefficients()) {
            throw new IllegalStateException(
                "smooth decomposition changed random coefficient count");
        }
        SparseLinearMixedModelResult sparse = prepared.fit(response,
            compiled.fixedDesign(), compiled.fixedColumns(),
            smooth.randomDesign());
        LinearMixedModelResult mixed = LinearMixedModelResult.fromSparse(
            sparse, response, compiled.fixedDesign(), observations,
            compiled.fixedColumns());
        double[] variances = sparse.varianceComponents();
        int smoothVarianceIndex = staticTerms.size();
        double smoothVariance = variances[smoothVarianceIndex];
        double residualVariance = variances[variances.length - 1];
        double[] beta = sparse.beta();
        double[] fixedCoefficients = java.util.Arrays.copyOfRange(
            beta, parametricColumns,
            parametricColumns + smooth.fixedColumns());
        RandomEffectEstimates random =
            sparse.randomEffects(SMOOTH_TERM);
        double edf = smooth.fixedColumns();
        for (double predictionErrorVariance
                : random.predictionErrorVariances()) {
            edf += 1.0 - predictionErrorVariance / smoothVariance;
        }
        edf = Math.max(smooth.fixedColumns(),
            Math.min(smooth.fixedColumns() + smooth.randomColumns(), edf));
        double[] contribution =
            PSplineMixedModelCompiler.contribution(
                smooth, fixedCoefficients, random.estimates());
        SmoothTermEstimate estimate = new SmoothTermEstimate(
            smooth.term(), smooth.fixedTransform(), smooth.fixedMeans(),
            fixedCoefficients, smooth.randomTransform(),
            smooth.randomMeans(), random.estimates(), contribution,
            residualVariance / smoothVariance, edf);
        Map<String, double[]> randomContributions =
            new LinkedHashMap<>();
        for (RandomEffectTerm term : staticTerms) {
            randomContributions.put(term.name(), multiply(
                term, sparse.randomEffects(term.name()).estimates()));
        }
        return new GammScanResult(
            mixed, estimate, randomContributions);
    }

    /** Uses a representative fit as the deterministic start for later genes. */
    public void warmStart(GammScanResult result) {
        requireOpen();
        prepared.warmStart(result.reml().varianceComponents());
    }

    public int observations() { return observations; }
    public int parametricColumns() { return parametricColumns; }
    public int basisDimension() { return basisDimension; }

    private static double[] multiply(
            RandomEffectTerm term, double[] coefficients) {
        double[] result = new double[term.observations()];
        if (term.sparse()) {
            int[] starts = term.rowPointers();
            int[] columns = term.columnIndices();
            double[] values = term.sparseValues();
            for (int row = 0; row < result.length; row++) {
                for (int index = starts[row];
                        index < starts[row + 1]; index++) {
                    result[row] += values[index]
                        * coefficients[columns[index]];
                }
            }
            return result;
        }
        double[] design = term.design();
        for (int row = 0; row < result.length; row++) {
            for (int column = 0;
                    column < coefficients.length; column++) {
                result[row] += design[
                    row * coefficients.length + column]
                    * coefficients[column];
            }
        }
        return result;
    }

    private void requireOpen() {
        if (closed) {
            throw new IllegalStateException(
                "prepared GAMM predictor scan is closed");
        }
    }

    @Override
    public void close() {
        if (closed) return;
        closed = true;
        prepared.close();
        context.close();
    }
}
