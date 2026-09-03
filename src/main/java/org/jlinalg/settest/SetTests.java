/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import jdistlib.T;
import jdistlib.accelerator.SymmetricEigenDecomposition;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.pipeline.VariantFilterResult;
import org.jlinalg.pipeline.VariantFilters;
import org.jlinalg.pipeline.VariantRecord;
import org.jlinalg.pipeline.VariantStatistics;

/** FP64 continuous-trait Burden, SKAT, and SKAT-O tests. */
public final class SetTests {
    private SetTests() { }

    public static SetTestResult burden(
            VariantSet set, LinearSetTestNullModel nullModel,
            SetTestOptions options) {
        PreparedVariantSet prepared = prepare(set, nullModel, options);
        double[] burden = new double[nullModel.observations()];
        for (int variant = 0; variant < prepared.dosagesView().length; variant++)
            for (int sample = 0; sample < burden.length; sample++)
                burden[sample] += prepared.weightsView()[variant]
                    * prepared.dosagesView()[variant][sample];
        double[] residualBurden = nullModel.residualize(
            new double[][] {burden})[0];
        double information = dot(residualBurden, residualBurden);
        if (!(information > 1e-14))
            throw new IllegalArgumentException(
                "burden is constant or collinear after covariate adjustment");
        double numerator = dot(residualBurden,
            nullModel.responseResidualView());
        double beta = numerator / information;
        int degrees = nullModel.residualDegreesOfFreedom() - 1;
        if (degrees < 1)
            throw new IllegalArgumentException(
                "burden test requires positive residual degrees of freedom");
        double residualSumSquares = Math.max(0,
            nullModel.residualSumSquares()
            - numerator * numerator / information);
        double standardError = Math.sqrt(
            residualSumSquares / degrees / information);
        double statistic = beta / standardError;
        double logP = (Math.log(2)
            + T.cumulative(-Math.abs(statistic), degrees, true, true))
            / Math.log(10);
        double pValue = Math.min(1, Math.exp(logP * Math.log(10)));
        return new SetTestResult(set.id(), "burden", set.variants().size(),
            prepared.includedVariants(), statistic, beta, standardError,
            degrees, pValue, logP, "student-t", new double[0],
            prepared.excludedVariants());
    }

    /** Related-sample burden score test using one retained REML projection. */
    public static SetTestResult burden(
            VariantSet set, RemlSetTestNullModel nullModel,
            SetTestOptions options) {
        PreparedVariantSet prepared = prepare(set, nullModel, options);
        double[] burden = new double[nullModel.observations()];
        for (int variant = 0; variant < prepared.dosagesView().length; variant++)
            for (int sample = 0; sample < burden.length; sample++)
                burden[sample] += prepared.weightsView()[variant]
                    * prepared.dosagesView()[variant][sample];
        SetTestScoreState state = nullModel.score(new double[][] {burden});
        double score = state.scoresView()[0];
        double information = state.informationView()[0];
        if (!(information > 1e-14) || !Double.isFinite(information))
            throw new IllegalArgumentException(
                "burden is constant or collinear after mixed-model adjustment");
        double beta = score / information;
        double standardError = Math.sqrt(1 / information);
        double statistic = score / Math.sqrt(information);
        double degrees = nullModel.degreesOfFreedom();
        double logP = (Math.log(2)
            + T.cumulative(-Math.abs(statistic), degrees, true, true))
            / Math.log(10);
        double pValue = Math.min(1, Math.exp(logP * Math.log(10)));
        return new SetTestResult(set.id(), "burden-reml-score",
            set.variants().size(), prepared.includedVariants(), statistic,
            beta, standardError, degrees, pValue, logP, "student-t",
            new double[0], prepared.excludedVariants());
    }

    /** Efficient burden score from a prepared set and any score null model. */
    public static SetTestResult burden(
            PreparedVariantSet prepared, SetTestScoreNullModel nullModel) {
        requireCompatible(prepared, nullModel);
        double[][] weighted = weighted(
            prepared.dosagesView(), prepared.weightsView());
        return burdenScore(prepared, nullModel, nullModel.score(weighted));
    }

    public static SetTestResult skat(
            VariantSet set, SetTestScoreNullModel nullModel,
            SetTestOptions options) {
        return skat(prepare(set, nullModel, options), nullModel);
    }

    public static SetTestResult skat(
            PreparedVariantSet prepared, SetTestScoreNullModel nullModel) {
        requireCompatible(prepared, nullModel);
        double[][] weighted = weighted(
            prepared.dosagesView(), prepared.weightsView());
        return kernelResult(prepared.id(), "skat",
            prepared.requestedVariants(), prepared.includedVariants(),
            nullModel.score(weighted), prepared.excludedVariants());
    }

    public static SkatOResult skatO(
            VariantSet set, SetTestScoreNullModel nullModel,
            SetTestOptions options) {
        return skatO(prepare(set, nullModel, options), nullModel, options);
    }

    public static SkatOResult skatO(
            PreparedVariantSet prepared, SetTestScoreNullModel nullModel,
            SetTestOptions options) {
        requireCompatible(prepared, nullModel);
        double[][] base = weighted(
            prepared.dosagesView(), prepared.weightsView());
        return skatO(prepared, nullModel.score(base), options);
    }

    /** Computes all three tests from one accelerated score projection. */
    public static SetTestSuiteResult scoreSuite(
            PreparedVariantSet prepared, SetTestScoreNullModel nullModel,
            SetTestOptions options) {
        requireCompatible(prepared, nullModel);
        double[][] base = weighted(
            prepared.dosagesView(), prepared.weightsView());
        SetTestScoreState state = nullModel.score(base);
        SetTestResult burden = burdenScore(prepared, nullModel, state);
        SetTestResult skat = kernelResult(prepared.id(), "skat",
            prepared.requestedVariants(), prepared.includedVariants(), state,
            prepared.excludedVariants());
        return new SetTestSuiteResult(
            burden, skat, skatO(prepared, state, options));
    }

    private static SkatOResult skatO(
            PreparedVariantSet prepared, SetTestScoreState baseState,
            SetTestOptions options) {
        double[] rhoGrid = options.skatORhoGrid();
        List<SkatOResult.Component> components = new ArrayList<>(rhoGrid.length);
        double minimumP = 1;
        for (double rho : rhoGrid) {
            SetTestResult result = kernelResult(prepared.id(),
                "skat-o[rho=" + rho + "]", prepared.requestedVariants(),
                prepared.includedVariants(), transform(baseState, rho),
                prepared.excludedVariants());
            components.add(new SkatOResult.Component(rho, result));
            minimumP = Math.min(minimumP, result.pValue());
        }
        double[] critical = new double[components.size()];
        for (int index = 0; index < critical.length; index++)
            critical[index] = QuadraticFormDistribution.critical(
                components.get(index).result().eigenvalues(), minimumP);
        Random random = new Random(options.randomSeed());
        ScoreSampler sampler = new ScoreSampler(
            baseState.informationView(), baseState.variants());
        int extreme = 0;
        double[] simulated = sampler.sample(
            random, options.skatOSimulations());
        for (int simulation = 0;
                simulation < options.skatOSimulations(); simulation++) {
            double squaredSum = 0;
            double sum = 0;
            for (int variant = 0; variant < baseState.variants(); variant++) {
                double score = simulated[
                    variant * options.skatOSimulations() + simulation];
                squaredSum += score * score;
                sum += score;
            }
            boolean exceeds = false;
            for (int component = 0; component < rhoGrid.length
                    && !exceeds; component++) {
                double rho = rhoGrid[component];
                double statistic = (1 - rho) * squaredSum
                    + rho * sum * sum;
                exceeds = statistic >= critical[component];
            }
            if (exceeds) extreme++;
        }
        double adjusted = (extreme + 1.0)
            / (options.skatOSimulations() + 1.0);
        return new SkatOResult(prepared.id(), prepared.requestedVariants(),
            prepared.includedVariants(), components, minimumP, adjusted,
            Math.log10(adjusted), options.skatOSimulations(),
            options.randomSeed(), prepared.excludedVariants());
    }

    private static SetTestResult kernelResult(
            String setId, String method, int requested,
            int included, SetTestScoreState scoreState,
            List<VariantFilterResult> excluded) {
        double statistic = 0;
        for (double score : scoreState.scoresView())
            statistic += score * score;
        double[] eigenvalues = eigenvalues(
            scoreState.informationView(), scoreState.variants());
        QuadraticFormDistribution.Tail tail =
            QuadraticFormDistribution.survival(statistic, eigenvalues);
        return new SetTestResult(setId, method, requested, included,
            statistic, Double.NaN, Double.NaN, Double.NaN,
            tail.pValue(), Math.log10(tail.pValue()), tail.method(),
            eigenvalues, excluded);
    }

    private static PreparedVariantSet prepare(
            VariantSet set, SetTestScoreNullModel nullModel,
            SetTestOptions options) {
        if (nullModel == null)
            throw new IllegalArgumentException("null model is required");
        return prepare(set, nullModel.observations(), options);
    }

    /** Performs filtering, effect-allele orientation, and imputation once. */
    public static PreparedVariantSet prepare(
            VariantSet set, int observations, SetTestOptions options) {
        if (set == null || options == null || observations < 1)
            throw new IllegalArgumentException(
                "variant set, observations, and options are required");
        List<double[]> dosages = new ArrayList<>();
        List<Double> weights = new ArrayList<>();
        List<VariantFilterResult> excluded = new ArrayList<>();
        for (WeightedVariant member : set.variants()) {
            VariantFilterResult filter = VariantFilters.evaluate(
                member.variant(), options.variantFilter());
            if (!filter.included()) {
                excluded.add(filter);
                continue;
            }
            dosages.add(oriented(member, filter.statistics(),
                options.missingPolicy(), observations));
            weights.add(member.weight());
        }
        if (dosages.isEmpty())
            throw new IllegalArgumentException(
                "variant set has no members after filtering: " + set.id());
        double[] numericWeights = new double[weights.size()];
        for (int index = 0; index < numericWeights.length; index++)
            numericWeights[index] = weights.get(index);
        return new PreparedVariantSet(set.id(), set.variants().size(),
            dosages.toArray(double[][]::new), numericWeights, excluded);
    }

    private static double[] oriented(
            WeightedVariant member, VariantStatistics statistics,
            SetTestMissingPolicy missingPolicy, int observations) {
        VariantRecord variant = member.variant();
        double[] source = variant.dosages();
        if (source.length != observations)
            throw new IllegalArgumentException(
                "variant dosage count does not match null-model observations");
        double altMean = statistics.dosageMean();
        double effectMean = member.effectAllele() == EffectAllele.ALTERNATE
            ? altMean : 2 - altMean;
        double[] result = new double[source.length];
        for (int sample = 0; sample < source.length; sample++) {
            if (!Double.isFinite(source[sample])) {
                result[sample] = switch (missingPolicy) {
                    case ERROR -> throw new IllegalArgumentException(
                        "variant contains a missing dosage: " + variant.id());
                    case MEAN_IMPUTE -> effectMean;
                    case ZERO -> 0;
                };
            } else {
                result[sample] = member.effectAllele() == EffectAllele.ALTERNATE
                    ? source[sample] : 2 - source[sample];
            }
        }
        return result;
    }

    private static double[][] weighted(
            double[][] rows, double[] weights) {
        double[][] result = new double[rows.length][];
        for (int row = 0; row < rows.length; row++) {
            result[row] = rows[row].clone();
            for (int column = 0; column < result[row].length; column++)
                result[row][column] *= weights[row];
        }
        return result;
    }

    private static SetTestResult burdenScore(
            PreparedVariantSet prepared, SetTestScoreNullModel nullModel,
            SetTestScoreState state) {
        double score = Arrays.stream(state.scoresView()).sum();
        double information = Arrays.stream(state.informationView()).sum();
        if (!(information > 1e-14) || !Double.isFinite(information))
            throw new IllegalArgumentException(
                "burden is constant or collinear after null-model adjustment");
        double beta = score / information;
        double standardError = Math.sqrt(1.0 / information);
        double statistic = score / Math.sqrt(information);
        double pValue = nullModel.burdenPValue(statistic);
        return new SetTestResult(prepared.id(), "burden-score",
            prepared.requestedVariants(), prepared.includedVariants(),
            statistic, beta, standardError, nullModel.degreesOfFreedom(),
            pValue, Math.log10(pValue), nullModel.burdenPValueMethod(),
            new double[0], prepared.excludedVariants());
    }

    private static SetTestScoreState transform(
            SetTestScoreState base, double rho) {
        int variants = base.variants();
        if (rho == 0) return base;
        double[] scores = base.scoresView();
        double[] information = base.informationView();
        if (rho == 1) {
            return new SetTestScoreState(
                new double[] {Arrays.stream(scores).sum()},
                new double[] {Arrays.stream(information).sum()}, 1);
        }
        int dimension = variants + 1;
        double[] transformedScores = new double[dimension];
        double[] transformedInformation = new double[dimension * dimension];
        double variantScale = Math.sqrt(1 - rho);
        double burdenScale = Math.sqrt(rho);
        double crossScale = variantScale * burdenScale;
        for (int left = 0; left < variants; left++) {
            transformedScores[left] = variantScale * scores[left];
            double rowSum = 0;
            for (int right = 0; right < variants; right++) {
                double value = information[left * variants + right];
                transformedInformation[left * dimension + right] =
                    (1 - rho) * value;
                rowSum += value;
            }
            transformedInformation[left * dimension + variants] =
                crossScale * rowSum;
            transformedInformation[variants * dimension + left] =
                crossScale * rowSum;
        }
        transformedScores[variants] = burdenScale
            * Arrays.stream(scores).sum();
        transformedInformation[dimension * dimension - 1] = rho
            * Arrays.stream(information).sum();
        return new SetTestScoreState(
            transformedScores, transformedInformation, dimension);
    }

    private static void requireCompatible(
            PreparedVariantSet prepared, SetTestScoreNullModel nullModel) {
        if (prepared == null || nullModel == null)
            throw new IllegalArgumentException(
                "prepared variant set and null model are required");
        if (prepared.observations() != nullModel.observations())
            throw new IllegalArgumentException(
                "prepared variant set does not match null-model observations");
    }

    private static double[] eigenvalues(double[] matrix, int dimension) {
        try (BackendContext context = BackendContext.select(
                org.jlinalg.compute.BackendPolicy.PREFERRED)) {
            SymmetricEigenDecomposition decomposition =
                context.backend().dsyev(symmetricCopy(matrix, dimension), dimension);
            double maximum = Arrays.stream(decomposition.eigenvalues())
                .map(Math::abs).max().orElse(0);
            double tolerance = 1e-12 * Math.max(1, maximum);
            return Arrays.stream(decomposition.eigenvalues())
                .filter(value -> value > tolerance)
                .sorted().toArray();
        }
    }

    private static double[] symmetricCopy(double[] matrix, int dimension) {
        double[] result = matrix.clone();
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < row; column++) {
                double value = 0.5 * (result[row * dimension + column]
                    + result[column * dimension + row]);
                result[row * dimension + column] = value;
                result[column * dimension + row] = value;
            }
        }
        return result;
    }

    private static double dot(double[] left, double[] right) {
        double result = 0;
        for (int index = 0; index < left.length; index++)
            result += left[index] * right[index];
        return result;
    }

    private static final class ScoreSampler {
        private final double[] leftVectors;
        private final double[] squareRootValues;
        private final int dimension;

        private ScoreSampler(double[] covariance, int dimension) {
            this.dimension = dimension;
            try (BackendContext context = BackendContext.select(
                    org.jlinalg.compute.BackendPolicy.PREFERRED)) {
                SymmetricEigenDecomposition decomposition =
                    context.backend().dsyev(
                        symmetricCopy(covariance, dimension), dimension);
                leftVectors = decomposition.eigenvectors();
                double[] eigenvalues = decomposition.eigenvalues();
                squareRootValues = new double[dimension];
                for (int index = 0; index < dimension; index++)
                    squareRootValues[index] = Math.sqrt(
                        Math.max(0, eigenvalues[index]));
            }
        }

        private double[] sample(Random random, int samples) {
            double[] gaussian = new double[dimension * samples];
            for (int row = 0; row < dimension; row++)
                for (int sample = 0; sample < samples; sample++)
                    gaussian[row * samples + sample] = random.nextGaussian()
                        * squareRootValues[row];
            try (BackendContext context = BackendContext.select(
                    org.jlinalg.compute.BackendPolicy.PREFERRED)) {
                return org.jlinalg.internal.MatrixOps.multiply(
                    context.backend(), leftVectors, dimension, dimension,
                    gaussian, samples);
            }
        }
    }
}
