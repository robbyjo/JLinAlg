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
import jdistlib.accelerator.SingularValueDecomposition;
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
        Prepared prepared = prepare(set, nullModel, options);
        double[] burden = new double[nullModel.observations()];
        for (int variant = 0; variant < prepared.dosages().length; variant++)
            for (int sample = 0; sample < burden.length; sample++)
                burden[sample] += prepared.weights()[variant]
                    * prepared.dosages()[variant][sample];
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
            prepared.dosages().length, statistic, beta, standardError,
            degrees, pValue, logP, "student-t", new double[0],
            prepared.excluded());
    }

    /** Related-sample burden score test using one retained REML projection. */
    public static SetTestResult burden(
            VariantSet set, RemlSetTestNullModel nullModel,
            SetTestOptions options) {
        Prepared prepared = prepare(set, nullModel, options);
        double[] burden = new double[nullModel.observations()];
        for (int variant = 0; variant < prepared.dosages().length; variant++)
            for (int sample = 0; sample < burden.length; sample++)
                burden[sample] += prepared.weights()[variant]
                    * prepared.dosages()[variant][sample];
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
            set.variants().size(), prepared.dosages().length, statistic,
            beta, standardError, degrees, pValue, logP, "student-t",
            new double[0], prepared.excluded());
    }

    public static SetTestResult skat(
            VariantSet set, GaussianSetTestNullModel nullModel,
            SetTestOptions options) {
        Prepared prepared = prepare(set, nullModel, options);
        double[][] weighted = weighted(
            prepared.dosages(), prepared.weights());
        return kernelResult(set.id(), "skat", set.variants().size(),
            prepared.dosages().length, weighted, nullModel,
            prepared.excluded());
    }

    public static SkatOResult skatO(
            VariantSet set, GaussianSetTestNullModel nullModel,
            SetTestOptions options) {
        Prepared prepared = prepare(set, nullModel, options);
        double[][] base = weighted(
            prepared.dosages(), prepared.weights());
        double[] burden = new double[nullModel.observations()];
        for (double[] variant : base)
            for (int sample = 0; sample < burden.length; sample++)
                burden[sample] += variant[sample];
        double[] rhoGrid = options.skatORhoGrid();
        List<SkatOResult.Component> components = new ArrayList<>(rhoGrid.length);
        double minimumP = 1;
        for (double rho : rhoGrid) {
            double[][] factor = factor(base, burden, rho);
            SetTestResult result = kernelResult(set.id(),
                "skat-o[rho=" + rho + "]", set.variants().size(),
                prepared.dosages().length, factor, nullModel,
                prepared.excluded());
            components.add(new SkatOResult.Component(rho, result));
            minimumP = Math.min(minimumP, result.pValue());
        }
        double[] critical = new double[components.size()];
        for (int index = 0; index < critical.length; index++)
            critical[index] = QuadraticFormDistribution.critical(
                components.get(index).result().eigenvalues(), minimumP);
        Random random = new Random(options.randomSeed());
        SetTestScoreState baseState = nullModel.score(base);
        ScoreSampler sampler = new ScoreSampler(
            baseState.informationView(), base.length);
        int extreme = 0;
        for (int simulation = 0;
                simulation < options.skatOSimulations(); simulation++) {
            double[] simulatedScores = sampler.sample(random);
            double squaredSum = 0;
            double sum = 0;
            for (double score : simulatedScores) {
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
        return new SkatOResult(set.id(), set.variants().size(),
            prepared.dosages().length, components, minimumP, adjusted,
            Math.log10(adjusted), options.skatOSimulations(),
            options.randomSeed(), prepared.excluded());
    }

    private static SetTestResult kernelResult(
            String setId, String method, int requested,
            int included, double[][] factor,
            GaussianSetTestNullModel nullModel,
            List<VariantFilterResult> excluded) {
        SetTestScoreState scoreState = nullModel.score(factor);
        double statistic = 0;
        for (double score : scoreState.scoresView())
            statistic += score * score;
        double[] eigenvalues = eigenvalues(
            scoreState.informationView(), factor.length);
        QuadraticFormDistribution.Tail tail =
            QuadraticFormDistribution.survival(statistic, eigenvalues);
        return new SetTestResult(setId, method, requested, included,
            statistic, Double.NaN, Double.NaN, Double.NaN,
            tail.pValue(), Math.log10(tail.pValue()), tail.method(),
            eigenvalues, excluded);
    }

    private static Prepared prepare(
            VariantSet set, GaussianSetTestNullModel nullModel,
            SetTestOptions options) {
        if (set == null || nullModel == null || options == null)
            throw new IllegalArgumentException(
                "variant set, null model, and options are required");
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
                options.missingPolicy(), nullModel.observations()));
            weights.add(member.weight());
        }
        if (dosages.isEmpty())
            throw new IllegalArgumentException(
                "variant set has no members after filtering: " + set.id());
        double[] numericWeights = new double[weights.size()];
        for (int index = 0; index < numericWeights.length; index++)
            numericWeights[index] = weights.get(index);
        return new Prepared(dosages.toArray(double[][]::new),
            numericWeights, excluded);
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

    private static double[][] factor(
            double[][] base, double[] burden, double rho) {
        int rows = (rho < 1 ? base.length : 0) + (rho > 0 ? 1 : 0);
        double[][] result = new double[rows][];
        int destination = 0;
        if (rho < 1) {
            double scale = Math.sqrt(1 - rho);
            for (double[] variant : base) {
                result[destination] = variant.clone();
                for (int sample = 0; sample < variant.length; sample++)
                    result[destination][sample] *= scale;
                destination++;
            }
        }
        if (rho > 0) {
            result[destination] = burden.clone();
            double scale = Math.sqrt(rho);
            for (int sample = 0; sample < burden.length; sample++)
                result[destination][sample] *= scale;
        }
        return result;
    }

    private static double[] eigenvalues(double[] matrix, int dimension) {
        try (BackendContext context = BackendContext.select(
                org.jlinalg.compute.BackendPolicy.CPU)) {
            SingularValueDecomposition decomposition =
                context.backend().dgesvd(matrix, dimension, dimension);
            return Arrays.stream(decomposition.singularValues())
                .filter(value -> value > 1e-12)
                .sorted().toArray();
        }
    }

    private static double dot(double[] left, double[] right) {
        double result = 0;
        for (int index = 0; index < left.length; index++)
            result += left[index] * right[index];
        return result;
    }

    private record Prepared(
        double[][] dosages, double[] weights,
        List<VariantFilterResult> excluded) { }

    private static final class ScoreSampler {
        private final double[] leftVectors;
        private final double[] squareRootValues;
        private final int dimension;

        private ScoreSampler(double[] covariance, int dimension) {
            this.dimension = dimension;
            try (BackendContext context = BackendContext.select(
                    org.jlinalg.compute.BackendPolicy.CPU)) {
                SingularValueDecomposition decomposition =
                    context.backend().dgesvd(
                        covariance, dimension, dimension);
                leftVectors = decomposition.leftSingularVectors();
                double[] singularValues = decomposition.singularValues();
                squareRootValues = new double[dimension];
                for (int index = 0;
                        index < Math.min(dimension, singularValues.length);
                        index++)
                    squareRootValues[index] = Math.sqrt(
                        Math.max(0, singularValues[index]));
            }
        }

        private double[] sample(Random random) {
            double[] gaussian = new double[dimension];
            for (int index = 0; index < dimension; index++)
                gaussian[index] = random.nextGaussian()
                    * squareRootValues[index];
            double[] result = new double[dimension];
            for (int row = 0; row < dimension; row++)
                for (int column = 0; column < dimension; column++)
                    result[row] += leftVectors[row * dimension + column]
                        * gaussian[column];
            return result;
        }
    }
}
