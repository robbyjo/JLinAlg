/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.ToDoubleFunction;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.Optimization;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceEstimation;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

/** Sparse grouped Gaussian likelihood with one estimated Cholesky factor per block. */
public final class SparseUnstructuredCorrelatedModel {
    private SparseUnstructuredCorrelatedModel() { }

    public static Result fit(double[] response, double[] fixed, int rows, int columns,
            List<String> groups, List<String> effectNames, double[][] design,
            RemlOptions options, BackendPolicy policy) {
        return fit(response, fixed, rows, columns, List.of(CorrelatedRandomEffectBlock.of(
            "unstructured", groups, effectNames, design)), options, policy);
    }

    /**
     * Fits b_g = L u_g with Var(u_g) = sigma^2 I. L contains all block parameters;
     * there is no redundant block scale. Only sigma^2 is profiled analytically.
     * Zero non-leading Cholesky diagonals permit singular covariance fits;
     * marginal random variances and residual variance obey the physical bounds.
     */
    public static Result fit(double[] response, double[] fixed, int rows, int columns,
            List<CorrelatedRandomEffectBlock> blocks, RemlOptions options, BackendPolicy policy) {
        return fit(response,fixed,rows,columns,blocks,List.of(),options,policy);
    }

    /** Jointly optimizes unstructured blocks and pedigree variance components. */
    public static Result fit(double[] response,double[] fixed,int rows,int columns,
            List<CorrelatedRandomEffectBlock> blocks,
            List<PedigreeRandomEffectTerm> pedigreeEffects,
            RemlOptions options,BackendPolicy policy) {
        try (Model model = new Model(response, fixed, rows, columns, blocks,
                pedigreeEffects, options, policy)) {
        Optimum optimum = model.optimize(model.initial(), model::objective);
        try (SparseLinearMixedModel.Prepared prepared = model.prepare(optimum.point())) {
            double[] scaleBounds = model.scaleBounds(optimum.point());
            SparseLinearMixedModelResult fitted = prepared.evaluateAt(response, fixed,
                columns, new double[model.termCount()], scaleBounds[0],
                scaleBounds[1], model.pevBlockSizes());
            double scale = fitted.varianceComponents()[model.termCount()];
            fitted=fitted.withVarianceComponents(model.reportedVariances(
                optimum.point(),scale));
            if (options.degreesOfFreedomMethod() != DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION) {
                if (!optimum.converged()) throw new IllegalStateException("finite-DF inference requires a converged fit");
                if (options.varianceEstimation() != VarianceEstimation.REML)
                    throw new IllegalArgumentException("finite-DF inference requires REML");
                try (BackendContext context = policy == BackendPolicy.PREFERRED
                        ? BackendContext.preferredSparse() : BackendContext.select(policy)) {
                    double[] absolute = model.absolute(optimum.point(), scale);
                    double[] marginal = new double[1 + pedigreeEffects.size()
                        + blocks.stream().mapToInt(CorrelatedRandomEffectBlock::effectCount).sum()];
                    int index = 0, start = 0;
                    for (CorrelatedRandomEffectBlock block : blocks) {
                        for (int r = 0; r < block.effectCount(); r++) marginal[index++] = absolute[start + r*(r+1)/2 + r];
                        start += block.effectCount() * (block.effectCount() + 1) / 2;
                    }
                    for(int extra=0;extra<pedigreeEffects.size();extra++)
                        marginal[index++]=absolute[model.blockParameters()+extra];
                    marginal[index] = scale;
                    SparseFiniteDf.requireInteriorVariances(marginal, options.minimumVariance(), options.maximumVariance());
                    SparseFiniteDf.Inference inference = SparseFiniteDf.compute(absolute, model.derivativeScales(absolute),
                        model::absolutePoint, columns, options.degreesOfFreedomMethod(), context.backend());
                    fitted = fitted.withInference(inference.covariance(), inference.degreesOfFreedom(),
                        options.degreesOfFreedomMethod(), inference.jointState());
                }
            }
            return new Result(fitted.withOptimization(optimum.evaluations(), optimum.converged()),
                blocks, model.factors(optimum.point()));
        }
        }
    }

    /** ML coefficient profile; all other fixed effects and covariance entries are refitted. */
    public static ProfileLikelihoodInterval profileFixedEffect(double[] response, double[] fixed,
            int rows, int columns, List<CorrelatedRandomEffectBlock> blocks,
            int coefficient, double confidence, double lower, double upper,
            RemlOptions options, BackendPolicy policy) {
        if (coefficient < 0 || coefficient >= columns)
            throw new IllegalArgumentException("coefficient index out of range");
        RemlOptions ml = options.toBuilder().varianceEstimation(VarianceEstimation.ML)
            .degreesOfFreedomMethod(DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION).build();
        try (Model model = new Model(response, fixed, rows, columns, blocks, ml, policy)) {
        Optimum optimum = model.optimize(model.initial(), model::objective);
        if (!optimum.converged()) throw new IllegalStateException("profile maximum did not converge");
        double[] reduced = new double[rows * (columns - 1)];
        for (int r = 0; r < rows; r++) for (int c = 0, k = 0; c < columns; c++)
            if (c != coefficient) reduced[r * (columns - 1) + k++] = fixed[r * columns + c];
        double estimate;
        try (SparseLinearMixedModel.Prepared prepared = model.prepare(optimum.point())) {
            estimate = prepared.evaluateAt(response, fixed, columns, new double[blocks.size()]).beta()[coefficient];
        }
        return MixedModelProfile.interval(value -> {
            double[] adjusted = response.clone();
            for (int r = 0; r < rows; r++) adjusted[r] -= fixed[r * columns + coefficient] * value;
            try (Model constrained = new Model(adjusted, reduced, rows, columns - 1, blocks, ml, policy)) {
            Optimum refit = constrained.optimize(optimum.point(), constrained::objective);
            if (!refit.converged()) throw new IllegalStateException("constrained profile refit did not converge");
            return -refit.value();
            }
        }, estimate, -optimum.value(), confidence, lower, upper, 12);
        }
    }

    private record Optimum(double[] point, double value, int evaluations, boolean converged) { }

    /** ML residual-SD profile with every block covariance refitted. */
    public static ProfileLikelihoodInterval profileResidualSd(double[] response, double[] fixed,
            int rows, int columns, List<CorrelatedRandomEffectBlock> blocks,
            double confidence, double lower, double upper, RemlOptions options, BackendPolicy policy) {
        if (!(lower > 0)) throw new IllegalArgumentException("residual SD lower bound must be positive");
        try (Model model = new Model(response, fixed, rows, columns, blocks, profileOptions(options), policy)) {
        Optimum optimum = model.optimize(model.initial(), model::objective);
        if (!optimum.converged()) throw new IllegalStateException("profile maximum did not converge");
        double sd = Math.sqrt(model.point(optimum.point(), Double.NaN).residualVariance());
        return MixedModelProfile.interval(value -> {
            Optimum refit = model.optimize(optimum.point(), theta ->
                -model.point(theta, value * value).logLikelihood());
            if (!refit.converged()) throw new IllegalStateException("residual SD profile did not converge");
            return -refit.value();
        }, sd, -optimum.value(), confidence, lower, upper, 12);
        }
    }

    /** ML marginal random-effect SD profile, including a possible zero-SD boundary. */
    public static ProfileLikelihoodInterval profileRandomSd(double[] response, double[] fixed,
            int rows, int columns, List<CorrelatedRandomEffectBlock> blocks, int block, int effect,
            double confidence, double lower, double upper, RemlOptions options, BackendPolicy policy) {
        if (block < 0 || block >= blocks.size() || effect < 0 || effect >= blocks.get(block).effectCount() || lower < 0)
            throw new IllegalArgumentException("random SD profile target is invalid");
        try (Model model = new Model(response, fixed, rows, columns, blocks, profileOptions(options), policy)) {
        Optimum optimum = model.optimize(model.initial(), model::objective);
        if (!optimum.converged()) throw new IllegalStateException("profile maximum did not converge");
        double scale = model.point(optimum.point(), Double.NaN).residualVariance();
        int start = 0;
        for (int i = 0; i < block; i++) start += blocks.get(i).effectCount() * (blocks.get(i).effectCount() + 1) / 2;
        final int rowStart = start + effect * (effect + 1) / 2, diagonal = rowStart + effect;
        double estimate = Math.sqrt(covariance(model.factors(optimum.point()).get(block), scale)[effect][effect]);
        boolean boundaryEstimate = estimate <= lower
            + 1e-7*Math.max(1,upper-lower);
        // Replace the constrained row diagonal with log(sigma); offdiagonals
        // become row ratios. This enforces a fixed row norm without PD rejection.
        double[] initial = optimum.point().clone(), lo = model.lower.clone(), hi = model.upper.clone();
        for (int j = rowStart; j < diagonal; j++) initial[j] /= Math.max(1e-12, initial[diagonal]);
        initial[diagonal] = .5 * Math.log(scale);
        lo[diagonal] = .5 * Math.log(options.minimumVariance());
        hi[diagonal] = .5 * Math.log(options.maximumVariance());
        java.util.function.DoubleUnaryOperator profile = value -> {
            Optimum refit = model.optimize(initial, free -> {
                double residual = Math.exp(2 * free[diagonal]);
                double norm = 1;
                for (int j = rowStart; j < diagonal; j++) norm += free[j] * free[j];
                double[] theta = free.clone();
                theta[diagonal] = value / Math.sqrt(residual * norm);
                for (int j = rowStart; j < diagonal; j++) theta[j] = free[j] * theta[diagonal];
                return -model.point(theta, residual).logLikelihood();
            }, lo, hi);
            if (!refit.converged()) throw new IllegalStateException("random SD profile did not converge");
            return -refit.value();
        };
        ProfileLikelihoodInterval interval = boundaryEstimate
            ? MixedModelProfile.boundaryInterval(profile, lower,
                -optimum.value(), confidence, lower, upper, 12)
            : MixedModelProfile.interval(profile, estimate,
                -optimum.value(), confidence, lower, upper, 12);
        if (!interval.lowerFound() && lower == 0)
            return new ProfileLikelihoodInterval(interval.estimate(), 0, interval.upper(), interval.cutoff(),
                false, interval.upperFound());
        return interval;
        }
    }

    private static RemlOptions profileOptions(RemlOptions options) {
        return options.toBuilder().varianceEstimation(VarianceEstimation.ML)
            .degreesOfFreedomMethod(DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION).build();
    }

    /**
     * Marginal-correlation profile for the first two coefficients of a block,
     * including +/-1 boundaries.
     */
    public static ProfileLikelihoodInterval profileCorrelation(double[] response, double[] fixed,
            int rows, int columns, List<CorrelatedRandomEffectBlock> blocks, int block,
            double confidence, RemlOptions options, BackendPolicy policy) {
        if (block < 0 || block >= blocks.size()
                || blocks.get(block).effectCount() < 2)
            throw new IllegalArgumentException(
                "correlation profiles require a block with at least two coefficients");
        try (Model model = new Model(response, fixed, rows, columns, blocks, profileOptions(options), policy)) {
        Optimum optimum = model.optimize(model.initial(), model::objective);
        if (!optimum.converged()) throw new IllegalStateException("profile maximum did not converge");
        int start = 0;
        for (int i = 0; i < block; i++) start += blocks.get(i).effectCount() * (blocks.get(i).effectCount() + 1) / 2;
        final int offDiagonal = start + 1;
        double[] full = optimum.point().clone();
        double slopeSd = Math.hypot(full[offDiagonal], full[offDiagonal + 1]);
        double estimate = full[offDiagonal] / slopeSd;
        if (!Double.isFinite(estimate)) throw new IllegalArgumentException("correlation is undefined at zero random SD");
        estimate = Math.max(-1, Math.min(1, estimate));
        full[offDiagonal + 1] = slopeSd;
        double[] initial = new double[full.length - 1], lo = new double[initial.length], hi = new double[initial.length];
        for (int i = 0, k = 0; i < full.length; i++) if (i != offDiagonal) {
            initial[k] = full[i]; lo[k] = model.lower[i]; hi[k++] = model.upper[i];
        }
        return MixedModelProfile.interval(value -> {
            Optimum refit = model.optimize(initial, free -> {
                double[] theta = new double[full.length];
                for (int i = 0, k = 0; i < theta.length; i++) if (i != offDiagonal) theta[i] = free[k++];
                double sd = theta[offDiagonal + 1];
                theta[offDiagonal] = value * sd;
                theta[offDiagonal + 1] = Math.sqrt(1 - value * value) * sd;
                return model.objective(theta);
            }, lo, hi);
            if (!refit.converged()) throw new IllegalStateException("correlation profile did not converge");
            return -refit.value();
        }, estimate, -optimum.value(), confidence, -1, 1, 12);
        }
    }

    /**
     * Profiles the marginal correlation of any two coefficients in an
     * unstructured block. Reordering the selected pair into the leading
     * Cholesky coordinates leaves the covariance family and likelihood
     * unchanged while retaining the exact +/-1 boundary representation.
     */
    public static ProfileLikelihoodInterval profileCorrelation(
            double[] response, double[] fixed, int rows, int columns,
            List<CorrelatedRandomEffectBlock> blocks, int block,
            int firstEffect, int secondEffect, double confidence,
            RemlOptions options, BackendPolicy policy) {
        if (blocks == null || block < 0 || block >= blocks.size()) {
            throw new IllegalArgumentException(
                "correlation profile block is out of range");
        }
        CorrelatedRandomEffectBlock selected = blocks.get(block);
        if (selected.observations() != rows) {
            throw new IllegalArgumentException(
                "correlated block rows must match model rows");
        }
        int effects = selected.effectCount();
        if (firstEffect < 0 || firstEffect >= effects
                || secondEffect < 0 || secondEffect >= effects
                || firstEffect == secondEffect) {
            throw new IllegalArgumentException(
                "correlation profile effects must be distinct and in range");
        }
        int[] order = new int[effects];
        order[0] = firstEffect;
        order[1] = secondEffect;
        for (int effect = 0, target = 2; effect < effects; effect++) {
            if (effect != firstEffect && effect != secondEffect) {
                order[target++] = effect;
            }
        }
        int[] groupIndices = selected.groupIndices();
        List<String> groups = new ArrayList<>(rows);
        for (int row = 0; row < rows; row++) {
            groups.add(selected.groupNames().get(groupIndices[row]));
        }
        double[] source = selected.effectDesign();
        double[][] design = new double[rows][effects];
        List<String> names = new ArrayList<>(effects);
        for (int target = 0; target < effects; target++) {
            int sourceEffect = order[target];
            names.add(selected.effectNames().get(sourceEffect));
            for (int row = 0; row < rows; row++) {
                design[row][target] = source[row * effects + sourceEffect];
            }
        }
        List<CorrelatedRandomEffectBlock> reordered =
            new ArrayList<>(blocks);
        reordered.set(block, CorrelatedRandomEffectBlock.of(
            selected.name(), groups, names, design));
        return profileCorrelation(response, fixed, rows, columns,
            reordered, block, confidence, options, policy);
    }

    private static final class Model implements AutoCloseable {
        final double[] response, fixed, lower, upper;
        final int rows, columns, parameters;
        final int blockParameters;
        final List<CorrelatedRandomEffectBlock> blocks;
        final List<PedigreeRandomEffectTerm> pedigreeEffects;
        final RemlOptions options;
        final BackendPolicy policy;
        final SparseLinearMixedModel.TransformedLikelihood likelihood;
        Model(double[] response, double[] fixed, int rows, int columns,
                List<CorrelatedRandomEffectBlock> blocks, RemlOptions options, BackendPolicy policy) {
            this(response,fixed,rows,columns,blocks,List.of(),options,policy);
        }
        Model(double[] response,double[] fixed,int rows,int columns,
                List<CorrelatedRandomEffectBlock> blocks,
                List<PedigreeRandomEffectTerm> pedigreeEffects,
                RemlOptions options,BackendPolicy policy) {
            if (columns > 0) MatrixOps.validateModelData(response, fixed, rows, columns);
            else if (columns != 0 || response == null || response.length != rows || fixed.length != 0)
                throw new IllegalArgumentException("invalid constrained design");
            if (blocks == null || blocks.isEmpty() || options == null || policy == null)
                throw new IllegalArgumentException("blocks and controls are required");
            this.response = response; this.fixed = fixed; this.rows = rows; this.columns = columns;
            this.blocks = List.copyOf(blocks); this.pedigreeEffects=List.copyOf(pedigreeEffects);
            this.options = options; this.policy = policy;
            int count = 0;
            java.util.HashSet<String> names = new java.util.HashSet<>();
            for (CorrelatedRandomEffectBlock block : blocks) {
                if (block.observations() != rows || !names.add(block.name()))
                    throw new IllegalArgumentException("block rows must match and names must be unique");
                count += block.effectCount() * (block.effectCount() + 1) / 2;
            }
            blockParameters=count;
            for(PedigreeRandomEffectTerm pedigree:this.pedigreeEffects)
                if(pedigree==null||pedigree.randomEffect().observations()!=rows
                        ||!names.add(pedigree.randomEffect().name()))
                    throw new IllegalArgumentException("pedigree term rows must match and names must be unique");
            parameters = count+this.pedigreeEffects.size();
            lower = new double[parameters]; upper = new double[parameters];
            double limit = Math.sqrt(options.maximumVariance() / options.minimumVariance());
            Arrays.fill(upper, limit);
            for (int b = 0, k = 0; b < blocks.size(); b++)
                for (int r = 0; r < blocks.get(b).effectCount(); r++) for (int c = 0; c <= r; c++)
                    lower[k++] = r == c ? (r == 0 ? 1 / limit : 0) : -limit;
            for(int k=blockParameters;k<parameters;k++){
                lower[k]=-.5*Math.log(limit*limit);upper[k]=.5*Math.log(limit*limit);
            }
            likelihood = new SparseLinearMixedModel.TransformedLikelihood(policy);
        }
        public void close() { likelihood.close(); }
        double[] initial() {
            double[] result = new double[parameters];
            double[] supplied = options.initialVariances();
            if (supplied != null && supplied.length != termCount() + 1)
                throw new IllegalArgumentException("initial variances require one block scale plus residual");
            int residualIndex=termCount();
            for (int b = 0, i = 0; b < blocks.size(); b++) {
                double sd = supplied == null ? 1 : Math.sqrt(supplied[b] / supplied[residualIndex]);
                int d = blocks.get(b).effectCount();
                for (int r = 0; r < d; r++) for (int c = 0; c <= r; c++) result[i++] = r == c ? sd : 0;
            }
            for(int extra=0;extra<pedigreeEffects.size();extra++)
                result[blockParameters+extra]=supplied==null?0:.5*Math.log(
                    supplied[blocks.size()+extra]/supplied[residualIndex]);
            return result;
        }
        List<double[][]> factors(double[] point) {
            List<double[][]> result = new ArrayList<>(); int k = 0;
            for (CorrelatedRandomEffectBlock block : blocks) {
                int d = block.effectCount(); double[][] l = new double[d][d];
                for (int r = 0; r < d; r++) for (int c = 0; c <= r; c++) l[r][c] = point[k++];
                result.add(l);
            }
            return result;
        }
        List<RandomEffectTerm> terms(double[] point) {
            List<double[][]> factors = factors(point);
            List<RandomEffectTerm> terms = new ArrayList<>();
            for (int b = 0; b < blocks.size(); b++) {
                CorrelatedRandomEffectBlock block = blocks.get(b);
                int d = block.effectCount(); int[] groups = block.groupIndices();
                double[] design = block.effectDesign(); double[][] l = factors.get(b);
                int[] starts = new int[rows + 1], indices = new int[rows * d];
                double[] values = new double[rows * d]; List<String> names = new ArrayList<>();
                for (String group : block.groupNames()) for (String effect : block.effectNames())
                    names.add(group + ":" + effect);
                for (int r = 0; r < rows; r++) {
                    starts[r] = r * d;
                    for (int c = 0; c < d; c++) {
                        indices[r * d + c] = groups[r] * d + c;
                        for (int j = c; j < d; j++) values[r * d + c] += design[r * d + j] * l[j][c];
                    }
                }
                starts[rows] = rows * d;
                terms.add(RandomEffectTerm.ofSparseCsr("latent:" + block.name(), rows,
                    block.groupCount() * d, starts, indices, values, names));
            }
            for(int extra=0;extra<pedigreeEffects.size();extra++)
                terms.add(scale(pedigreeEffects.get(extra).randomEffect(),
                    Math.exp(point[blockParameters+extra])));
            return terms;
        }
        List<SparsePrecisionMatrix> precisions(){
            List<SparsePrecisionMatrix> result=new ArrayList<>();
            for(CorrelatedRandomEffectBlock block:blocks)
                result.add(SparsePrecisionMatrix.identity(block.groupCount()*block.effectCount()));
            for(PedigreeRandomEffectTerm pedigree:pedigreeEffects)result.add(pedigree.precision());
            return result;
        }
        private RandomEffectTerm scale(RandomEffectTerm term,double multiplier){
            if(term.sparse()){
                double[] values=term.sparseValues();for(int i=0;i<values.length;i++)values[i]*=multiplier;
                return RandomEffectTerm.ofSparseCsr(term.name(),term.observations(),term.coefficients(),
                    term.rowPointers(),term.columnIndices(),values,term.coefficientNames());
            }
            double[] values=term.design();for(int i=0;i<values.length;i++)values[i]*=multiplier;
            return RandomEffectTerm.of(term.name(),values,term.observations(),term.coefficients(),term.coefficientNames());
        }
        SparseLinearMixedModel.Prepared prepare(double[] point) {
            return SparseLinearMixedModel.prepareWithPrecision(rows, terms(point),precisions(), options.toBuilder()
                .degreesOfFreedomMethod(DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION).build(), policy);
        }
        int[] pevBlockSizes(){
            int[] result=new int[blocks.size()+pedigreeEffects.size()];
            for(int i=0;i<blocks.size();i++)result[i]=blocks.get(i).effectCount();
            for(int i=blocks.size();i<result.length;i++)result[i]=1;
            return result;
        }
        int termCount(){return blocks.size()+pedigreeEffects.size();}
        int blockParameters(){return blockParameters;}
        double[] reportedVariances(double[] theta,double scale){
            double[] result=new double[termCount()+1];
            Arrays.fill(result,0,blocks.size(),scale);
            for(int extra=0;extra<pedigreeEffects.size();extra++)
                result[blocks.size()+extra]=scale*Math.exp(
                    2*theta[blockParameters+extra]);
            result[termCount()]=scale;return result;
        }
        SparseFiniteDf.Point point(double[] theta, double scale) {
            double[] bounds = scaleBounds(theta);
            return likelihood.point(response, fixed, rows, columns, terms(theta),precisions(), options.varianceEstimation(), scale,
                bounds[0], bounds[1]);
        }
        double[] scaleBounds(double[] theta) {
            double lo = options.minimumVariance(), hi = options.maximumVariance();
            for (double[][] factor : factors(theta)) for (int r = 0; r < factor.length; r++) {
                double ratio = 0;
                for (int c = 0; c <= r; c++) ratio += factor[r][c] * factor[r][c];
                if (!(ratio > 0)) throw new IllegalArgumentException("random variance below physical lower bound");
                lo = Math.max(lo, options.minimumVariance() / ratio);
                hi = Math.min(hi, options.maximumVariance() / ratio);
            }
            for(int extra=0;extra<pedigreeEffects.size();extra++){
                double ratio=Math.exp(2*theta[blockParameters+extra]);
                lo=Math.max(lo,options.minimumVariance()/ratio);
                hi=Math.min(hi,options.maximumVariance()/ratio);
            }
            if (lo > hi * (1 + 1e-12)) throw new IllegalArgumentException("infeasible covariance bounds");
            return new double[] {Math.min(lo, hi), hi};
        }
        double objective(double[] theta) {
            try {
                double value = -point(theta, Double.NaN).logLikelihood();
                return Double.isFinite(value) ? value : 1e100;
            } catch (IllegalArgumentException | IllegalStateException failure) { return 1e100; }
        }
        Optimum optimize(double[] initial, ToDoubleFunction<double[]> objective) {
            Optimum best = optimize(initial, objective, lower, upper);
            // At a zero Cholesky diagonal, coordinate-wise stationarity does not
            // exclude a better full-rank solution: off-diagonal coordinates may
            // have to move together with that diagonal. Check distinct interior
            // starts before accepting a singular solution (including profiles).
            boolean singular = false;
            for (double[][] factor : factors(best.point())) if (factor.length > 1) {
                double largest = 0;
                for (int r = 0; r < factor.length; r++)
                    for (int c = 0; c <= r; c++) largest = Math.max(largest, Math.abs(factor[r][c]));
                for (int r = 0; r < factor.length; r++)
                    singular |= factor[r][r] < 1e-5 * Math.max(1, largest);
            }
            int evaluations = best.evaluations();
            if (singular) for (double sd : new double[] {.3, 3}) {
                double[] restart = new double[parameters]; int k = 0;
                for (CorrelatedRandomEffectBlock block : blocks)
                    for (int r = 0; r < block.effectCount(); r++) for (int c = 0; c <= r; c++) {
                        restart[k] = Math.max(lower[k], Math.min(upper[k], r == c ? sd : 0));
                        k++;
                    }
                Optimum alternative = optimize(restart, objective, lower, upper);
                evaluations += alternative.evaluations();
                if (alternative.value() < best.value()) best = alternative;
            }
            return new Optimum(best.point(), best.value(), evaluations, best.converged());
        }
        Optimum optimize(double[] initial, ToDoubleFunction<double[]> objective, double[] lower, double[] upper) {
            int parameters = initial.length;
            int[] calls = {0};
            ToDoubleFunction<double[]> counted = x -> {
                calls[0]++;
                try {
                    double value = objective.applyAsDouble(x);
                    return Double.isFinite(value) ? value : 1e100;
                } catch (IllegalArgumentException | IllegalStateException failure) { return 1e100; }
            };
            double[] point;
            int limit = Math.max(200, options.maximumIterations() * 20);
            try {
                if (parameters == 1) {
                    double distance = 1, left = Math.max(lower[0], initial[0] - distance);
                    double right = Math.min(upper[0], initial[0] + distance);
                    double center = counted.applyAsDouble(initial);
                    while (distance < 1e6 && ((left > lower[0] && counted.applyAsDouble(new double[] {left}) < center)
                            || (right < upper[0] && counted.applyAsDouble(new double[] {right}) < center))) {
                        distance *= 2;
                        left = Math.max(lower[0], initial[0] - distance);
                        right = Math.min(upper[0], initial[0] + distance);
                    }
                    double value = Optimization.optimize(x -> counted.applyAsDouble(new double[] {x}),
                        left, right, Math.max(1e-8, options.relativeTolerance()), limit);
                    point = new double[] {value};
                    if (counted.applyAsDouble(new double[] {left}) < counted.applyAsDouble(point)) point[0] = left;
                    if (counted.applyAsDouble(new double[] {right}) < counted.applyAsDouble(point)) point[0] = right;
                } else {
                    point = Bobyqa.bobyqa(initial.clone(), lower, upper, counted::applyAsDouble,
                        2 * parameters + 1, .2, Math.max(1e-8, options.relativeTolerance()), limit, true).mX;
                }
            } catch (ArithmeticException | ArrayIndexOutOfBoundsException failure) { point = initial.clone(); }
            double best = counted.applyAsDouble(point), step = .01;
            while (calls[0] < limit && step > 1e-6) {
                boolean improved = false;
                for (int i = 0; i < parameters; i++) for (int sign : new int[] {-1, 1}) {
                    double[] candidate = point.clone();
                    candidate[i] = Math.max(lower[i], Math.min(upper[i], point[i] + sign * step));
                    double value = counted.applyAsDouble(candidate);
                    if (value < best - 1e-12) { best = value; point = candidate; improved = true; }
                }
                if (!improved) step *= .5;
            }
            boolean stationary = best < 1e100;
            for (int i = 0; i < parameters; i++) for (int sign : new int[] {-1, 1}) {
                double[] candidate = point.clone(); double delta = 1e-4 * Math.max(1, Math.abs(point[i]));
                candidate[i] = Math.max(lower[i], Math.min(upper[i], point[i] + sign * delta));
                if (counted.applyAsDouble(candidate) < best - delta * Math.max(1e-4, options.scoreTolerance() * 10))
                    stationary = false;
            }
            return new Optimum(point, best, calls[0], stationary);
        }
        double[] absolute(double[] theta, double scale) {
            double[] result = new double[parameters + 1]; int k = 0;
            for (double[][] l : factors(theta)) {
                double[][] covariance = covariance(l, scale);
                for (int r = 0; r < l.length; r++) for (int c = 0; c <= r; c++) result[k++] = covariance[r][c];
            }
            for(int extra=0;extra<pedigreeEffects.size();extra++)
                result[k++]=scale*Math.exp(2*theta[blockParameters+extra]);
            result[k] = scale; return result;
        }
        SparseFiniteDf.Point absolutePoint(double[] absolute) {
            double scale = absolute[parameters];
            if (!(scale > 0)) throw new IllegalArgumentException("nonpositive scale");
            double[] theta = new double[parameters]; int k = 0;
            for (CorrelatedRandomEffectBlock block : blocks) {
                int d = block.effectCount(); double[][] covariance = new double[d][d];
                for (int r = 0; r < d; r++) for (int c = 0; c <= r; c++)
                    covariance[r][c] = covariance[c][r] = absolute[k++] / scale;
                double[][] l = cholesky(covariance); int start = k - d * (d + 1) / 2;
                for (int r = 0; r < d; r++) for (int c = 0; c <= r; c++) theta[start++] = l[r][c];
            }
            for(int extra=0;extra<pedigreeEffects.size();extra++){
                if(!(absolute[k]>0))throw new IllegalArgumentException("pedigree variance must be positive");
                theta[blockParameters+extra]=.5*Math.log(absolute[k++]/scale);
            }
            return point(theta, scale);
        }
        double[] derivativeScales(double[] absolute) {
            double[] scales = new double[absolute.length];
            int start = 0, k = 0;
            for (CorrelatedRandomEffectBlock block : blocks) {
                int d = block.effectCount();
                for (int r = 0; r < d; r++) for (int c = 0; c <= r; c++)
                    scales[k++] = Math.sqrt(absolute[start + r*(r+1)/2 + r]
                        * absolute[start + c*(c+1)/2 + c]);
                start = k;
            }
            for(int extra=0;extra<pedigreeEffects.size();extra++)scales[k]=absolute[k++];
            scales[parameters] = absolute[parameters];
            return scales;
        }
    }

    private static double[][] cholesky(double[][] covariance) {
        int d = covariance.length; double[][] l = new double[d][d];
        for (int r = 0; r < d; r++) for (int c = 0; c <= r; c++) {
            double value = covariance[r][c];
            for (int k = 0; k < c; k++) value -= l[r][k] * l[c][k];
            if (r == c) {
                if (!(value > 0)) throw new IllegalArgumentException("covariance must be interior positive definite");
                l[r][c] = Math.sqrt(value);
            } else l[r][c] = value / l[c][c];
        }
        return l;
    }
    private static double[][] covariance(double[][] l, double scale) {
        double[][] result = new double[l.length][l.length];
        for (int r = 0; r < l.length; r++) for (int c = 0; c < l.length; c++)
            for (int k = 0; k <= Math.min(r, c); k++) result[r][c] += scale * l[r][k] * l[c][k];
        return result;
    }

    public static final class Result {
        private final SparseLinearMixedModelResult fit;
        private final List<CorrelatedRandomEffectEstimates> randomEffects;
        private final List<double[][]> shapes;
        private Result(SparseLinearMixedModelResult fit, List<CorrelatedRandomEffectBlock> blocks,
                List<double[][]> factors) {
            this.fit = fit; shapes = new ArrayList<>();
            List<CorrelatedRandomEffectEstimates> estimates = new ArrayList<>();
            double scale = fit.varianceComponents()[fit.varianceComponents().length-1];
            for (int b = 0; b < blocks.size(); b++) {
                CorrelatedRandomEffectBlock block = blocks.get(b); int d = block.effectCount();
                double[][] l = factors.get(b); shapes.add(covariance(l, 1));
                double[] cov = new double[d * d]; double[][] matrix = covariance(l, scale);
                for (int r = 0; r < d; r++) System.arraycopy(matrix[r], 0, cov, r * d, d);
                RandomEffectEstimates latentResult=fit.randomEffects().get(b);
                double[] latent = latentResult.estimates(), modes = new double[latent.length];
                double[] latentPev=latentResult.predictionErrorCovariances();
                double[] pev=new double[block.groupCount()*d*d];
                for (int g = 0; g < block.groupCount(); g++) for (int r = 0; r < d; r++)
                    for (int c = 0; c <= r; c++) modes[g * d + r] += l[r][c] * latent[g * d + c];
                for(int g=0;g<block.groupCount();g++)for(int r=0;r<d;r++)for(int c=0;c<d;c++)
                    for(int a=0;a<d;a++)for(int j=0;j<d;j++)
                        pev[g*d*d+r*d+c]+=l[r][a]*latentPev[g*d*d+a*d+j]*l[c][j];
                estimates.add(new CorrelatedRandomEffectEstimates(block.name(), block.groupNames(),
                    block.effectNames(), cov, modes,pev));
            }
            randomEffects = List.copyOf(estimates);
        }
        /** Equation diagnostics: random coefficients here are latent u, not b=L u. */
        public SparseLinearMixedModelResult fit() { return fit; }
        public List<CorrelatedRandomEffectEstimates> randomEffects() { return randomEffects; }
        /** First block covariance relative to the residual variance (legacy accessor). */
        public double[][] covarianceShape() {
            double[][] value = shapes.get(0), result = new double[value.length][];
            for (int i = 0; i < value.length; i++) result[i] = value[i].clone();
            return result;
        }
        public int evaluations() { return fit.functionEvaluations(); }
        public boolean converged() { return fit.converged(); }
        public CorrelatedLinearMixedModelResult correlatedFit() {
            return new CorrelatedLinearMixedModelResult(fit.associationStatistics(), fit.fixedEffectCovariance(),
                randomEffects, fit.varianceComponents()[fit.varianceComponents().length-1], fit.fittedValues(), fit.residuals(),
                fit.logLikelihood(), fit.varianceEstimation(), evaluations(), converged(), fit.backend());
        }
        public Result withObservationScale(double[] weights, double[] offsets) {
            return new Result(fit.withObservationScale(weights, offsets), randomEffects, shapes, true);
        }
        private Result(SparseLinearMixedModelResult fit, List<CorrelatedRandomEffectEstimates> effects,
                List<double[][]> shapes, boolean ignored) {
            this.fit = fit; this.randomEffects = effects; this.shapes = shapes;
        }
    }
}
