/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTriangle;
import jdistlib.accelerator.PreparedSparseCholesky;
import jdistlib.accelerator.SparseOrdering;
import jdistlib.matrix.CsrMatrix;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.LaplaceTunableFamily;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;

/**
 * First-order Laplace GLMM using sparse grouped designs and coefficient-space
 * precision matrices. Observation-scale covariance matrices are never formed.
 */
public final class SparseGlmmLaplace {

    private SparseGlmmLaplace() { }

    /** Fits a grouped sparse GLMM with default controls and backend selection. */
    public static GlmmLaplaceResult fit(
            double[] response, double[][] fixedEffects,
            GlmFamily family, List<RandomEffectTerm> randomEffects) {
        if (response == null || fixedEffects == null
                || fixedEffects.length == 0 || fixedEffects[0] == null)
            throw new IllegalArgumentException(
                "response and fixed effects are required");
        return fit(response, MatrixOps.rowMajor(fixedEffects, response.length),
            response.length, fixedEffects[0].length, family, randomEffects,
            null, null, GlmmLaplaceOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    public static GlmmLaplaceResult fit(
            double[] response, double[] fixedEffects, int rows, int columns,
            GlmFamily family, List<RandomEffectTerm> randomEffects,
            double[] priorWeights, double[] offset,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        return fitWithPrecision(response, fixedEffects, rows, columns, family,
            randomEffects, null, priorWeights, offset, options, backendPolicy);
    }

    public static GlmmLaplaceResult fitWithPrecision(
            double[] response, double[] fixedEffects, int rows, int columns,
            GlmFamily family, List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            double[] priorWeights, double[] offset,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        try (Prepared prepared = prepareWithPrecision(rows, family,
                randomEffects, precisionBases, options, backendPolicy)) {
            return prepared.fit(response, fixedEffects, columns,
                priorWeights, offset);
        }
    }

    /** Prepares scan-owned sparse structure and one numerical factor per worker. */
    public static Prepared prepare(
            int rows, GlmFamily family,
            List<RandomEffectTerm> randomEffects,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        return prepareWithPrecision(rows, family, randomEffects, null,
            options, backendPolicy);
    }

    /** Prepares a scan with caller-supplied random-coefficient precisions. */
    public static Prepared prepareWithPrecision(
            int rows, GlmFamily family,
            List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
        validate(rows, family, randomEffects, options, backendPolicy);
        return new Prepared(rows, family, randomEffects,
            precisions(randomEffects, precisionBases), options, backendPolicy);
    }

    /** Reusable scan state. Close it after all worker tasks have completed. */
    public static final class Prepared implements AutoCloseable {
        private final int rows;
        private final GlmFamily family;
        private final List<RandomEffectTerm> terms;
        private final GlmmLaplaceOptions options;
        private final CombinedDesign design;
        private final SparsePattern pattern;
        private final PrecisionData precisionData;
        private final BackendContext context;
        private final ComputeBackend backend;
        private final ConcurrentLinkedQueue<PreparedSparseCholesky> factors =
            new ConcurrentLinkedQueue<>();
        private final ThreadLocal<PreparedSparseCholesky> localFactor;
        private volatile double[] sharedLogVariances;
        private volatile boolean closed;

        private Prepared(
                int rows, GlmFamily family,
                List<RandomEffectTerm> terms,
                List<SparsePrecisionMatrix> precisions,
                GlmmLaplaceOptions options, BackendPolicy backendPolicy) {
            this.rows = rows;
            this.family = family;
            this.terms = List.copyOf(terms);
            this.options = options;
            BackendContext selectedContext =
                BackendContext.select(backendPolicy);
            ComputeBackend selectedBackend = selectedContext.backend();
            CombinedDesign combinedDesign;
            PrecisionData preparedPrecision;
            SparsePattern preparedPattern;
            try {
                combinedDesign = combine(this.terms, rows);
                preparedPrecision = precisionData(
                    this.terms, precisions, combinedDesign, selectedBackend);
                preparedPattern = pattern(
                    combinedDesign, preparedPrecision);
            } catch (RuntimeException | Error failure) {
                selectedContext.close();
                throw failure;
            }
            context = selectedContext;
            backend = selectedBackend;
            design = combinedDesign;
            precisionData = preparedPrecision;
            pattern = preparedPattern;
            double[] initial = initial(options, terms.size());
            sharedLogVariances = initial;
            localFactor = ThreadLocal.withInitial(() -> {
                PreparedSparseCholesky factor = backend.prepareDcsrpotrf(
                    pattern.matrix(new double[rows], design, precisionData,
                        initial), MatrixTriangle.LOWER,
                    SparseOrdering.MINIMUM_DEGREE);
                factors.add(factor);
                return factor;
            });
        }

        public GlmmLaplaceResult fit(
                double[] response, double[] fixedEffects, int columns) {
            return fit(response, fixedEffects, columns, null, null);
        }

        public GlmmLaplaceResult fit(
                double[] response, double[] fixedEffects, int columns,
                double[] priorWeights, double[] offset) {
            // Tunable families hold mutable likelihood parameters. A prepared
            // instance may otherwise fit responses concurrently using local factors.
            if (family instanceof LaplaceTunableFamily) {
                synchronized (family) { return fitMarginal(response,fixedEffects,columns,priorWeights,offset); }
            }
            return fitMarginal(response,fixedEffects,columns,priorWeights,offset);
        }

        private GlmmLaplaceResult fitMarginal(double[] response, double[] fixedEffects, int columns,
                double[] priorWeights, double[] offset) {
            if(closed)throw new IllegalStateException("prepared sparse GLMM is closed");
            MatrixOps.validateModelData(response,fixedEffects,rows,columns);
            double[] weights=weights(priorWeights,rows),offsets=offsets(offset,rows);
            for(int row=0;row<rows;row++)family.validateResponse(response[row],weights[row]);
            double[] scales=new double[columns],fixed=fixedEffects.clone();
            for(int j=0;j<columns;j++) {
                for(int i=0;i<rows;i++)scales[j]=Math.max(scales[j],Math.abs(fixed[i*columns+j]));
                if(scales[j]==0)throw new IllegalArgumentException("fixed design has a zero column");
                for(int i=0;i<rows;i++)fixed[i*columns+j]/=scales[j];
            }
            var tunable=family instanceof LaplaceTunableFamily value ? value : null;
            double[] familyStart=tunable==null?new double[0]:tunable.laplaceParameters();
            int varianceCount=terms.size(),dimension=columns+varianceCount+familyStart.length;
            double[] initial=new double[dimension],lower=new double[dimension],upper=new double[dimension];
            Arrays.fill(lower,Double.NEGATIVE_INFINITY);Arrays.fill(upper,Double.POSITIVE_INFINITY);
            double[] betaStart=new double[columns];initializeIntercept(response,fixed,columns,offsets,betaStart);
            System.arraycopy(betaStart,0,initial,0,columns);
            for(int i=0;i<varianceCount;i++) {
                initial[columns+i]=clamp(sharedLogVariances[i]);
                lower[columns+i]=Math.log(options.minimumVariance());upper[columns+i]=Math.log(options.maximumVariance());
            }
            for(int i=0;i<familyStart.length;i++) {
                initial[columns+varianceCount+i]=familyStart[i];
                lower[columns+varianceCount+i]=tunable.minimumLaplaceParameter(i);
                upper[columns+varianceCount+i]=tunable.maximumLaplaceParameter(i);
            }
            PreparedSparseCholesky factor=localFactor.get();
            LaplaceOptimization.Objective objective=point->{
                if(tunable!=null)tunable.setLaplaceParameters(Arrays.copyOfRange(point,columns+varianceCount,dimension));
                double[] beta=Arrays.copyOf(point,columns),logV=Arrays.copyOfRange(point,columns,columns+varianceCount);
                try {
                    Mode fitted=conditionalMode(response,fixed,columns,weights,offsets,logV,factor,beta);
                    return fitted.converged() ? -fitted.laplaceLogLikelihood() : Double.POSITIVE_INFINITY;
                }catch(IllegalArgumentException | IllegalStateException invalid){return Double.POSITIVE_INFINITY;}
            };
            var optimum=LaplaceOptimization.fit(objective,initial,lower,upper,
                options.maximumOuterIterations(),Math.max(1e-6,10*options.relativeTolerance()),columns,options.initialLogVarianceStep());
            // An extremely small supplied variance can flatten a log-variance
            // score far from the optimum. Compare a finite-scale restart.
            boolean smallStart=false;
            for(int i=0;i<varianceCount;i++)smallStart|=initial[columns+i]<-10;
            if(smallStart) {
                double[] restart=optimum.point().clone();
                for(int i=0;i<varianceCount;i++)restart[columns+i]=Math.max(lower[columns+i],Math.min(upper[columns+i],0));
                var alternative=LaplaceOptimization.fit(objective,restart,lower,upper,
                    options.maximumOuterIterations(),Math.max(1e-6,10*options.relativeTolerance()),columns,options.initialLogVarianceStep());
                if(alternative.value()<optimum.value())optimum=alternative;
            }
            double[] point=optimum.point(),covariance=LaplaceOptimization.fixedCovariance(objective,point,lower,upper,columns);
            if(tunable!=null)tunable.setLaplaceParameters(Arrays.copyOfRange(point,columns+varianceCount,dimension));
            double[] logV=Arrays.copyOfRange(point,columns,columns+varianceCount);
            Mode fitted=conditionalMode(response,fixed,columns,weights,offsets,logV,factor,Arrays.copyOf(point,columns));
            double[] beta=fitted.beta().clone();
            for(int i=0;i<columns;i++) {
                beta[i]/=scales[i];
                for(int j=0;j<columns;j++)covariance[i*columns+j]=
                    unscaleCovariance(covariance[i*columns+j],scales[i],scales[j]);
            }
            boolean converged=optimum.converged()&&fitted.converged();
            for(double value:covariance)converged&=Double.isFinite(value);
            if(!converged)Arrays.fill(covariance,Double.NaN);
            return result(new Mode(beta,fitted.random(),covariance,fitted.linear(),fitted.means(),
                fitted.laplaceLogLikelihood(),fitted.iterations(),fitted.converged()),columns,logV,optimum.iterations(),converged);
        }


        /** Sets deterministic variance starts shared by subsequent workers. */
        public void warmStart(double... varianceComponents) {
            if (closed) throw new IllegalStateException(
                "prepared sparse GLMM is closed");
            if (varianceComponents == null
                    || varianceComponents.length != terms.size())
                throw new IllegalArgumentException(
                    "one variance is required per random-effect term");
            double[] values = new double[varianceComponents.length];
            for (int index = 0; index < values.length; index++) {
                double value = varianceComponents[index];
                if (!(value > 0.0) || !Double.isFinite(value))
                    throw new IllegalArgumentException(
                        "warm-start variances must be finite and positive");
                values[index] = Math.log(value);
            }
            sharedLogVariances = values;
        }

        private double clamp(double value) {
            return Math.max(Math.log(options.minimumVariance()),
                Math.min(Math.log(options.maximumVariance()), value));
        }

        /** Random modes at FIXED beta; beta is optimized in the marginal objective. */
        private Mode conditionalMode(double[] response,double[] fixed,int fixedColumns,
                double[] priorWeights,double[] offsets,double[] logVariances,
                PreparedSparseCholesky factor,double[] beta) {
            double[] random=new double[design.columns()],means=new double[rows],linear=new double[rows];
            double[] weights=new double[rows];boolean converged=false;int iteration=0;
            double previous=penalized(response,fixed,fixedColumns,beta,random,priorWeights,offsets,logVariances);
            for(;iteration<options.maximumModeIterations();iteration++) {
                double[] z=working(response,fixed,fixedColumns,beta,random,priorWeights,offsets,means,linear,weights);
                factor.refactor(pattern.matrix(weights,design,precisionData,logVariances));
                double[] right=new double[random.length];
                for(int row=0;row<rows;row++) {
                    double value=weights[row]*(z[row]-fixedValue(fixed,row,fixedColumns,beta));
                    for(int j=design.rowStarts()[row];j<design.rowStarts()[row+1];j++)right[design.columnIndices()[j]]+=design.values()[j]*value;
                }
                factor.solveInPlace(right,1);
                double change=relativeChange(random,right);
                if(change<=Math.min(options.relativeTolerance(),1e-9)) {random=right;converged=true;break;}
                boolean accepted=false;double[] candidate=new double[random.length];
                for(double alpha=1;alpha>=1e-10;alpha*=.5) {
                    for(int i=0;i<random.length;i++)candidate[i]=random[i]+alpha*(right[i]-random[i]);
                    double value=penalized(response,fixed,fixedColumns,beta,candidate,priorWeights,offsets,logVariances);
                    if(Double.isFinite(value)&&value>=previous-1e-12*(1+Math.abs(previous))) {
                        previous=value;random=candidate.clone();accepted=true;break;
                    }
                }
                if(!accepted)break;
            }
            double conditional=penalized(response,fixed,fixedColumns,beta,random,priorWeights,offsets,logVariances);
            for(int row=0;row<rows;row++) {
                double eta=offsets[row]+fixedValue(fixed,row,fixedColumns,beta)+design.rowProduct(row,random);
                linear[row]=eta;means[row]=family.inverseLink(eta);
                weights[row]=observedWeight(response[row],eta,means[row],priorWeights[row]);
                if(!Double.isFinite(weights[row]))throw new IllegalArgumentException("nonfinite observed family curvature");
            }
            // The Laplace determinant is the observed random Hessian at the
            // accepted mode, not an earlier Fisher-scoring working matrix.
            factor.refactor(pattern.matrix(weights,design,precisionData,logVariances));
            double likelihood=conditional+.5*precisionData.logDeterminant(logVariances)-.5*factor.logDeterminant();
            return new Mode(beta,random,null,linear,means,likelihood,iteration+1,converged&&Double.isFinite(likelihood));
        }

        private double observedWeight(double response,double eta,double mean,double weight) {
            if(family instanceof LaplaceFamilyDerivatives exact)
                return exact.linearPredictorInformation(response,eta,weight);
            if(family==org.jlinalg.glm.GlmFamilies.binomial()||family==org.jlinalg.glm.GlmFamilies.poisson())
                return family.workingWeight(response,eta,mean,weight);
            double h=1e-5*Math.max(1,Math.abs(eta));
            return -(familyScore(response,eta+h,weight)-familyScore(response,eta-h,weight))/(2*h);
        }
        private double familyScore(double response,double eta,double weight) {
            if(family instanceof LaplaceFamilyDerivatives exact)return exact.linearPredictorScore(response,eta,weight);
            double mean=family.inverseLink(eta);
            return family.workingWeight(response,eta,mean,weight)
                *(family.workingResponse(response,eta,mean,weight,0)-eta);
        }
        private double penalized(double[] response,double[] fixed,int fixedColumns,double[] beta,
                double[] random,double[] weights,double[] offsets,double[] logVariances) {
            double value=-.5*precisionData.quadratic(random,logVariances);
            for(int row=0;row<rows;row++) {
                double eta=offsets[row]+fixedValue(fixed,row,fixedColumns,beta)+design.rowProduct(row,random);
                value+=family.logLikelihoodAtPredictor(response[row],eta,family.inverseLink(eta),weights[row],1);
            }
            return value;
        }


        private double[] working(
                double[] response, double[] fixed, int fixedColumns,
                double[] beta, double[] random, double[] priorWeights,
                double[] offsets, double[] means, double[] linear,
                double[] workingWeights) {
            double[] workingResponse = new double[rows];
            for (int row = 0; row < rows; row++) {
                double eta = offsets[row]
                    + fixedValue(fixed, row, fixedColumns, beta)
                    + design.rowProduct(row, random);
                linear[row] = eta;
                means[row] = family.inverseLink(eta);
            }
            for (int row = 0; row < rows; row++) {
                double eta = linear[row];
                double mean = means[row];
                double weight = family.workingWeight(response[row], eta,
                    mean, priorWeights[row]);
                if (!(weight > 0.0) || !Double.isFinite(weight))
                    throw new IllegalArgumentException("family produced unrepresentable Laplace working precision");
                workingWeights[row] = weight;
                workingResponse[row] = family.workingResponse(response[row],
                    eta, mean, priorWeights[row], offsets[row]);
                if (!Double.isFinite(workingResponse[row]))
                    throw new IllegalArgumentException("family produced nonfinite Laplace working response");
            }
            return workingResponse;
        }


        private GlmmLaplaceResult result(
                Mode mode, int fixedColumns, double[] logVariances,
                int outerIterations, boolean converged) {
            double[] standardErrors = new double[fixedColumns];
            for (int column = 0; column < fixedColumns; column++)
                standardErrors[column] = Math.sqrt(Math.max(0.0,
                    mode.fixedCovariance()[column * fixedColumns + column]));
            double[] variances = new double[logVariances.length];
            for (int index = 0; index < variances.length; index++)
                variances[index] = Math.exp(logVariances[index]);
            List<String> names = new ArrayList<>(terms.size());
            Map<String, double[]> predictors = new LinkedHashMap<>();
            Map<String, double[]> coefficientsByTerm = new LinkedHashMap<>();
            for (int term = 0; term < terms.size(); term++) {
                RandomEffectTerm value = terms.get(term);
                names.add(value.name());
                double[] coefficients = Arrays.copyOfRange(mode.random(),
                    design.termStarts()[term],
                    design.termStarts()[term] + value.coefficients());
                coefficientsByTerm.put(value.name(), coefficients);
                predictors.put(value.name(), multiply(value, coefficients));
            }
            return new GlmmLaplaceResult(family.name(), names, variances,
                AssociationStatistics.normal(mode.beta(), standardErrors),
                mode.fixedCovariance(), coefficientsByTerm, predictors,
                mode.linear(), mode.means(),
                mode.laplaceLogLikelihood(), outerIterations,
                mode.iterations(), converged);
        }

        private void initializeIntercept(
                double[] response, double[] fixed, int columns,
                double[] offsets, double[] beta) {
            int intercept = -1;
            for (int column = 0; column < columns; column++) {
                boolean allOne = true;
                for (int row = 0; row < rows; row++)
                    if (Math.abs(fixed[row * columns + column] - 1.0)
                            > 1e-12) {
                        allOne = false;
                        break;
                    }
                if (allOne) {
                    intercept = column;
                    break;
                }
            }
            if (intercept < 0) return;
            double total = 0.0;
            for (double value : response) total += value;
            double mean = Math.max(1e-7,
                Math.min(1.0 - 1e-7, total / rows));
            double offsetMean = 0.0;
            for (double value : offsets) offsetMean += value;
            beta[intercept] = family.link(mean) - offsetMean / rows;
        }

        @Override public void close() {
            if (closed) return;
            closed = true;
            localFactor.remove();
            for (PreparedSparseCholesky factor : factors) factor.close();
            factors.clear();
            context.close();
        }
    }

    private static SparsePattern pattern(
            CombinedDesign design, PrecisionData precision) {
        int dimension = design.columns();
        TreeMap<Long, Boolean> keys = new TreeMap<>();
        for (int row = 0; row < design.rows(); row++)
            for (int left = design.rowStarts()[row];
                    left < design.rowStarts()[row + 1]; left++)
                for (int right = design.rowStarts()[row];
                        right <= left; right++) {
                    int r = Math.max(design.columnIndices()[left],
                        design.columnIndices()[right]);
                    int c = Math.min(design.columnIndices()[left],
                        design.columnIndices()[right]);
                    keys.put((long) r * dimension + c, Boolean.TRUE);
                }
        for (int index = 0; index < precision.rows().length; index++)
            keys.put((long) precision.rows()[index] * dimension
                + precision.columns()[index], Boolean.TRUE);
        int[] rowStarts = new int[dimension + 1];
        int[] columns = new int[keys.size()];
        int position = 0;
        int currentRow = 0;
        for (long key : keys.keySet()) {
            int row = (int) (key / dimension);
            while (currentRow < row) rowStarts[++currentRow] = position;
            columns[position++] = (int) (key % dimension);
        }
        while (currentRow < dimension) rowStarts[++currentRow] = position;
        return new SparsePattern(dimension, rowStarts, columns);
    }

    private static PrecisionData precisionData(
            List<RandomEffectTerm> terms,
            List<SparsePrecisionMatrix> precisions,
            CombinedDesign design, ComputeBackend backend) {
        List<Integer> rows = new ArrayList<>();
        List<Integer> columns = new ArrayList<>();
        List<Double> values = new ArrayList<>();
        List<Integer> termIndices = new ArrayList<>();
        double[] logDeterminants = new double[terms.size()];
        CovarianceIdentification identification = new CovarianceIdentification(design.rows(), terms.size());
        for (int term = 0; term < terms.size(); term++) {
            SparsePrecisionMatrix precision = precisions.get(term);
            int[] starts = precision.rowStarts();
            int[] cols = precision.columnIndices();
            double[] numeric = precision.values();
            int lowerCount = 0;
            for (int row = 0; row < precision.dimension(); row++)
                for (int index = starts[row]; index < starts[row + 1]; index++)
                    if (cols[index] <= row) lowerCount++;
            int[] lowerRows = new int[precision.dimension() + 1];
            int[] lowerColumns = new int[lowerCount];
            double[] lowerValues = new double[lowerCount];
            int lowerPosition = 0;
            for (int row = 0; row < precision.dimension(); row++) {
                lowerRows[row] = lowerPosition + 1;
                for (int index = starts[row]; index < starts[row + 1]; index++) {
                    int column = cols[index];
                    if (column <= row) {
                        lowerColumns[lowerPosition] = column + 1;
                        lowerValues[lowerPosition++] = numeric[index];
                        rows.add(design.termStarts()[term] + row);
                        columns.add(design.termStarts()[term] + column);
                        values.add(numeric[index]);
                        termIndices.add(term);
                    }
                }
            }
            lowerRows[precision.dimension()] = lowerPosition + 1;
            CsrMatrix lower = new CsrMatrix(precision.dimension(),
                precision.dimension(), lowerValues, lowerColumns, lowerRows);
            var factor = backend.dcsrpotrf(lower,
                MatrixTriangle.LOWER,
                SparseOrdering.MINIMUM_DEGREE);
            logDeterminants[term] = factor.logDeterminant();
            identification.add(terms.get(term), factor::solve);
        }
        return new PrecisionData(toIntArray(rows), toIntArray(columns),
            toDoubleArray(values), toIntArray(termIndices), logDeterminants,
            terms.stream().mapToInt(RandomEffectTerm::coefficients).toArray());
    }

    private static CombinedDesign combine(
            List<RandomEffectTerm> terms, int rows) {
        int[] termStarts = new int[terms.size()];
        int columns = 0;
        int nonzeros = 0;
        for (int term = 0; term < terms.size(); term++) {
            termStarts[term] = columns;
            columns += terms.get(term).coefficients();
            nonzeros += terms.get(term).nonzeroCount();
        }
        int[] rowStarts = new int[rows + 1];
        int[] columnIndices = new int[nonzeros];
        double[] values = new double[nonzeros];
        TermData[] data = new TermData[terms.size()];
        for (int term = 0; term < terms.size(); term++)
            data[term] = termData(terms.get(term));
        int position = 0;
        for (int row = 0; row < rows; row++) {
            rowStarts[row] = position;
            for (int term = 0; term < terms.size(); term++) {
                TermData value = data[term];
                for (int index = value.rowStarts()[row];
                        index < value.rowStarts()[row + 1]; index++) {
                    double entry = value.values()[index];
                    if (entry != 0.0) {
                        columnIndices[position] = termStarts[term]
                            + value.columnIndices()[index];
                        values[position++] = entry;
                    }
                }
            }
        }
        rowStarts[rows] = position;
        return new CombinedDesign(rows, columns, rowStarts,
            Arrays.copyOf(columnIndices, position),
            Arrays.copyOf(values, position), termStarts);
    }

    static double unscaleCovariance(double value, double first, double second) {
        if (value == 0.0 || !Double.isFinite(value)) return value;
        int a=Math.getExponent(value),b=Math.getExponent(first),c=Math.getExponent(second);
        return Math.scalb(Math.scalb(value,-a)/Math.scalb(first,-b)/Math.scalb(second,-c),a-b-c);
    }

    private static TermData termData(RandomEffectTerm term) {
        if (term.sparse())
            return new TermData(term.rowPointers(), term.columnIndices(),
                term.sparseValues());
        double[] dense = term.design();
        int count = 0;
        for (double value : dense) if (value != 0.0) count++;
        int[] starts = new int[term.observations() + 1];
        int[] columns = new int[count];
        double[] values = new double[count];
        int position = 0;
        for (int row = 0; row < term.observations(); row++) {
            starts[row] = position;
            for (int column = 0; column < term.coefficients(); column++) {
                double value = dense[row * term.coefficients() + column];
                if (value != 0.0) {
                    columns[position] = column;
                    values[position++] = value;
                }
            }
        }
        starts[term.observations()] = position;
        return new TermData(starts, columns, values);
    }

    private static double[] multiply(
            RandomEffectTerm term, double[] coefficients) {
        double[] result = new double[term.observations()];
        if (term.sparse()) {
            int[] starts = term.rowPointers();
            int[] columns = term.columnIndices();
            double[] values = term.sparseValues();
            for (int row = 0; row < result.length; row++)
                for (int index = starts[row]; index < starts[row + 1]; index++)
                    result[row] += values[index] * coefficients[columns[index]];
        } else {
            double[] design = term.design();
            for (int row = 0; row < result.length; row++)
                for (int column = 0; column < coefficients.length; column++)
                    result[row] += design[row * coefficients.length + column]
                        * coefficients[column];
        }
        return result;
    }

    private static List<SparsePrecisionMatrix> precisions(
            List<RandomEffectTerm> terms,
            List<SparsePrecisionMatrix> supplied) {
        if (supplied == null) {
            List<SparsePrecisionMatrix> result = new ArrayList<>(terms.size());
            for (RandomEffectTerm term : terms)
                result.add(SparsePrecisionMatrix.identity(term.coefficients()));
            return List.copyOf(result);
        }
        if (supplied.size() != terms.size())
            throw new IllegalArgumentException(
                "one precision matrix is required per random-effect term");
        for (int index = 0; index < terms.size(); index++)
            if (supplied.get(index) == null
                    || supplied.get(index).dimension()
                        != terms.get(index).coefficients())
                throw new IllegalArgumentException(
                    "precision dimensions must match random coefficients");
        return List.copyOf(supplied);
    }

    private static void validate(
            int rows, GlmFamily family,
            List<RandomEffectTerm> terms,
            GlmmLaplaceOptions options, BackendPolicy policy) {
        if (rows < 1 || family == null || terms == null || terms.isEmpty()
                || options == null || policy == null)
            throw new IllegalArgumentException(
                "rows, family, random effects, controls, and backend are required");
        if(!family.fixedDispersion())throw new IllegalArgumentException(
            "Laplace requires a fixed-dispersion or internally parameterized likelihood; use Gaussian REML for estimated residual dispersion");
        java.util.HashSet<String> names = new java.util.HashSet<>();
        for (RandomEffectTerm term : terms)
            if (term == null || term.observations() != rows
                    || !names.add(term.name()))
                throw new IllegalArgumentException(
                    "random effects need unique names and matching rows");
        double[] supplied = options.initialVariances();
        if (supplied != null && supplied.length != terms.size())
            throw new IllegalArgumentException(
                "one initial variance is required per random-effect term");
    }

    private static double[] initial(
            GlmmLaplaceOptions options, int components) {
        double[] supplied = options.initialVariances();
        double[] result = new double[components];
        if (supplied == null) return result;
        for (int index = 0; index < components; index++)
            result[index] = Math.log(supplied[index]);
        return result;
    }

    private static double[] weights(double[] supplied, int rows) {
        if (supplied == null) {
            double[] result = new double[rows];
            Arrays.fill(result, 1.0);
            return result;
        }
        if (supplied.length != rows)
            throw new IllegalArgumentException(
                "prior-weight length must match observations");
        double[] result = supplied.clone();
        for (double value : result)
            if (!(value > 0.0) || !Double.isFinite(value))
                throw new IllegalArgumentException(
                    "prior weights must be finite and positive");
        return result;
    }

    private static double[] offsets(double[] supplied, int rows) {
        if (supplied == null) return new double[rows];
        if (supplied.length != rows)
            throw new IllegalArgumentException(
                "offset length must match observations");
        return MatrixOps.finiteCopy(supplied, "offset");
    }

    private static double fixedValue(
            double[] fixed, int row, int columns, double[] beta) {
        double value = 0.0;
        for (int column = 0; column < columns; column++)
            value += fixed[row * columns + column] * beta[column];
        return value;
    }

    private static double relativeChange(double[] previous, double[] next) {
        double maximum = 0.0;
        for (int index = 0; index < previous.length; index++)
            maximum = Math.max(maximum, Math.abs(next[index] - previous[index])
                / (1.0 + Math.abs(previous[index])));
        return maximum;
    }

    private static int[] toIntArray(List<Integer> values) {
        int[] result = new int[values.size()];
        for (int index = 0; index < result.length; index++)
            result[index] = values.get(index);
        return result;
    }

    private static double[] toDoubleArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int index = 0; index < result.length; index++)
            result[index] = values.get(index);
        return result;
    }

    private record TermData(
            int[] rowStarts, int[] columnIndices, double[] values) { }

    private record CombinedDesign(
            int rows, int columns, int[] rowStarts,
            int[] columnIndices, double[] values, int[] termStarts) {
        double rowProduct(int row, double[] coefficients) {
            double result = 0.0;
            for (int index = rowStarts[row]; index < rowStarts[row + 1]; index++)
                result += values[index] * coefficients[columnIndices[index]];
            return result;
        }
    }

    private record SparsePattern(
            int dimension, int[] rowStarts, int[] columnIndices) {
        CsrMatrix matrix(
                double[] weights, CombinedDesign design,
                PrecisionData precision, double[] logVariances) {
            double[] numeric = new double[columnIndices.length];
            for (int observation = 0;
                    observation < design.rows(); observation++) {
                double weight = weights[observation];
                for (int left = design.rowStarts()[observation];
                        left < design.rowStarts()[observation + 1]; left++)
                    for (int right = design.rowStarts()[observation];
                            right <= left; right++) {
                        int row = Math.max(design.columnIndices()[left],
                            design.columnIndices()[right]);
                        int column = Math.min(design.columnIndices()[left],
                            design.columnIndices()[right]);
                        int position = Arrays.binarySearch(columnIndices,
                            rowStarts[row], rowStarts[row + 1], column);
                        numeric[position] += weight * design.values()[left]
                            * design.values()[right];
                    }
            }
            for (int index = 0; index < precision.values().length; index++) {
                int row = precision.rows()[index];
                int position = Arrays.binarySearch(columnIndices,
                    rowStarts[row], rowStarts[row + 1],
                    precision.columns()[index]);
                numeric[position] += precision.values()[index]
                    / Math.exp(logVariances[precision.terms()[index]]);
            }
            int[] csrRows = rowStarts.clone();
            int[] csrColumns = columnIndices.clone();
            for (int index = 0; index < csrRows.length; index++) csrRows[index]++;
            for (int index = 0; index < csrColumns.length; index++)
                csrColumns[index]++;
            return new CsrMatrix(dimension, dimension, numeric,
                csrColumns, csrRows);
        }
    }

    private record PrecisionData(
            int[] rows, int[] columns, double[] values, int[] terms,
            double[] logDeterminants, int[] dimensions) {
        double logDeterminant(double[] logVariances) {
            double result = 0.0;
            for (int term = 0; term < logDeterminants.length; term++)
                result += logDeterminants[term]
                    - dimensions[term] * logVariances[term];
            return result;
        }

        double quadratic(double[] random, double[] logVariances) {
            double result = 0.0;
            for (int index = 0; index < values.length; index++) {
                int row = rows[index];
                int column = columns[index];
                double contribution = values[index] * random[row]
                    * random[column]
                    / Math.exp(logVariances[terms[index]]);
                result += row == column ? contribution : 2.0 * contribution;
            }
            return result;
        }
    }

    private record Mode(
            double[] beta, double[] random, double[] fixedCovariance,
            double[] linear, double[] means, double laplaceLogLikelihood,
            int iterations, boolean converged) { }
}
