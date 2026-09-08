/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.SymmetricEigenDecomposition;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.reml.VarianceComponent;

/** Laplace GLMM with dense positive-semidefinite covariance bases.
 * Bases are factorized once; the common coefficient-space engine optimizes the
 * complete marginal objective over fixed effects and covariance scales. */
public final class GlmmLaplace {
    private GlmmLaplace() { }
    public static GlmmLaplaceResult fit(double[] response,double[][] fixedEffects,
            GlmFamily family,List<VarianceComponent> components) {
        if(response==null||fixedEffects==null||fixedEffects.length==0||fixedEffects[0]==null)
            throw new IllegalArgumentException("response and fixed design are required");
        return fit(response,MatrixOps.rowMajor(fixedEffects,response.length),response.length,
            fixedEffects[0].length,family,components,null,null,GlmmLaplaceOptions.defaults(),BackendPolicy.PREFERRED);
    }
    public static GlmmLaplaceResult fit(double[] response,double[] fixedEffects,int observations,int fixedColumns,
            GlmFamily family,List<VarianceComponent> components,double[] priorWeights,double[] offset,
            GlmmLaplaceOptions options,BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response,fixedEffects,observations,fixedColumns);
        if(family==null||components==null||components.isEmpty()||options==null||backendPolicy==null)
            throw new IllegalArgumentException("family, covariance components, controls, and backend are required");
        List<RandomEffectTerm> terms=new ArrayList<>();
        try(BackendContext context=BackendContext.select(backendPolicy)) {
            List<ComponentFactor> decomposed=factors(components,observations,context.backend());
            for(int i=0;i<components.size();i++) {
                ComponentFactor factor=decomposed.get(i);List<String> labels=new ArrayList<>();
                for(int j=0;j<factor.columns();j++)labels.add("eigen"+j);
                terms.add(RandomEffectTerm.of(components.get(i).name(),factor.values(),observations,factor.columns(),labels));
            }
        }
        return SparseGlmmLaplace.fit(response,fixedEffects,observations,fixedColumns,family,terms,
            priorWeights,offset,options,backendPolicy);
    }

    private static List<ComponentFactor> factors(
            List<VarianceComponent> components,
            int observations,
            ComputeBackend backend) {
        List<ComponentFactor> result = new ArrayList<>(components.size());
        java.util.HashSet<String> names = new java.util.HashSet<>();
        for (VarianceComponent component : components) {
            if (component == null || component.dimension() != observations
                    || !names.add(component.name())) {
                throw new IllegalArgumentException(
                    "components must have unique names and matching observations");
            }
            SymmetricEigenDecomposition eigen = backend.dsyev(
                component.covariance(), observations);
            double[] values = eigen.eigenvalues();
            double[] vectors = eigen.eigenvectors();
            double maximum = Arrays.stream(values).map(Math::abs).max().orElse(0.0);
            double tolerance = 1e-10 * maximum;
            int rank = 0;
            for (double value : values) {
                if (value < -tolerance) {
                    throw new IllegalArgumentException(
                        "covariance component is not positive semidefinite: "
                            + component.name());
                }
                if (value > tolerance) rank++;
            }
            if (rank == 0) {
                throw new IllegalArgumentException(
                    "covariance component has zero rank: " + component.name());
            }
            double[] factor = new double[observations * rank];
            int destination = 0;
            for (int source = 0; source < observations; source++) {
                if (values[source] <= tolerance) continue;
                double scale = Math.sqrt(values[source]);
                for (int row = 0; row < observations; row++) {
                    factor[row * rank + destination] =
                        vectors[row * observations + source] * scale;
                }
                destination++;
            }
            result.add(new ComponentFactor(factor, observations, rank));
        }
        return List.copyOf(result);
    }


    private record ComponentFactor(double[] values,int observations,int columns) { }
}
