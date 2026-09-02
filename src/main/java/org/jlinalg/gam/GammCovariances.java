/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;
import org.jlinalg.genetics.GenomicRelationshipMatrix;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.reml.VarianceComponent;

/** Covariance adapters for GAMM pedigree, GRM, and grouped effects. */
public final class GammCovariances {
    private GammCovariances() { }

    /** Creates an observation-aligned additive pedigree covariance. */
    public static VarianceComponent pedigree(
            String name,
            Pedigree pedigree,
            List<String> observationIndividualIds) {
        if (pedigree == null || observationIndividualIds == null
                || observationIndividualIds.isEmpty()) {
            throw new IllegalArgumentException(
                "pedigree and observation individual IDs are required");
        }
        int observations = observationIndividualIds.size();
        double[] covariance = new double[observations * observations];
        for (int row = 0; row < observations; row++) {
            String first = observationIndividualIds.get(row);
            if (first == null) {
                throw new IllegalArgumentException(
                    "observation individual IDs must not be null");
            }
            for (int column = 0; column <= row; column++) {
                String second = observationIndividualIds.get(column);
                if (second == null) {
                    throw new IllegalArgumentException(
                        "observation individual IDs must not be null");
                }
                double value = pedigree.relationship(first, second);
                covariance[row * observations + column] = value;
                covariance[column * observations + row] = value;
            }
        }
        return new VarianceComponent(name, observations, covariance);
    }

    /** Creates an observation-aligned GRM or cryptic-relatedness covariance. */
    public static VarianceComponent genomicRelationship(
            String name,
            GenomicRelationshipMatrix relationship,
            List<String> observationSampleIds) {
        if (relationship == null) {
            throw new IllegalArgumentException(
                "genomic relationship matrix is required");
        }
        return relationship.varianceComponent(name, observationSampleIds);
    }

    /** Creates a grouped random-intercept covariance. */
    public static VarianceComponent randomIntercept(
            String name, List<?> groups) {
        return VarianceComponent.randomIntercept(name, groups);
    }

    /** Creates an independent grouped random-slope covariance. */
    public static VarianceComponent randomSlope(
            String name, List<?> groups, double[] covariate) {
        return VarianceComponent.randomSlope(name, groups, covariate);
    }
}
