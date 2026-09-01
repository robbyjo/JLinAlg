/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.pipeline.VariantRecord;
import org.jlinalg.reml.VarianceComponent;

/**
 * Dense additive genomic relationship matrix aligned to immutable sample IDs.
 *
 * <p>The genotype builder mean-imputes missing calls and averages products of
 * per-variant genotypes standardized by {@code sqrt(2 p (1-p))}. Callers
 * should normally supply LD-pruned common variants for cryptic-relatedness
 * estimation.</p>
 */
public final class GenomicRelationshipMatrix {
    private final List<String> sampleIds;
    private final Map<String, Integer> sampleIndex;
    private final double[] relationship;
    private final int variantsConsidered;
    private final int variantsUsed;
    private final BackendProvenance computationBackend;

    /** Wraps a caller-supplied row-major relationship matrix. */
    public GenomicRelationshipMatrix(
            List<String> sampleIds, double[] relationship) {
        this(sampleIds, relationship, 0, 0, null);
    }

    private GenomicRelationshipMatrix(
            List<String> sampleIds, double[] relationship,
            int variantsConsidered, int variantsUsed,
            BackendProvenance computationBackend) {
        if (sampleIds == null || sampleIds.isEmpty())
            throw new IllegalArgumentException("sample IDs are required");
        this.sampleIds = List.copyOf(sampleIds);
        sampleIndex = index(this.sampleIds);
        int samples = this.sampleIds.size();
        if (relationship == null
                || relationship.length != samples * samples)
            throw new IllegalArgumentException(
                "relationship matrix dimensions are invalid");
        this.relationship = MatrixOps.finiteCopy(
            relationship, "genomic relationship matrix");
        validateRelationship(this.relationship, samples);
        this.variantsConsidered = variantsConsidered;
        this.variantsUsed = variantsUsed;
        this.computationBackend = computationBackend;
    }

    /**
     * Builds a GRM from variant-by-sample alternate-allele dosages.
     * Non-finite missing values must be represented by {@code NaN}.
     */
    public static GenomicRelationshipMatrix fromVariantDosages(
            double[][] variantDosages,
            List<String> sampleIds,
            GenomicRelationshipOptions options,
            BackendPolicy backendPolicy) {
        if (variantDosages == null || variantDosages.length == 0
                || sampleIds == null || sampleIds.isEmpty()
                || options == null || backendPolicy == null)
            throw new IllegalArgumentException(
                "variant dosages, sample IDs, options, and backend are required");
        int samples = sampleIds.size();
        List<double[]> standardized = new ArrayList<>();
        for (double[] dosage : variantDosages) {
            if (dosage == null || dosage.length != samples)
                throw new IllegalArgumentException(
                    "every variant must contain one dosage per sample");
            int called = 0;
            double sum = 0;
            for (double value : dosage) {
                if (Double.isNaN(value)) continue;
                if (!Double.isFinite(value) || value < 0 || value > 2)
                    throw new IllegalArgumentException(
                        "called dosages must be finite and lie in [0,2]");
                called++;
                sum += value;
            }
            if (called == 0
                    || called / (double) samples < options.minimumCallRate())
                continue;
            double mean = sum / called;
            double frequency = mean / 2;
            double maf = Math.min(frequency, 1 - frequency);
            double alleleVariance = 2 * frequency * (1 - frequency);
            if (maf < options.minimumMinorAlleleFrequency()
                    || !(alleleVariance > 1e-15))
                continue;
            double scale = Math.sqrt(alleleVariance);
            double[] values = new double[samples];
            for (int sample = 0; sample < samples; sample++)
                values[sample] = Double.isNaN(dosage[sample])
                    ? 0 : (dosage[sample] - mean) / scale;
            standardized.add(values);
        }
        if (standardized.isEmpty())
            throw new IllegalArgumentException(
                "no variants remain for GRM construction");
        int variants = standardized.size();
        double[] sampleByVariant = new double[samples * variants];
        for (int variant = 0; variant < variants; variant++)
            for (int sample = 0; sample < samples; sample++)
                sampleByVariant[sample * variants + variant] =
                    standardized.get(variant)[sample];
        double[] matrix = new double[samples * samples];
        BackendProvenance provenance;
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            context.backend().dgemm(
                MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
                samples, samples, variants, 1.0 / variants,
                sampleByVariant, sampleByVariant, 0, matrix);
            provenance = context.provenance();
        }
        for (int row = 0; row < samples; row++)
            for (int column = 0; column < row; column++) {
                double value = 0.5 * (matrix[row * samples + column]
                    + matrix[column * samples + row]);
                matrix[row * samples + column] = value;
                matrix[column * samples + row] = value;
            }
        return new GenomicRelationshipMatrix(sampleIds, matrix,
            variantDosages.length, variants, provenance);
    }

    /** Convenience overload for already aligned pipeline variant records. */
    public static GenomicRelationshipMatrix fromVariants(
            List<VariantRecord> variants,
            List<String> sampleIds,
            GenomicRelationshipOptions options,
            BackendPolicy backendPolicy) {
        if (variants == null || variants.isEmpty())
            throw new IllegalArgumentException("variants are required");
        double[][] dosages = new double[variants.size()][];
        for (int index = 0; index < dosages.length; index++) {
            VariantRecord variant = variants.get(index);
            if (variant == null)
                throw new IllegalArgumentException(
                    "variants must not contain null records");
            dosages[index] = variant.dosages();
        }
        return fromVariantDosages(
            dosages, sampleIds, options, backendPolicy);
    }

    public int samples() { return sampleIds.size(); }
    public List<String> sampleIds() { return sampleIds; }
    public double[] relationshipMatrix() { return relationship.clone(); }
    public int variantsConsidered() { return variantsConsidered; }
    public int variantsUsed() { return variantsUsed; }
    public int variantsExcluded() {
        return variantsConsidered - variantsUsed;
    }
    public Optional<BackendProvenance> computationBackend() {
        return Optional.ofNullable(computationBackend);
    }

    public int indexOf(String sampleId) {
        Integer value = sampleIndex.get(sampleId);
        if (value == null)
            throw new IllegalArgumentException(
                "unknown GRM sample ID: " + sampleId);
        return value;
    }

    public double relationship(String firstId, String secondId) {
        return relationship[indexOf(firstId) * samples() + indexOf(secondId)];
    }

    /** Kinship is one half of the additive relationship coefficient. */
    public double kinshipCoefficient(String firstId, String secondId) {
        return 0.5 * relationship(firstId, secondId);
    }

    /** Returns positive off-diagonal relationships at or above a threshold. */
    public List<RelatednessPair> relatedPairs(double minimumRelationship) {
        if (!Double.isFinite(minimumRelationship)
                || minimumRelationship < 0)
            throw new IllegalArgumentException(
                "minimum relationship must be finite and nonnegative");
        List<RelatednessPair> result = new ArrayList<>();
        for (int row = 0; row < samples(); row++)
            for (int column = 0; column < row; column++) {
                double value = relationship[row * samples() + column];
                if (value >= minimumRelationship)
                    result.add(new RelatednessPair(
                        sampleIds.get(column), sampleIds.get(row),
                        value, 0.5 * value));
            }
        result.sort((left, right) ->
            Double.compare(right.relationship(), left.relationship()));
        return List.copyOf(result);
    }

    /** Returns pairs at or above a conventional kinship-coefficient cutoff. */
    public List<RelatednessPair> relatedPairsByKinship(
            double minimumKinshipCoefficient) {
        if (!Double.isFinite(minimumKinshipCoefficient)
                || minimumKinshipCoefficient < 0)
            throw new IllegalArgumentException(
                "minimum kinship coefficient must be finite and nonnegative");
        return relatedPairs(2 * minimumKinshipCoefficient);
    }

    /** Adapts the GRM to REML, GLMM, GWAS, and related set-test null models. */
    public VarianceComponent varianceComponent(String name) {
        return new VarianceComponent(name, samples(), relationship);
    }

    /**
     * Expands {@code Z K Z'} for repeated or reordered observation IDs.
     */
    public VarianceComponent varianceComponent(
            String name, List<String> observationSampleIds) {
        if (observationSampleIds == null
                || observationSampleIds.isEmpty())
            throw new IllegalArgumentException(
                "observation sample IDs are required");
        int observations = observationSampleIds.size();
        int[] aligned = new int[observations];
        for (int row = 0; row < observations; row++) {
            String id = observationSampleIds.get(row);
            if (id == null)
                throw new IllegalArgumentException(
                    "observation sample IDs must not be null");
            aligned[row] = indexOf(id);
        }
        double[] covariance = new double[observations * observations];
        for (int row = 0; row < observations; row++)
            for (int column = 0; column < observations; column++)
                covariance[row * observations + column] =
                    relationship[aligned[row] * samples() + aligned[column]];
        return new VarianceComponent(name, observations, covariance);
    }

    private static Map<String, Integer> index(List<String> ids) {
        Map<String, Integer> result = new HashMap<>();
        for (int index = 0; index < ids.size(); index++) {
            String id = ids.get(index);
            if (id == null || id.isBlank() || result.put(id, index) != null)
                throw new IllegalArgumentException(
                    "sample IDs must be nonblank and unique");
        }
        return Map.copyOf(result);
    }

    private static void validateRelationship(
            double[] matrix, int dimension) {
        double maximum = 0;
        for (double value : matrix)
            maximum = Math.max(maximum, Math.abs(value));
        double tolerance = 1e-10 * Math.max(1, maximum);
        for (int row = 0; row < dimension; row++) {
            if (!(matrix[row * dimension + row] > 0))
                throw new IllegalArgumentException(
                    "GRM diagonal entries must be positive");
            for (int column = 0; column < row; column++)
                if (Math.abs(matrix[row * dimension + column]
                        - matrix[column * dimension + row]) > tolerance)
                    throw new IllegalArgumentException(
                        "genomic relationship matrix must be symmetric");
        }
    }
}
