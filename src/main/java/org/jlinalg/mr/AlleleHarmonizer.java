/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Aligns two sets of biallelic SNP summary associations by variant and allele. */
public final class AlleleHarmonizer {
    private AlleleHarmonizer() {
    }

    /** Harmonizes with default frequency tolerance. */
    public static HarmonizationResult harmonize(
            List<SummaryAssociation> exposure,
            List<SummaryAssociation> outcome) {
        return harmonize(exposure, outcome, HarmonizationOptions.defaults());
    }

    /** Harmonizes outcome effects to exposure effect alleles. */
    public static HarmonizationResult harmonize(
            List<SummaryAssociation> exposure,
            List<SummaryAssociation> outcome,
            HarmonizationOptions options) {
        if (exposure == null || outcome == null || options == null) {
            throw new IllegalArgumentException(
                "exposure, outcome, and options are required");
        }
        DuplicateIndex exposureIndex = index(exposure, "exposure");
        DuplicateIndex outcomeIndex = index(outcome, "outcome");
        List<HarmonizedInstrument> retained = new ArrayList<>();
        List<HarmonizationExclusion> exclusions = new ArrayList<>();
        Set<String> handled = new HashSet<>();

        for (SummaryAssociation exposureAssociation : exposure) {
            if (exposureAssociation == null) {
                throw new IllegalArgumentException(
                    "exposure associations must not contain null");
            }
            String variant = exposureAssociation.variantId();
            if (!handled.add(variant)) {
                continue;
            }
            if (exposureIndex.duplicates().contains(variant)) {
                exclusions.add(exclusion(variant,
                    HarmonizationExclusionReason.DUPLICATE_EXPOSURE,
                    "variant occurs more than once in exposure associations"));
                continue;
            }
            if (outcomeIndex.duplicates().contains(variant)) {
                exclusions.add(exclusion(variant,
                    HarmonizationExclusionReason.DUPLICATE_OUTCOME,
                    "variant occurs more than once in outcome associations"));
                continue;
            }
            SummaryAssociation outcomeAssociation = outcomeIndex.values().get(variant);
            if (outcomeAssociation == null) {
                exclusions.add(exclusion(variant,
                    HarmonizationExclusionReason.MISSING_OUTCOME,
                    "variant has no outcome association"));
                continue;
            }
            if (exposureAssociation.effect() == 0.0) {
                exclusions.add(exclusion(variant,
                    HarmonizationExclusionReason.ZERO_EXPOSURE_EFFECT,
                    "Wald ratios are undefined for a zero exposure effect"));
                continue;
            }
            Alignment alignment = alignment(
                exposureAssociation, outcomeAssociation, options);
            if (alignment.exclusionReason() != null) {
                exclusions.add(exclusion(variant, alignment.exclusionReason(),
                    alignment.detail()));
                continue;
            }
            double alignedFrequency = outcomeAssociation.hasEffectAlleleFrequency()
                ? (alignment.flipEffect()
                    ? 1.0 - outcomeAssociation.effectAlleleFrequency()
                    : outcomeAssociation.effectAlleleFrequency())
                : Double.NaN;
            retained.add(new HarmonizedInstrument(
                variant,
                exposureAssociation.effectAllele(),
                exposureAssociation.otherAllele(),
                exposureAssociation.effect(),
                exposureAssociation.standardError(),
                alignment.flipEffect()
                    ? -outcomeAssociation.effect() : outcomeAssociation.effect(),
                outcomeAssociation.standardError(),
                exposureAssociation.effectAlleleFrequency(),
                alignedFrequency,
                alignment.flipEffect(),
                alignment.strandComplemented()));
        }
        return new HarmonizationResult(retained, exclusions);
    }

    private static Alignment alignment(
            SummaryAssociation exposure,
            SummaryAssociation outcome,
            HarmonizationOptions options) {
        boolean palindromic = isPalindromic(
            exposure.effectAllele(), exposure.otherAllele());
        if (palindromic) {
            return palindromicAlignment(exposure, outcome, options);
        }

        String exposureEffect = exposure.effectAllele();
        String exposureOther = exposure.otherAllele();
        String outcomeEffect = outcome.effectAllele();
        String outcomeOther = outcome.otherAllele();
        if (exposureEffect.equals(outcomeEffect)
                && exposureOther.equals(outcomeOther)) {
            return Alignment.keep(false, false);
        }
        if (exposureEffect.equals(outcomeOther)
                && exposureOther.equals(outcomeEffect)) {
            return Alignment.keep(true, false);
        }
        String complementEffect = complement(outcomeEffect);
        String complementOther = complement(outcomeOther);
        if (exposureEffect.equals(complementEffect)
                && exposureOther.equals(complementOther)) {
            return Alignment.keep(false, true);
        }
        if (exposureEffect.equals(complementOther)
                && exposureOther.equals(complementEffect)) {
            return Alignment.keep(true, true);
        }
        return Alignment.exclude(HarmonizationExclusionReason.ALLELE_MISMATCH,
            "exposure and outcome allele pairs do not match under swapping or strand complement");
    }

    private static Alignment palindromicAlignment(
            SummaryAssociation exposure,
            SummaryAssociation outcome,
            HarmonizationOptions options) {
        Set<String> exposureAlleles = Set.of(
            exposure.effectAllele(), exposure.otherAllele());
        Set<String> outcomeAlleles = Set.of(
            outcome.effectAllele(), outcome.otherAllele());
        if (!exposureAlleles.equals(outcomeAlleles)) {
            return Alignment.exclude(HarmonizationExclusionReason.ALLELE_MISMATCH,
                "palindromic exposure and outcome allele pairs differ");
        }
        if (!exposure.hasEffectAlleleFrequency()
                || !outcome.hasEffectAlleleFrequency()) {
            return Alignment.exclude(
                HarmonizationExclusionReason.PALINDROMIC_AMBIGUOUS,
                "palindromic SNP requires exposure and outcome effect-allele frequencies");
        }
        double sameDistance = Math.abs(exposure.effectAlleleFrequency()
            - outcome.effectAlleleFrequency());
        double flippedDistance = Math.abs(exposure.effectAlleleFrequency()
            - (1.0 - outcome.effectAlleleFrequency()));
        double tolerance = options.alleleFrequencyTolerance();
        boolean samePlausible = sameDistance <= tolerance;
        boolean flippedPlausible = flippedDistance <= tolerance;
        if (samePlausible == flippedPlausible) {
            HarmonizationExclusionReason reason = samePlausible
                ? HarmonizationExclusionReason.PALINDROMIC_AMBIGUOUS
                : HarmonizationExclusionReason.FREQUENCY_MISMATCH;
            String detail = samePlausible
                ? "allele frequencies cannot distinguish palindromic orientation"
                : "palindromic allele frequencies disagree beyond tolerance";
            return Alignment.exclude(reason, detail);
        }
        return Alignment.keep(flippedPlausible, false);
    }

    private static boolean isPalindromic(String first, String second) {
        return complement(first).equals(second);
    }

    private static String complement(String allele) {
        return switch (allele) {
            case "A" -> "T";
            case "T" -> "A";
            case "C" -> "G";
            case "G" -> "C";
            default -> throw new IllegalArgumentException("invalid DNA allele: " + allele);
        };
    }

    private static DuplicateIndex index(
            List<SummaryAssociation> associations, String name) {
        Map<String, SummaryAssociation> values = new LinkedHashMap<>();
        Set<String> duplicates = new HashSet<>();
        for (SummaryAssociation association : associations) {
            if (association == null) {
                throw new IllegalArgumentException(
                    name + " associations must not contain null");
            }
            if (values.putIfAbsent(association.variantId(), association) != null) {
                duplicates.add(association.variantId());
            }
        }
        return new DuplicateIndex(Map.copyOf(values), Set.copyOf(duplicates));
    }

    private static HarmonizationExclusion exclusion(
            String variant,
            HarmonizationExclusionReason reason,
            String detail) {
        return new HarmonizationExclusion(variant, reason, detail);
    }

    private record DuplicateIndex(
        Map<String, SummaryAssociation> values, Set<String> duplicates) {
    }

    private record Alignment(
            boolean flipEffect,
            boolean strandComplemented,
            HarmonizationExclusionReason exclusionReason,
            String detail) {
        private static Alignment keep(boolean flip, boolean complemented) {
            return new Alignment(flip, complemented, null, "");
        }

        private static Alignment exclude(
                HarmonizationExclusionReason reason, String detail) {
            return new Alignment(false, false, reason, detail);
        }
    }
}
