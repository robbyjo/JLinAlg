/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.enrichment;

import java.util.List;

/** One term in the complete tested family; log probabilities survive ordinary underflow. */
public record EnrichmentResult(String id, String name, int universeSize,
        int selectedSize, int setSize, int overlap, double weightedOverlap,
        double expectedOverlap, double foldEnrichment, double oddsRatio,
        double selectionOdds, double pValue, double logPValue,
        double adjustedPValue, double logAdjustedPValue, List<String> hitGenes) {
    public EnrichmentResult { hitGenes = List.copyOf(hitGenes); }
}
