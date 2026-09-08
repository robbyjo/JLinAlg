/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import static org.junit.jupiter.api.Assertions.*;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class GeneticCovarianceAuditTest {
    @Test void zeroVarianceSampleIsAValidMeanImputedGrmRow() {
        var grm=GenomicRelationshipMatrix.fromVariantDosages(new double[][]{{0,1,2}},List.of("a","b","c"),GenomicRelationshipOptions.defaults(),BackendPolicy.CPU);
        assertArrayEquals(new double[]{2,0,-2,0,0,0,-2,0,2},grm.relationshipMatrix(),1e-15);
        assertDoesNotThrow(()->grm.varianceComponent("g"));
        assertDoesNotThrow(()->new GenomicRelationshipMatrix(grm.sampleIds(),grm.relationshipMatrix()));
    }
    @Test void rejectsIndefiniteOrNonzeroCovarianceOnAZeroVarianceRow() {
        assertThrows(IllegalArgumentException.class,()->new GenomicRelationshipMatrix(List.of("a","b"),new double[]{1,2,2,1}));
        assertThrows(IllegalArgumentException.class,()->new GenomicRelationshipMatrix(List.of("a","b"),new double[]{1,.1,.1,0}));
        assertThrows(IllegalArgumentException.class,()->GeneticCovarianceValidation.requirePositiveSemidefinite(new double[]{1,.9,.9,.9,1,-.9,.9,-.9,1},3));
    }
    @Test void psdCheckIsRelativeToEachVariableScale() {
        assertDoesNotThrow(()->GeneticCovarianceValidation.requirePositiveSemidefinite(new double[]{1e-200,.5,.5,1e200},2));
        assertThrows(IllegalArgumentException.class,()->GeneticCovarianceValidation.requirePositiveSemidefinite(new double[]{1e-200,2,2,1e200},2));
    }
}
