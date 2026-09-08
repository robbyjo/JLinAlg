/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;

final class SemBackendPolicyTest {
    private static final double[][] DATA = {{-1,-1},{-1,1},{1,-1},{1,1}};
    private static final double[][] MISSING = {
        {-1,-1},{-1,1},{1,-1},{1,1},{0,Double.NaN},{Double.NaN,0}};
    private static final double[] COVARIANCE = {1,0,0,1};

    @ParameterizedTest
    @EnumSource(value = BackendPolicy.class, names = {"CPU", "PREFERRED", "AUTO"},
        mode = EnumSource.Mode.EXCLUDE)
    void explicitUnsupportedPoliciesAreRejectedAcrossPublicFitPaths(BackendPolicy policy) {
        SemModel model = SemModel.builder("x","y").build();
        SemModel means = model.toBuilder().meanStructure().build();
        List<org.junit.jupiter.api.function.Executable> calls = List.of(
            () -> Sem.fit(DATA, model, SemOptions.defaults(), policy),
            () -> Sem.fitCovariance(COVARIANCE, 4, model, SemOptions.defaults(), policy),
            () -> Sem.fitMoments(COVARIANCE, new double[2], 4, means, SemOptions.defaults(), policy),
            () -> SemFiml.fit(MISSING, means, 10000, 1e-8, policy));
        for (var call : calls) {
            UnsupportedOperationException exception = assertThrows(UnsupportedOperationException.class, call);
            assertTrue(exception.getMessage().contains(policy.name()));
            assertTrue(exception.getMessage().contains("CPU, PREFERRED, or AUTO"));
        }
    }

    @ParameterizedTest
    @EnumSource(value = BackendPolicy.class, names = {"CPU", "PREFERRED", "AUTO"})
    void cpuAndFallbackPoliciesRetainValidFitsAndTruthfulProvenance(BackendPolicy policy) {
        SemModel model = SemModel.builder("x","y").build();
        SemModel means = model.toBuilder().meanStructure().build();
        List<SemFitResult> fits = List.of(
            Sem.fit(DATA, model, SemOptions.defaults(), policy),
            Sem.fitCovariance(COVARIANCE, 4, model, SemOptions.defaults(), policy),
            Sem.fitMoments(COVARIANCE, new double[2], 4, means, SemOptions.defaults(), policy),
            SemFiml.fit(MISSING, means, 10000, 1e-8, policy).fit());
        for (SemFitResult fit : fits) {
            assertTrue(fit.converged(), "score=" + fit.scoreNorm());
            assertTrue(fit.informationAvailable());
            assertEquals(policy, fit.backend().requested());
            assertEquals("CPU", fit.backend().selectedBackend());
        }
    }
}
