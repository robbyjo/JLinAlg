/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;

/** Deterministic PET-PEESE and median-reflection trim-and-fill corrections. */
public final class MetaBiasCorrections {
    private MetaBiasCorrections() { }

    public static Result pet(List<MetaStudy> studies) { return regression(studies, false); }
    public static Result peese(List<MetaStudy> studies) { return regression(studies, true); }

    private static Result regression(List<MetaStudy> studies, boolean variancePredictor) {
        if (studies == null || studies.size() < 3) throw new IllegalArgumentException("bias correction needs at least three studies");
        double s00 = 0.0, s01 = 0.0, s11 = 0.0, t0 = 0.0, t1 = 0.0;
        for (MetaStudy study : studies) { double se = study.standardError(), predictor = variancePredictor ? se * se : se, weight = 1.0 / (se * se); s00 += weight; s01 += weight * predictor; s11 += weight * predictor * predictor; t0 += weight * study.effectSize(); t1 += weight * predictor * study.effectSize(); }
        double determinant = s00 * s11 - s01 * s01; if (!(determinant > 0.0)) throw new IllegalArgumentException("bias correction predictor is singular");
        double intercept = (t0 * s11 - t1 * s01) / determinant, slope = (s00 * t1 - s01 * t0) / determinant;
        return new Result(variancePredictor ? "PEESE" : "PET", intercept, slope, studies.size());
    }

    /** A transparent deterministic reflection correction around the study median. */
    public static TrimFillResult trimAndFill(List<MetaStudy> studies) {
        if (studies == null || studies.size() < 3) throw new IllegalArgumentException("trim-and-fill needs at least three studies");
        double[] effects = studies.stream().mapToDouble(MetaStudy::effectSize).sorted().toArray(); double median = effects[effects.length / 2];
        int positive = 0, negative = 0; for (double effect : effects) { if (effect > median) positive++; else if (effect < median) negative++; }
        int imputed = Math.abs(positive - negative); double sum = 0.0; for (double effect : effects) sum += effect;
        if (positive > negative) for (int i = 0; i < imputed; i++) sum += 2.0 * median - effects[effects.length - 1 - i];
        else for (int i = 0; i < imputed; i++) sum += 2.0 * median - effects[i];
        return new TrimFillResult(sum / (effects.length + imputed), imputed, "median-reflection");
    }

    public record Result(String method, double intercept, double slope, int studyCount) { }
    public record TrimFillResult(double adjustedEffect, int imputedStudyCount, String method) { }
}
