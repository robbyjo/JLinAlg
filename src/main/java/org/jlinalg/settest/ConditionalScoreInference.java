/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Local summary-score conditioning under the exported null curvature.
 *
 * <p>The Schur complement is evaluated separately within each independent
 * cohort, then aligned scores and information are summed. This is deliberately
 * not a nonlinear cohort likelihood refit.</p>
 */
public final class ConditionalScoreInference {
    private ConditionalScoreInference() { }

    public record Result(List<String> targetVariantKeys,
            List<String> conditioningVariantKeys, SetTestScoreState state,
            int[] cohortCounts, String[] directions, List<String> cohorts) {
        public Result {
            targetVariantKeys = List.copyOf(targetVariantKeys);
            conditioningVariantKeys = List.copyOf(conditioningVariantKeys);
            cohortCounts = cohortCounts.clone();
            directions = directions.clone();
            cohorts = List.copyOf(cohorts);
        }
        @Override public int[] cohortCounts() { return cohortCounts.clone(); }
        @Override public String[] directions() { return directions.clone(); }
    }

    /** Condition targets on the supplied variants using complete cohort score
     * blocks. All cohorts must contain every target, condition and pairwise
     * covariance and must share a compatible model contract. */
    public static Result condition(List<ConditionalScoreStudy> studies,
            List<String> targetVariantKeys,
            List<String> conditioningVariantKeys) {
        if (studies == null || studies.isEmpty())
            throw new IllegalArgumentException("at least one conditional score cohort is required");
        List<String> targets = validated(targetVariantKeys, "target");
        List<String> conditions = validated(conditioningVariantKeys,
            "conditioning");
        if (conditions.isEmpty())
            throw new IllegalArgumentException("at least one conditioning variant is required");
        Set<String> identities = new HashSet<>();
        for (String value : targets)
            identities.add(ScoreVariantKey.parse(value).unorientedIdentity());
        for (String value : conditions)
            if (!identities.add(ScoreVariantKey.parse(value).unorientedIdentity()))
                throw new IllegalArgumentException(
                    "target and conditioning variants overlap: " + value);

        ConditionalScoreStudy first = studies.get(0);
        Set<String> fitted = new HashSet<>();
        for (String value : first.contract().fittedConditioningVariants())
            fitted.add(ScoreVariantKey.parse(value).unorientedIdentity());
        for (String value : targets) if (fitted.contains(
                ScoreVariantKey.parse(value).unorientedIdentity()))
            throw new IllegalArgumentException("target variant is already in the fitted null conditioning set: " + value);
        for (String value : conditions) if (fitted.contains(
                ScoreVariantKey.parse(value).unorientedIdentity()))
            throw new IllegalArgumentException("conditioning variant is already in the fitted null conditioning set: " + value);
        Set<String> cohortNames = new HashSet<>(), nullModels = new HashSet<>();
        List<String> cohorts = new ArrayList<>();
        double[][] scores = new double[studies.size()][];
        double[][] information = new double[studies.size()][];
        int t = targets.size(), c = conditions.size(), n = t + c;
        int[] ti = new int[t], ci = new int[c];
        for (int i = 0; i < t; i++) ti[i] = i;
        for (int i = 0; i < c; i++) ci[i] = t + i;

        for (int h = 0; h < studies.size(); h++) {
            ConditionalScoreStudy study = studies.get(h);
            if (!cohortNames.add(study.cohort()))
                throw new IllegalArgumentException(
                    "duplicate conditional score cohort: " + study.cohort());
            if (!nullModels.add(study.nullModelId()))
                throw new IllegalArgumentException(
                    "duplicate null_model_id would double-count a cohort: "
                    + study.nullModelId());
            first.contract().requireCompatible(study.contract());
            cohorts.add(study.cohort());

            Map<String,Integer> source = new HashMap<>();
            List<String> keys = study.variantKeys();
            for (int i = 0; i < keys.size(); i++)
                source.put(ScoreVariantKey.parse(keys.get(i)).unorientedIdentity(), i);
            List<String> requested = new ArrayList<>(targets); requested.addAll(conditions);
            int[] order = new int[n], sign = new int[n];
            for (int i = 0; i < n; i++) {
                ScoreVariantKey wanted = ScoreVariantKey.parse(requested.get(i));
                Integer found = source.get(wanted.unorientedIdentity());
                if (found == null)
                    throw new IllegalArgumentException("cohort " + study.cohort()
                        + " lacks requested score variant " + wanted);
                ScoreVariantKey available = ScoreVariantKey.parse(keys.get(found));
                order[i] = found;
                sign[i] = wanted.alignmentSign(available);
                if (sign[i] < 0 && !study.contract().permitsRefAltSwap())
                    throw new IllegalArgumentException("cohort " + study.cohort()
                        + " requires a REF/ALT score swap for " + wanted
                        + ", but its null model has no intercept and is not "
                        + "translation invariant; request the exported orientation");
            }
            double[] su = study.state().scores(), sv = study.state().information();
            int columns = study.state().variants();
            double[] u = new double[n], v = new double[Math.multiplyExact(n, n)];
            for (int i = 0; i < n; i++) {
                u[i] = sign[i] * su[order[i]];
                for (int j = 0; j < n; j++)
                    v[i*n+j] = sign[i] * sign[j]
                        * sv[order[i]*columns+order[j]];
            }
            SetTestScoreState adjusted = SummaryScoreModels.condition(
                new SetTestScoreState(u, v, n), ti, ci);
            scores[h] = adjusted.scores();
            information[h] = adjusted.information();
        }
        ScoreMetaAnalysis.Pooled pooled = ScoreMetaAnalysis.pool(scores,
            information, studies.size());
        if (pooled == null || pooled.state().variants() != t)
            throw new IllegalArgumentException(
                "not every target has positive conditional information in every cohort");
        int[] indices = pooled.indices();
        for (int i = 0; i < indices.length; i++) if (indices[i] != i)
            throw new IllegalArgumentException("conditional target coverage is incomplete");
        return new Result(targets, conditions, pooled.state(),
            pooled.cohortCounts(), pooled.directions(), cohorts);
    }

    private static List<String> validated(List<String> values, String label) {
        if (values == null || values.isEmpty()) {
            if (label.equals("target"))
                throw new IllegalArgumentException("at least one target variant is required");
            return List.of();
        }
        List<String> result = List.copyOf(values); Set<String> seen = new HashSet<>();
        for (String value : result) {
            ScoreVariantKey key = ScoreVariantKey.parse(value);
            if (!seen.add(key.unorientedIdentity()))
                throw new IllegalArgumentException("duplicate " + label
                    + " variant: " + value);
        }
        return result;
    }
}
