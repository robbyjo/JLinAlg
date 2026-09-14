/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

/** One validated cohort score block imported from a conditional-GWAS export. */
public final class ConditionalScoreStudy {
    /** Fields that must describe the same null-score estimand before cohorts
     * can be pooled. Cohort sample sizes and Gaussian dispersion may differ. */
    public record ModelContract(String exportSchema, String genomeBuild,
            String model, String family, String formula,
            List<String> nullCovariateColumns, String phenotypeTransform,
            String caseValue, String controlValue, String analysisSample,
            String genotypeCoding, String missingGenotypes,
            String scoreCovarianceScale, String varianceMethod,
            String relatednessAdjustment, String ties,
            String tailCalibration, List<String> fittedConditioningVariants) {
        public ModelContract {
            for (String value : new String[]{exportSchema, genomeBuild, model,
                    family, formula, phenotypeTransform, genotypeCoding,
                    missingGenotypes, analysisSample, scoreCovarianceScale,
                    varianceMethod, relatednessAdjustment, tailCalibration})
                if (value == null || value.isBlank())
                    throw new IllegalArgumentException(
                        "conditional score model fields must be nonblank");
            ties = ties == null ? "" : ties;
            caseValue = caseValue == null ? "" : caseValue;
            controlValue = controlValue == null ? "" : controlValue;
            if (caseValue.isEmpty() != controlValue.isEmpty()
                    || !caseValue.isEmpty() && caseValue.equals(controlValue))
                throw new IllegalArgumentException(
                    "case/control values must both be absent or distinct");
            nullCovariateColumns = List.copyOf(Objects.requireNonNull(
                nullCovariateColumns, "nullCovariateColumns"));
            if (!exportSchema.equals("jlinalg-conditional-score-v1"))
                throw new IllegalArgumentException(
                    "unsupported conditional score export schema: " + exportSchema);
            Set<String> covariates = new HashSet<>();
            for (String value : nullCovariateColumns)
                if (value == null || value.isBlank() || !covariates.add(value))
                    throw new IllegalArgumentException(
                        "null covariate columns must be unique and nonblank");
            fittedConditioningVariants = List.copyOf(Objects.requireNonNull(
                fittedConditioningVariants, "fittedConditioningVariants"));
            Set<String> conditions = new HashSet<>();
            for (String value : fittedConditioningVariants) {
                ScoreVariantKey key = ScoreVariantKey.parse(value);
                if (!conditions.add(key.unorientedIdentity()))
                    throw new IllegalArgumentException(
                        "duplicate fitted conditioning variant: " + value);
            }
        }

        void requireCompatible(ModelContract other) {
            require("export schema", exportSchema, other.exportSchema);
            require("genome build", genomeBuild, other.genomeBuild);
            require("model", model, other.model);
            require("family", family, other.family);
            require("formula", formula, other.formula);
            require("null covariate columns", covariateSignature(this),
                covariateSignature(other));
            require("phenotype transform", phenotypeTransform,
                other.phenotypeTransform);
            require("case value", caseValue, other.caseValue);
            require("control value", controlValue, other.controlValue);
            require("analysis-sample policy", analysisSample,
                other.analysisSample);
            require("genotype coding", genotypeCoding, other.genotypeCoding);
            require("missing-genotype policy", missingGenotypes,
                other.missingGenotypes);
            require("score covariance scale", scoreCovarianceScale,
                other.scoreCovarianceScale);
            require("variance method", varianceMethod, other.varianceMethod);
            require("relatedness adjustment", relatednessAdjustment,
                other.relatednessAdjustment);
            require("ties method", ties, other.ties);
            require("tail calibration", tailCalibration, other.tailCalibration);
            require("fitted conditioning set", conditionSet(this),
                conditionSet(other));
        }

        /**
         * Whether replacing an ALT dosage {@code g} by {@code 2-g} leaves the
         * efficient score space translation invariant. Cox partial likelihood
         * centers genotype values within risk sets; the other exported models
         * require an explicitly fitted intercept.
         */
        public boolean permitsRefAltSwap() {
            return model.equals("cox")
                || nullCovariateColumns.contains("(Intercept)");
        }

        private static Set<String> conditionSet(ModelContract contract) {
            Set<String> result = new HashSet<>();
            for (String value : contract.fittedConditioningVariants) {
                ScoreVariantKey key = ScoreVariantKey.parse(value);
                result.add(contract.permitsRefAltSwap()
                    ? key.unorientedIdentity() : key.toString());
            }
            return result;
        }

        private static List<String> covariateSignature(ModelContract value) {
            Set<String> conditions = conditionSet(value);
            java.util.ArrayList<String> result = new java.util.ArrayList<>();
            for (String column : value.nullCovariateColumns) {
                try {
                    ScoreVariantKey key = ScoreVariantKey.parse(column);
                    String identity = value.permitsRefAltSwap()
                        ? key.unorientedIdentity() : key.toString();
                    if (conditions.contains(identity)) continue;
                } catch (IllegalArgumentException ignored) {
                    // Ordinary covariate name, not an allele key.
                }
                result.add("column:" + column);
            }
            return List.copyOf(result);
        }

        private static void require(String field, Object first, Object second) {
            if (!first.equals(second))
                throw new IllegalArgumentException(
                    "incompatible conditional score " + field + ": "
                    + first + " versus " + second);
        }
    }

    private final String cohort;
    private final String nullModelId;
    private final String conditioningSetId;
    private final ModelContract contract;
    private final List<String> variantKeys;
    private final SetTestScoreState state;

    public ConditionalScoreStudy(String cohort, String nullModelId,
            String conditioningSetId, ModelContract contract,
            List<String> variantKeys, SetTestScoreState state) {
        if (cohort == null || cohort.isBlank() || !cohort.equals(cohort.trim()))
            throw new IllegalArgumentException("cohort name must be nonblank and unpadded");
        if (nullModelId == null || nullModelId.isBlank()
                || conditioningSetId == null || conditioningSetId.isBlank())
            throw new IllegalArgumentException(
                "null-model and conditioning-set IDs are required");
        this.cohort = cohort;
        this.nullModelId = nullModelId;
        this.conditioningSetId = conditioningSetId;
        this.contract = Objects.requireNonNull(contract, "contract");
        this.variantKeys = List.copyOf(Objects.requireNonNull(
            variantKeys, "variantKeys"));
        this.state = Objects.requireNonNull(state, "state");
        if (this.variantKeys.size() != state.variants())
            throw new IllegalArgumentException(
                "one score variant key is required per score-state column");
        Set<String> variants = new HashSet<>();
        for (String value : this.variantKeys) {
            ScoreVariantKey key = ScoreVariantKey.parse(value);
            if (!variants.add(key.unorientedIdentity()))
                throw new IllegalArgumentException(
                    "duplicate or oppositely coded score variant: " + value);
        }
        SummarySetTests.validate(state);
    }

    public String cohort() { return cohort; }
    public String nullModelId() { return nullModelId; }
    public String conditioningSetId() { return conditioningSetId; }
    public ModelContract contract() { return contract; }
    public List<String> variantKeys() { return variantKeys; }
    public SetTestScoreState state() { return state; }
}
