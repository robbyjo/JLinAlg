/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.settest.ConditionalScoreInference;
import org.jlinalg.settest.ConditionalScoreStudy;
import org.jlinalg.settest.ScoreVariantKey;
import org.jlinalg.settest.SetTestResult;
import org.jlinalg.settest.SummarySetTests;

/** Summary-only conditioning of versioned conditional-GWAS score exports. */
final class ConditionalScoreCli {
    private ConditionalScoreCli() { }

    static int run(String[] arguments, PrintStream output,
            PrintStream errorOutput) {
        try {
            if (Arrays.asList(arguments).contains("--help")
                    || Arrays.asList(arguments).contains("-h")) {
                output.println(help()); return 0;
            }
            Map<String,String> options = XwasFiles.options(arguments,
                "--cohorts", "--targets", "--condition-on", "--out");
            Path cohortPath = XwasFiles.path(options, "--cohorts");
            Path outputPath = XwasFiles.path(options, "--out");
            List<String> targets = XwasFiles.names(
                XwasFiles.required(options, "--targets"));
            List<String> conditions = XwasFiles.names(
                XwasFiles.required(options, "--condition-on"));
            validateRequested(targets, conditions);
            List<String> requested = new ArrayList<>(targets);
            requested.addAll(conditions);

            DelimitedData cohorts = DelimitedData.read(cohortPath);
            int nameColumn = cohorts.column("cohort");
            int summaryColumn = cohorts.column("summary");
            int covarianceColumn = cohorts.column("covariance");
            int manifestColumn = cohorts.column("manifest");
            Set<String> names = new HashSet<>();
            List<ConditionalScoreStudy> studies = new ArrayList<>();
            Path directory = cohortPath.getParent();
            for (String[] row : cohorts.rows()) {
                String name = XwasFiles.id(row[nameColumn]);
                if (!names.add(name))
                    throw new IllegalArgumentException(
                        "duplicate conditional score cohort: " + name);
                studies.add(ConditionalScoreImporter.read(name,
                    resolve(directory, row[summaryColumn], "summary"),
                    resolve(directory, row[covarianceColumn], "covariance"),
                    resolve(directory, row[manifestColumn], "manifest"),
                    requested));
            }
            ConditionalScoreInference.Result fit =
                ConditionalScoreInference.condition(studies, targets, conditions);
            publish(outputPath, fit, studies.get(0).contract());
            output.println("conditional-score complete: " + targets.size()
                + " target(s), " + conditions.size() + " condition(s), "
                + studies.size() + " independent cohort(s)");
            return 0;
        } catch (IOException | RuntimeException failure) {
            errorOutput.println("jlinalg: " + failure.getMessage());
            return 2;
        }
    }

    private static void validateRequested(List<String> targets,
            List<String> conditions) {
        Set<String> identities = new HashSet<>();
        for (String value : targets) {
            ScoreVariantKey key = ScoreVariantKey.parse(value);
            if (!identities.add(key.unorientedIdentity()))
                throw new IllegalArgumentException("duplicate target variant: " + value);
        }
        for (String value : conditions) {
            ScoreVariantKey key = ScoreVariantKey.parse(value);
            if (!identities.add(key.unorientedIdentity()))
                throw new IllegalArgumentException(
                    "target and conditioning variants overlap: " + value);
        }
    }

    private static Path resolve(Path directory, String value, String field) {
        if (value == null || value.isBlank() || !value.equals(value.trim()))
            throw new IllegalArgumentException("cohort " + field
                + " path must be nonblank and unpadded");
        Path path = Path.of(value);
        return (path.isAbsolute() ? path : directory.resolve(path))
            .toAbsolutePath().normalize();
    }

    private static void publish(Path output,
            ConditionalScoreInference.Result fit,
            ConditionalScoreStudy.ModelContract contract) throws IOException {
        double[] scores = fit.state().scores(), covariance = fit.state().information();
        int n = scores.length;
        String conditioned = String.join(",", fit.conditioningVariantKeys());
        StringBuilder results = new StringBuilder(
            "score_variant_key\tscore_u\tscore_variance\tbeta_score\tse_score\t"
            + "z_score\tp_score_normal\tcohorts\tdirection\tconditioned_on\t"
            + "inference_scope\n");
        for (int i = 0; i < n; i++) {
            SetTestResult result = SummarySetTests.singleVariant(
                fit.targetVariantKeys().get(i), scores[i], covariance[i*n+i]);
            results.append(fit.targetVariantKeys().get(i)).append('\t')
                .append(scores[i]).append('\t').append(covariance[i*n+i])
                .append('\t').append(result.beta()).append('\t')
                .append(result.standardError()).append('\t')
                .append(result.statistic()).append('\t').append(result.pValue())
                .append('\t').append(fit.cohortCounts()[i]).append('\t')
                .append(fit.directions()[i]).append('\t').append(conditioned)
                .append("\tlocal_schur_one_step\n");
        }
        StringBuilder matrix = new StringBuilder(
            "variant_i\tvariant_j\tscore_covariance\n");
        for (int i = 0; i < n; i++) for (int j = i; j < n; j++)
            matrix.append(fit.targetVariantKeys().get(i)).append('\t')
                .append(fit.targetVariantKeys().get(j)).append('\t')
                .append(covariance[i*n+j]).append('\n');
        String metadata = "key\tvalue\n"
            + "source_schema\t" + field(contract.exportSchema()) + "\n"
            + "genome_build\t" + field(contract.genomeBuild()) + "\n"
            + "model\t" + field(contract.model()) + "\n"
            + "family\t" + field(contract.family()) + "\n"
            + "formula\t" + field(contract.formula()) + "\n"
            + "case_value\t" + field(contract.caseValue()) + "\n"
            + "control_value\t" + field(contract.controlValue()) + "\n"
            + "analysis_sample\t" + field(contract.analysisSample()) + "\n"
            + "fitted_conditioning_variants\t" + field(String.join(",",
                contract.fittedConditioningVariants())) + "\n"
            + "cohorts\t" + field(String.join(",", fit.cohorts())) + "\n"
            + "targets\t" + field(String.join(",", fit.targetVariantKeys())) + "\n"
            + "conditioned_on\t" + field(conditioned) + "\n"
            + "method\tcohort-wise Schur complement followed by independent fixed score pooling\n"
            + "inference_scope\tlocal null-score curvature; one-step normal inference\n"
            + "nonlinear_cohort_refit\tfalse\n"
            + "cross_block_covariance\tunknown and rejected; never filled with zero or external LD\n"
            + "cohort_overlap\tnot modeled; input cohorts must be independent\n"
            + "allele_alignment\texact forward-strand REF/ALT match; swaps require a Cox or explicit-intercept null; no complement inference\n";
        Map<Path,String> files = new LinkedHashMap<>();
        files.put(output, results.toString());
        files.put(Path.of(output + ".score-cov.tsv"), matrix.toString());
        files.put(Path.of(output + ".metadata.tsv"), metadata);
        XwasFiles.publish(files);
    }

    private static String field(String value) {
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

    static String help() {
        return """
            Usage: jlinalg conditional-score --cohorts FILE
              --targets CHR:POS:REF:ALT[,KEY...] --condition-on KEY[,KEY...]
              --out FILE.tsv
            Cohort manifest columns: cohort, summary, covariance, manifest.
            Paths may be absolute or relative to the cohort manifest.
            Each trio must be a completed jlinalg-conditional-score-v1 export.
            Exact forward-strand REF/ALT swaps are aligned only for Cox or
            explicit-intercept null models; complements are not inferred.
            Every requested target, condition and covariance must occur in one complete
            covariance block per cohort. Cross-block covariance is unknown, not zero.
            Null model/family/formula/covariate/coding contracts must match exactly.
            The Schur complement is applied within each independent cohort before pooling.
            Output is local one-step normal score inference, not a nonlinear cohort refit
            and not the separate mr-estimate Gaussian conditional interface.
            Produces FILE.tsv, FILE.tsv.score-cov.tsv and FILE.tsv.metadata.tsv.
            Existing outputs are never replaced.
            """;
    }
}
