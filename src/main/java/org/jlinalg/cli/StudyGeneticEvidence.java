/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyNeed;
import org.jlinalg.study.StudyOperationExecutor;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.jlinalg.study.StudyOperationExecutor.Result;
import org.jlinalg.study.StudyOperationExecutor.State;

/** Reviewed regional evidence -> allele-aligned SuSiE, coloc and cis-MR. */
final class StudyGeneticEvidence implements StudyOperationExecutor {
    record Variant(String id, String chromosome, long position, String effect, String other,
            double beta, double se, double p) { }
    record Aligned(Variant exposure, Variant outcome, int outcomeSign, int ldSign) { }
    private static final class NotEstimableException extends IOException {
        private static final long serialVersionUID = 1L;
        private final String interpretation;
        NotEstimableException(String message, String interpretation) {
            super(message);
            this.interpretation = interpretation;
        }
    }

    @Override public Result execute(Context context) throws IOException {
        Files.createDirectories(context.outputDirectory());
        List<Object> results = new ArrayList<>();
        List<StudyNeed> needs = new ArrayList<>();
        Set<String> ids = new LinkedHashSet<>();
        List<?> configuredPairs = list(StudyBioinformatics.config(context).get("genetic_pairs"));
        Map<String, Object> automaticSources = Map.of();
        String automaticSourceKind = "";
        if (configuredPairs.isEmpty() && StudyMrSources.configured(context)) {
            try {
                StudyMrSources.Prepared prepared = StudyMrSources.prepare(context);
                configuredPairs = prepared.pairs(); automaticSources = prepared.audit(); automaticSourceKind = "mr";
            } catch (IllegalArgumentException | IOException exception) {
                Map<String, Object> details = new LinkedHashMap<>();
                details.put("automatic_mr_sources", Map.of("status", "waiting_for_input", "reason", String.valueOf(exception.getMessage())));
                StudyNeed need = StudyNeed.review("mr-sources:" + context.operation().get("hypothesis"),
                    "Can you resolve the automatic MR source prerequisite?",
                    String.valueOf(exception.getMessage()) + ". Review the concise analysis.bioinformatics.mr declaration and automatic-mr-sources.yaml when present; do not interpret an unavailable public source as a null result.", true);
                details.put("needs", List.of(need.document()));
                write(context.outputDirectory().resolve("genetic-evidence.yaml"), details);
                return new Result(State.WAITING_FOR_INPUT, "Automatic MR sources require input: " + exception.getMessage(),
                    List.of("genetic-evidence.yaml"), details);
            }
        } else if (configuredPairs.isEmpty() && StudyColocalizationSources.configured(context)) {
            try {
                StudyColocalizationSources.Prepared prepared = StudyColocalizationSources.prepare(context);
                configuredPairs = prepared.pairs(); automaticSources = prepared.audit(); automaticSourceKind = "colocalization";
            } catch (IllegalArgumentException | IOException exception) {
                Map<String, Object> details = new LinkedHashMap<>();
                details.put("automatic_colocalization_sources", Map.of("status", "waiting_for_input", "reason", String.valueOf(exception.getMessage())));
                StudyNeed need = StudyNeed.review("colocalization-sources:" + context.operation().get("hypothesis"),
                    "Can you resolve the automatic colocalization source prerequisite?",
                    String.valueOf(exception.getMessage()) + ". Review related-phenotype advice, source-discovery.yaml, and automatic-colocalization-sources.yaml when present. Missing or incomplete database coverage is not negative evidence.", true);
                details.put("needs", List.of(need.document()));
                write(context.outputDirectory().resolve("genetic-evidence.yaml"), details);
                return new Result(State.WAITING_FOR_INPUT, "Automatic colocalization sources require input: " + exception.getMessage(),
                    List.of("genetic-evidence.yaml"), details);
            }
        }
        for (Object value : configuredPairs) {
            Map<String, Object> pair = map(value);
            String id = required(pair, "id");
            if (!id.matches("[A-Za-z0-9][A-Za-z0-9_-]{0,63}") || !ids.add(id)) throw new IllegalArgumentException("Genetic pair IDs must be safe and unique");
            Map<String, Object> record = new LinkedHashMap<>();
            record.put("pair_id", id);
            record.put("feature_id", text(pair.get("feature_id")));
            record.put("sources", pair);
            results.add(record);
            Path directory = context.outputDirectory().resolve(id);
            Files.createDirectories(directory);
            try {
                executePair(context, pair, directory, record);
                record.put("status", "complete");
            } catch (NotEstimableException exception) {
                record.put("status", "not_estimable");
                record.put("reason", String.valueOf(exception.getMessage()));
                record.put("interpretation", exception.interpretation);
            } catch (IllegalArgumentException | IOException exception) {
                record.put("status", "waiting_for_input");
                record.put("reason", String.valueOf(exception.getMessage()));
                needs.add(StudyNeed.review("genetic-pair:" + context.operation().get("hypothesis") + ":" + id,
                    "Can you resolve the scientific or input prerequisite for genetic evidence pair " + id + "?",
                    String.valueOf(exception.getMessage()) + ". Review genetic-evidence.yaml and amend the pair/source declarations; do not treat a blocked analysis as negative evidence.", true));
            }
        }
        Map<String, Object> details = new LinkedHashMap<>();
        details.put("pairs", results);
        if (!automaticSources.isEmpty()) details.put("automatic_" + automaticSourceKind + "_sources", automaticSources);
        List<Object> caveats = results.stream().map(StudyPublicEvidence::map)
            .filter(row -> row.containsKey("required_manuscript_caveat"))
            .map(row -> (Object) Map.of("pair_id", row.get("pair_id"), "caveat", row.get("required_manuscript_caveat"))).toList();
        if (!caveats.isEmpty()) details.put("required_manuscript_caveats", caveats);
        List<Object> methodPairs = new ArrayList<>(); Set<String> methodReferences = new LinkedHashSet<>();
        for (Object raw : results) {
            var row = map(raw); var method = new LinkedHashMap<String, Object>();
            for (String key : List.of("pair_id", "status", "reason", "sources", "effective_parameters", "citation_methods", "directory",
                    "evidence_grade", "required_manuscript_caveat"))
                if (row.containsKey(key)) method.put(key, row.get(key));
            methodPairs.add(method); if ("complete".equals(row.get("status"))) for (Object ref : list(row.get("citation_methods"))) methodReferences.add(text(ref));
        }
        details.put("method_report", StudyMethods.record(text(context.operation().get("type")), Map.of("pairs", methodPairs),
            List.copyOf(methodReferences), "See the scientist-approved pair rationale and instrument-validity rationale; no additional scientific justification inferred.", List.of()));
        details.put("needs", needs.stream().map(StudyNeed::document).toList());
        details.put("interpretation", "Exploratory regional genetic evidence. Association, colocalization and MR have different assumptions; neither proves a safe or effective drug target. Selection of the omics hit and instruments can introduce selection bias.");
        write(context.outputDirectory().resolve("genetic-evidence.yaml"), details);
        long completed = results.stream().map(StudyPublicEvidence::map).filter(row -> "complete".equals(row.get("status"))).count();
        long notEstimable = results.stream().map(StudyPublicEvidence::map).filter(row -> "not_estimable".equals(row.get("status"))).count();
        return new Result(needs.isEmpty() ? State.COMPLETE : State.WAITING_FOR_INPUT,
            context.operation().get("type") + ": " + completed + " completed pairs, " + notEstimable
                + " completed but not estimable, " + needs.size() + " require input",
            List.of("genetic-evidence.yaml"), details);
    }

    private static void executePair(Context context, Map<String, Object> pair, Path directory, Map<String, Object> record) throws IOException {
        String feature = required(pair, "feature_id");
        String operation = text(context.operation().get("type"));
        boolean mrOnly = operation.equals("mendelian_randomization");
        boolean partialColocalization = operation.equals("colocalization") && Boolean.TRUE.equals(pair.get("partial_colocalization"));
        int featureLimit = mrOnly ? StudyBioinformatics.mrFeatureLimit(context)
            : StudyBioinformatics.featureLimit(context);
        StudyBioinformatics.Feature finding = StudyBioinformatics.features(context,
            featureLimit).stream().filter(item -> item.id().equals(feature)).findFirst()
            .orElseThrow(() -> new IllegalArgumentException("Genetic pair must identify an eligible unambiguously annotated finding from this hypothesis"));
        required(pair, "rationale");
        if (!Boolean.TRUE.equals(pair.get("approved"))) throw new IllegalArgumentException("Set approved:true only after scientist review of phenotype, molecular feature, tissue, ancestry, source citations and instrument assumptions");
        Map<String, Object> exposure = map(pair.get("exposure")), outcome = map(pair.get("outcome")), ld = map(pair.get("ld"));
        boolean automaticLd = !ld.containsKey("file_role") && !ld.containsKey("url");
        if (automaticLd) {
            String population = text(ld.getOrDefault("population", map(map(context.study().get("resources")).get("ld")).getOrDefault("population", "auto")));
            String database = text(ld.getOrDefault("database", map(map(context.study().get("resources")).get("ld")).getOrDefault("database", "1000g-phase3")));
            if (database.equals(ThousandGenomesGrch38Ld.DATABASE)) {
                String ancestry = ThousandGenomesGrch38Ld.ancestry(population);
                ld.put("database", database); ld.put("population", population);
                ld.put("genome_build", "GRCh38"); ld.put("ancestry", ancestry);
            } else {
                LdDatabaseSpec reference = LdDatabaseCatalog.find(database);
                if (reference == null) throw new IllegalArgumentException("Unknown automatic LD database: " + database);
                LdPanelSpec panel = reference.panels().stream().filter(item -> item.id().equals(population)).findFirst()
                    .orElseThrow(() -> new IllegalArgumentException("Choose a scientifically appropriate 1000G population (AFR/AMR/EAS/EUR/SAS); ancestry is not guessed"));
                ld.put("database", database); ld.put("population", population); ld.put("genome_build", reference.genomeBuild()); ld.put("ancestry", panel.ancestry());
            }
        }
        String build = required(exposure, "genome_build");
        if (!Set.of("GRCh37", "GRCh38").contains(build) || !build.equals(required(outcome, "genome_build")) || !build.equals(required(ld, "genome_build")))
            throw new IllegalArgumentException("Exposure, outcome and LD need the same explicit GRCh37/GRCh38 build; no implicit liftover");
        String ancestry = required(exposure, "ancestry");
        if (!ancestry.equals(required(outcome, "ancestry")) || !ancestry.equals(required(ld, "ancestry")))
            throw new IllegalArgumentException("Exposure, outcome and LD ancestry must match; heterogeneous ancestry requires a separate reviewed analysis");
        if (!feature.equals(required(exposure, "feature_id"))) throw new IllegalArgumentException("QTL feature_id must match the measured feature exactly; gene eQTL cannot silently stand in for a transcript/splicing QTL");
        required(exposure, "assay"); required(exposure, "tissue"); required(exposure, "trait"); required(outcome, "trait");
        String outcomeTrait = required(outcome, "trait");
        Set<String> declaredTraits = new LinkedHashSet<>();
        declaredTraits.add(text(StudyBioinformatics.hypothesis(context).get("outcome")));
        list(StudyBioinformatics.config(context).get("phenotype_terms")).forEach(term -> declaredTraits.add(text(term)));
        list(map(StudyBioinformatics.config(context).get("colocalization")).get("related_phenotypes")).forEach(term -> declaredTraits.add(text(term)));
        if (!declaredTraits.contains(outcomeTrait) && !Boolean.TRUE.equals(pair.get("automatic_colocalization_source")))
            throw new IllegalArgumentException("Outcome trait is not declared for this hypothesis; add a reviewed phenotype term or an amendment");
        if (!mrOnly && !partialColocalization && (!Boolean.TRUE.equals(exposure.get("complete_region")) || !Boolean.TRUE.equals(outcome.get("complete_region"))))
            throw new IllegalArgumentException("Colocalization/fine-mapping need full regional statistics for both sources, not significant-only or lead-SNP lists");
        if (partialColocalization && (!Boolean.TRUE.equals(exposure.get("complete_region")) || !Boolean.TRUE.equals(outcome.get("complete_region")))) {
            String caveat = "Exploratory partial-data colocalization was run because at least one regional source was incomplete. Posterior probabilities are conditional on the observed, selected variants and must not be interpreted as conventional complete-region evidence for or against a shared causal variant. Missing variants, significance filtering, unequal coverage, and winner's curse can materially bias the result; absence of colocalization is not evidence of absence.";
            record.put("evidence_grade", "exploratory_partial_region_major_limitation");
            record.put("required_manuscript_caveat", caveat);
            record.put("regional_completeness", Map.of("exposure_complete_region", Boolean.TRUE.equals(exposure.get("complete_region")),
                "outcome_complete_region", Boolean.TRUE.equals(outcome.get("complete_region")), "policy", "run_and_caveat"));
        }
        if (mrOnly && !Boolean.TRUE.equals(exposure.get("complete_region"))) {
            double exported = number(exposure.get("export_p_threshold"));
            double instruments = number(pair.getOrDefault("instrument_p", 5e-8));
            if (!(exported >= instruments && exported <= 1))
                throw new IllegalArgumentException("Filtered MR exposure must declare an export threshold at least as inclusive as the instrument threshold");
            record.put("exposure_selection_limitation", "Significance-filtered exposure; selected instruments may have winner's-curse bias. This input is not valid for regional colocalization.");
        }
        if (mrOnly && !Boolean.TRUE.equals(outcome.get("complete_region"))
                && !"exposure_instrument_lookup".equals(outcome.get("selection_basis")))
            throw new IllegalArgumentException("MR outcome statistics must be unfiltered regional results or an exposure_instrument_lookup, never selected by outcome significance");
        String overlap = required(pair, "sample_overlap");
        if (!overlap.equals("none") && !(partialColocalization && overlap.equals("unknown")))
            throw new IllegalArgumentException("Automatic pair inference requires nonoverlapping samples, except exploratory partial colocalization may retain unknown overlap as a major limitation");
        required(pair, "sample_overlap_basis");
        if (partialColocalization && overlap.equals("unknown")) {
            String overlapCaveat = " Sample overlap between QTL and GWAS sources was not established and can inflate evidence for sharing.";
            record.put("required_manuscript_caveat", text(record.get("required_manuscript_caveat")) + overlapCaveat);
        }
        int nx = StudyBioinformatics.integer(exposure, "sample_size", -1, 3, Integer.MAX_VALUE);
        int ny = StudyBioinformatics.integer(outcome, "sample_size", -1, 3, Integer.MAX_VALUE);
        Map<String, Object> region = map(pair.get("region"));
        String chromosome = required(region, "chromosome").replaceFirst("^chr", "");
        int start = StudyBioinformatics.integer(region, "start", -1, 1, Integer.MAX_VALUE);
        int end = StudyBioinformatics.integer(region, "end", -1, start, Integer.MAX_VALUE);
        if (end - (long) start > 10000000) throw new IllegalArgumentException("Regional window exceeds 10 Mb");
        if (!finding.status().equals("located") || !build.equals(finding.build()) || !chromosome.equals(finding.chromosome().replaceFirst("^chr", ""))
                || Long.parseLong(finding.start()) < start || Long.parseLong(finding.end()) > end)
            throw new IllegalArgumentException("Regional evidence must contain the measured feature on its verified genome build and chromosome");
        Path xPath = source(context, exposure, directory, "exposure-source.tsv", required(pair, "rationale") + " Molecular exposure statistics for " + feature + ".");
        Path yPath;
        if (outcome.containsKey("gwas_catalog_accession")) {
            yPath = directory.resolve("outcome-source.tsv");
            GwasCatalogRegionalOutcome.Result acquired;
            try {
                acquired = GwasCatalogRegionalOutcome.extract(context,
                    required(outcome, "gwas_catalog_accession"), ancestry, chromosome, start, end, yPath);
            } catch (IOException exception) {
                if (partialColocalization && String.valueOf(exception.getMessage()).startsWith(
                        "Fewer than two unambiguous effect estimates"))
                    throw new NotEstimableException(String.valueOf(exception.getMessage()),
                        "The prespecified public GWAS source was acquired and queried, but the region did not contain enough unambiguous effect estimates for colocalization. The source is exhausted, so this is a completed source-limitation result rather than a request for investigator input or evidence against sharing.");
                throw exception;
            }
            record.put("outcome_acquisition", acquired.provenance());
        } else yPath = source(context, outcome, directory, "outcome-source.tsv", required(pair, "rationale") + " Phenotype-matched outcome statistics for " + outcomeTrait + ".");
        Map<String, Variant> x = variants(xPath, chromosome, start, end, feature, null);
        Set<String> outcomeFilter = mrOnly || partialColocalization ? x.keySet() : null;
        Map<String, Variant> y = variants(yPath, chromosome, start, end, null, outcomeFilter);
        Path ldPath, ldAlleles;
        if (automaticLd) {
            Path[] generated;
            try {
                generated = automaticLd(context, ld, x, y, directory, mrOnly ? 1 : 2);
            } catch (IOException exception) {
                if (mrOnly && Boolean.TRUE.equals(pair.get("automatic_mr_source"))
                        && String.valueOf(exception.getMessage()).startsWith(
                            "Fewer than 1 requested variant is usable"))
                    throw new NotEstimableException(String.valueOf(exception.getMessage()),
                        "The prespecified GRCh38 LD panel was queried locally, but none of the outcome-matched cis-irQTL candidates was usable. The LD source is exhausted, so this transcript is completed but not estimable rather than waiting for investigator input.");
                if (partialColocalization && String.valueOf(exception.getMessage()).startsWith(
                        "Fewer than 2 requested variants are usable"))
                    throw new NotEstimableException(String.valueOf(exception.getMessage()),
                        "The prespecified GRCh38 LD panel was queried locally, but fewer than two requested variants were usable. The LD source is exhausted, so this pair is completed but not estimable rather than waiting for investigator input.");
                throw exception;
            }
            ldPath = generated[0]; ldAlleles = generated[1];
        } else {
            ldPath = source(context, ld, directory, "ld-source.tsv", required(pair, "rationale") + " Signed linkage disequilibrium for the harmonized regional variants.");
            Map<String, Object> alleleSource = new LinkedHashMap<>(ld);
            alleleSource.remove("url");
            alleleSource.put("file_role", required(ld, "alleles_file_role"));
            ldAlleles = source(context, alleleSource, directory, "ld-alleles.tsv", required(pair, "rationale") + " Verify effect-allele orientation of the LD matrix.");
        }
        DelimitedData alleles = DelimitedData.read(ldAlleles);
        Map<String, String[]> orientation = new LinkedHashMap<>();
        for (String[] row : alleles.rows()) {
            String id = row[alleles.column("variant_id")];
            if (orientation.put(id, new String[] {row[alleles.column("effect_allele")].toUpperCase(java.util.Locale.ROOT), row[alleles.column("other_allele")].toUpperCase(java.util.Locale.ROOT)}) != null)
                throw new IOException("Duplicate LD allele identifier: " + id);
        }
        List<Aligned> aligned = new ArrayList<>();
        List<Object> excluded = new ArrayList<>();
        for (Variant v : x.values()) {
            Variant w = y.get(v.id());
            String[] a = orientation.get(v.id());
            String reason = w == null ? "absent_from_outcome" : a == null ? "absent_from_ld_alleles"
                : !v.chromosome().equals(w.chromosome()) || v.position() != w.position() ? "coordinate_mismatch"
                : palindrome(v.effect(), v.other()) ? "palindromic_alleles_excluded" : "";
            int sign = w == null ? 0 : sign(v.effect(), v.other(), w.effect(), w.other());
            int ldSign = a == null ? 0 : sign(v.effect(), v.other(), a[0], a[1]);
            if (reason.isEmpty() && (sign == 0 || ldSign == 0)) reason = "incompatible_alleles";
            if (!reason.isEmpty()) { excluded.add(Map.of("variant_id", v.id(), "reason", reason)); continue; }
            aligned.add(new Aligned(v, w, sign, ldSign));
        }
        for (String id : y.keySet()) if (!x.containsKey(id)) excluded.add(Map.of("variant_id", id, "reason", "absent_from_exposure"));
        aligned.sort(Comparator.comparingLong((Aligned a) -> a.exposure().position()).thenComparing(a -> a.exposure().id()));
        record.put("exclusions", excluded);
        record.put("exposure_variants", x.size()); record.put("outcome_variants", y.size()); record.put("aligned_variants", aligned.size());
        int minimumAligned = mrOnly ? 1 : 2;
        if (aligned.size() < minimumAligned) {
            String reason = "Fewer than " + minimumAligned + " usable aligned variant"
                + (minimumAligned == 1 ? "" : "s") + " in the region";
            if (mrOnly && Boolean.TRUE.equals(pair.get("automatic_mr_source")))
                throw new NotEstimableException(reason,
                    "The prespecified exposure, Pan-UKB outcome and GRCh38 LD sources were exhausted, but no cis-irQTL could be aligned unambiguously. This transcript is completed but not estimable rather than waiting for investigator input.");
            if (partialColocalization) throw new NotEstimableException(reason,
                "The prespecified exposure, public outcome, and local LD sources were exhausted, but fewer than two variants could be aligned. The pair is completed but not estimable rather than waiting for investigator input or providing negative colocalization evidence.");
            throw new IllegalArgumentException(reason);
        }
        double exposureCoverage = aligned.size() / (double) x.size();
        double outcomeCoverage = aligned.size() / (double) y.size();
        double coverage = aligned.size() / (double) (mrOnly ? x.size() : Math.max(x.size(), y.size()));
        record.put("exposure_variant_match_fraction", exposureCoverage); record.put("outcome_variant_match_fraction", outcomeCoverage);
        record.put("overlap_denominator", mrOnly ? "exposure variant set" : "larger observed regional variant set");
        record.put("regional_overlap_fraction", coverage);
        String overlapKey = partialColocalization ? "minimum_exposure_overlap" : "minimum_overlap";
        if (!mrOnly || pair.containsKey(overlapKey)) {
            double minCoverage = number(pair.getOrDefault(overlapKey, partialColocalization ? 0.5 : 0.8));
            if (!(minCoverage > 0 && minCoverage <= 1)) throw new IllegalArgumentException("Declared overlap threshold must be in (0,1]");
            if ((partialColocalization ? exposureCoverage : coverage) < minCoverage)
                throw new IllegalArgumentException(partialColocalization
                    ? "Too few exported QTL variants matched the GWAS and LD sources for the declared exploratory threshold"
                    : "Regional overlap is below the declared minimum; inspect exclusions and supply matched full-region sources");
        } else if (exposureCoverage < 1) {
            record.put("instrument_lookup_limitation", "Only " + aligned.size() + " of " + x.size()
                + " exposure candidates had usable aligned outcome and LD evidence; MR used the retained prespecified instruments.");
        }
        // Automatic MR LD is a correlation matrix computed directly from the
        // reference genotypes and is used only for pairwise r^2 clumping.  A
        // cubic Cholesky PSD check adds no information there and can dominate
        // runtime for large cis-instrument families.  Keep that check for
        // user-supplied MR matrices and every fine-mapping/colocalization use.
        boolean verifyPositiveSemidefinite = !mrOnly || !automaticLd;
        double[][] matrix = matrix(ldPath, aligned, verifyPositiveSemidefinite);
        record.put("ld_validation", Map.of(
            "square_symmetric_bounded_unit_diagonal", true,
            "positive_semidefinite", verifyPositiveSemidefinite
                ? "verified_by_cholesky"
                : "guaranteed_by_direct_genotype_correlation; cubic recheck omitted for MR clumping",
            "usage", mrOnly ? "pairwise_r_squared_clumping" : "fine_mapping"));
        Path matrixPath = directory.resolve("aligned-ld.tsv");
        writeMatrix(matrixPath, matrix);
        StringBuilder audit = new StringBuilder("variant_id\tchromosome\tposition\teffect_allele\tother_allele\tleft_beta\tright_beta\tright_sign\tld_sign\n");
        for (Aligned a : aligned) audit.append(a.exposure().id()).append('\t').append(a.exposure().chromosome()).append('\t').append(a.exposure().position()).append('\t')
            .append(a.exposure().effect()).append('\t').append(a.exposure().other()).append('\t').append(a.exposure().beta()).append('\t')
            .append(a.outcome().beta() * a.outcomeSign()).append('\t').append(a.outcomeSign()).append('\t').append(a.ldSign()).append('\n');
        Files.writeString(directory.resolve("harmonization.tsv"), audit);
        if (operation.equals("mendelian_randomization")) mr(pair, aligned, matrix, directory, record);
        else {
            fine(aligned, false, nx, matrixPath, directory, "exposure", record, partialColocalization);
            fine(aligned, true, ny, matrixPath, directory, "outcome", record, partialColocalization);
            if (operation.equals("colocalization")) {
                if (!hasDataRows(directory.resolve("exposure.effects.tsv")))
                    throw new NotEstimableException("Exposure SuSiE converged but produced no credible signal for colocalization",
                        "The prespecified fine-mapping step completed, but no credible exposure signal was available for signal-pair colocalization; this is not evidence against a shared causal variant.");
                if (!hasDataRows(directory.resolve("outcome.effects.tsv")))
                    throw new NotEstimableException("Outcome SuSiE converged but produced no credible signal for colocalization",
                        "The prespecified fine-mapping step completed, but no credible outcome signal was available for signal-pair colocalization; this is not evidence against a shared causal variant.");
                // Make the previously implicit coloc defaults explicit in both commands and reporting.
                List<String> priors = List.of("--trait1-prior", "1e-4", "--trait2-prior", "1e-4", "--minimum-overlap", "0.5");
                List<String> primary = new ArrayList<>(List.of("--trait1", directory.resolve("exposure.effects.tsv").toString(), "--trait2", directory.resolve("outcome.effects.tsv").toString(),
                    "--shared-prior", "5e-6", "--out", directory.resolve("colocalization.tsv").toString())); primary.addAll(priors);
                List<String> sensitivity = new ArrayList<>(List.of("--trait1", directory.resolve("exposure.effects.tsv").toString(), "--trait2", directory.resolve("outcome.effects.tsv").toString(),
                    "--shared-prior", "1e-6", "--out", directory.resolve("colocalization-prior-sensitivity.tsv").toString())); sensitivity.addAll(priors);
                run("coloc", primary, directory.resolve("coloc-command.log"));
                run("coloc", sensitivity, directory.resolve("coloc-sensitivity-command.log"));
                record.put("colocalization", preview(directory.resolve("colocalization.tsv")));
                record.put("prior_sensitivity", preview(directory.resolve("colocalization-prior-sensitivity.tsv")));
                var effective = new LinkedHashMap<>(map(record.get("effective_parameters")));
                effective.put("colocalization", Map.of("trait1_prior", 1e-4, "trait2_prior", 1e-4, "shared_prior", 5e-6,
                    "shared_prior_sensitivity", 1e-6, "minimum_component_overlap", 0.5, "posterior_decision_threshold", "Not applied; posterior probabilities are reported",
                    "signal_model", "SuSiE single-effect components", "partial_region_exploratory", partialColocalization,
                    "primary_command", primary, "sensitivity_command", sensitivity));
                record.put("effective_parameters", effective);
            }
        }
        record.put("directory", directory.getFileName().toString());
        record.put("citation_methods", operation.equals("mendelian_randomization")
            ? list(map(record.get("effective_parameters")).get("estimators")).size() > 2
                ? List.of("https://doi.org/10.1002/gepi.21758", "https://doi.org/10.1093/ije/dyv080", "https://doi.org/10.1002/gepi.21965")
                : List.of("https://doi.org/10.1002/gepi.21758")
            : operation.equals("colocalization") ? outcome.containsKey("gwas_catalog_accession")
                ? List.of("https://doi.org/10.1371/journal.pgen.1009440", "https://doi.org/10.1371/journal.pgen.1004383", GwasCatalogRegionalOutcome.DOCUMENTATION)
                : List.of("https://doi.org/10.1371/journal.pgen.1009440", "https://doi.org/10.1371/journal.pgen.1004383")
            : List.of("https://doi.org/10.1371/journal.pgen.1009440"));
    }

    private static void fine(List<Aligned> aligned, boolean outcome, int n, Path ld,
            Path directory, String name, Map<String, Object> record,
            boolean exploratoryConvergenceFallback) throws IOException {
        StringBuilder summary = new StringBuilder("variant_id\tz\n");
        for (Aligned a : aligned) summary.append(a.exposure().id()).append('\t').append(outcome ? a.outcomeSign() * a.outcome().beta() / a.outcome().se() : a.exposure().beta() / a.exposure().se()).append('\n');
        Path input = directory.resolve(name + "-z.tsv");
        Files.writeString(input, summary);
        Path output = directory.resolve(name + "-fine-mapping.tsv");
        Path effectsOutput = directory.resolve(name + ".effects.tsv");
        int requestedEffects = Math.min(10, aligned.size());
        int selectedEffects = 0;
        List<Object> attempts = new ArrayList<>();
        for (int effects = requestedEffects; effects >= 1; effects--) {
            Path attemptOutput = directory.resolve(name + "-fine-mapping-L" + effects + ".tsv");
            Path attemptEffects = directory.resolve(name + "-L" + effects + ".effects.tsv");
            Path attemptCommand = directory.resolve(name + "-command-L" + effects + ".log");
            run("susie", List.of("--summary", input.toString(), "--ld", ld.toString(), "--sample-size", Integer.toString(n),
                "--z-column", "z", "--effects", Integer.toString(effects), "--max-iterations", "2000", "--backend", "cpu",
                "--out", attemptOutput.toString(), "--effects-output", attemptEffects.toString()), attemptCommand);
            boolean converged = Files.readString(Path.of(attemptOutput + ".log")).contains("converged=true");
            attempts.add(Map.of("effects", effects, "converged", converged,
                "command_log", attemptCommand.getFileName().toString(),
                "backend_log", attemptOutput.getFileName() + ".log"));
            if (converged) {
                selectedEffects = effects;
                Files.copy(attemptOutput, output, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                Files.copy(attemptEffects, effectsOutput, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                Files.copy(Path.of(attemptOutput + ".log"), Path.of(output + ".log"), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                Files.copy(attemptCommand, directory.resolve(name + "-command.log"), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                break;
            }
            if (!exploratoryConvergenceFallback) break;
        }
        if (selectedEffects == 0) throw new IOException("SuSiE did not converge for " + name);
        record.put(name + "_fine_mapping", preview(output));
        var effective = new LinkedHashMap<>(map(record.get("effective_parameters")));
        Map<String, Object> fineParameters = new LinkedHashMap<>();
        fineParameters.put("sample_size", n);
        fineParameters.put("requested_maximum_effects", requestedEffects);
        fineParameters.put("selected_effects", selectedEffects);
        fineParameters.put("maximum_iterations_per_attempt", 2000);
        fineParameters.put("backend", "cpu");
        fineParameters.put("convergence_required", true);
        fineParameters.put("exploratory_largest_convergent_effect_count_fallback", exploratoryConvergenceFallback);
        fineParameters.put("fallback_used", selectedEffects < requestedEffects);
        fineParameters.put("attempts", attempts);
        fineParameters.put("command_log", name + "-command.log");
        fineParameters.put("effective_backend_log", name + "-fine-mapping.tsv.log");
        effective.put(name + "_fine_mapping", fineParameters);
        record.put("effective_parameters", effective);
    }

    private static boolean hasDataRows(Path path) throws IOException {
        try (var reader = Files.newBufferedReader(path)) {
            if (reader.readLine() == null) throw new IOException("Missing table header: " + path);
            return reader.readLine() != null;
        }
    }

    private static Path[] automaticLd(Context context, Map<String, Object> ld,
            Map<String, Variant> x, Map<String, Variant> y, Path directory,
            int minimumVariants) throws IOException {
        List<String> requested = x.keySet().stream().filter(y::containsKey).sorted().toList();
        if (requested.isEmpty()) throw new IOException("No overlapping SNPs for LD extraction");
        if (ThousandGenomesGrch38Ld.DATABASE.equals(ld.get("database"))) {
            List<Variant> variants = requested.stream().map(x::get).toList();
            ThousandGenomesGrch38Ld.Result result = ThousandGenomesGrch38Ld.extract(
                context, required(ld, "population"), variants, directory,
                minimumVariants);
            return new Path[] {result.matrix(), result.alleles()};
        }
        StudyResourceManager.LdReference installed;
        try (PrintStream log = new PrintStream(Files.newOutputStream(directory.resolve("ld-download.log")), true, StandardCharsets.UTF_8)) {
            try { installed = StudyResourceManager.ld(context.project(), context.study(), log); }
            catch (InterruptedException exception) { Thread.currentThread().interrupt(); throw new IOException("LD acquisition interrupted", exception); }
        }
        if (installed == null) throw new IOException("1000G reference unavailable under current download/network policy; install it locally or resume online");
        Path prefix = org.jlinalg.genetics.LdReferenceLayout.panelPrefix(installed.directory(), required(ld, "population"));
        var signed = org.jlinalg.genetics.PlinkBedLdClumper.signedCorrelation(prefix, requested);
        double[][] values = signed.correlation();
        Path matrix = directory.resolve("ld-source.tsv"), alleles = directory.resolve("ld-alleles.tsv");
        StringBuilder rows = new StringBuilder("variant_id\t" + String.join("\t", requested) + "\n");
        StringBuilder orientation = new StringBuilder("variant_id\teffect_allele\tother_allele\n");
        for (int i = 0; i < requested.size(); i++) {
            Variant v = x.get(requested.get(i));
            if (!signed.chromosomes().get(i).replaceFirst("^chr", "").equals(v.chromosome()) || signed.positions().get(i) != v.position())
                throw new IOException("LD reference position mismatch for " + v.id());
            rows.append(v.id());
            for (double correlation : values[i]) rows.append('\t').append(correlation);
            rows.append('\n');
            orientation.append(v.id()).append('\t').append(signed.effectAlleles().get(i)).append('\t').append(signed.otherAlleles().get(i)).append('\n');
        }
        Files.writeString(matrix, rows); Files.writeString(alleles, orientation);
        Map<String, Object> hashes = new LinkedHashMap<>();
        for (String suffix : List.of(".bed", ".bim", ".fam")) hashes.put(prefix + suffix, FollowupSupport.hash(Path.of(prefix + suffix)));
        write(directory.resolve("ld-reference.yaml"), Map.of("database", installed.database(), "population", ld.get("population"),
            "founders", signed.founders(), "source_sha256", hashes, "method", "Pearson A1 dosage correlation among founders; missing genotypes mean-imputed after a 95% call-rate check",
            "citation", "https://doi.org/10.1038/nature15393", "genome_build", "GRCh37"));
        return new Path[] {matrix, alleles};
    }

    private static void mr(Map<String, Object> pair, List<Aligned> aligned, double[][] ld, Path directory, Map<String, Object> record) throws IOException {
        required(pair, "instrument_validity_rationale");
        double threshold = number(pair.getOrDefault("instrument_p", 5e-8));
        double minF = number(pair.getOrDefault("minimum_f", 10));
        double r2 = number(pair.getOrDefault("clump_r2", 0.001));
        if (!(threshold > 0 && threshold <= 5e-8) || minF < 10 || !(r2 > 0 && r2 <= 0.01)) throw new IllegalArgumentException("MR thresholds must be conservative: p<=5e-8, F>=10, clump r2<=0.01");
        List<Integer> order = java.util.stream.IntStream.range(0, aligned.size()).boxed()
            .sorted(Comparator.comparingDouble((Integer i) -> aligned.get(i).exposure().p()).thenComparing(i -> aligned.get(i).exposure().id())).toList();
        List<Integer> selected = new ArrayList<>();
        List<Object> decisions = new ArrayList<>();
        for (int i : order) {
            Variant v = aligned.get(i).exposure();
            double f = Math.pow(v.beta() / v.se(), 2);
            String reason = v.p() > threshold ? "above_instrument_p" : f < minF ? "weak_instrument" : "retained";
            if (reason.equals("retained") && selected.stream().anyMatch(j -> ld[i][j] * ld[i][j] > r2)) reason = "ld_clumped";
            decisions.add(Map.of("variant_id", v.id(), "F", f, "exposure_p", v.p(), "decision", reason));
            if (reason.equals("retained")) selected.add(i);
        }
        record.put("instrument_selection", decisions);
        if (selected.isEmpty()) throw new NotEstimableException(
            "No independent strong instrument remains after the prespecified p-value, strength and LD filters",
            "The prespecified cis-irQTL family was exhausted after significance, F-statistic and LD-clumping filters. This transcript is completed but not estimable; thresholds were not relaxed post hoc.");
        Path input = directory.resolve("mr-instruments.tsv");
        StringBuilder table = new StringBuilder("variant_id\tbeta_exposure\tse_exposure\tbeta_outcome\tse_outcome\n");
        for (int i : selected) {
            Aligned a = aligned.get(i);
            table.append(a.exposure().id()).append('\t').append(a.exposure().beta()).append('\t').append(a.exposure().se()).append('\t')
                .append(a.outcome().beta() * a.outcomeSign()).append('\t').append(a.outcome().se()).append('\n');
        }
        Files.writeString(input, table);
        List<String> methods = selected.size() == 1
            ? new ArrayList<>(List.of("wald-ratio"))
            : new ArrayList<>(List.of("ivw-fixed", "ivw-random"));
        if (selected.size() >= 3) methods.addAll(List.of("egger", "weighted-median"));
        else if (selected.size() == 1) record.put("sensitivity_limitation",
            "Only one independent strong instrument remained; a Wald ratio is reported. Heterogeneity, MR-Egger, weighted median, pleiotropy and leave-one-instrument-out diagnostics are not estimable.");
        else record.put("sensitivity_limitation", "Fewer than three instruments; Egger and weighted median not attempted");
        List<Object> estimates = new ArrayList<>();
        for (String method : methods) {
            Path output = directory.resolve("mr-" + method + ".tsv");
            run("mr-estimate", List.of("--input", input.toString(), "--method", method, "--seed", "2026", "--output", output.toString()), directory.resolve(method + "-command.log"));
            estimates.addAll(preview(output));
        }
        record.put("mr_estimates", estimates);
        record.put("effective_parameters", Map.of("instrument_p", threshold, "minimum_f", minF, "clump_r2", r2,
            "clumping_scope", "Supplied aligned regional LD matrix; no separate sliding window", "estimators", methods,
            "seed", 2026, "leave_one_instrument_out", selected.size() >= 3, "parameter_source", "Exact adapter options; estimator-specific settings are retained in command/output logs"));
        // Leave-one-instrument-out estimates use the same harmonized, clumped family.
        if (selected.size() >= 3) for (int omitted = 0; omitted < selected.size(); omitted++) {
            List<String> lines = Files.readAllLines(input);
            lines.remove(omitted + 1);
            Path leave = directory.resolve("mr-loo-" + omitted + ".tsv");
            Files.write(leave, lines);
            run("mr-estimate", List.of("--input", leave.toString(), "--method", "ivw-fixed", "--output", directory.resolve("mr-loo-" + omitted + "-estimate.tsv").toString()), directory.resolve("loo-" + omitted + "-command.log"));
        }
    }

    static Path source(Context context, Map<String, Object> source, Path directory, String name, String rationale) throws IOException {
        required(source, "citation"); required(source, "source_id");
        if (rationale == null || rationale.isBlank()) throw new IllegalArgumentException("Genetic evidence acquisition requires a scientific rationale");
        Map<String, Object> intent = new LinkedHashMap<>();
        intent.put("scientific_rationale", rationale); intent.put("source", source);
        intent.put("hypothesis_id", context.operation().get("hypothesis")); intent.put("operation", context.operation().get("type"));
        intent.put("started_at", java.time.Instant.now().toString()); intent.put("request_attempted", false);
        write(directory.resolve(name + ".intent.yaml"), intent);
        if (source.containsKey("prepared_path")) {
            Path input = Path.of(required(source, "prepared_path")).toAbsolutePath().normalize();
            Path trusted = context.outputDirectory().toAbsolutePath().normalize();
            if (!input.startsWith(trusted) || !Files.isRegularFile(input))
                throw new IOException("Prepared MR sources must be regular files inside this operation output directory");
            Path copy = directory.resolve(name);
            Files.copy(input, copy);
            write(directory.resolve(name + ".source.yaml"), Map.of("source", source,
                "source_path", input.toString(), "sha256", FollowupSupport.hash(copy),
                "retrieval", "automatically prepared and audited public/local source"));
            return copy;
        }
        if (source.containsKey("url")) {
            required(source, "license");
            String sha = required(source, "sha256").toLowerCase(java.util.Locale.ROOT);
            if (!sha.matches("[a-f0-9]{64}")) throw new IllegalArgumentException("Remote regional sources require a pinned SHA-256");
            java.net.URI uri = java.net.URI.create(required(source, "url"));
            StudyPublicEvidence.validate(uri);
            Path cached = StudyResourceManager.resourceRoot(context.project(), context.study()).resolve("regional-evidence").resolve(sha + ".tsv");
            if (!Files.isRegularFile(cached)) {
                String blocked = StudyPublicEvidence.blocked(context.study());
                if (!blocked.isEmpty()) throw new IOException("Regional acquisition not sent: " + blocked);
                Files.createDirectories(cached.getParent());
                Path stage = Files.createTempFile(cached.getParent(), ".download-", ".part");
                Map<String, Object> acquisition = new LinkedHashMap<>();
                acquisition.put("url", uri.toString()); acquisition.put("started_at", java.time.Instant.now().toString());
                acquisition.put("expected_sha256", sha); acquisition.put("request_attempted", true);
                acquisition.put("scientific_rationale", rationale); acquisition.put("citation", source.get("citation"));
                write(directory.resolve(name + ".acquisition.yaml"), acquisition);
                try {
                    PublicResourceDownloader.downloadStrict(uri, stage, 512L * 1024 * 1024);
                    String actual = FollowupSupport.hash(stage);
                    acquisition.put("actual_sha256", actual);
                    if (!actual.equals(sha)) throw new IOException("Remote regional source checksum mismatch");
                    Files.move(stage, cached);
                    acquisition.put("status", "downloaded_and_verified");
                } catch (InterruptedException exception) {
                    Thread.currentThread().interrupt();
                    throw new IOException("Regional acquisition interrupted", exception);
                } finally {
                    acquisition.put("finished_at", java.time.Instant.now().toString());
                    write(directory.resolve(name + ".acquisition.yaml"), acquisition);
                    Files.deleteIfExists(stage);
                }
            }
            if (!FollowupSupport.hash(cached).equals(sha)) throw new IOException("Cached regional source checksum mismatch");
            Path copy = directory.resolve(name);
            Files.copy(cached, copy);
            write(directory.resolve(name + ".source.yaml"), Map.of("source", source, "sha256", sha, "retrieval", "checksum-pinned public source cache"));
            return copy;
        }
        String datasetId = required(source, "dataset_id"), role = required(source, "file_role");
        Map<String, Object> dataset = StudyBioinformatics.datasets(context).stream().filter(item -> datasetId.equals(item.get("id"))).findFirst()
            .orElseThrow(() -> new IllegalArgumentException("Evidence dataset must be an explicit hypothesis input: " + datasetId));
        String path = required(map(dataset.get("files")), role);
        Path input = context.project().resolve(path).normalize();
        if (!Files.isRegularFile(input)) throw new IOException("Source file is missing for " + datasetId + ":" + role);
        if (Files.size(input) > 512L * 1024 * 1024) throw new IOException("Supply a regional subset smaller than 512 MiB with its acquisition provenance");
        Path copy = directory.resolve(input.getFileName().toString().toLowerCase(java.util.Locale.ROOT).endsWith(".csv") ? name.replace(".tsv", ".csv") : name);
        Files.copy(input, copy);
        write(directory.resolve(name + ".source.yaml"), Map.of("source", source, "source_path", input.toString(), "sha256", FollowupSupport.hash(copy)));
        return copy;
    }

    private static Map<String, Variant> variants(Path path, String chromosome, int start,
            int end, String feature, Set<String> allowedIds) throws IOException {
        DelimitedData data = DelimitedData.read(path);
        boolean probability = data.header().contains("p_value");
        boolean logarithm = data.header().contains("log10_p");
        if (!probability && !logarithm) throw new IOException("Regional adapter requires p_value or log10_p");
        Map<String, Variant> result = new LinkedHashMap<>();
        for (String[] row : data.rows()) {
            if (feature != null && data.header().contains("feature_id") && !feature.equals(row[data.column("feature_id")] )) continue;
            String chr = row[data.column("chromosome")].replaceFirst("^chr", "");
            long pos;
            try { pos = Long.parseLong(row[data.column("position")]); } catch (NumberFormatException exception) { throw new IOException("Variant positions must be 1-based integers", exception); }
            if (!chr.equals(chromosome) || pos < start || pos > end) continue;
            String id = row[data.column("variant_id")];
            String[] idParts = id.split(":", -1);
            if (idParts.length > 1 && idParts[idParts.length - 1].matches("rs[1-9][0-9]*")) id = idParts[idParts.length - 1];
            if (allowedIds != null && !allowedIds.contains(id)) continue;
            String ea = row[data.column("effect_allele")].toUpperCase(java.util.Locale.ROOT), oa = row[data.column("other_allele")].toUpperCase(java.util.Locale.ROOT);
            double b = number(row[data.column("beta")]), se = number(row[data.column("se")]);
            double p;
            if (probability) p = number(row[data.column("p_value")]);
            else {
                double logp = number(row[data.column("log10_p")]);
                if (logp > 0) throw new IOException("log10_p must be log10(p) and therefore <=0");
                p = Math.pow(10, logp);
            }
            if (id.isBlank() || !ea.matches("[ACGT]") || !oa.matches("[ACGT]") || ea.equals(oa) || se <= 0 || p < 0 || p > 1)
                throw new IOException("Regional adapter requires unique biallelic SNPs, finite beta, positive SE and p in [0,1]");
            if (result.put(id, new Variant(id, chr, pos, ea, oa, b, se, p)) != null) throw new IOException("Duplicate variant " + id);
            if (result.size() > 5000) throw new IOException("Region exceeds 5000-variant in-memory LD budget");
        }
        if (result.isEmpty()) throw new IOException("Source has no variants in the declared region");
        return result;
    }

    static int sign(String ea, String oa, String otherEa, String otherOa) {
        if (ea.equals(otherEa) && oa.equals(otherOa)) return 1;
        if (ea.equals(otherOa) && oa.equals(otherEa)) return -1;
        String complementEa = complement(otherEa), complementOa = complement(otherOa);
        if (ea.equals(complementEa) && oa.equals(complementOa)) return 1;
        if (ea.equals(complementOa) && oa.equals(complementEa)) return -1;
        return 0;
    }
    private static String complement(String allele) { return switch (allele) { case "A" -> "T"; case "T" -> "A"; case "C" -> "G"; case "G" -> "C"; default -> ""; }; }
    private static boolean palindrome(String a, String b) { return complement(a).equals(b); }
    static double[][] matrix(Path path, List<Aligned> aligned,
            boolean verifyPositiveSemidefinite) throws IOException {
        DelimitedData table = DelimitedData.read(path);
        int n = table.header().size() - 1;
        if (n > 5000 || table.rows().size() != n) throw new IOException("LD matrix must be square and contain at most 5000 variants");
        Map<String, String[]> rows = new LinkedHashMap<>();
        for (String[] row : table.rows()) if (rows.put(row[0], row) != null) throw new IOException("Duplicate LD matrix row");
        Set<String> columns = new LinkedHashSet<>(table.header().subList(1, table.header().size()));
        if (columns.size() != n || !rows.keySet().equals(columns)) throw new IOException("LD row and column IDs must match uniquely");
        double[][] result = new double[aligned.size()][aligned.size()];
        for (int i = 0; i < aligned.size(); i++) for (int j = 0; j < aligned.size(); j++) {
            Aligned a = aligned.get(i), b = aligned.get(j);
            if (!rows.containsKey(a.exposure().id()) || !columns.contains(b.exposure().id())) throw new IOException("Missing aligned variant in LD matrix");
            result[i][j] = number(rows.get(a.exposure().id())[table.column(b.exposure().id())]) * a.ldSign() * b.ldSign();
        }
        for (int i = 0; i < result.length; i++) for (int j = 0; j < result.length; j++) {
            if (Math.abs(result[i][j]) > 1 || Math.abs(result[i][j] - result[j][i]) > 1e-10 || (i == j && Math.abs(result[i][i] - 1) > 1e-10)) throw new IOException("LD must be symmetric, bounded and have unit diagonal");
        }
        if (!verifyPositiveSemidefinite) return result;
        // Reference LD can legitimately be rank deficient (more SNPs than founders).
        // Check positive semidefiniteness without silently replacing or ridging it.
        double[][] l = new double[result.length][result.length];
        for (int i = 0; i < result.length; i++) for (int j = 0; j <= i; j++) {
            double v = result[i][j];
            for (int k = 0; k < j; k++) v -= l[i][k] * l[j][k];
            if (i == j) { if (v < -1e-8) throw new IOException("LD matrix is not positive semidefinite; review variants/reference"); l[i][j] = Math.sqrt(Math.max(0, v)); }
            else if (l[j][j] > 1e-8) l[i][j] = v / l[j][j];
            else if (Math.abs(v) > 1e-8) throw new IOException("LD matrix is not positive semidefinite; review variants/reference");
        }
        return result;
    }
    private static void writeMatrix(Path path, double[][] matrix) throws IOException {
        try (var writer = Files.newBufferedWriter(path)) { for (double[] row : matrix) {
            for (int j = 0; j < row.length; j++) { if (j > 0) writer.write('\t'); writer.write(Double.toString(row[j])); } writer.write('\n');
        } }
    }
    private static void run(String command, List<String> args, Path log) throws IOException {
        ByteArrayOutputStream messages = new ByteArrayOutputStream();
        int exit;
        try (PrintStream stream = new PrintStream(messages, true, StandardCharsets.UTF_8)) {
            stream.println("command=" + command); stream.println("arguments=" + args);
            exit = command.equals("mr-estimate") ? MrEstimatorCli.run(args.toArray(String[]::new), stream, stream)
                : FineMappingCli.run(command, args.toArray(String[]::new), stream, stream);
        }
        Files.writeString(log, messages.toString(StandardCharsets.UTF_8));
        if (exit != 0) throw new IOException(command + " could not complete; inspect " + log.getFileName());
    }
    private static List<Object> preview(Path path) throws IOException {
        DelimitedData table = DelimitedData.read(path);
        List<Object> rows = new ArrayList<>();
        for (String[] values : table.rows().stream().limit(20).toList()) {
            Map<String, Object> row = new LinkedHashMap<>();
            for (int i = 0; i < table.header().size(); i++) row.put(table.header().get(i), values[i]);
            rows.add(row);
        }
        return rows;
    }
    static String required(Map<String, Object> value, String key) {
        String result = text(value.get(key));
        if (result.isBlank()) throw new IllegalArgumentException("Missing required declaration: " + key);
        return result;
    }
    private static double number(Object value) {
        double result = Double.parseDouble(text(value));
        if (!Double.isFinite(result)) throw new IllegalArgumentException("Nonfinite numeric input");
        return result;
    }
}
