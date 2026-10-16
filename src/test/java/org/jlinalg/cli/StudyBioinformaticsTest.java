package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.jlinalg.study.StudyOperationExecutor.State;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.LoaderOptions;
import org.yaml.snakeyaml.constructor.SafeConstructor;

class StudyBioinformaticsTest {
    @TempDir Path project;

    @Test void sharedProviderMetadataIsWrittenWithoutYamlAliases() throws Exception {
        Map<String, Object> cached = Map.of("genes", List.of(Map.of("name", "PCBP2")));
        Map<String, Object> document = Map.of("evidence", List.of(
            Map.of("feature_id", "isoform-1", "metadata", cached),
            Map.of("feature_id", "isoform-2", "metadata", cached)));
        Path path = project.resolve("shared-evidence.yaml");
        StudyPublicEvidence.write(path, document);
        LoaderOptions options = new LoaderOptions();
        options.setMaxAliasesForCollections(0);
        assertEquals(document, new Yaml(new SafeConstructor(options)).load(Files.readString(path)));
    }

    @Test void publicInvestigationRecordsProteinModelsReferenceSequenceAndQueries() throws Exception {
        List<String> requests = new ArrayList<>();
        var client = new StudyPublicEvidence(uri -> {
            requests.add(uri.toString());
            return new StudyPublicEvidence.Reply(200, switch (uri.getHost()) {
                case "rest.uniprot.org" -> "{\"results\":[{\"primaryAccession\":\"P04637\"}]}";
                case "alphafold.ebi.ac.uk" -> "[{\"entryId\":\"AF-P04637-F1\",\"globalMetricValue\":75.06}]";
                case "api.genome.ucsc.edu" -> "{\"dna\":\"ACGT\"}";
                default -> throw new AssertionError(uri);
            });
        });
        Context context = context("bioinformatics_investigation", Map.of());
        var result = new StudyBioinformatics(client).execute(context);
        assertEquals(State.COMPLETE, result.state());
        assertEquals(3, requests.size());
        assertTrue(requests.get(2).endsWith("start=99;end=103"));
        String evidence = Files.readString(context.outputDirectory().resolve("bioinformatics-investigation.yaml"));
        assertTrue(evidence.contains("gene_level_match"));
        assertTrue(evidence.contains("models_available"));
        assertTrue(evidence.contains("sequence_sha256"));
        assertTrue(Files.readString(context.outputDirectory().resolve("queries/0001/request.yaml")).contains("response_sha256"));
        assertTrue(requests.stream().noneMatch(query -> query.contains("PRIVATE_PATIENT")));
    }

    @Test void denyAllAndOfflineNeverInvokeTransport() throws Exception {
        for (Map<String, Object> settings : List.of(
                Map.<String, Object>of("project", Map.of("data_policy", Map.of("remote_transmission", "deny_all"))),
                Map.<String, Object>of("resources", Map.of("network", "offline")),
                Map.<String, Object>of("project", Map.of("data_policy", Map.of("llm", Map.of("require_human_approval", true)))))) {
            var client = new StudyPublicEvidence(uri -> { throw new AssertionError("Unexpected network access"); });
            var result = new StudyBioinformatics(client).execute(context("bioinformatics_investigation", settings));
            assertEquals(State.WAITING_FOR_INPUT, result.state());
            assertEquals(true, result.details().get("retryable"));
        }
    }

    @Test void retiredServiceIsNotReportedAsNoQtlHits() throws Exception {
        var client = new StudyPublicEvidence(uri -> {
            if (uri.getPath().contains("eqtl")) return new StudyPublicEvidence.Reply(410, "Gone");
            if (uri.getHost().equals("gtexportal.org")) return new StudyPublicEvidence.Reply(200, "{\"data\":[]}");
            return new StudyPublicEvidence.Reply(200, uri.getPath().contains("associations")
                ? "{\"_embedded\":{\"associations\":[]}}" : "{\"_embedded\":{\"studies\":[]}}");
        });
        Context context = context("source_discovery", Map.of());
        var result = new StudyBioinformatics(client).execute(context);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        String evidence = Files.readString(context.outputDirectory().resolve("source-discovery.yaml"));
        assertTrue(evidence.contains("endpoint_retired"));
        assertTrue(evidence.contains("no_matches"));
        assertTrue(evidence.contains("genetic-source-selection"));
        assertTrue(evidence.contains("not evidence of no scientific matches"));
    }

    @Test void literatureNoMatchesIsACompletedBoundedSearch() throws Exception {
        Context context = context("literature_search", Map.of());
        var result = new StudyBioinformatics(new StudyPublicEvidence(uri -> new StudyPublicEvidence.Reply(200,
            "{\"hitCount\":0,\"resultList\":{\"result\":[]}}"))).execute(context);
        assertEquals(State.COMPLETE, result.state());
        assertTrue(Files.readString(context.outputDirectory().resolve("literature-search.yaml")).contains("no_matches"));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("reference-catalog.html")));
        assertTrue(result.details().toString().contains("bounded_coverage_complete_scientific_review_required"));
    }

    @Test void fieldSearchDoesNotDependOnSignificantGenesAndCatalogRetainsProvenance() throws Exception {
        Context context = context("literature_search", Map.of());
        Path annotation = context.dependencyDirectories().get("harmonize").resolve("feature-annotation.tsv");
        Files.writeString(annotation, Files.readString(annotation).replace("located", "ambiguous"));
        List<String> queries = new ArrayList<>();
        var result = new StudyBioinformatics(new StudyPublicEvidence(uri -> {
            queries.add(java.net.URLDecoder.decode(uri.getQuery(), java.nio.charset.StandardCharsets.UTF_8));
            return new StudyPublicEvidence.Reply(200, """
                {"hitCount":1,"resultList":{"result":[{"doi":"10.1000/fixture","pmid":"123","title":"Fixture paper","authorString":"A Author","pubYear":"2025","abstractText":"A source abstract for field coverage with enough words to excerpt safely and reproducibly."}]}}
                """);
        })).execute(context);
        assertEquals(State.COMPLETE, result.state()); assertEquals(4, queries.size());
        assertTrue(queries.stream().allMatch(q -> q.contains("transcriptomics") && q.contains("body mass index") && !q.contains("PRIVATE_PATIENT")));
        var catalog = StudyPublicEvidence.map(new Yaml().load(Files.readString(context.outputDirectory().resolve("reference-catalog.yaml"))));
        assertEquals(1, StudyPublicEvidence.list(catalog.get("references")).size());
        var ref = StudyPublicEvidence.map(StudyPublicEvidence.list(catalog.get("references")).get(0));
        assertEquals(4, StudyPublicEvidence.list(ref.get("retrieval_context")).size());
        assertTrue(ref.toString().contains("response_sha256"));
        assertTrue(ref.toString().contains("started_at"));
        assertTrue(ref.toString().contains("unreviewed_candidate"));
    }

    @Test void unknownFieldRequiresInputAndPolicyBlocksFieldSearchesAsWellAsGeneSearches() throws Exception {
        var unknown = context("literature_search", Map.of("datasets", List.of(Map.of("id", "d1", "type", "unknown"))));
        var result = new StudyBioinformatics(new StudyPublicEvidence(uri -> new StudyPublicEvidence.Reply(200,
            "{\"hitCount\":0,\"resultList\":{\"result\":[]}}"))).execute(unknown);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(result.details().toString().contains("literature_terms"));
        for (Map<String, Object> setting : List.of(Map.<String, Object>of("resources", Map.of("network", "offline")),
                Map.<String, Object>of("project", Map.of("data_policy", Map.of("remote_transmission", "deny_all"))))) {
            var blocked = new StudyBioinformatics(new StudyPublicEvidence(uri -> { throw new AssertionError("No transport permitted"); }))
                .execute(context("literature_search", setting));
            assertEquals(State.WAITING_FOR_INPUT, blocked.state());
            assertTrue(blocked.details().toString().contains("incomplete_connection_or_policy"));
        }
    }

    @Test void explicitFieldTermsAreUsedAndTruncatedSearchesAreVisible() throws Exception {
        Context context = context("literature_search", Map.of("hypotheses", List.of(Map.of("id", "h1", "outcome", "sleep",
            "analysis", Map.of("bioinformatics", Map.of("literature_terms", List.of("isoform ratio"), "max_pages", 1))))));
        List<String> queries = new ArrayList<>();
        var result = new StudyBioinformatics(new StudyPublicEvidence(uri -> {
            queries.add(java.net.URLDecoder.decode(uri.getQuery(), java.nio.charset.StandardCharsets.UTF_8));
            return new StudyPublicEvidence.Reply(200, "{\"hitCount\":21,\"nextCursorMark\":\"next\",\"resultList\":{\"result\":[{\"doi\":\"10.1000/test\"}]}}");
        })).execute(context);
        assertTrue(queries.get(0).contains("isoform ratio")); assertFalse(queries.get(0).contains("transcriptomics"));
        assertTrue(result.details().toString().contains("bounded_results_truncated"));
        assertTrue(result.details().toString().contains("additional_search_or_input_required"));
    }

    @Test void cohortAndMethodSearchesHaveSeparateRationalesAndDoNotSendPrivateDescriptions() throws Exception {
        Context ctx = context("literature_search", Map.of("datasets", List.of(Map.of("id", "d1", "type", "transcriptomics",
            "cohorts", List.of("fhs", "cardia", "whi", "jhs"), "analysis_description", Map.of("references", List.of("10.1000/fixture"), "model", "PRIVATE_MODEL_NOT_SENT"))),
            "cohorts", List.of(Map.of("id", "fhs", "description", "PRIVATE_COHORT_NOT_SENT", "references", List.of("pmid:123")))));
        List<String> queries = new ArrayList<>();
        var result = new StudyBioinformatics(new StudyPublicEvidence(uri -> {
            queries.add(java.net.URLDecoder.decode(uri.getQuery(), java.nio.charset.StandardCharsets.UTF_8));
            return new StudyPublicEvidence.Reply(200, "{\"hitCount\":1,\"resultList\":{\"result\":[{\"pmid\":\"123\",\"doi\":\"10.1000/fixture\",\"title\":\"Fixture source\"}]}}");
        })).execute(ctx);
        assertEquals(State.COMPLETE, result.state());
        for (String name : List.of("Framingham Heart Study", "Coronary Artery Risk Development in Young Adults", "Women's Health Initiative", "Jackson Heart Study"))
            assertTrue(queries.stream().anyMatch(q -> q.contains(name)));
        assertTrue(queries.stream().anyMatch(q -> q.contains("EXT_ID:123 AND SRC:MED")));
        assertTrue(queries.stream().anyMatch(q -> q.contains("DOI:\"10.1000/fixture\"")));
        assertTrue(queries.stream().noneMatch(q -> q.contains("PRIVATE_")));
        String catalog = Files.readString(ctx.outputDirectory().resolve("reference-catalog.yaml"));
        assertTrue(catalog.contains("methods_searches")); assertTrue(catalog.contains("cohort_id: fhs"));
        assertTrue(catalog.contains("bounded_coverage_complete_scientific_review_required"));
    }

    @Test void unknownCohortNeedsPublicNameAndDenyAllBlocksCohortLookups() throws Exception {
        var privateCohort = Map.<String, Object>of("datasets", List.of(Map.of("id", "d1", "cohort", "PRIVATE_COHORT", "type", "transcriptomics")));
        var result = new StudyBioinformatics(new StudyPublicEvidence(uri -> {
            assertFalse(uri.toString().contains("PRIVATE_COHORT"));
            return new StudyPublicEvidence.Reply(200, "{\"hitCount\":0,\"resultList\":{\"result\":[]}}");
        })).execute(context("literature_search", privateCohort));
        assertEquals(State.WAITING_FOR_INPUT, result.state()); assertTrue(result.details().toString().contains("literature_name"));
        var denied = Map.<String, Object>of("datasets", List.of(Map.of("id", "d1", "cohort", "fhs", "type", "transcriptomics")),
            "project", Map.of("data_policy", Map.of("remote_transmission", "deny_all")));
        assertEquals(State.WAITING_FOR_INPUT, new StudyBioinformatics(new StudyPublicEvidence(uri -> { throw new AssertionError("No cohort request allowed"); }))
            .execute(context("literature_search", denied)).state());
    }

    @Test void changedProviderSchemaCannotBecomeNegativeEvidence() throws Exception {
        assertThrows(java.io.IOException.class, () -> new StudyBioinformatics(new StudyPublicEvidence(uri ->
            new StudyPublicEvidence.Reply(200, "{\"unexpected\":[]}"))).execute(context("literature_search", Map.of())));
    }

    @Test void catalogHalEmptyPageIsNotAMalformedOrPositivePage() throws Exception {
        Object empty = SimpleJson.parse("{\"_links\":{\"self\":{\"href\":\"https://www.ebi.ac.uk/gwas/rest/api/v2/associations\"}},\"page\":{\"size\":100,\"totalElements\":0,\"totalPages\":0,\"number\":0}}");
        assertTrue(StudyBioinformatics.catalogItems(empty, "associations").isEmpty());
        for (String malformed : List.of("{}", "{\"page\":{\"totalElements\":1,\"totalPages\":1,\"number\":0}}",
                "{\"_embedded\":{\"unexpected\":[]},\"page\":{\"totalElements\":0,\"totalPages\":0,\"number\":0}}"))
            assertThrows(java.io.IOException.class, () -> StudyBioinformatics.catalogItems(SimpleJson.parse(malformed), "associations"));
    }

    @Test void humanGeneticResourcesAreNotSubstitutedForOtherSpecies() throws Exception {
        var client = new StudyPublicEvidence(uri -> { throw new AssertionError("Unexpected human source query"); });
        var result = new StudyBioinformatics(client).execute(context("source_discovery",
            Map.of("datasets", List.of(Map.of("id", "d1", "species", "mouse", "tissue", "brain")))));
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(client.requests().isEmpty());
        assertTrue(result.details().toString().contains("species-matched"));
    }

    @Test void coordinateWindowChecksBuildAndQueriesOnlyMatchingSnpPhenotype() throws Exception {
        List<String> requests = new ArrayList<>();
        var client = new StudyPublicEvidence(uri -> {
            requests.add(uri.toString());
            String body = uri.getPath().endsWith("metadata") ? "{\"gene_build\":\"GRCh38.p14\"}"
                : uri.getPath().contains("single-nucleotide") ? "{\"_embedded\":{\"snps\":[{\"rs_id\":\"rs123\",\"locations\":[{\"chromosome_name\":\"1\",\"chromosome_position\":102}]}]}}"
                : uri.getPath().contains("associations") ? "{\"_embedded\":{\"associations\":[]}}"
                : uri.getHost().equals("gtexportal.org") ? "{\"data\":[]}"
                : uri.getPath().contains("eqtl") ? "[]" : "{\"_embedded\":{\"studies\":[]}}";
            return new StudyPublicEvidence.Reply(200, body);
        });
        Context context = context("source_discovery", Map.of("hypotheses", List.of(Map.of(
            "id", "h1", "outcome", "body mass index", "analysis", Map.of("bioinformatics",
                Map.of("remote_feature_queries", true))))));
        new StudyBioinformatics(client).execute(context);
        assertTrue(requests.stream().anyMatch(url -> url.contains("bp_start=1&bp_end=500103")));
        assertTrue(requests.stream().anyMatch(url -> url.contains("rs_id=rs123&efo_trait=body+mass+index")));
        assertTrue(Files.readString(context.outputDirectory().resolve("source-discovery.yaml")).contains("catalogue_metadata_audit"));
    }

    @Test void defaultDiscoveryDoesNotTransmitSelectedFeatureCoordinates() throws Exception {
        List<String> requests = new ArrayList<>();
        var client = new StudyPublicEvidence(uri -> {
            requests.add(uri.toString());
            String body = uri.getPath().contains("associations") ? "{\"_embedded\":{\"associations\":[]}}"
                : uri.getHost().equals("gtexportal.org") ? "{\"data\":[]}"
                : uri.getPath().contains("eqtl") ? "[]" : "{\"_embedded\":{\"studies\":[]}}";
            return new StudyPublicEvidence.Reply(200, body);
        });
        Context context = context("source_discovery", Map.of("hypotheses", List.of(Map.of(
            "id", "h1", "outcome", "body mass index", "analysis", Map.of("bioinformatics",
                Map.of("colocalization", Map.of("exposure_ancestry", "European")))))));
        var result = new StudyBioinformatics(client).execute(context);
        assertEquals(State.COMPLETE, result.state());
        assertTrue(requests.stream().noneMatch(url -> url.contains("bp_start=") || url.contains("rs_id=")
            || url.contains("TP53") || url.contains("ENSG00000141510")));
        String discovery = Files.readString(context.outputDirectory().resolve("source-discovery.yaml"));
        assertTrue(discovery.contains("feature_specific_public_queries"));
        assertTrue(discovery.contains("not_requested"));
        assertFalse(discovery.contains("genetic-source-selection"));
    }

    @Test void sourceDiscoverySearchesAdvisorPhenotypesWithoutTreatingThemAsApprovedBiology() throws Exception {
        Context base = context("source_discovery", Map.of());
        Path advice = Files.createTempDirectory(project, "phenotype-advice-");
        StudyPublicEvidence.write(advice.resolve("decision.yaml"), Map.of("purpose", "related_phenotypes", "status", "completed_unreviewed",
            "phenotypes", List.of(Map.of("name", "glycated hemoglobin", "rationale", "Fixture", "literature", List.of()))));
        Map<String, Path> dependencies = new LinkedHashMap<>(base.dependencyDirectories()); dependencies.put("phenotype-advice", advice);
        Context context = new Context(base.project(), base.outputDirectory(), base.operation(), base.study(), base.plan(), dependencies);
        List<String> requests = new ArrayList<>();
        var client = new StudyPublicEvidence(uri -> {
            requests.add(java.net.URLDecoder.decode(uri.toString(), java.nio.charset.StandardCharsets.UTF_8));
            String body = uri.getPath().endsWith("metadata") ? "{\"gene_build\":\"GRCh38.p14\"}"
                : uri.getPath().contains("single-nucleotide") ? "{\"_embedded\":{\"snps\":[]}}"
                : uri.getPath().contains("associations") ? "{\"_embedded\":{\"associations\":[]}}"
                : uri.getHost().equals("gtexportal.org") ? "{\"data\":[]}"
                : uri.getPath().contains("eqtl") ? "[]" : "{\"_embedded\":{\"studies\":[]}}";
            return new StudyPublicEvidence.Reply(200, body);
        });
        var result = new StudyBioinformatics(client).execute(context);
        assertEquals(State.WAITING_FOR_INPUT, result.state());
        assertTrue(requests.stream().anyMatch(url -> url.contains("efo_trait=glycated hemoglobin")));
        String discovery = Files.readString(context.outputDirectory().resolve("source-discovery.yaml"));
        assertTrue(discovery.contains("llm_advisor_proposal_unreviewed"));
    }

    @Test void fullTestedFamilyControlsSelectionAndAmbiguityIsExcluded() throws Exception {
        Context context = context("literature_search", Map.of());
        var features = StudyBioinformatics.features(context);
        assertEquals(1, features.size());
        assertEquals(1e-7, features.get(0).fdr(), 1e-15);
        Path annotation = context.dependencyDirectories().get("harmonize").resolve("feature-annotation.tsv");
        Files.writeString(annotation, Files.readString(annotation).replace("located", "ambiguous"));
        assertTrue(StudyBioinformatics.features(context).isEmpty());
    }

    @Test void mrFeatureScopeIsIndependentOfContextualFeatureScope() throws Exception {
        Context context = context("literature_search", Map.of("hypotheses", List.of(Map.of(
            "id", "h1", "outcome", "body mass index", "analysis", Map.of(
                "bioinformatics", Map.of("max_features", 10, "mr", Map.of("max_features", 1000)))))));
        assertEquals(10, StudyBioinformatics.featureLimit(context));
        assertEquals(1000, StudyBioinformatics.mrFeatureLimit(context));
        assertEquals(1000, StudyBioinformatics.geneticFeatureLimit(context));
    }

    @Test void mrDefaultIncludesAllSignificantTranscriptsProteinsAndProbesBeyondOldLimits() throws Exception {
        Context context = context("mendelian_randomization", Map.of("hypotheses", List.of(Map.of("id", "h1", "analysis",
            Map.of("bioinformatics", Map.of("max_features", 1))))));
        significantFeatures(context, 10005);
        var findings = StudyBioinformatics.mrFeatures(context);
        assertEquals(10005, findings.size());
        assertEquals(1, StudyBioinformatics.features(context).size());
        assertEquals(Integer.MAX_VALUE, StudyBioinformatics.mrFeatureLimit(context));
        assertTrue(findings.stream().allMatch(feature -> feature.fdr() < .05));
        for (String kind : List.of("transcript", "protein", "probe"))
            assertTrue(findings.stream().anyMatch(feature -> feature.id().startsWith(kind)));
        assertEquals(findings.stream().map(StudyBioinformatics.Feature::id).sorted().toList(), findings.stream().map(StudyBioinformatics.Feature::id).toList());
        assertEquals("all_significant_features", StudyBioinformatics.mrFeatureSelection(context).get("scope"));
        assertEquals("unlimited", StudyBioinformatics.mrFeatureSelection(context).get("max_features"));
    }

    @Test void mrExplicitCapStillAppliesAfterFullFamilySelection() throws Exception {
        Context context = context("mendelian_randomization", Map.of("hypotheses", List.of(Map.of("id", "h1", "analysis",
            Map.of("bioinformatics", Map.of("max_features", 1, "mr", Map.of("max_features", 3)))))));
        significantFeatures(context, 30);
        assertEquals(3, StudyBioinformatics.mrFeatures(context).size());
        assertEquals(1, StudyBioinformatics.features(context).size());
        assertEquals(3, StudyBioinformatics.geneticFeatures(context).size());
        assertEquals("significant_features_with_explicit_cap", StudyBioinformatics.mrFeatureSelection(context).get("scope"));
        Path primary = context.dependencyDirectories().get("harmonize").resolve("standard-results/d1.tsv");
        Files.writeString(primary, "feature_id\tp_value\ntranscript-00000\t0.04\n" +
            java.util.stream.IntStream.range(0,99).mapToObj(i -> "background-"+i+"\t0.5\n").collect(java.util.stream.Collectors.joining()));
        assertTrue(StudyBioinformatics.mrFeatures(context).isEmpty(), "BH must not be recomputed within the capped or annotated subset");
    }

    @Test void mrUsesStrictConfiguredFdrThresholdAndFullTestedFamily() throws Exception {
        Context context = context("mendelian_randomization", Map.of());
        Path primary = context.dependencyDirectories().get("harmonize").resolve("standard-results/d1.tsv");
        Files.writeString(primary, "feature_id\tp_value\nbest\t0.05\n");
        assertTrue(StudyBioinformatics.mrFeatures(context).isEmpty());
        assertEquals(1, StudyBioinformatics.features(context).size(), "Contextual boundary behavior remains unchanged");
        Files.writeString(primary, "feature_id\tp_value\nbest\t0.049\n");
        assertEquals(1, StudyBioinformatics.mrFeatures(context).size());
        Context custom = context("mendelian_randomization", Map.of("hypotheses", List.of(Map.of("id", "h1", "analysis", Map.of("bioinformatics", Map.of("fdr", .1))))));
        Path customPrimary = custom.dependencyDirectories().get("harmonize").resolve("standard-results/d1.tsv");
        Files.writeString(customPrimary, Files.readString(customPrimary).replace("1e-9", "0.0008"));
        assertEquals(.08, StudyBioinformatics.mrFeatures(custom).get(0).fdr(), 1e-12);
        assertEquals(.1, StudyBioinformatics.mrFeatureSelection(custom).get("fdr_threshold"));
    }

    @Test void mrPrefersFullPooledMetaFamilyOverPrimaryOrIndividualCohortRows() throws Exception {
        Context base = context("mendelian_randomization", Map.of());
        Path meta = Files.createTempDirectory(project, "meta-");
        Files.writeString(meta.resolve("manifest.yaml"), "operation_type: meta_analysis\nstate: complete\n");
        String background = java.util.stream.IntStream.range(0,99).mapToObj(i -> "background-"+i+"\tpooled\tok\t0.5\n").collect(java.util.stream.Collectors.joining());
        Files.writeString(meta.resolve("meta-analysis.tsv"), "feature_id\tterm\tstatus\tp_value\nbest\tpooled\tok\t0.0004\nbest\tcohort\tok\t0.9\n"+background);
        var dependencies = new LinkedHashMap<>(base.dependencyDirectories()); dependencies.put("meta", meta);
        Context context = new Context(base.project(),base.outputDirectory(),base.operation(),base.study(),base.plan(),dependencies);
        assertEquals(.04, StudyBioinformatics.mrFeatures(context).get(0).fdr(), 1e-12);
        Files.writeString(meta.resolve("meta-analysis.tsv"), Files.readString(meta.resolve("meta-analysis.tsv")).replace("0.0004", "0.0008"));
        assertTrue(StudyBioinformatics.mrFeatures(context).isEmpty());
    }

    @Test void qtlReviewExpandsForDeclaredOrPlannedMrButNotColocalizationAlone() throws Exception {
        Context base = context("irqtl_candidate_review", Map.of()); significantFeatures(base,30);
        assertEquals(10, StudyBioinformatics.geneticFeatures(base).size());
        for (String goal : List.of("mendelian_randomization", "colocalization")) {
            Context planned = new Context(base.project(),base.outputDirectory(),base.operation(),base.study(),
                Map.of("operations",List.of(Map.of("hypothesis","h1","type",goal))),base.dependencyDirectories());
            assertEquals(goal.equals("mendelian_randomization") ? 30 : 10,StudyBioinformatics.geneticFeatures(planned).size());
        }
        Map<String,Object> study = new LinkedHashMap<>(base.study());
        study.put("hypotheses",List.of(Map.of("id","h1","analysis",Map.of("goals",List.of("mendelian_randomization")))));
        Context declared = new Context(base.project(),base.outputDirectory(),base.operation(),study,base.plan(),base.dependencyDirectories());
        assertEquals(30,StudyBioinformatics.geneticFeatures(declared).size());
    }

    @Test void localQtlCandidateReviewIncludesSignificantTranscriptsBeyondContextualCap() throws Exception {
        Context base=context("irqtl_candidate_review",Map.of("hypotheses",List.of(Map.of("id","h1","analysis",Map.of("goals",List.of("mendelian_randomization"))))));
        significantFeatures(base,30);
        Path source=project.resolve("local-cis.csv");
        Files.writeString(source,"SNP,SNP_Chr,SNP_Pos,EAF,N,Tx,Is_Cis,RSq,Fx,T,log10P\nchr1:123:G:A:rs1,1,123,0.2,1000,transcript-00027,1,.02,.1,20,-10\n");
        Map<String,Object> study=new LinkedHashMap<>(base.study());
        study.put("datasets",List.of(Map.of("id","d1","type","transcriptomics","files",Map.of("irqtl",source.toString()),
            "columns",Map.of("irqtl",Map.of("genome_build","GRCh38","effect_allele","ALT","complete_region",false,"export_p_threshold",1e-7,"rationale","Synthetic cis-only review")))));
        Context context=new Context(base.project(),base.outputDirectory(),base.operation(),study,base.plan(),base.dependencyDirectories());
        var result=new StudyIrQtlReview().execute(context);
        assertEquals(State.COMPLETE,result.state());
        assertEquals(30,((List<?>)result.details().get("selected_features")).size());
        assertTrue(Files.readString(context.outputDirectory().resolve("irqtl-1-candidates.tsv")).contains("transcript-00027"));
        assertEquals("all_significant_features",((Map<?,?>)result.details().get("mr_feature_selection")).get("scope"));
    }

    private void significantFeatures(Context context, int count) throws Exception {
        var primary = new StringBuilder("feature_id\tp_value\n");
        var annotation = new StringBuilder("feature_id\tgene_symbol\tgene_id\tchromosome\tstart\tend\tgenome_build\tmapping_status\n");
        for (int i=0;i<count;i++) {
            String id=List.of("transcript","protein","probe").get(i%3)+String.format(java.util.Locale.ROOT,"-%05d",i);
            primary.append(id).append("\t1e-9\n");
            annotation.append(id).append("\tTEST\tENSG1\t1\t100\t300\tGRCh38\tlocated\n");
        }
        Path directory=context.dependencyDirectories().get("harmonize");
        Files.writeString(directory.resolve("standard-results/d1.tsv"),primary);
        Files.writeString(directory.resolve("feature-annotation.tsv"),annotation);
    }

    @Test void sequenceBudgetDoesNotMasqueradeAsMissingAnnotation() throws Exception {
        Context base = context("bioinformatics_investigation", Map.of());
        Path annotation = base.dependencyDirectories().get("harmonize")
            .resolve("feature-annotation.tsv");
        Files.writeString(annotation, Files.readString(annotation)
            .replace("\t100\t103\t", "\t100\t200000\t"));
        var client = new StudyPublicEvidence(uri -> new StudyPublicEvidence.Reply(200,
            uri.getHost().equals("rest.uniprot.org")
                ? "{\"results\":[{\"primaryAccession\":\"P04637\"}]}"
                : "[{\"entryId\":\"AF-P04637-F1\",\"globalMetricValue\":75.06}]"));
        var result = new StudyBioinformatics(client).execute(base);
        assertEquals(State.COMPLETE, result.state());
        assertTrue(result.details().toString().contains("interval_too_large"));
        assertTrue(result.details().toString().contains("annotation is present"));
        assertFalse(result.details().toString().contains("bioinformatics:h1:annotation"));
    }

    @Test void untrustedUrlsAreRejectedBeforeTransportAndDuplicateQueriesAreCached() throws Exception {
        AtomicInteger calls = new AtomicInteger();
        var client = new StudyPublicEvidence(uri -> { calls.incrementAndGet(); return new StudyPublicEvidence.Reply(200, "{}"); });
        Context context = context("literature_search", Map.of());
        for (String url : List.of("http://www.ebi.ac.uk/data", "https://www.ebi.ac.uk.evil.test/data", "https://user:pass@www.ebi.ac.uk/data", "https://127.0.0.1/data", "https://www.ebi.ac.uk:8443/data"))
            assertThrows(IllegalArgumentException.class, () -> client.get(context, "test", url, "citation", "Check a public source"));
        assertThrows(IllegalArgumentException.class, () -> client.get(context, "test", "https://www.ebi.ac.uk/data", "citation", " "));
        assertEquals(0, calls.get());
        client.get(context, "test", "https://www.ebi.ac.uk/data", "citation", "First feature rationale");
        client.get(context, "test", "https://www.ebi.ac.uk/data", "citation", "Second feature rationale");
        assertEquals(1, calls.get());
        assertEquals(2, client.requests().size());
        String audit = Files.readString(context.outputDirectory().resolve("queries/0002/request.yaml"));
        Map<?, ?> saved = new Yaml().load(audit);
        assertEquals("Second feature rationale", saved.get("scientific_rationale"));
        assertTrue(saved.containsKey("cached_response_audit"));
        assertEquals(false, saved.get("request_attempted"));
    }

    private Context context(String operation, Map<String, Object> overrides) throws Exception {
        Path upstream = Files.createTempDirectory(project, "upstream-");
        Files.createDirectories(upstream.resolve("standard-results"));
        Files.writeString(upstream.resolve("manifest.yaml"), "operation_type: validate_harmonize\nstate: complete\n");
        StringBuilder results = new StringBuilder("feature_id\tp_value\nbest\t1e-9\n");
        for (int i = 0; i < 99; i++) results.append("background-").append(i).append("\t0.5\n");
        Files.writeString(upstream.resolve("standard-results/d1.tsv"), results);
        Files.writeString(upstream.resolve("feature-annotation.tsv"), "feature_id\tgene_symbol\tgene_id\tchromosome\tstart\tend\tgenome_build\tmapping_status\nbest\tTP53\tENSG00000141510\t1\t100\t103\tGRCh38\tlocated\n");
        Map<String, Object> study = new LinkedHashMap<>();
        study.put("datasets", List.of(Map.of("id", "d1", "type", "transcriptomics", "species", "human", "tissue", "Whole Blood", "metadata", "PRIVATE_PATIENT_SECRET")));
        study.put("hypotheses", List.of(Map.of("id", "h1", "outcome", "body mass index", "analysis", Map.of())));
        study.putAll(overrides);
        Path output = Files.createTempDirectory(project, "output-");
        return new Context(project, output, Map.of("id", "test", "hypothesis", "h1", "type", operation, "datasets", List.of("d1")), study, Map.of(), Map.of("harmonize", upstream));
    }
}
