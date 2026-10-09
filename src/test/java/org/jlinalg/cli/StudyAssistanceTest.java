/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import org.jlinalg.study.StudyOperationExecutor;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyAssistanceTest {
    @TempDir Path root;

    @Test void deniedAndOfflineNeverReadCredentialsOrSendRequests() throws Exception {
        StudyAssistance client = new StudyAssistance((uri, provider, key, body, timeout) -> {
            fail("network must not be called"); return null;
        }, key -> { fail("credentials must not be read"); return null; });
        var denied = context("denied", "advisor_review", Map.of("data_policy", Map.of("llm", Map.of("remote_transmission", "deny_all"))), Map.of());
        var result = client.execute(denied);
        assertEquals("policy_denied", result.details().get("status"));
        assertEquals(false, result.details().get("request_sent"));
        assertFalse(Files.exists(denied.outputDirectory().resolve("llm/interaction-001/request.json")));
        var offline = context("offline", "advisor_review", Map.of(), Map.of("network", "offline"));
        assertEquals("offline_policy", client.execute(offline).details().get("status"));
        var review = context("approval", "advisor_review", Map.of("data_policy", Map.of("llm", Map.of("require_human_approval", true))), Map.of());
        assertEquals("additional_disclosure_review_required", client.execute(review).details().get("status"));
    }

    @Test void missingCredentialsAndDeniedEndpointHaveAuditableDispositions() throws Exception {
        var context = context("credentials", "advisor_review", Map.of(), Map.of());
        var missing = new StudyAssistance((uri, provider, key, body, timeout) -> { fail("No key"); return null; }, key -> null).execute(context);
        assertEquals("credentials_missing", missing.details().get("status"));
        assertTrue(Files.readString(context.outputDirectory().resolve("llm/interaction-001/request.json")).contains("proposals"));
        var blocked = context("blocked", "advisor_review", Map.of(), Map.of());
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> new StudyAssistance.Reply(403, "{\"error\":\"test-secret forbidden\"}"), key -> "test-secret").execute(blocked);
        assertEquals("authorization_or_endpoint_denied", result.details().get("status"));
        String response = Files.readString(blocked.outputDirectory().resolve("llm/interaction-001/response.json"));
        assertTrue(response.contains("REDACTED_API_KEY"));
        assertFalse(response.contains("test-secret"));
        assertTrue(Files.readString(blocked.outputDirectory().resolve("llm/interaction-001/metadata.yaml")).contains("finished_at"));
    }

    @Test void proposalsCannotRunCodeAndPayloadOmitsPatientFieldsAndDictionaryText() throws Exception {
        AtomicInteger requests = new AtomicInteger();
        var context = context("advice", "advisor_review", Map.of("description", "PATIENT_SECRET"), Map.of());
        String advice = "{\"proposals\":[{\"action\":\"leave_one_cohort_out\",\"rationale\":\"Check robustness\",\"required_inputs\":[]},{\"action\":\"curl sensitive-data\",\"rationale\":\"untrusted\",\"required_inputs\":[]}]}";
        var client = new StudyAssistance((uri, provider, key, body, timeout) -> {
            requests.incrementAndGet();
            assertEquals("https://api.openai.com/v1/responses", uri.toString());
            assertFalse(body.contains("PATIENT_SECRET"));
            assertFalse(body.contains("PRIVATE_FEATURE"));
            assertFalse(body.contains("sample_id"));
            assertFalse(body.contains("description"));
            assertTrue(body.contains("BMI"));
            assertTrue(body.contains("GENE1"));
            assertTrue(body.contains("study_contract_context"));
            assertTrue(body.contains("ALL observed results"));
            assertTrue(body.contains("exact metadata field/dataset IDs"));
            assertTrue(body.contains("\"store\":false"));
            return openai(advice, List.of());
        }, key -> "test-secret");
        var result = client.execute(context);
        assertEquals(StudyOperationExecutor.State.COMPLETE, result.state());
        assertEquals(1, requests.get());
        String decision = Files.readString(context.outputDirectory().resolve("decision.yaml"));
        assertTrue(decision.contains("scheduled_by_approved_plan_subject_to_prerequisites"));
        assertTrue(decision.contains("unsupported_requires_scientist_review"));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("source-map.yaml")));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("study-context-review.yaml")));
    }

    @Test void malformedOrIncompleteProviderOutputIsNotSuccessfulAdvice() throws Exception {
        for (String answer : List.of("not JSON", "{}", "{\"proposals\":[{\"action\":\"leave_one_cohort_out\"}]}")) {
            var context = context("invalid-" + Math.abs(answer.hashCode()), "advisor_review", Map.of(), Map.of());
            assertEquals(StudyOperationExecutor.State.WAITING_FOR_INPUT,
                new StudyAssistance((uri, provider, key, body, timeout) -> openai(answer, List.of()), key -> "test-secret").execute(context).state());
        }
        var incomplete = context("incomplete", "advisor_review", Map.of(), Map.of());
        assertEquals(StudyOperationExecutor.State.WAITING_FOR_INPUT,
            new StudyAssistance((uri, provider, key, body, timeout) -> new StudyAssistance.Reply(200,
                "{\"status\":\"incomplete\",\"output\":[]}"), key -> "test-secret").execute(incomplete).state());
    }

    @Test void literatureUsesRetrievalCitationsAndDoesNotAcceptInventedReferences() throws Exception {
        var cited = context("cited", "literature_support", Map.of(), Map.of());
        var citation = Map.of("type", "url_citation", "url", "https://pubmed.ncbi.nlm.nih.gov/123/", "title", "Fixture primary paper");
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> {
            assertTrue(body.contains("web_search"));
            assertTrue(body.contains("\"tool_choice\":\"required\""));
            return openai("Limited relevant evidence [1]", List.of(citation));
        }, key -> "test-secret").execute(cited);
        assertEquals(StudyOperationExecutor.State.COMPLETE, result.state());
        assertTrue(Files.readString(cited.outputDirectory().resolve("decision.yaml")).contains("https://pubmed.ncbi.nlm.nih.gov/123/"));
        var uncited = context("uncited", "literature_support", Map.of(), Map.of());
        var noSources = new StudyAssistance((uri, provider, key, body, timeout) -> openai("Invented reference PMID:999999", List.of()), key -> "test-secret").execute(uncited);
        assertEquals("no_citable_sources_returned", noSources.details().get("status"));
        assertEquals(StudyOperationExecutor.State.WAITING_FOR_INPUT, noSources.state());
    }

    @Test void relatedPhenotypeAdviceIsBoundedCitedAndSearchOnly() throws Exception {
        var context = context("phenotypes", "related_phenotype_advice", Map.of(), Map.of());
        String answer = "{\"phenotypes\":[{\"name\":\"glycated hemoglobin\",\"rationale\":\"BMI is related to glycemic regulation.\",\"literature\":[{\"title\":\"Fixture review\",\"year\":\"2025\",\"doi\":\"10.1000/fixture\",\"pmid\":\"\",\"url\":\"\"}]}]}";
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> {
            assertTrue(body.contains("twenty phenotypes"));
            assertTrue(body.contains("web_search"));
            assertTrue(body.contains("related_phenotypes"));
            return openai(answer, List.of(Map.of("type", "url_citation", "url", "https://doi.org/10.1000/fixture", "title", "Fixture review")));
        }, key -> "test-secret").execute(context);
        assertEquals(StudyOperationExecutor.State.COMPLETE, result.state(), result.details().toString());
        String decision = Files.readString(context.outputDirectory().resolve("decision.yaml"));
        assertTrue(decision.contains("glycated hemoglobin"));
        assertTrue(decision.contains("advisor_proposal_unreviewed"));
        Map<String, Object> saved = StudyPublicEvidence.map(new org.yaml.snakeyaml.Yaml().load(decision));
        assertTrue(saved.get("search_authorization").toString().toLowerCase(java.util.Locale.ROOT).contains("do not amend the hypothesis"));
    }

    @Test void manuscriptAdviceReceivesCatalogedPublicEvidenceAndStructuredIntroductionPolicy() throws Exception {
        var base = context("manuscript-catalog", "manuscript_draft", Map.of("description", "PATIENT_SECRET"), Map.of());
        Path literature = root.resolve("literature-input"); Files.createDirectory(literature);
        var query = Map.of("search_strategy", "field_background", "query", "PRIVATE_QUERY_AUDIT", "status", "candidate_matches", "articles", List.of(
            Map.of("pmid", "123", "title", "Public field paper", "abstractText", "PRIVATE_FULL_ABSTRACT", "excerpt", "Literal source excerpt.")));
        StudyPublicEvidence.write(literature.resolve("manifest.yaml"), Map.of("operation_type", "literature_search", "state", "complete", "output_sha256", Map.of(), "details", Map.of("evidence", List.of(query))));
        var dependencies = new LinkedHashMap<>(base.dependencyDirectories()); dependencies.put("literature", literature);
        var ctx = new StudyOperationExecutor.Context(base.project(), base.outputDirectory(), base.operation(), base.study(), base.plan(), dependencies);
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> {
            assertTrue(body.contains("four or five connected prose paragraphs"));
            assertTrue(body.contains("Always start Methods with a description of every participating cohort"));
            assertTrue(body.contains("results_discussion_policy"));
            assertTrue(body.contains("most impactful finding and why it matters"));
            assertTrue(body.contains("future_potential, next_steps and limitations"));
            assertTrue(body.contains("Public field paper")); assertTrue(body.contains("reference_candidates"));
            assertTrue(body.contains("Cite every journal-derived claim wherever it occurs"));
            assertTrue(body.contains("Review the study_contract_context against ALL observed results"));
            assertTrue(body.contains("pmid:123"));
            assertFalse(body.contains("PRIVATE_QUERY_AUDIT")); assertFalse(body.contains("PRIVATE_FULL_ABSTRACT"));
            assertFalse(body.contains("PATIENT_SECRET")); assertFalse(body.contains("PRIVATE_FEATURE"));
            return openai("Introduction requires scientist review.", List.of());
        }, key -> "fixture-key").execute(ctx);
        assertEquals(StudyOperationExecutor.State.COMPLETE, result.state());
    }

    @Test void geminiAndManuscriptAreAuditedAndDraftHtmlEscapesProviderContent() throws Exception {
        var base = context("gemini", "advisor_review", Map.of(), Map.of());
        Map<String, Object> study = new LinkedHashMap<>(base.study());
        Map<String, Object> config = new LinkedHashMap<>(cast(study.get("assistance")));
        config.put("provider", "gemini");
        config.put("model", "gemini-fixture");
        study.put("assistance", config);
        var context = new StudyOperationExecutor.Context(base.project(), base.outputDirectory(), base.operation(), study, base.plan(), base.dependencyDirectories());
        var response = Map.of("modelVersion", "gemini-fixture-version", "candidates", List.of(Map.of("finishReason", "STOP", "content",
            Map.of("parts", List.of(Map.of("text", "{\"proposals\":[]}"))))));
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> {
            assertTrue(uri.toString().endsWith("gemini-fixture:generateContent"));
            assertFalse(uri.toString().contains("test-secret"));
            assertTrue(body.contains("responseJsonSchema"));
            return new StudyAssistance.Reply(200, StudyAssistance.json(response));
        }, key -> "test-secret").execute(context);
        assertEquals(StudyOperationExecutor.State.COMPLETE, result.state());
        var paper = context("paper", "manuscript_draft", Map.of(), Map.of());
        assertEquals(StudyOperationExecutor.State.COMPLETE,
            new StudyAssistance((uri, provider, key, body, timeout) -> openai("Background\n<script>alert(1)</script>\n[VERIFY]", List.of()), key -> "test-secret").execute(paper).state());
        String html = Files.readString(paper.outputDirectory().resolve("manuscript.html"));
        assertTrue(html.contains("&lt;script&gt;"));
        assertFalse(html.contains("<script>"));
        assertTrue(html.contains("source-map.yaml"));
        assertTrue(html.contains("Not for submission"));
    }

    @Test void networkFailureIsNotMisreportedAsNoScientificEvidence() throws Exception {
        var context = context("network-failure", "literature_support", Map.of(), Map.of());
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> { throw new IOException("DNS fixture"); }, key -> "test-secret").execute(context);
        assertEquals("connection_or_response_error", result.details().get("status"));
        assertEquals(true, result.details().get("request_attempted"));
        assertEquals("unconfirmed", result.details().get("request_sent"));
        assertFalse(result.details().containsKey("search_status"));
        assertTrue(Files.readString(context.outputDirectory().resolve("llm/interaction-001/metadata.yaml")).contains("java.io.IOException"));
    }

    @Test void editorialPolicyDenialPrecedesCredentialsAndSuccessfulEditsRemainProposals() throws Exception {
        var base = context("editorial-denied", "manuscript_editorial", Map.of("data_policy", Map.of("llm", Map.of("remote_transmission", "deny_all"))), Map.of());
        Map<String, Object> study = new LinkedHashMap<>(base.study()), config = new LinkedHashMap<>(cast(study.get("assistance")));
        config.put("editorial", true); study.put("assistance", config);
        var ctx = new StudyOperationExecutor.Context(base.project(), base.outputDirectory(), base.operation(), study, base.plan(), base.dependencyDirectories());
        var client = new StudyAssistance((uri, provider, key, body, timeout) -> { fail("denied must not send"); return null; }, ignored -> { fail("denied must not read credentials"); return null; });
        var result = client.execute(ctx, Map.of("manuscript_sections", List.of(Map.of("text", "Aggregate draft"))));
        assertEquals("policy_denied", result.details().get("status"));
        assertFalse(Files.exists(base.outputDirectory().resolve("llm/interaction-001/request.json")));
    }

    @Test void keylessLocalEditorialUsesOnlyLoopbackAndNeverReadsAnApiKey() throws Exception {
        var base = context("local-editorial", "manuscript_editorial", Map.of("data_policy", Map.of("llm", Map.of("remote_transmission", "deny_all"))), Map.of());
        Map<String, Object> study = new LinkedHashMap<>(base.study());
        study.put("assistance", Map.of("enabled", true, "editorial", true, "provider", "local", "model", "local/model",
            "endpoint", "http://127.0.0.1:1234/v1/chat/completions"));
        var ctx = new StudyOperationExecutor.Context(base.project(), base.outputDirectory(), base.operation(), study, base.plan(), base.dependencyDirectories());
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> {
            assertEquals("127.0.0.1", uri.getHost()); assertNull(key); assertEquals("local", provider);
            assertTrue(body.contains("messages")); assertFalse(body.contains("tools"));
            return new StudyAssistance.Reply(200, StudyAssistance.json(Map.of("model", "local/model", "choices", List.of(Map.of("finish_reason", "stop", "message", Map.of("content", "Section main: clarify heterogeneity."))))));
        }, ignored -> { fail("Local editorial must not read credentials"); return null; }).execute(ctx, Map.of("title", "Aggregate draft"));
        assertEquals("proposal_only_not_applied", result.details().get("editorial_status"));
        assertTrue(Files.readString(base.outputDirectory().resolve("llm/interaction-001/metadata.yaml")).contains("loopback_only_user_managed_model_server"));
        for (String bad : List.of("http://example.org:1234/v1/chat/completions", "http://localhost:1234/v1/chat/completions",
                "http://127.0.0.1:1234/v1/chat/completions?redirect=remote", "http://user@127.0.0.1:1234/v1/chat/completions", "https://127.0.0.1:1234/v1/chat/completions"))
            assertThrows(IllegalArgumentException.class, () -> StudyAssistance.localEndpoint(bad));
    }

    @Test void llamaCppIsKeylessLocalEditorialAndRejectsTruncation() throws Exception {
        for (String finish : List.of("stop", "length")) {
            var base = context("llama-" + finish, "manuscript_editorial",
                Map.of("data_policy", Map.of("remote_transmission", "deny_all")), Map.of());
            Map<String, Object> study = new LinkedHashMap<>(base.study());
            study.put("assistance", Map.of("enabled", true, "editorial", true, "provider", "llama_cpp",
                "model", "editorial-local", "endpoint", "http://127.0.0.1:8080/v1/chat/completions"));
            var ctx = new StudyOperationExecutor.Context(base.project(), base.outputDirectory(), base.operation(), study, base.plan(), base.dependencyDirectories());
            var result = new StudyAssistance((uri, provider, key, body, timeout) -> {
                assertEquals("llama_cpp", provider); assertNull(key);
                assertEquals("http://127.0.0.1:8080/v1/chat/completions", uri.toString());
                assertTrue(body.contains("\"stream\":false")); assertTrue(body.contains("\"max_tokens\":4000"));
                assertFalse(body.contains("tools"));
                return new StudyAssistance.Reply(200, StudyAssistance.json(Map.of("model", "editorial-local",
                    "choices", List.of(Map.of("finish_reason", finish, "message", Map.of("content", "Clarify limitations."))))));
            }, ignored -> { fail("llama.cpp must not read credentials"); return null; }).execute(ctx, Map.of("title", "Aggregate draft"));
            if (finish.equals("stop")) assertEquals("proposal_only_not_applied", result.details().get("editorial_status"));
            else assertEquals(StudyOperationExecutor.State.WAITING_FOR_INPUT, result.state());
            String metadata = Files.readString(base.outputDirectory().resolve("llm/interaction-001/metadata.yaml"));
            assertTrue(metadata.contains("loopback_only_user_managed_model_server"));
            assertTrue(metadata.contains("llama.cpp"));
        }
    }

    @Test void typedQuestionsAndCodeAreAuditedWithoutGrantingExecution() throws Exception {
        var context = context("typed-advice", "advisor_review", Map.of(), Map.of());
        String advice = StudyAssistance.json(Map.of("proposals", List.of(
            Map.of("action", "scientist_question", "rationale", "Interpret tissue relevance", "required_inputs", List.of(), "content", "Should brain sQTL evidence complement blood?", "language", ""),
            Map.of("action", "related_phenotype", "rationale", "Review shared biology", "required_inputs", List.of(), "content", "Sleep duration", "language", ""),
            Map.of("action", "code_proposal", "rationale", "Review sensitivity calculation", "required_inputs", List.of("aggregate table"), "content", "print('review only')", "language", "python"))));
        var result = new StudyAssistance((uri, provider, key, body, timeout) -> openai(advice, List.of()), key -> "test-secret").execute(context);
        assertEquals(StudyOperationExecutor.State.COMPLETE, result.state());
        assertEquals(1, ((List<?>) result.details().get("needs")).size());
        String decision = Files.readString(context.outputDirectory().resolve("decision.yaml"));
        assertTrue(decision.contains("llm_suggested_review"));
        assertTrue(decision.contains("code_sha256"));
        assertTrue(decision.contains("test-model-version"));
        assertTrue(Files.isRegularFile(context.outputDirectory().resolve("llm/interaction-001/proposed-code-3.txt")));
        assertFalse(Files.exists(context.outputDirectory().resolve("proposed-code.py")));
    }

    @Test void compatibleApiUsesExplicitHttpsEndpointAndEnvCredentialsWithoutFallback() throws Exception {
        for (String auth : List.of("api_key", "oauth")) {
            var base = context("compatible-" + auth, "advisor_review", Map.of(), Map.of());
            Map<String,Object> study = new LinkedHashMap<>(base.study());
            Map<String,Object> config = new LinkedHashMap<>(cast(study.get("assistance")));
            config.put("provider", "openai_compatible"); config.put("endpoint", "https://gateway.example.org/v1/chat/completions");
            config.put("auth", auth); config.put(auth.equals("oauth") ? "oauth_token_env" : "api_key_env", "TEST_CREDENTIAL");
            config.put("model", "vendor/model"); study.put("assistance", config);
            var ctx = new StudyOperationExecutor.Context(base.project(), base.outputDirectory(), base.operation(), study, base.plan(), base.dependencyDirectories());
            var result = new StudyAssistance((uri, provider, key, body, timeout) -> {
                assertEquals("https://gateway.example.org/v1/chat/completions", uri.toString());
                assertEquals("openai_compatible", provider); assertEquals("PRIVATE_TOKEN", key);
                assertTrue(body.contains("messages")); assertFalse(body.contains("PRIVATE_TOKEN"));
                assertFalse(body.contains("web_search"));
                return new StudyAssistance.Reply(200, StudyAssistance.json(Map.of("model", "served-model", "choices", List.of(
                    Map.of("finish_reason", "stop", "message", Map.of("content", "{\"proposals\":[]}", "debug", "PRIVATE_TOKEN"))))));
            }, name -> { assertEquals("TEST_CREDENTIAL", name); return "PRIVATE_TOKEN"; }).execute(ctx);
            assertEquals(StudyOperationExecutor.State.COMPLETE, result.state());
            try (var files = Files.walk(ctx.outputDirectory())) {
                for (Path file : files.filter(Files::isRegularFile).toList()) assertFalse(Files.readString(file).contains("PRIVATE_TOKEN"));
            }
            assertTrue(Files.readString(ctx.outputDirectory().resolve("llm/interaction-001/metadata.yaml")).contains("authentication: " + auth));
        }
        for (String bad : List.of("http://gateway.example.org/v1/chat/completions", "https://user:secret@gateway.example.org/v1/chat/completions",
                "https://gateway.example.org/v1/chat/completions?key=secret", "https://gateway.example.org/v1/chat/completions#fragment", "https://gateway.example.org/v1/responses"))
            assertThrows(IllegalArgumentException.class, () -> StudyAssistance.compatibleEndpoint(bad));
    }

    @Test void compatibleEndpointValidationPrecedesCredentialRead() throws Exception {
        var base = context("bad-compatible", "advisor_review", Map.of(), Map.of());
        Map<String,Object> study = new LinkedHashMap<>(base.study()), config = new LinkedHashMap<>(cast(study.get("assistance")));
        config.put("provider", "openai_compatible"); config.put("endpoint", "http://insecure.example/v1/chat/completions"); study.put("assistance", config);
        var ctx = new StudyOperationExecutor.Context(base.project(), base.outputDirectory(), base.operation(), study, base.plan(), base.dependencyDirectories());
        var result = new StudyAssistance((u,p,k,b,t) -> { fail("No network"); return null; }, k -> { fail("No credentials"); return null; }).execute(ctx);
        assertEquals("invalid_provider_response_or_configuration", result.details().get("status"));
    }

    private StudyOperationExecutor.Context context(String name, String type, Map<String, Object> project, Map<String, Object> resources) throws IOException {
        Path output = root.resolve(name);
        Files.createDirectories(output);
        Path input = root.resolve("input");
        Files.createDirectories(input);
        Files.writeString(input.resolve("manifest.yaml"), "operation_type: meta_analysis\nstate: complete\noutput_sha256: {}\n");
        Files.writeString(input.resolve("meta-analysis.tsv"), "feature_id\tbeta\tse\tp_value\nPRIVATE_FEATURE\t1\t.1\t.001\n");
        Files.writeString(input.resolve("feature-annotation.tsv"), "feature_id\tgene_symbol\tmapping_status\nPRIVATE_FEATURE\tGENE1\tlocated\n");
        Map<String, Object> study = new LinkedHashMap<>();
        study.put("project", project);
        study.put("resources", resources);
        study.put("datasets", List.of(Map.of("files", Map.of("matrix", "PATIENT_SECRET"), "metadata", Map.of("dictionary", "PATIENT_SECRET"))));
        study.put("hypotheses", List.of(Map.of("id", "H1", "outcome", "BMI", "description", "PATIENT_SECRET")));
        study.put("assistance", Map.of("enabled", true, "provider", "openai", "model", "test-model",
            "next_steps", true, "related_phenotypes", true, "literature", true, "manuscript", true, "allow_analyses", List.of("leave_one_cohort_out")));
        return new StudyOperationExecutor.Context(root, output, Map.of("type", type, "hypothesis", "H1"), study, Map.of(), Map.of("meta", input));
    }

    private static StudyAssistance.Reply openai(String answer, List<?> annotations) {
        return new StudyAssistance.Reply(200, StudyAssistance.json(Map.of("status", "completed", "model", "test-model-version", "id", "response-fixture",
            "output", List.of(Map.of("type", "message", "content", List.of(Map.of("type", "output_text", "text", answer, "annotations", annotations)))))));
    }
    @SuppressWarnings("unchecked") private static Map<String, Object> cast(Object value) { return (Map<String, Object>) value; }
}
