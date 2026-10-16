/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyProject;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class StudyManuscriptRevisionTest {
    @TempDir Path directory;

    @Test void explicitLocalOnlyRevisionRetainsFullTextWithoutWeakeningEditorialLimits() throws Exception {
        Path p=project(), r=request(); String ledger=FollowupSupport.hash(p.resolve("ledger/events.jsonl"));
        String longText="Complete local manuscript context ".repeat(4500);
        Files.writeString(r,Files.readString(r).replace("<script>alert('not executable')</script> Editorial context.",longText));
        var oversized=assertThrows(java.io.IOException.class,()->StudyManuscriptRevision.create(p,r,directory.resolve("too-large"),StudyProject.Actor.human("test")));
        assertTrue(oversized.getMessage().contains("120000"));assertFalse(Files.exists(directory.resolve("too-large")));
        Files.writeString(r,Files.readString(r)+"\neditorial_payload: {mode: omit, rationale: Keep the complete manuscript local without an external editorial handoff.}\n");
        Path out=directory.resolve("local-only");StudyManuscriptRevision.create(p,r,out,StudyProject.Actor.human("test"));
        StudyManuscriptRevision.verifyBundle(out);
        assertTrue(Files.readString(out.resolve("artifacts/manuscript.html")).contains(longText));
        assertTrue(Files.readString(out.resolve("artifacts/editorial-payload-status.yaml")).contains("omitted_by_explicit_request"));
        assertFalse(Files.exists(out.resolve("artifacts/editorial-payload.json")));
        var missing=assertThrows(java.io.IOException.class,()->StudyEditorialReview.request(out,directory.resolve("not-sent"),"unused",Map.of(),StudyProject.Actor.human("test"),false,
            new StudyAssistance((uri,provider,key,body,timeout)->{fail("No omitted payload may be transmitted");return null;},ignored->{fail("No credential access");return null;})));
        assertTrue(missing.getMessage().contains("no editorial disclosure payload"));assertFalse(Files.exists(directory.resolve("not-sent")));
        assertEquals(ledger,FollowupSupport.hash(p.resolve("ledger/events.jsonl")));
    }

    @Test void editorialPayloadOptOutRequiresAnExplicitReason() throws Exception {
        Path p=project(),r=request();Files.writeString(r,Files.readString(r)+"\neditorial_payload: {mode: omit}\n");
        assertThrows(java.io.IOException.class,()->StudyManuscriptRevision.create(p,r,directory.resolve("missing-reason"),StudyProject.Actor.human("test")));
        assertFalse(Files.exists(directory.resolve("missing-reason")));
    }

    @Test void importedReferencesAreHashPinnedOfflineAndStrictCitationsPersistInPayloadAndPaper() throws Exception {
        Path p = project(), refs = directory.resolve("refs.yaml");
        Files.writeString(refs, "schema_version: 1\nreferences:\n  - {id: 'pmid:123', rationale: 'Cohort design reference.'}\n");
        // This source denies downloads: no transport call is permitted.
        var denied = StudyCitationSources.fetch(p, refs, directory.resolve("denied"), StudyProject.Actor.human("test"), false,
            new StudyPublicEvidence(uri -> { throw new AssertionError("Denied policy must prevent network"); }));
        assertTrue(Files.readString(denied.resolve("metadata/manifest.yaml")).contains("needs_input"));
        Path policy = directory.resolve("public-policy"); Files.createDirectory(policy);
        Files.writeString(policy.resolve("resolved-study.yaml"), "resources: {downloads: automatic}\n");
        var acquired = StudyCitationSources.fetch(policy, refs, directory.resolve("refs"), StudyProject.Actor.human("test"), false,
            new StudyPublicEvidence(uri -> new StudyPublicEvidence.Reply(200, "{\"resultList\":{\"result\":[{\"source\":\"MED\",\"id\":\"123\",\"title\":\"Cohort study\",\"pubYear\":\"2020\"}]}}")));
        Path manifest = acquired.resolve("metadata/manifest.yaml"), r = request();
        Files.writeString(r, Files.readString(r).replace("<script>alert('not executable')</script> Editorial context.", "Context[[cite:pmid:123]].")
            .replace("section: introduction", "section: introduction\n    citation_basis: literature").replace("section: results", "section: results\n    citation_basis: study_results")
            + "\ncitation_policy: strict\nreference_bundles:\n  - manifest: " + manifest.toString().replace('\\', '/') + "\n    sha256: " + FollowupSupport.hash(manifest) + "\n");
        Path out = directory.resolve("cited"); StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("test")); StudyManuscriptRevision.verifyBundle(out);
        String html = Files.readString(out.resolve("artifacts/manuscript.html"));
        assertTrue(html.contains("href=\"#reference-1\"")); assertTrue(html.lastIndexOf("<h2>References</h2>") > html.indexOf("Reproducibility and use"));
        assertTrue(Files.readString(out.resolve("artifacts/editorial-payload.json")).contains("Cohort study"));
        assertTrue(Files.exists(out.resolve("metadata/reference-bundles/1/queries/0001/response.txt")));
        Files.writeString(acquired.resolve("artifacts/queries.yaml"), "tampered");
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, directory.resolve("tampered"), StudyProject.Actor.human("test")));
    }

    private Path project() throws Exception {
        return project("project", false);
    }
    private Path project(String name, boolean characteristics) throws Exception {
        Path p = directory.resolve(name); Files.createDirectory(p);
        Files.writeString(p.resolve("a.tsv"), "feature_id\tbeta\tse\tp_value\nGENE1\t1\t.1\t.001\n");
        Files.writeString(p.resolve("b.tsv"), "feature_id\tbeta\tse\tp_value\nGENE1\t2\t.2\t.001\n");
        Files.writeString(p.resolve("study.yaml"), """
            schema_version: 1
            project: {id: manuscript-fixture, title: Local revision, description: Test reporting without network.}
            resources: {downloads: deny}
            datasets:
              - {id: a, cohort: a, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: a.tsv}}
              - {id: b, cohort: b, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: b.tsv}}
            hypotheses:
              - id: H1
                description: Isoform ratios associate with the phenotype.
                analysis: {mode: manual, goals: [meta_analysis]}
            """);
        if (characteristics) {
            Files.writeString(p.resolve("phenotype.csv"), "id,Age,Sex,BMI\na,20,1,25\nb,40,2,31\n");
            String study = Files.readString(p.resolve("study.yaml")).replace("files: {results: a.tsv}",
                "files: {results: a.tsv, phenotype: phenotype.csv}, columns: {characteristics: {id: id, main_phenotype: BMI, male_value: 1, female_value: 2}}");
            Files.writeString(p.resolve("study.yaml"), study);
        }
        ByteArrayOutputStream err = new ByteArrayOutputStream();
        assertEquals(0, StudyCli.run(new String[]{"auto", p.resolve("study.yaml").toString(), "--offline"},
            new PrintStream(new ByteArrayOutputStream()), new PrintStream(err)), err.toString());
        return p;
    }
    private Path request() throws Exception {
        Path r = directory.resolve("revision.yaml");
        Files.writeString(r, """
            schema_version: 1
            id: revision-1
            title: Molecular associations
            rationale: Integrate the verified current findings.
            sections:
              - id: intro
                section: introduction
                text: "<script>alert('not executable')</script> Editorial context."
                rationale: Preserve text literally.
                sources: []
              - id: main
                section: results
                title: Pooled associations
                text: "Tested {{H1.meta_analysis#/tested_features}} features."
                sources: [H1.meta_analysis]
                rationale: Summarize the primary completed analysis.
            comments:
              - {id: review-1, comment: Verify the count., response: Bound to native evidence., section: main, status: addressed}
            """);
        return r;
    }
    @Test void participantRefreshIsNativeHashPinnedAndDoesNotReplaceAssociationClaims() throws Exception {
        Path p = project("original", true), refresh = project("refresh", true), r = request();
        Files.writeString(r, Files.readString(r) + "\nparticipant_project:\n  path: " + refresh.toString().replace('\\', '/')
            + "\n  ledger_sha256: " + FollowupSupport.hash(refresh.resolve("ledger/events.jsonl")) + "\n  rationale: Updated local demographics.\n");
        Path out = directory.resolve("with-refresh");
        StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("scientist"));
        StudyManuscriptRevision.verifyBundle(out);
        assertTrue(Files.isRegularFile(out.resolve("artifacts/reference-catalog.html")));
        assertTrue(Files.readString(out.resolve("artifacts/introduction-review.yaml")).contains("older prose is preserved"));
        assertTrue(Files.readString(out.resolve("artifacts/methods-inventory.yaml")).contains("effective_defaults"));
        assertTrue(Files.readString(out.resolve("artifacts/methods-template.yaml")).contains("methods_role: cohorts"));
        String html = Files.readString(out.resolve("artifacts/manuscript.html"));
        String methods = html.substring(html.indexOf("<h2>Methods</h2>"), html.indexOf("<h2>Results</h2>"));
        assertTrue(methods.indexOf("describe cohort design") < methods.indexOf("Execution and parameter inventory"));
        assertTrue(Files.exists(out.resolve("metadata/participant-study.yaml")));
        assertTrue(Files.readString(out.resolve("artifacts/manuscript.html")).contains("Tested 1 features."));
        assertTrue(Files.readString(out.resolve("artifacts/participant-characteristics.csv")).contains("14.142135623730951"));
        assertTrue(Files.readString(out.resolve("artifacts/supplement.html")).contains("Supplemental Table 1"));
        Files.writeString(r, Files.readString(r).replace(FollowupSupport.hash(refresh.resolve("ledger/events.jsonl")), "0".repeat(64)));
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, directory.resolve("bad-refresh"), StudyProject.Actor.human("scientist")));
        assertFalse(Files.exists(directory.resolve("bad-refresh")));
    }
    @Test void nativeRevisionUsesVerifiedEvidenceWithoutMutatingSourceOrExecutingHtml() throws Exception {
        Path p = project(), r = request(), out = directory.resolve("revision-1");
        String ledger = Files.readString(p.resolve("ledger/events.jsonl"));
        String oldIndex = Files.readString(p.resolve("index/index.json"));
        StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("scientist"));
        String html = Files.readString(out.resolve("artifacts/manuscript.html"));
        assertTrue(html.contains("Tested 1 features."));
        assertTrue(html.contains("&lt;script&gt;alert"));
        assertFalse(html.contains("<script>alert"));
        assertTrue(html.indexOf("Table 1.") < html.indexOf("Pooled associations"));
        assertEquals(ledger, Files.readString(p.resolve("ledger/events.jsonl")));
        assertEquals(oldIndex, Files.readString(p.resolve("index/index.json")));
        StudyManuscriptRevision.verifyBundle(out);
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("scientist")));
    }
    @Test void orderedNarrativeLinksContextAndImpactReasonsInAuditedBundleWithoutFiveResultCap() throws Exception {
        Path p = project(), r = directory.resolve("narrative.yaml"), out = directory.resolve("narrative-bundle");
        var sections = new ArrayList<Object>();
        for (int i = 1; i <= 6; i++) {
            var result = new java.util.LinkedHashMap<String, Object>(Map.of("id", "result-" + i, "section", "results", "results_role", "meta_analysis",
                "title", "Synthetic result " + i, "text", "Tested {{H1.meta_analysis#/tested_features}} features.", "sources", List.of("H1.meta_analysis"), "rationale", "Test the writing contract."));
            if (i == 1) result.put("supporting_evidence", List.of(Map.of("source", "H1.meta_analysis", "relationship",
                "Context from the same {{H1.meta_analysis#/tested_features}}-feature analysis is not independent evidence. <script>literal only</script>")));
            sections.add(result);
        }
        sections.add(Map.of("id", "impact", "section", "discussion", "discussion_role", "principal_finding", "impact_rank", 1,
            "impact_reason", "The pooled estimate is the only completed inference in this fixture.", "result_ids", List.of("result-1"), "text", "This illustrative finding motivates follow-up, subject to uncertainty.",
            "sources", List.of("H1.meta_analysis"), "rationale", "Do not invent three findings."));
        for (String role : StudyResultsDiscussion.DISCUSSION.subList(1, 5)) sections.add(Map.of("id", role.replace('_', '-'), "section", role.equals("overall_conclusion") ? "conclusion" : "discussion",
            "discussion_role", role, "text", "Narrative paragraph for " + role + ".", "sources", List.of(), "rationale", "Fixture prose, not scientific claims."));
        Files.writeString(r, new org.yaml.snakeyaml.Yaml().dump(Map.of("schema_version", 1, "id", "narrative", "title", "Synthetic narrative", "rationale", "Verify the contract.",
            "discussion_selection_note", "Only one substantive finding is available in this fixture.", "sections", sections)));
        StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("scientist"));
        StudyManuscriptRevision.verifyBundle(out);
        String html = Files.readString(out.resolve("artifacts/manuscript.html"));
        assertTrue(html.indexOf("Table 1.") < html.indexOf("Synthetic result 1"));
        assertTrue(html.indexOf("1-feature analysis") < html.indexOf("Synthetic result 2"));
        assertTrue(html.contains("Synthetic result 6")); assertTrue(html.contains("&lt;script&gt;literal only&lt;/script&gt;"));
        assertFalse(html.contains("<script>literal only")); assertFalse(html.contains("Fewer than three main points"));
        assertTrue(html.contains("href=\"#result-1\""));
        assertTrue(html.indexOf("Narrative paragraph for future_potential") < html.indexOf("Narrative paragraph for next_steps"));
        assertTrue(html.indexOf("Narrative paragraph for limitations") < html.indexOf("Narrative paragraph for overall_conclusion"));
        String claims = Files.readString(out.resolve("artifacts/publication-claims.yaml"));
        assertTrue(claims.contains("supporting_source: H1.meta_analysis")); assertTrue(claims.contains("impact_rank: 1"));
        String payload = Files.readString(out.resolve("artifacts/editorial-payload.json"));
        assertTrue(payload.contains("supporting_evidence")); assertTrue(payload.contains("impact_reason"));
        assertFalse(payload.contains("{{H1.meta_analysis")); assertFalse(payload.contains(directory.toString()));
        assertTrue(Files.readString(out.resolve("artifacts/results-discussion-review.yaml")).contains("not independent validation"));
    }
    @Test void explicitParagraphReorderingPreservesProseAttributionAndRejectsIncompleteOrders() throws Exception {
        Path p = project(), original = directory.resolve("original-order");
        StudyManuscriptRevision.create(p, request(), original, StudyProject.Actor.human("original-author"));
        Path r = directory.resolve("reorder.yaml");
        var edit = new java.util.LinkedHashMap<String, Object>(Map.of("schema_version", 1, "id", "reordered", "title", "Reordered narrative", "rationale", "Explicit ordering decision.",
            "previous_revision", Map.of("manifest", original.resolve("metadata/manifest.yaml").toString(), "sha256", FollowupSupport.hash(original.resolve("metadata/manifest.yaml"))),
            "sections", List.of(Map.of("id", "before-main", "section", "results", "text", "An earlier result paragraph.", "sources", List.of("H1.meta_analysis"), "rationale", "Place before main.")),
            "paragraph_order", List.of("intro", "before-main", "main")));
        Files.writeString(r, new org.yaml.snakeyaml.Yaml().dump(edit));
        Path out = directory.resolve("new-order");
        StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("editor"));
        String html = Files.readString(out.resolve("artifacts/manuscript.html"));
        assertTrue(html.indexOf("An earlier result paragraph") < html.indexOf("Tested 1 features"));
        var attribution = StudyPublicEvidence.map(new org.yaml.snakeyaml.Yaml().load(Files.readString(out.resolve("metadata/section-attribution.yaml"))));
        assertTrue(attribution.get("main").toString().contains("original-author"));
        assertTrue(attribution.get("before-main").toString().contains("editor"));
        assertTrue(Files.readString(out.resolve("metadata/change-log.yaml")).contains("paragraphs_reordered"));
        StudyManuscriptRevision.verifyBundle(out);
        edit.put("paragraph_order", List.of("intro", "main")); Files.writeString(r, new org.yaml.snakeyaml.Yaml().dump(edit));
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, directory.resolve("bad-order"), StudyProject.Actor.human("editor")));
        assertFalse(Files.exists(directory.resolve("bad-order")));
    }
    @Test void editorialCliDefaultsToLocalBrowserPacketWithoutCredentialsOrModel() throws Exception {
        Path p = project(), revision = directory.resolve("browser-revision"), out = directory.resolve("browser-editorial");
        StudyManuscriptRevision.create(p, request(), revision, StudyProject.Actor.human("scientist"));
        String hash = FollowupSupport.hash(revision.resolve("artifacts/editorial-payload.json"));
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        assertEquals(0, StudyCli.run(new String[]{"manuscript-editorial", "--project", revision.toString(), "--out", out.toString(),
            "--approved-payload-sha256", hash, "--actor-id", "scientist", "--offline"}, new PrintStream(new ByteArrayOutputStream()), new PrintStream(error)), error.toString());
        StudyManuscriptRevision.verifyBundle(out);
        assertTrue(Files.isRegularFile(out.resolve("artifacts/handoff/prompt.txt")));
        assertTrue(Files.readString(out.resolve("artifacts/index.html")).contains("handoff/index.html"));
        assertTrue(Files.readString(out.resolve("artifacts/decision.yaml")).contains("waiting_for_browser_response"));
        assertFalse(Files.readString(out.resolve("artifacts/decision.yaml")).contains("assistance.response_files"));
        assertEquals(2, StudyCli.run(new String[]{"manuscript-editorial", "--model", "ambiguous-remote-model"},
            new PrintStream(new ByteArrayOutputStream()), new PrintStream(error)));
    }
    @Test void llamaCppEditorialCliIsExplicitLocalAndOfflineDoesNotContactServer() throws Exception {
        Path p = project(), revision = directory.resolve("llama-revision"), out = directory.resolve("llama-review");
        StudyManuscriptRevision.create(p, request(), revision, StudyProject.Actor.human("scientist"));
        String hash = FollowupSupport.hash(revision.resolve("artifacts/editorial-payload.json"));
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        assertEquals(0, StudyCli.run(new String[]{"manuscript-editorial", "--project", revision.toString(), "--out", out.toString(),
            "--provider", "llama_cpp", "--model", "editorial-local", "--endpoint", "http://127.0.0.1:8080/v1/chat/completions",
            "--approved-payload-sha256", hash, "--actor-id", "scientist", "--offline"},
            new PrintStream(new ByteArrayOutputStream()), new PrintStream(error)), error.toString());
        StudyManuscriptRevision.verifyBundle(out);
        assertTrue(Files.readString(out.resolve("artifacts/decision.yaml")).contains("offline_policy"));
        assertFalse(Files.exists(out.resolve("artifacts/llm/interaction-001/request.json")));
        for (String[] invalid : List.of(
                new String[]{"--api-key-env", "UNUSED"}, new String[]{"--auth", "oauth"},
                new String[]{"--endpoint", "http://example.org:8080/v1/chat/completions"})) {
            assertEquals(2, StudyCli.run(new String[]{"manuscript-editorial", "--provider", "llama_cpp", "--model", "editorial-local",
                invalid[0], invalid[1]}, new PrintStream(new ByteArrayOutputStream()), new PrintStream(error)));
        }
    }

    @Test void subsequentRevisionPreservesUntouchedProseAndAttributionAndRecordsReplacement() throws Exception {
        Path p = project(), first = directory.resolve("first");
        StudyManuscriptRevision.create(p, request(), first, StudyProject.Actor.human("author-one"));
        Path next = directory.resolve("next.yaml");
        String manifest = first.resolve("metadata/manifest.yaml").toString().replace('\\', '/');
        Files.writeString(next, """
            schema_version: 1
            id: revision-2
            title: Revised molecular associations
            rationale: Clarify interpretation without changing analysis.
            previous_revision:
              manifest: '%s'
              sha256: %s
            sections:
              - id: main
                section: results
                title: Updated pooled association
                text: "Current tested count {{H1.meta_analysis#/tested_features}}; association is not causality."
                sources: [H1.meta_analysis]
                rationale: Preserve the distinction between association and causality.
            """.formatted(manifest, FollowupSupport.hash(first.resolve("metadata/manifest.yaml"))));
        Path second = directory.resolve("second");
        StudyManuscriptRevision.create(p, next, second, StudyProject.Actor.human("author-two"));
        String html = Files.readString(second.resolve("artifacts/manuscript.html"));
        assertTrue(html.contains("Editorial context.")); assertTrue(html.contains("Current tested count 1"));
        String authors = Files.readString(second.resolve("metadata/section-attribution.yaml"));
        assertTrue(authors.contains("author-one")); assertTrue(authors.contains("author-two"));
        assertTrue(Files.readString(second.resolve("metadata/change-log.yaml")).contains("replaced"));
        assertTrue(Files.readString(second.resolve("metadata/resolved-revision.yaml")).contains("review-1"));
        Files.writeString(first.resolve("artifacts/manuscript.html"), "tampered");
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, next, directory.resolve("third"), StudyProject.Actor.human("author-two")));
    }
    @Test void invalidBindingFailsBeforePublishingAndSourceTamperingIsRejected() throws Exception {
        Path p = project(), r = request();
        Files.writeString(r, Files.readString(r).replace("#/tested_features", "#/missing"));
        Path out = directory.resolve("bad");
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("a")));
        assertFalse(Files.exists(out));
        try (var files = Files.walk(p.resolve("runs"))) {
            Path manifest = files.filter(x -> x.getFileName().toString().equals("manifest.yaml")).findFirst().orElseThrow();
            Files.writeString(manifest, Files.readString(manifest) + "\ntampered: true\n");
        }
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("a")));
    }
    @Test void pointersRejectMissingComplexOrMalformedValuesAndTrackExactScalar() throws Exception {
        Map<String, Map<String, Object>> evidence = Map.of("H.meta", Map.of("manifest", "manifest.yaml", "key_findings", Map.of("n", 3, "rows", List.of(Map.of("beta", "2.5")))));
        List<Object> claims = new ArrayList<>();
        assertEquals("n=3 beta=2.5", StudyManuscriptRevision.expand("n={{H.meta#/n}} beta={{H.meta#/rows/0/beta}}", evidence, claims, "main"));
        assertEquals(2, claims.size());
        for (String invalid : List.of("{{H.meta#/absent}}", "{{H.meta#/rows}}", "{{H.meta#/rows/99/beta}}", "{{oops}}"))
            assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.expand(invalid, evidence, claims, "main"));
    }
    @Test void yamlTypesAndUnknownFieldsCannotSilentlyDiscardEditorialInstructions() throws Exception {
        Path p = project(), r = request();
        Files.writeString(r, Files.readString(r) + "\nexecute_code: false\n");
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, directory.resolve("bad"), StudyProject.Actor.human("a")));
        Files.writeString(r, "schema_version: 1\nid: r\ntitle: t\nrationale: x\nsections: not-a-list\n");
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, directory.resolve("bad"), StudyProject.Actor.human("a")));
    }

    @Test void editorialPassRequiresExactDisclosureHashAndPreservesManuscript() throws Exception {
        Path p = project(), revision = directory.resolve("review-source");
        StudyManuscriptRevision.create(p, request(), revision, StudyProject.Actor.human("scientist"));
        String before = FollowupSupport.hash(revision.resolve("artifacts/manuscript.html"));
        Map<String, Object> config = Map.of("provider", "openai", "model", "fixture-model");
        var client = new StudyAssistance((uri, provider, key, body, timeout) -> {
            assertTrue(body.contains("existing manuscript")); assertTrue(body.contains("Tested 1 features."));
            assertTrue(body.contains("four or five connected prose paragraphs"));
            assertTrue(body.contains("importance_and_gap"));
            assertTrue(body.contains("Always start Methods with a description of every participating cohort"));
            assertTrue(body.contains("methods_policy"));
            assertTrue(body.contains("results_discussion_policy"));
            assertTrue(body.contains("impact_reason"));
            assertTrue(body.contains("same order as Methods"));
            assertTrue(body.contains("literature_search_coverage"));
            assertFalse(body.contains("web_search")); assertFalse(body.contains(directory.toString()));
            return new StudyAssistance.Reply(200, StudyAssistance.json(Map.of("status", "completed", "model", "returned-fixture-model",
                "output", List.of(Map.of("type", "message", "content", List.of(Map.of("type", "output_text", "text", "Section main: clarify uncertainty. <script>not executed</script>")))))));
        }, ignored -> "fixture-key");
        assertThrows(java.io.IOException.class, () -> StudyEditorialReview.request(revision, directory.resolve("wrong"), "wrong", config, StudyProject.Actor.human("scientist"), false, client));
        assertFalse(Files.exists(directory.resolve("wrong")));
        String payloadHash = FollowupSupport.hash(revision.resolve("artifacts/editorial-payload.json"));
        Path review = StudyEditorialReview.request(revision, directory.resolve("review"), payloadHash, config, StudyProject.Actor.human("scientist"), false, client);
        StudyManuscriptRevision.verifyBundle(review);
        assertEquals(before, FollowupSupport.hash(revision.resolve("artifacts/manuscript.html")));
        assertTrue(Files.readString(review.resolve("artifacts/editorial-review.html")).contains("&lt;script&gt;"));
        assertTrue(Files.readString(review.resolve("artifacts/decision.yaml")).contains("proposal_only_not_applied"));
        assertTrue(Files.readString(review.resolve("artifacts/llm/interaction-001/metadata.yaml")).contains("returned-fixture-model"));
        Path offline = StudyEditorialReview.request(revision, directory.resolve("offline"), payloadHash, config, StudyProject.Actor.human("scientist"), true,
            new StudyAssistance((uri, provider, key, body, timeout) -> { fail("offline must not send"); return null; }, ignored -> { fail("offline must not read keys"); return null; }));
        assertTrue(Files.readString(offline.resolve("artifacts/decision.yaml")).contains("offline_policy"));
    }

    @Test void typedIntroductionRendersInOrderWithoutMiniHeadingsAndHypothesisIsLast() throws Exception {
        Path p = project(), r = directory.resolve("intro.yaml"), out = directory.resolve("intro-revision");
        var request = new java.util.LinkedHashMap<>(StudyIntroduction.template("Fixture introduction", List.of("Declared hypothesis.")));
        List<Object> sections = new ArrayList<>();
        for (Object raw : StudyPublicEvidence.list(request.get("sections"))) {
            var section = StudyPublicEvidence.map(raw);
            section.put("text", "Fixture " + section.get("introduction_role") + ".");
            sections.add(section);
        }
        request.put("sections", sections); StudyPublicEvidence.write(r, request);
        StudyManuscriptRevision.create(p, r, out, StudyProject.Actor.human("scientist"));
        StudyManuscriptRevision.verifyBundle(out);
        String html = Files.readString(out.resolve("artifacts/manuscript.html"));
        String intro = html.substring(html.indexOf("<h2>Introduction</h2>"), html.indexOf("<h2>Methods</h2>"));
        assertFalse(intro.contains("<h3>")); assertFalse(intro.contains("Editorial rationale"));
        int previous = -1;
        for (String role : StudyIntroduction.ROLES) {
            int position = intro.indexOf("Fixture " + role); assertTrue(position > previous); previous = position;
        }
        assertTrue(Files.readString(out.resolve("artifacts/introduction-review.yaml")).contains("cite cataloged literature"));
        java.util.Collections.swap(sections, 0, 1); StudyPublicEvidence.write(r, request);
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, directory.resolve("bad-order"), StudyProject.Actor.human("scientist")));
        assertFalse(Files.exists(directory.resolve("bad-order")));
        java.util.Collections.swap(sections, 0, 1);
        var first = StudyPublicEvidence.map(sections.get(0)); first.put("references", List.of("10.1000/invented")); sections.set(0, first);
        StudyPublicEvidence.write(r, request);
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, r, directory.resolve("invented-reference"), StudyProject.Actor.human("scientist")));
        assertFalse(Files.exists(directory.resolve("invented-reference")));
    }
    @Test void methodsTemplateIsAcceptedAndReversedRolesFailBeforeBundleCreation() throws Exception {
        Path p = project(), first = directory.resolve("methods-source"), out = directory.resolve("methods-updated");
        StudyManuscriptRevision.create(p, request(), first, StudyProject.Actor.human("scientist"));
        Path template = first.resolve("artifacts/methods-template.yaml");
        StudyManuscriptRevision.create(p, template, out, StudyProject.Actor.human("scientist"));
        StudyManuscriptRevision.verifyBundle(out);
        String html = Files.readString(out.resolve("artifacts/manuscript.html"));
        String methods = html.substring(html.indexOf("<h2>Methods</h2>"), html.indexOf("<h2>Results</h2>"));
        assertTrue(methods.indexOf("describe cohort design") < methods.indexOf("describe association"));
        assertTrue(methods.contains("Parameters and recorded justifications"));
        var parsed = StudyPublicEvidence.map(new org.yaml.snakeyaml.Yaml().load(Files.readString(template)));
        List<Object> sections = new ArrayList<>(StudyPublicEvidence.list(parsed.get("sections")));
        java.util.Collections.reverse(sections); parsed.put("sections", sections);
        Path reversed = directory.resolve("reversed-methods.yaml"); StudyPublicEvidence.write(reversed, parsed);
        assertThrows(java.io.IOException.class, () -> StudyManuscriptRevision.create(p, reversed, directory.resolve("bad-methods"), StudyProject.Actor.human("scientist")));
        assertFalse(Files.exists(directory.resolve("bad-methods")));
    }
}
