/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor;
import org.jlinalg.study.StudyProject;
import org.jlinalg.study.StudyValidator;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.yaml.snakeyaml.Yaml;

class StudyBrowserHandoffTest {
    @TempDir Path root;
    private final StudyAssistance offlineClient = new StudyAssistance((u,p,k,b,t) -> {
        fail("Browser handoff must not call a network endpoint"); return null;
    }, key -> { fail("Browser handoff must not read credentials"); return null; });

    private StudyOperationExecutor.Context context(String name, Map<String, Object> policy, Map<String, Object> extra) throws Exception {
        Path out = root.resolve(name); Files.createDirectory(out);
        Map<String, Object> config = new LinkedHashMap<>(Map.of("enabled", true, "next_steps", true)); config.putAll(extra);
        return new StudyOperationExecutor.Context(root, out, Map.of("type", "advisor_review", "hypothesis", "H1"),
            Map.of("project", Map.of("data_policy", policy), "resources", Map.of("network", "offline"), "assistance", config,
                "hypotheses", List.of(Map.of("id", "H1", "outcome", "BMI")),
                "datasets", List.of(Map.of("files", Map.of("phenotypes", "NEVER_READ_PATIENTS.csv")))), Map.of(), Map.of());
    }
    private Map<String, Object> reply(Path packet, String answer) throws Exception {
        Map<String, Object> r = load(packet.resolve("response-template.yaml"));
        r.put("provider", "Google AI Mode"); r.put("model", "not_disclosed"); r.put("url", "https://google.com/ai");
        r.put("submitted_at", Instant.now().toString()); r.put("received_at", Instant.now().toString());
        r.put("exact_prompt_submitted", true); r.put("disclosure_approved", true); r.put("no_patient_data", true);
        r.put("response", answer); r.put("rationale", "Review scientific suggestions against native results.");
        return r;
    }
    private Path save(String name, Map<String, Object> r) throws Exception {
        Path p = root.resolve(name); Files.writeString(p, new Yaml().dump(r)); return p;
    }
    @SuppressWarnings("unchecked") private Map<String, Object> load(Path p) throws Exception { return new Yaml().load(Files.readString(p)); }

    @Test void defaultProviderProducesOfflineScriptFreePacketAndNeedsExplicitHumanResponse() throws Exception {
        var ctx = context("start", Map.of("llm", Map.of("require_human_approval", true)), Map.of());
        var result = offlineClient.execute(ctx);
        assertEquals("waiting_for_browser_response", result.details().get("status"));
        assertEquals(false, result.details().get("request_sent"));
        String prompt = Files.readString(ctx.outputDirectory().resolve("handoff/prompt.txt"));
        assertFalse(prompt.contains("NEVER_READ_PATIENTS")); assertTrue(prompt.contains("BMI"));
        String html = Files.readString(ctx.outputDirectory().resolve("handoff/index.html"));
        assertFalse(html.contains("<script")); assertFalse(html.contains("<form")); assertFalse(html.contains("src=\"https:"));
        assertTrue(html.contains("default-src 'none'")); assertTrue(html.contains("not observe"));
        assertThrows(IOException.class, () -> StudyBrowserHandoff.importReply(ctx.outputDirectory().resolve("handoff"),
            ctx.outputDirectory().resolve("handoff/response-template.yaml"), root.resolve("invalid"), StudyProject.Actor.human("scientist")));
    }

    @Test void importedResponseIsImmutableAuditedAndResumesWithoutNetworkOrImplicitExecution() throws Exception {
        var first = context("first", Map.of(), Map.of("response_files", List.of("reply/receipt.yaml")));
        offlineClient.execute(first); Path packet = first.outputDirectory().resolve("handoff");
        String original = FollowupSupport.hash(packet.resolve("packet.yaml"));
        Path response = save("response.yaml", reply(packet, "{\"proposals\":[{\"action\":\"run_shell\",\"rationale\":\"untrusted\",\"required_inputs\":[]}]}"));
        Path imported = StudyBrowserHandoff.importReply(packet, response, root.resolve("reply"), StudyProject.Actor.human("scientist"));
        StudyBrowserHandoff.verifyReceipt(imported.resolve("receipt.yaml"));
        assertEquals(original, FollowupSupport.hash(packet.resolve("packet.yaml")));
        assertThrows(IOException.class, () -> StudyBrowserHandoff.importReply(packet, response, imported, StudyProject.Actor.human("scientist")));
        var next = context("next", Map.of(), Map.of("response_files", List.of("reply/receipt.yaml")));
        var result = offlineClient.execute(next);
        assertEquals(StudyOperationExecutor.State.COMPLETE, result.state());
        assertEquals(false, result.details().get("request_attempted"));
        assertEquals("human_attested_not_observed", result.details().get("browser_submission"));
        assertTrue(Files.readString(next.outputDirectory().resolve("decision.yaml")).contains("unsupported_requires_scientist_review"));
        assertTrue(Files.isRegularFile(next.outputDirectory().resolve("browser-response/receipt.yaml")));
        Files.writeString(imported.resolve("response.yaml"), "tampered");
        var corrupted = context("corrupted", Map.of(), Map.of("response_files", List.of("reply/receipt.yaml")));
        assertEquals("browser_response_invalid", offlineClient.execute(corrupted).details().get("status"));
    }

    @Test void policyDenialAndUnsupportedDisclosureChecksCannotBeBypassedWithHandoff() throws Exception {
        int n = 0;
        for (Map<String, Object> policy : List.of(Map.<String,Object>of("remote_transmission", "deny_all"),
                Map.<String,Object>of("llm", Map.of("minimum_cell_count", 10)), Map.<String,Object>of("llm", Map.of("eligible_remote_artifacts", List.of("public"))))) {
            var ctx = context("denied" + n++, policy, Map.of());
            assertEquals("policy_denied", offlineClient.execute(ctx).details().get("status"));
            Path packet = ctx.outputDirectory().resolve("handoff");
            Path response = save("denied" + n + ".yaml", reply(packet, "{}"));
            assertThrows(IOException.class, () -> StudyBrowserHandoff.importReply(packet, response, root.resolve("not-created"), StudyProject.Actor.human("scientist")));
            assertFalse(Files.exists(root.resolve("not-created")));
        }
    }

    @Test void importRejectsWrongHashMissingAttestationsUnsafeUrlsAndInvalidTimes() throws Exception {
        var ctx = context("validation", Map.of(), Map.of()); offlineClient.execute(ctx);
        Path packet = ctx.outputDirectory().resolve("handoff");
        int i = 0;
        for (Map.Entry<String, Object> mutation : Map.<String,Object>of("prompt_sha256", "bad", "request_id", "bad",
                "exact_prompt_submitted", false, "disclosure_approved", false, "no_patient_data", false,
                "url", "https://google.com/ai?token=SECRET", "submitted_at", "2100-01-01T00:00:00Z", "model", "").entrySet()) {
            Map<String,Object> r = reply(packet, "{}"); r.put(mutation.getKey(), mutation.getValue());
            Path response = save("bad-" + i++ + ".yaml", r);
            assertThrows(IOException.class, () -> StudyBrowserHandoff.importReply(packet, response, root.resolve("failed"), StudyProject.Actor.human("scientist")));
        }
        Path response = save("valid.yaml", reply(packet, "{}"));
        assertThrows(IOException.class, () -> StudyBrowserHandoff.importReply(packet, response, root.resolve("failed"), new StudyProject.Actor("llm", "assistant", "model")));
        assertThrows(IOException.class, () -> StudyBrowserHandoff.importReply(packet, response, root.resolve("failed"), StudyProject.Actor.human("unspecified")));
        Files.writeString(packet.resolve("prompt.txt"), "tampered");
        assertThrows(IOException.class, () -> StudyBrowserHandoff.importReply(packet, response, root.resolve("failed"), StudyProject.Actor.human("scientist")));
    }

    @Test void browserHtmlEscapesResponseAndDoesNotInventModelIdentity() throws Exception {
        var ctx = context("escape", Map.of(), Map.of()); offlineClient.execute(ctx);
        Path packet = ctx.outputDirectory().resolve("handoff");
        Path response = save("escape-response.yaml", reply(packet, "<script>alert(1)</script>"));
        Path receipt = StudyBrowserHandoff.importReply(packet, response, root.resolve("escaped"), StudyProject.Actor.human("scientist"));
        String html = Files.readString(receipt.resolve("index.html"));
        assertTrue(html.contains("&lt;script&gt;")); assertFalse(html.contains("<script>"));
        assertTrue(html.contains("not_disclosed"));
        assertEquals(0, cli("handoff-verify", receipt.resolve("receipt.yaml").toString()));
    }

    @Test void openingIsExplicitAndOnlyAllowsFixedBrowserNamesAndVerifiedLocalPackets() throws Exception {
        var ctx = context("viewer", Map.of(), Map.of()); offlineClient.execute(ctx);
        Path packet = ctx.outputDirectory().resolve("handoff");
        String report = StudyBrowserLauncher.open(packet, "none");
        assertTrue(report.contains(packet.resolve("index.html").toString()));
        assertThrows(IllegalArgumentException.class, () -> StudyBrowserLauncher.open(packet, "chrome --remote-debugging-port=1234"));
        Files.writeString(packet.resolve("index.html"), "<script>unsafe</script>");
        assertThrows(IOException.class, () -> StudyBrowserLauncher.open(packet, "none"));
        assertThrows(IllegalArgumentException.class, () -> StudyBrowserLauncher.candidates("shell", "linux", Map.of()));
    }

    @Test void unchangedPipelineCommandConsumesReceiptAndReusesUpstreamOperations() throws Exception {
        Path p = root.resolve("pipeline"); Files.createDirectory(p);
        Files.writeString(p.resolve("a.tsv"), "feature_id\tbeta\tse\tp_value\nGENE1\t1\t.1\t.001\n");
        Files.writeString(p.resolve("b.tsv"), "feature_id\tbeta\tse\tp_value\nGENE1\t2\t.2\t.001\n");
        Files.writeString(p.resolve("study.yaml"), """
            schema_version: 1
            project: {id: browser-fixture, title: Browser fixture, description: Offline audited handoff test.}
            resources: {downloads: deny}
            assistance: {enabled: true, next_steps: true, response_files: [responses/review/receipt.yaml]}
            datasets:
              - {id: a, cohort: a, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: a.tsv}}
              - {id: b, cohort: b, type: transcriptomics, stage: association_results, feature_id_type: gene_symbol, files: {results: b.tsv}}
            hypotheses:
              - id: H1
                description: Molecular feature associates with the phenotype.
                outcome: BMI
                analysis: {mode: manual, goals: [meta_analysis]}
            """);
        assertTrue(StudyValidator.validate(p.resolve("study.yaml"), false).valid());
        assertEquals(0, cli("auto", p.resolve("study.yaml").toString(), "--offline"));
        Path packet;
        try (var files = Files.walk(p.resolve("runs"))) { packet = files.filter(f -> f.getFileName().toString().equals("packet.yaml")).findFirst().orElseThrow().getParent(); }
        Path response = save("pipeline-response.yaml", reply(packet, "{\"proposals\":[]}"));
        assertEquals(0, cli("handoff-import", response.toString(), "--handoff", packet.toString(), "--out", p.resolve("responses/review").toString(), "--actor-id", "scientist"));
        assertEquals(0, cli("auto", p.resolve("study.yaml").toString(), "--offline"));
        StudyProject.verify(p);
        assertTrue(Files.readString(p.resolve("ledger/events.jsonl")).contains("operation_reused"));
        assertTrue(Files.readString(p.resolve("SESSION_HISTORY.md")).contains("Imported human-attested browser response"));
        Map<String,Object> index = cast(SimpleJson.parse(Files.readString(p.resolve("index/index.json"))));
        var current = ((List<?>) index.get("entries")).stream().map(StudyBrowserHandoffTest::cast)
            .filter(e -> Boolean.TRUE.equals(e.get("current")) && "advisor_review".equals(e.get("operation_type"))).findFirst().orElseThrow();
        assertEquals("complete", current.get("state"));
        assertTrue(Files.readString(p.resolve("index/index.html")).contains("handoff/index.html"));
        // A modified companion must invalidate cache reuse even if receipt.yaml is unchanged.
        Files.writeString(p.resolve("responses/review/response.yaml"), "corrupt");
        assertEquals(0, cli("auto", p.resolve("study.yaml").toString(), "--offline"));
        assertTrue(Files.readString(p.resolve("index/index.json")).contains("browser_response_invalid"));
    }

    private int cli(String... args) {
        ByteArrayOutputStream errors = new ByteArrayOutputStream();
        int code = StudyCli.run(args, new PrintStream(new ByteArrayOutputStream()), new PrintStream(errors));
        assertEquals(0, code, errors.toString()); return code;
    }
    @SuppressWarnings("unchecked") private static Map<String,Object> cast(Object value) { return (Map<String,Object>) value; }
}
