/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** Deterministic rule planner over typed datasets and registered operations. */
public final class StudyPlanner {
    private StudyPlanner() { }

    public static StudyPlan compile(StudyContract resolved) {
        Map<String, Object> study = resolved.document();
        String resolvedHash = StudyProject.resolvedHash(study);
        Map<String, Object> project = StudyIO.map(study.get("project"), "project");
        Map<String, Map<String, Object>> datasets = byId(
            StudyIO.list(study.get("datasets"), "datasets"));
        List<Object> operations = new ArrayList<>();
        List<Object> decisions = pendingAssumptions(study);

        for (Object value : StudyIO.list(study.get("hypotheses"), "hypotheses")) {
            Map<String, Object> hypothesis = StudyIO.map(value, "hypothesis");
            if ("retired".equals(hypothesis.get("status"))) continue;
            String hypothesisId = StudyIO.requireString(
                hypothesis, "id", "hypothesis");
            Map<String, Object> inputs = StudyIO.map(
                hypothesis.get("inputs"), "hypothesis.inputs");
            List<String> selected = selected(inputs);
            Map<String, Object> analysis = StudyIO.map(
                hypothesis.get("analysis"), "hypothesis.analysis");
            Set<String> goals = goals(analysis);
            Map<String, Object> bioinformatics = StudyIO.optionalMap(
                analysis.get("bioinformatics"), "hypothesis.analysis.bioinformatics");
            boolean automaticColocalization = !StudyIO.optionalMap(
                bioinformatics.get("colocalization"), "hypothesis.analysis.bioinformatics.colocalization").isEmpty();
            if (automaticColocalization) {
                goals.add("source_discovery");
                goals.add("colocalization");
            }
            Map<String, Object> assist = StudyIO.optionalMap(study.get("assistance"), "assistance");
            if (Boolean.TRUE.equals(assist.get("enabled")) && Boolean.TRUE.equals(assist.get("literature")))
                goals.add("literature_support");
            Counter counter = new Counter(hypothesisId);

            Map<String, Object> harmonize = operation(counter, "validate_harmonize",
                hypothesisId, selected, List.of());
            operations.add(harmonize);
            List<String> resultOperations = new ArrayList<>();
            for (String datasetId : selected) {
                Map<String, Object> dataset = datasets.get(datasetId);
                if (dataset == null || "retired".equals(dataset.get("status"))) continue;
                String stage = String.valueOf(dataset.get("stage"));
                if (goals.contains("association") && Set.of(
                        "raw_quantified", "quantified_matrix",
                        "processed_matrix").contains(stage)) {
                    Map<String, Object> association = operation(counter,
                        "association", hypothesisId, List.of(datasetId),
                        List.of(String.valueOf(harmonize.get("id"))));
                    association.put("dataset_id", datasetId);
                    operations.add(association);
                    resultOperations.add(String.valueOf(association.get("id")));
                } else if (stage.endsWith("_results")) {
                    resultOperations.add("dataset:" + datasetId);
                }
            }
            String effectComparison = null;
            if (selected.size() > 1) {
                Map<String, Object> operation = operation(counter,
                    "effect_comparison", hypothesisId, selected,
                    dependencies(harmonize, resultOperations));
                operations.add(operation);
                effectComparison = String.valueOf(operation.get("id"));
            }
            String comparison = null;
            if (goals.contains("discovery_replication")) {
                Map<String, Object> operation = operation(counter,
                    "discovery_replication", hypothesisId, selected,
                    dependencies(harmonize, resultOperations));
                operation.put("roles", StudyIO.copyMap(inputs));
                operations.add(operation);
                comparison = String.valueOf(operation.get("id"));
            }
            String meta = null;
            if (goals.contains("meta_analysis")) {
                Map<String, Object> operation = operation(counter,
                    "meta_analysis", hypothesisId, selected,
                    dependencies(harmonize, resultOperations));
                operation.put("configuration",
                    operationConfiguration(analysis, "meta_analysis"));
                operations.add(operation);
                meta = String.valueOf(operation.get("id"));
            }
            List<String> followupDependencies = new ArrayList<>();
            if (meta != null) followupDependencies.add(meta);
            else if (comparison != null) followupDependencies.add(comparison);
            else {
                for (String result : resultOperations)
                    if (!result.startsWith("dataset:"))
                        followupDependencies.add(result);
                if (followupDependencies.isEmpty())
                    followupDependencies.add(String.valueOf(harmonize.get("id")));
            }
            List<String> evidenceOperations = new ArrayList<>();
            String irqtlReview = null;
            if (selected.stream().anyMatch(id -> StudyIO.optionalMap(datasets.get(id).get("files"), "files").containsKey("irqtl"))) {
                List<String> inputsForIrqtl = new ArrayList<>(followupDependencies);
                if (!inputsForIrqtl.contains(String.valueOf(harmonize.get("id")))) inputsForIrqtl.add(String.valueOf(harmonize.get("id")));
                Map<String, Object> irqtl = operation(counter, "irqtl_candidate_review", hypothesisId, selected, inputsForIrqtl);
                irqtl.put("allow_incomplete_dependencies", true);
                operations.add(irqtl);
                irqtlReview = String.valueOf(irqtl.get("id"));
                evidenceOperations.add(irqtlReview);
            }
            String sourceDiscovery = null;
            String phenotypeAdvice = null;
            if (automaticColocalization && Boolean.TRUE.equals(assist.get("enabled"))
                    && Boolean.TRUE.equals(assist.get("related_phenotypes"))) {
                List<String> inputsForAdvice = new ArrayList<>(followupDependencies);
                if (!inputsForAdvice.contains(String.valueOf(harmonize.get("id"))))
                    inputsForAdvice.add(String.valueOf(harmonize.get("id")));
                Map<String, Object> advice = operation(counter, "related_phenotype_advice", hypothesisId, selected, inputsForAdvice);
                advice.put("allow_incomplete_dependencies", true);
                operations.add(advice);
                phenotypeAdvice = String.valueOf(advice.get("id"));
                evidenceOperations.add(phenotypeAdvice);
            }
            for (String goal : List.of("source_discovery", "bioinformatics_investigation", "literature_search")) {
                if (!goals.contains(goal)) continue;
                List<String> annotated = new ArrayList<>(followupDependencies);
                if (!annotated.contains(String.valueOf(harmonize.get("id")))) annotated.add(String.valueOf(harmonize.get("id")));
                for (String result : resultOperations)
                    if (!result.startsWith("dataset:") && !annotated.contains(result)) annotated.add(result);
                if (goal.equals("source_discovery") && phenotypeAdvice != null && !annotated.contains(phenotypeAdvice))
                    annotated.add(phenotypeAdvice);
                Map<String, Object> followup = operation(counter, goal, hypothesisId, selected, annotated);
                followup.put("allow_incomplete_dependencies", true);
                operations.add(followup);
                evidenceOperations.add(String.valueOf(followup.get("id")));
                if (goal.equals("source_discovery")) sourceDiscovery = String.valueOf(followup.get("id"));
            }
            // Annotation assessment accompanies harmonized data. Preserve the single
            // upstream placeholder contract for explicitly configured CLI adapters.
            for (String goal : List.of("gene_ontology", "pathway_enrichment",
                    "variant_scoring", "fine_mapping", "colocalization",
                    "mendelian_randomization", "literature_support")) {
                if (!goals.contains(goal) && !(goals.contains("genetic_support")
                        && Set.of("variant_scoring", "fine_mapping",
                            "colocalization", "mendelian_randomization")
                            .contains(goal))) continue;
                Map<String, Object> operation = operation(counter, goal,
                    hypothesisId, selected, followupDependencies);
                operation.put("configuration", operationConfiguration(analysis, goal));
                if (!operationConfiguration(analysis, goal).containsKey("command")) {
                    List<String> annotated = new ArrayList<>(followupDependencies);
                    if (!annotated.contains(String.valueOf(harmonize.get("id"))))
                        annotated.add(String.valueOf(harmonize.get("id")));
                    for (String result : resultOperations)
                        if (!result.startsWith("dataset:") && !annotated.contains(result))
                            annotated.add(result);
                    if (irqtlReview != null && Set.of("fine_mapping", "colocalization", "mendelian_randomization").contains(goal)
                            && !annotated.contains(irqtlReview)) {
                        annotated.add(irqtlReview);
                        operation.put("allow_incomplete_dependencies", true);
                    }
                    operation.put("depends_on", annotated);
                    if (sourceDiscovery != null && Set.of("fine_mapping", "colocalization", "mendelian_randomization").contains(goal)) {
                        annotated.add(sourceDiscovery);
                        operation.put("allow_incomplete_dependencies", true);
                    }
                }
                operations.add(operation);
                evidenceOperations.add(String.valueOf(operation.get("id")));
            }
            List<String> reviewDependencies = new ArrayList<>(followupDependencies);
            boolean characteristicsRequested = !"manual".equals(analysis.get("mode")) || selected.stream().anyMatch(id ->
                StudyIO.optionalMap(datasets.get(id).get("columns"), "columns").containsKey("characteristics"));
            if (characteristicsRequested && selected.stream().anyMatch(id -> StudyIO.optionalMap(datasets.get(id).get("files"), "files").containsKey("phenotype"))) {
                Map<String, Object> characteristics = operation(counter, "participant_characteristics", hypothesisId, selected, List.of());
                operations.add(characteristics);
                reviewDependencies.add(String.valueOf(characteristics.get("id")));
            }
            if (!reviewDependencies.contains(String.valueOf(harmonize.get("id"))))
                reviewDependencies.add(String.valueOf(harmonize.get("id")));
            for (String result : resultOperations)
                if (!result.startsWith("dataset:") && !reviewDependencies.contains(result))
                    reviewDependencies.add(result);
            if (comparison != null && !reviewDependencies.contains(comparison))
                reviewDependencies.add(comparison);
            if (effectComparison != null
                    && !reviewDependencies.contains(effectComparison))
                reviewDependencies.add(effectComparison);
            if (meta != null && !reviewDependencies.contains(meta))
                reviewDependencies.add(meta);
            reviewDependencies.addAll(evidenceOperations);
            if (goals.contains("candidate_prioritization")) {
                Map<String, Object> priority = operation(counter,
                    "candidate_prioritization", hypothesisId, selected,
                    reviewDependencies);
                priority.put("allow_incomplete_dependencies", true);
                operations.add(priority);
                reviewDependencies.add(String.valueOf(priority.get("id")));
            }
            Map<String, Object> assistance = StudyIO.optionalMap(
                study.get("assistance"), "assistance");
            Set<String> selectedCohorts = new LinkedHashSet<>();
            for (String datasetId : selected) {
                Map<String, Object> dataset = datasets.get(datasetId);
                if (dataset != null && dataset.get("cohort") instanceof String cohort
                        && !cohort.isBlank()) selectedCohorts.add(cohort);
            }
            boolean automaticLeaveOneOut = meta != null
                && selectedCohorts.size() >= 3;
            boolean requestedSensitivity = !StudyIO.optionalList(
                assistance.get("allow_analyses"), "assistance.allow_analyses")
                .isEmpty();
            if (automaticLeaveOneOut || requestedSensitivity) {
                Map<String, Object> sensitivity = operation(counter,
                    "sensitivity_analysis", hypothesisId, selected,
                    reviewDependencies);
                sensitivity.put("allow_incomplete_dependencies", true);
                Map<String, Object> configuration = new LinkedHashMap<>(
                    operationConfiguration(analysis, "sensitivity_analysis"));
                configuration.put("automatic_leave_one_cohort_out",
                    automaticLeaveOneOut);
                sensitivity.put("configuration", configuration);
                operations.add(sensitivity);
                reviewDependencies.add(String.valueOf(sensitivity.get("id")));
            }
            Map<String, Object> diagnostics = operation(counter,
                "diagnostic_review", hypothesisId, selected,
                reviewDependencies.isEmpty()
                    ? List.of(String.valueOf(harmonize.get("id")))
                    : reviewDependencies);
            diagnostics.put("allow_incomplete_dependencies", true);
            operations.add(diagnostics);
            reviewDependencies.add(String.valueOf(diagnostics.get("id")));
            if (Boolean.TRUE.equals(assistance.get("enabled")) && Boolean.TRUE.equals(assistance.get("next_steps"))) {
                Map<String, Object> advisor = operation(counter, "advisor_review", hypothesisId, selected, reviewDependencies);
                advisor.put("allow_incomplete_dependencies", true);
                operations.add(advisor);
                reviewDependencies.add(String.valueOf(advisor.get("id")));
            }
            if (Boolean.TRUE.equals(assistance.get("enabled")) && Boolean.TRUE.equals(assistance.get("manuscript"))) {
                Map<String, Object> draft = operation(counter, "manuscript_draft", hypothesisId, selected, reviewDependencies);
                draft.put("allow_incomplete_dependencies", true);
                operations.add(draft);
                reviewDependencies.add(String.valueOf(draft.get("id")));
            }
            Map<String, Object> report = operation(counter, "report",
                hypothesisId, selected, reviewDependencies);
            report.put("allow_incomplete_dependencies", true);
            operations.add(report);
        }
        if (operations.isEmpty())
            throw new IllegalArgumentException("No active hypotheses can be planned");

        Map<String, Object> plan = new LinkedHashMap<>();
        plan.put("schema_version", 1);
        plan.put("project_id", project.get("id"));
        plan.put("resolved_study_hash", resolvedHash);
        plan.put("status", "proposed");
        plan.put("operations", operations);
        plan.put("pending_decisions", decisions);
        String identity = StudyIO.hash(plan);
        plan.put("plan_id", "plan-" + identity.substring(0, 12));
        plan.put("seed", Long.parseUnsignedLong(identity.substring(0, 15), 16));
        plan.put("plan_hash", hashWithoutPlanHash(plan));
        return new StudyPlan(plan);
    }

    public static StudyPlan create(Path project, StudyProject.Actor actor)
            throws IOException {
        return StudyProjectLock.with(project, () -> createLocked(project, actor));
    }

    private static StudyPlan createLocked(Path project, StudyProject.Actor actor) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        StudyPlan plan = compile(StudyProject.resolved(root));
        Path path = root.resolve("plans").resolve(plan.id())
            .resolve("analysis-plan.yaml");
        if (!Files.exists(path)) StudyIO.writeYaml(path, plan.document());
        Map<String, Object> status = StudyProject.status(root);
        Object active = status.get("active_plan");
        if (active instanceof Map<?, ?> map && plan.id().equals(map.get("plan_id"))
                && plan.hash().equals(map.get("plan_hash")))
            return plan;
        Map<String, Object> payload = new LinkedHashMap<>();
        payload.put("plan_id", plan.id());
        payload.put("plan_hash", plan.hash());
        payload.put("resolved_study_hash",
            plan.document().get("resolved_study_hash"));
        payload.put("state", "waiting_for_plan_approval");
        payload.put("phase", "planning");
        payload.put("summary", "Compiled deterministic analysis plan " + plan.id());
        StudyProject.append(root, "plan_created", actor, payload);
        StudyQuestions.reconcile(root, plan, null);
        StudyEvidenceIndex.rebuild(root);
        return plan;
    }

    public static void approve(Path project, String planId,
            StudyProject.Actor approver, boolean acceptAssumptions)
            throws IOException {
        StudyProjectLock.with(project, () -> { approveLocked(project, planId, approver, acceptAssumptions); return null; });
    }

    private static void approveLocked(Path project, String planId,
            StudyProject.Actor approver, boolean acceptAssumptions) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        Path path = root.resolve("plans").resolve(planId)
            .resolve("analysis-plan.yaml");
        Map<String, Object> plan = StudyIO.load(path);
        if (!planId.equals(plan.get("plan_id")))
            throw new IllegalArgumentException("Plan ID does not match its path");
        String calculated = hashWithoutPlanHash(plan);
        if (!calculated.equals(plan.get("plan_hash")))
            throw new IllegalArgumentException("Plan hash is invalid");
        Map<String, Object> status = StudyProject.status(root);
        Map<String, Object> active = StudyIO.map(
            status.get("active_plan"), "active_plan");
        if (!planId.equals(active.get("plan_id"))
                || !calculated.equals(active.get("plan_hash"))
                || !"proposed".equals(active.get("status")))
            throw new IllegalArgumentException(
                "Only the current proposed plan can be approved");
        String currentStudyHash = StudyProject.resolvedHash(
            StudyProject.resolved(root).document());
        if (!currentStudyHash.equals(plan.get("resolved_study_hash")))
            throw new IllegalArgumentException(
                "Plan was compiled from a superseded resolved study");
        List<Object> decisions = StudyIO.optionalList(
            plan.get("pending_decisions"), "pending_decisions");
        List<Object> approvableAssumptions = decisions.stream()
            .filter(item -> !Boolean.TRUE.equals(StudyIO.map(item, "decision").get("requires_explicit_review")))
            .toList();
        if (!approvableAssumptions.isEmpty() && !acceptAssumptions)
            throw new IllegalArgumentException(
                "Plan has pending assumptions; review them and use --accept-assumptions to approve explicitly");
        Map<String, Object> payload = new LinkedHashMap<>();
        payload.put("plan_id", planId);
        payload.put("plan_hash", plan.get("plan_hash"));
        payload.put("accepted_assumptions", acceptAssumptions
            ? approvableAssumptions.stream().map(item -> StudyIO.map(item, "decision").get("id")).toList()
            : List.of());
        List<Object> gatedAssumptions = decisions.stream()
            .filter(item -> Boolean.TRUE.equals(StudyIO.map(item, "decision").get("requires_explicit_review")))
            .map(item -> StudyIO.map(item, "decision").get("id")).toList();
        if (!gatedAssumptions.isEmpty()) payload.put("assumptions_deferred_to_scientist_questions", gatedAssumptions);
        payload.put("state", "ready");
        payload.put("phase", "execution");
        payload.put("summary", "Approved analysis plan " + planId);
        StudyProject.append(root, "plan_approved", approver, payload);
    }

    public static StudyPlan loadApproved(Path project) throws IOException {
        Path root = project.toAbsolutePath().normalize();
        Map<String, Object> status = StudyProject.status(root);
        Map<String, Object> active = StudyIO.map(
            status.get("active_plan"), "active_plan");
        if (!"approved".equals(active.get("status")))
            throw new IllegalArgumentException("The active plan is not approved");
        String id = String.valueOf(active.get("plan_id"));
        Map<String, Object> document = StudyIO.load(root.resolve("plans")
            .resolve(id).resolve("analysis-plan.yaml"));
        if (!hashWithoutPlanHash(document).equals(document.get("plan_hash")))
            throw new IllegalArgumentException("Approved plan hash is invalid");
        if (!String.valueOf(active.get("plan_hash")).equals(document.get("plan_hash")))
            throw new IllegalArgumentException("Approved plan differs from the stored plan");
        String currentStudyHash = StudyProject.resolvedHash(
            StudyProject.resolved(root).document());
        if (!currentStudyHash.equals(document.get("resolved_study_hash")))
            throw new IllegalArgumentException(
                "Approved plan was compiled from a superseded resolved study");
        return new StudyPlan(document);
    }

    static String hashWithoutPlanHash(Map<String, Object> plan) {
        Map<String, Object> copy = StudyIO.copyMap(plan);
        copy.remove("plan_hash");
        return StudyIO.hash(copy);
    }

    private static Map<String, Object> operation(Counter counter, String type,
            String hypothesis, List<String> datasets, List<String> dependencies) {
        StudyOperationRegistry.Descriptor descriptor =
            StudyOperationRegistry.descriptor(type);
        Map<String, Object> operation = new LinkedHashMap<>();
        operation.put("id", counter.next(type));
        operation.put("type", type);
        operation.put("hypothesis", hypothesis);
        operation.put("datasets", List.copyOf(datasets));
        operation.put("depends_on", List.copyOf(new LinkedHashSet<>(dependencies)));
        operation.put("executor", descriptor.executor());
        operation.put("outputs", descriptor.outputs());
        operation.put("diagnostics", descriptor.diagnostics());
        operation.put("status", "planned");
        return operation;
    }

    private static List<String> dependencies(Map<String, Object> harmonize,
            List<String> results) {
        List<String> dependencies = new ArrayList<>();
        dependencies.add(String.valueOf(harmonize.get("id")));
        for (String result : results)
            if (!result.startsWith("dataset:")) dependencies.add(result);
        return dependencies;
    }

    private static Map<String, Object> operationConfiguration(
            Map<String, Object> analysis, String type) {
        Map<String, Object> operations = StudyIO.optionalMap(
            analysis.get("operations"), "analysis.operations");
        return StudyIO.optionalMap(operations.get(type),
            "analysis.operations." + type);
    }

    private static Set<String> goals(Map<String, Object> analysis) {
        Set<String> result = new LinkedHashSet<>();
        for (Object value : StudyIO.list(analysis.get("goals"), "analysis.goals"))
            result.add(String.valueOf(value).toLowerCase(Locale.ROOT)
                .replace('-', '_'));
        return result;
    }

    private static List<String> selected(Map<String, Object> inputs) {
        List<String> result = new ArrayList<>();
        for (Object value : inputs.values()) {
            if (value instanceof String text) {
                if (!result.contains(text)) result.add(text);
            } else if (value instanceof Collection<?> values)
                for (Object item : values) {
                    String id = String.valueOf(item);
                    if (!result.contains(id)) result.add(id);
                }
        }
        return result;
    }

    private static Map<String, Map<String, Object>> byId(List<Object> values) {
        Map<String, Map<String, Object>> result = new LinkedHashMap<>();
        for (Object value : values) {
            Map<String, Object> item = StudyIO.map(value, "dataset");
            result.put(StudyIO.requireString(item, "id", "dataset"), item);
        }
        return result;
    }

    private static List<Object> pendingAssumptions(Map<String, Object> study) {
        List<Object> result = new ArrayList<>();
        for (Object value : StudyIO.optionalList(
                study.get("assumptions"), "assumptions")) {
            Map<String, Object> assumption = StudyIO.map(value, "assumption");
            if ("pending".equals(assumption.get("review_status"))) {
                Map<String, Object> decision = new LinkedHashMap<>();
                decision.put("id", assumption.get("id"));
                decision.put("type", "assumption_approval");
                decision.put("statement", assumption.get("statement"));
                decision.put("requires_explicit_review",
                    "potentially_dependent".equals(assumption.get("type")));
                result.add(decision);
            }
        }
        return result;
    }

    private static final class Counter {
        private final String hypothesis;
        private int value;

        private Counter(String hypothesis) {
            this.hypothesis = hypothesis;
        }

        String next(String type) {
            return hypothesis + "-op-" + String.format("%03d", ++value)
                + "-" + type;
        }
    }
}
