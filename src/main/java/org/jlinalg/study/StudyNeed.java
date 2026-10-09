/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** A typed scientific need. Its response is never an executable instruction. */
public record StudyNeed(String rule, String key, String prompt, String why,
        String responseType, List<String> options, Map<String, Object> target,
        boolean blocking) {
    public StudyNeed {
        for (String value : List.of(rule, key, prompt, why, responseType))
            if (value.isBlank() || value.length() > 8000)
                throw new IllegalArgumentException("Need fields must be nonblank and bounded");
        if (!Set.of("text", "single_choice", "local_file").contains(responseType))
            throw new IllegalArgumentException("Unsupported question response type: " + responseType);
        options = List.copyOf(options);
        target = StudyIO.copyMap(target);
        if (responseType.equals("single_choice") && options.isEmpty())
            throw new IllegalArgumentException("Choice questions require options");
    }

    @Override public Map<String, Object> target() { return StudyIO.copyMap(target); }

    public Map<String, Object> document() {
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("rule", rule);
        result.put("rule_version", 1);
        result.put("semantic_key", key);
        result.put("prompt", prompt);
        result.put("why", why);
        result.put("response", Map.of("type", responseType, "options", options,
            "allow_unknown", true, "allow_defer", true, "allow_decline", !blocking));
        result.put("target", target());
        result.put("blocking", blocking);
        result.put("consequences", "An answer records the scientist's statement; it does not authorize execution. "
            + (target.isEmpty() ? "Review and amend the plan or supply the missing inputs before retrying."
                : "A supported answer creates a proposed amendment, requiring separate scientist approval."));
        return result;
    }

    public static StudyNeed dataset(String id, String field, String prompt, String why,
            String responseType, List<String> options) {
        return new StudyNeed("dataset." + field, "dataset:" + id + ":" + field,
            prompt, why, responseType, options,
            Map.of("collection", "datasets", "id", id, "field", field), true);
    }

    public static StudyNeed review(String key, String prompt, String why, boolean blocking) {
        return new StudyNeed("scientific_review", key, prompt, why, "text", List.of(), Map.of(), blocking);
    }

    static StudyNeed from(Map<String, Object> record) {
        Map<String, Object> response = StudyIO.map(record.get("response"), "need.response");
        return new StudyNeed(StudyIO.requireString(record, "rule", "need"),
            StudyIO.requireString(record, "semantic_key", "need"),
            StudyIO.requireString(record, "prompt", "need"),
            StudyIO.requireString(record, "why", "need"),
            StudyIO.requireString(response, "type", "need.response"),
            StudyIO.list(response.get("options"), "options").stream().map(String::valueOf).toList(),
            StudyIO.map(record.get("target"), "need.target"), Boolean.TRUE.equals(record.get("blocking")));
    }
}
