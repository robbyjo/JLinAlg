/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.util.List;

/** Structural and material validation findings for a study contract. */
public record StudyValidation(List<Issue> issues) {
    public StudyValidation {
        issues = List.copyOf(issues);
    }

    public enum Severity { ERROR, WARNING }

    public record Issue(Severity severity, String path, String message) {
        public Issue {
            if (severity == null || path == null || message == null)
                throw new IllegalArgumentException("validation issue fields are required");
        }
    }

    public boolean valid() {
        return issues.stream().noneMatch(issue -> issue.severity() == Severity.ERROR);
    }

    public List<Issue> errors() {
        return issues.stream()
            .filter(issue -> issue.severity() == Severity.ERROR).toList();
    }

    public List<Issue> warnings() {
        return issues.stream()
            .filter(issue -> issue.severity() == Severity.WARNING).toList();
    }

    public String report() {
        if (issues.isEmpty()) return "valid: no validation findings";
        StringBuilder result = new StringBuilder();
        for (Issue issue : issues)
            result.append(issue.severity().name().toLowerCase())
                .append('\t').append(issue.path()).append('\t')
                .append(issue.message()).append('\n');
        return result.toString();
    }
}
