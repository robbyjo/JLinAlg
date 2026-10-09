/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.nio.file.Path;
import java.util.List;
import java.util.Map;

/** Adapter from a typed study operation to an auditable concrete execution. */
@FunctionalInterface
public interface StudyOperationExecutor {
    Result execute(Context context) throws Exception;

    record Context(Path project, Path outputDirectory,
            Map<String, Object> operation, Map<String, Object> study,
            Map<String, Object> plan, Map<String, Path> dependencyDirectories) {
        public Context {
            operation = Map.copyOf(operation);
            study = Map.copyOf(study);
            plan = Map.copyOf(plan);
            dependencyDirectories = Map.copyOf(dependencyDirectories);
        }
    }

    record Result(State state, String summary, List<String> outputs,
            Map<String, Object> details) {
        public Result {
            if (state == null || summary == null)
                throw new IllegalArgumentException("operation state and summary are required");
            outputs = outputs == null ? List.of() : List.copyOf(outputs);
            details = details == null ? Map.of() : Map.copyOf(details);
        }

        public static Result complete(String summary, List<String> outputs) {
            return new Result(State.COMPLETE, summary, outputs, Map.of());
        }

        public static Result waiting(String summary, Map<String, Object> details) {
            return new Result(State.WAITING_FOR_INPUT, summary, List.of(), details);
        }

        public static Result waiting(String summary, Map<String, Object> details, List<StudyNeed> needs) {
            Map<String, Object> typed = new java.util.LinkedHashMap<>(details);
            typed.put("needs", needs.stream().map(StudyNeed::document).toList());
            return waiting(summary, typed);
        }
    }

    enum State {
        COMPLETE("complete"),
        NOT_APPLICABLE("not_applicable"),
        WAITING_FOR_INPUT("waiting_for_input"),
        FAILED("failed");

        private final String value;

        State(String value) {
            this.value = value;
        }

        public String value() {
            return value;
        }
    }
}
