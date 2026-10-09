/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.util.Map;

/** Immutable deterministic analysis plan. */
public final class StudyPlan {
    private final Map<String, Object> document;

    StudyPlan(Map<String, Object> document) {
        this.document = StudyIO.copyMap(document);
    }

    public Map<String, Object> document() {
        return StudyIO.copyMap(document);
    }

    public String id() {
        return String.valueOf(document.get("plan_id"));
    }

    public String hash() {
        return String.valueOf(document.get("plan_hash"));
    }
}
