/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.study;

import java.nio.file.Path;
import java.util.Map;

/** Validated human or resolved study document with a canonical identity. */
public final class StudyContract {
    private final Path source;
    private final Map<String, Object> document;
    private final String hash;

    StudyContract(Path source, Map<String, Object> document) {
        this.source = source == null ? null : source.toAbsolutePath().normalize();
        this.document = StudyIO.copyMap(document);
        this.hash = StudyIO.hash(document);
    }

    public Path source() {
        return source;
    }

    public String projectId() {
        return StudyIO.requireString(
            StudyIO.map(document.get("project"), "project"), "id", "project");
    }

    public String hash() {
        return hash;
    }

    /** Returns a defensive deep copy suitable for inspection or extension. */
    public Map<String, Object> document() {
        return StudyIO.copyMap(document);
    }
}
