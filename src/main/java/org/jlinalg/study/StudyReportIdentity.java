/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.study;

import java.util.LinkedHashMap;
import java.util.Map;

/** Privacy-preserving identity for a reported participant selection. */
public final class StudyReportIdentity {
    private StudyReportIdentity() { }
    public static String selection(Map<String, Object> characteristics) {
        Map<String, Object> selected = new LinkedHashMap<>();
        for (String key : new String[]{"id", "filters", "membership_role", "membership_id_column"})
            if (characteristics.containsKey(key)) selected.put(key, characteristics.get(key));
        return StudyIO.hash(selected);
    }
}
