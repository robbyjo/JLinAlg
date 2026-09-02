/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Ridge-penalized one-hot basis equivalent to mgcv's {@code bs="re"}. */
public final class RandomEffectSmoothTerm {
    private RandomEffectSmoothTerm() { }

    public static QuadraticSmoothTerm of(String name, List<?> groups) {
        if (groups == null || groups.isEmpty()) {
            throw new IllegalArgumentException("groups are required");
        }
        Map<Object, Integer> levels = new LinkedHashMap<>();
        for (Object group : groups) {
            if (group == null) throw new IllegalArgumentException("groups must not be null");
            levels.computeIfAbsent(group, ignored -> levels.size());
        }
        int columns = levels.size();
        double[] design = new double[groups.size() * columns];
        for (int row = 0; row < groups.size(); row++) {
            design[row * columns + levels.get(groups.get(row))] = 1.0;
        }
        double[] penalty = new double[columns * columns];
        for (int column = 0; column < columns; column++) {
            penalty[column * columns + column] = 1.0;
        }
        return new QuadraticSmoothTerm(name, groups.size(), columns,
            design, List.of(penalty));
    }
}
