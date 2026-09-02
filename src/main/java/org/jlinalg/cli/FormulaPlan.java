/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** Parsed CLI formula extensions around the library's R-style compiler. */
final class FormulaPlan {
    private static final String OMICS = "<omics>";
    private static final Pattern SURV = Pattern.compile(
        "(?i)^Surv\\(([^,()]+),([^,()]+)(?:,([^,()]+))?\\)$");
    private final String original;
    private final String withoutOmics;
    private final String response;
    private final Survival survival;
    private final boolean omics;
    private final boolean random;

    private FormulaPlan(
            String original, String withoutOmics, String response,
            Survival survival, boolean omics, boolean random) {
        this.original = original;
        this.withoutOmics = withoutOmics;
        this.response = response;
        this.survival = survival;
        this.omics = omics;
        this.random = random;
    }

    static FormulaPlan parse(String formula) {
        int separator = formula.indexOf('~');
        if (separator < 1 || formula.indexOf('~', separator + 1) >= 0)
            throw new IllegalArgumentException(
                "formula must contain exactly one '~'");
        String compact = formula.replaceAll("\\s+", "");
        int first = compact.indexOf(OMICS);
        boolean hasOmics = first >= 0;
        if (hasOmics && compact.indexOf(OMICS, first + OMICS.length()) >= 0)
            throw new IllegalArgumentException(
                "formula may contain <omics> only once");
        if (compact.contains(OMICS + ":") || compact.contains(":" + OMICS)
                || compact.contains(OMICS + "*") || compact.contains("*" + OMICS))
            throw new IllegalArgumentException(
                "the initial CLI supports <omics> as a main effect only");
        String cleaned = compact.replace(OMICS, "")
            .replaceAll("\\++", "+")
            .replace("~+", "~")
            .replaceAll("\\+$", "");
        if (cleaned.endsWith("~")) cleaned += "1";
        String left = compact.substring(0, compact.indexOf('~'));
        Matcher matcher = SURV.matcher(left);
        Survival survival = null;
        String response = left;
        if (matcher.matches()) {
            String firstTime = matcher.group(1).trim();
            String second = matcher.group(2).trim();
            String third = matcher.group(3);
            survival = third == null
                ? new Survival(null, firstTime, second)
                : new Survival(firstTime, second, third.trim());
            response = survival.stop();
        }
        return new FormulaPlan(formula, cleaned, response, survival,
            hasOmics, compact.matches(".*\\([^()]*\\|\\|?[^()]+\\).*"));
    }

    String original() { return original; }
    String withoutOmics() { return withoutOmics; }
    String response() { return response; }
    Survival survival() { return survival; }
    boolean hasOmics() { return omics; }
    boolean hasRandomEffects() { return random; }
    boolean isCox() { return survival != null; }

    record Survival(String start, String stop, String event) { }
}
