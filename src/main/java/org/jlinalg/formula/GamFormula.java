/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.gam.PSplineTerm;

/** Compiler for fixed terms plus univariate {@code s()} P-spline terms. */
public final class GamFormula {
    private GamFormula() { }

    /** Compiles with treatment contrasts and no weight column. */
    public static CompiledGamFormula compile(String formula, ModelTable table) {
        return compile(formula, table, FormulaOptions.defaults());
    }

    /**
     * Compiles expressions such as {@code y ~ sex + s(age, k=12) + offset(e)}.
     * Supported smooth options are {@code k}, {@code degree}, and {@code m}.
     */
    public static CompiledGamFormula compile(
            String formula, ModelTable table, FormulaOptions options) {
        if (formula == null || table == null || options == null) {
            throw new IllegalArgumentException(
                "formula, table, and options are required");
        }
        String compact = formula.replaceAll("\\s+", "");
        String[] sides = compact.split("~", -1);
        if (sides.length != 2 || sides[0].isEmpty() || sides[1].isEmpty()) {
            throw new IllegalArgumentException(
                "formula must have the form response ~ terms");
        }
        List<String> tokens = splitTopLevel(sides[1]);
        List<SmoothSpec> smoothSpecs = new ArrayList<>();
        List<String> fixedTokens = new ArrayList<>();
        for (String token : tokens) {
            if (token.startsWith("s(") && token.endsWith(")")) {
                smoothSpecs.add(parseSmooth(token));
            } else {
                fixedTokens.add(token);
            }
        }
        if (smoothSpecs.isEmpty()) {
            throw new IllegalArgumentException(
                "a GAM formula requires at least one s() term");
        }
        String fixedRight = fixedTokens.isEmpty()
            ? "1" : String.join("+", fixedTokens);
        CompiledFormula fixed = Formula.compile(
            sides[0] + "~" + fixedRight, table, options);
        List<PSplineTerm> smoothTerms = new ArrayList<>(smoothSpecs.size());
        for (SmoothSpec spec : smoothSpecs) {
            smoothTerms.add(PSplineTerm.of(
                "s(" + spec.column() + ")", table.numeric(spec.column()),
                spec.basisDimension(), spec.degree(), spec.differenceOrder()));
        }
        return new CompiledGamFormula(fixed, smoothTerms);
    }

    private static SmoothSpec parseSmooth(String token) {
        String inside = token.substring(2, token.length() - 1);
        List<String> arguments = splitTopLevel(inside);
        if (arguments.isEmpty() || arguments.get(0).isEmpty()) {
            throw new IllegalArgumentException("smooth covariate is required");
        }
        String column = arguments.get(0);
        int basisDimension = 10;
        int degree = 3;
        int differenceOrder = 2;
        for (int index = 1; index < arguments.size(); index++) {
            String argument = arguments.get(index);
            String[] pair = argument.split("=", -1);
            if (pair.length != 2 || pair[1].isEmpty()) {
                throw new IllegalArgumentException(
                    "invalid smooth option: " + argument);
            }
            switch (pair[0]) {
                case "k" -> basisDimension = parsePositive(pair[1], "k");
                case "degree" -> degree = parsePositive(pair[1], "degree");
                case "m" -> differenceOrder = parsePositive(pair[1], "m");
                case "bs" -> {
                    String value = pair[1].replace("\"", "")
                        .replace("'", "");
                    if (!"ps".equals(value)) {
                        throw new IllegalArgumentException(
                            "only bs='ps' is currently implemented");
                    }
                }
                default -> throw new IllegalArgumentException(
                    "unknown smooth option: " + pair[0]);
            }
        }
        return new SmoothSpec(
            column, basisDimension, degree, differenceOrder);
    }

    private static int parsePositive(String value, String name) {
        try {
            int parsed = Integer.parseInt(value);
            if (parsed < 1) throw new NumberFormatException();
            return parsed;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                name + " must be a positive integer", exception);
        }
    }

    private static List<String> splitTopLevel(String expression) {
        List<String> result = new ArrayList<>();
        int depth = 0;
        int start = 0;
        for (int index = 0; index < expression.length(); index++) {
            char character = expression.charAt(index);
            if (character == '(') depth++;
            else if (character == ')') depth--;
            else if ((character == '+' || character == ',') && depth == 0) {
                result.add(expression.substring(start, index));
                start = index + 1;
            }
            if (depth < 0) {
                throw new IllegalArgumentException("unbalanced formula parentheses");
            }
        }
        if (depth != 0) {
            throw new IllegalArgumentException("unbalanced formula parentheses");
        }
        result.add(expression.substring(start));
        return result;
    }

    private record SmoothSpec(
            String column,
            int basisDimension,
            int degree,
            int differenceOrder) { }
}
