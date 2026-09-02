/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.gam.PSplineTerm;
import org.jlinalg.gam.QuadraticSmoothTerm;
import org.jlinalg.gam.TensorInteractionPSplineTerm;
import org.jlinalg.gam.TensorProductPSplineTerm;

/** R-like compiler for multi-penalty additive formulas. */
public final class AdvancedGamFormula {
    private AdvancedGamFormula() { }

    public static CompiledQuadraticGamFormula compile(
            String formula, ModelTable table) {
        return compile(formula, table, FormulaOptions.defaults());
    }

    /** Compiles fixed terms plus {@code s(x)}, {@code te(x,z)}, and {@code ti(x,z)}. */
    public static CompiledQuadraticGamFormula compile(
            String formula, ModelTable table, FormulaOptions options) {
        if (formula == null || table == null || options == null) {
            throw new IllegalArgumentException("formula, table, and options are required");
        }
        String compact = formula.replaceAll("\\s+", "");
        String[] sides = compact.split("~", -1);
        if (sides.length != 2 || sides[0].isEmpty() || sides[1].isEmpty()) {
            throw new IllegalArgumentException("formula must have the form response ~ terms");
        }
        List<String> fixedTokens = new ArrayList<>();
        List<QuadraticSmoothTerm> smooths = new ArrayList<>();
        for (String token : splitTopLevel(sides[1], '+')) {
            if (token.startsWith("s(") && token.endsWith(")")) {
                SmoothSpec spec = parseSmooth(token);
                PSplineTerm term = PSplineTerm.of("s(" + spec.first() + ")",
                    table.numeric(spec.first()), spec.firstDimension(),
                    spec.degree(), spec.differenceOrder());
                smooths.add(QuadraticSmoothTerm.from(term));
            } else if ((token.startsWith("te(") || token.startsWith("ti("))
                    && token.endsWith(")")) {
                SmoothSpec spec = parseTensor(token);
                double[] first = table.numeric(spec.first());
                double[] second = table.numeric(spec.second());
                String name = token.substring(0, token.indexOf('('))
                    + "(" + spec.first() + "," + spec.second() + ")";
                smooths.add(token.startsWith("te(")
                    ? TensorProductPSplineTerm.of(name, first, second,
                        spec.firstDimension(), spec.secondDimension())
                    : TensorInteractionPSplineTerm.of(name, first, second,
                        spec.firstDimension(), spec.secondDimension()));
            } else {
                fixedTokens.add(token);
            }
        }
        if (smooths.isEmpty()) {
            throw new IllegalArgumentException("an additive formula needs a smooth term");
        }
        String fixedRight = fixedTokens.isEmpty() ? "1" : String.join("+", fixedTokens);
        CompiledFormula fixed = Formula.compile(
            sides[0] + "~" + fixedRight, table, options);
        return new CompiledQuadraticGamFormula(fixed, smooths);
    }

    private static SmoothSpec parseSmooth(String token) {
        List<String> arguments = splitTopLevel(
            token.substring(2, token.length() - 1), ',');
        if (arguments.isEmpty() || arguments.get(0).isEmpty()) {
            throw new IllegalArgumentException("smooth covariate is required");
        }
        int dimension = 10;
        int degree = 3;
        int order = 2;
        for (int index = 1; index < arguments.size(); index++) {
            String[] pair = pair(arguments.get(index));
            switch (pair[0]) {
                case "k" -> dimension = positive(pair[1], "k");
                case "degree" -> degree = positive(pair[1], "degree");
                case "m" -> order = positive(pair[1], "m");
                case "bs" -> requirePs(pair[1]);
                default -> throw new IllegalArgumentException(
                    "unknown smooth option: " + pair[0]);
            }
        }
        return new SmoothSpec(arguments.get(0), null, dimension, 0, degree, order);
    }

    private static SmoothSpec parseTensor(String token) {
        int prefix = token.indexOf('(') + 1;
        List<String> arguments = splitTopLevel(
            token.substring(prefix, token.length() - 1), ',');
        if (arguments.size() < 2 || arguments.get(0).isEmpty()
                || arguments.get(1).isEmpty()) {
            throw new IllegalArgumentException("tensor smooth requires two covariates");
        }
        int firstDimension = 7;
        int secondDimension = 7;
        for (int index = 2; index < arguments.size(); index++) {
            String[] pair = pair(arguments.get(index));
            switch (pair[0]) {
                case "k", "kx" -> firstDimension = positive(pair[1], pair[0]);
                case "kz" -> secondDimension = positive(pair[1], "kz");
                case "bs" -> requirePs(pair[1]);
                default -> throw new IllegalArgumentException(
                    "unknown tensor option: " + pair[0]);
            }
        }
        if (arguments.stream().noneMatch(value -> value.startsWith("kz="))) {
            secondDimension = firstDimension;
        }
        return new SmoothSpec(arguments.get(0), arguments.get(1),
            firstDimension, secondDimension, 3, 2);
    }

    private static String[] pair(String option) {
        String[] result = option.split("=", -1);
        if (result.length != 2 || result[1].isEmpty()) {
            throw new IllegalArgumentException("invalid smooth option: " + option);
        }
        return result;
    }

    private static void requirePs(String value) {
        String unquoted = value.replace("\"", "").replace("'", "");
        if (!"ps".equals(unquoted)) {
            throw new IllegalArgumentException("only bs='ps' is implemented here");
        }
    }

    private static int positive(String value, String name) {
        try {
            int result = Integer.parseInt(value);
            if (result < 1) throw new NumberFormatException();
            return result;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(name + " must be positive", exception);
        }
    }

    private static List<String> splitTopLevel(String expression, char separator) {
        List<String> result = new ArrayList<>();
        int depth = 0;
        int start = 0;
        for (int index = 0; index < expression.length(); index++) {
            char character = expression.charAt(index);
            if (character == '(') depth++;
            else if (character == ')') depth--;
            else if (character == separator && depth == 0) {
                result.add(expression.substring(start, index));
                start = index + 1;
            }
            if (depth < 0) throw new IllegalArgumentException("unbalanced parentheses");
        }
        if (depth != 0) throw new IllegalArgumentException("unbalanced parentheses");
        result.add(expression.substring(start));
        return result;
    }

    private record SmoothSpec(
            String first,
            String second,
            int firstDimension,
            int secondDimension,
            int degree,
            int differenceOrder) { }
}
