/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.CorrelatedRandomEffectBlock;

/** Compiler for lme4-style independent random-intercept and slope terms. */
public final class MixedFormula {
    private static final Pattern RANDOM_TERM =
        Pattern.compile("\\(([^|]+)(\\|\\|?)([^()]+)\\)");

    private MixedFormula() { }

    /**
     * Compiles terms such as {@code y ~ x + (1|family) + (0+x|family)}.
     * The compiled sparse random designs are reused by every subsequent fit.
     */
    public static CompiledMixedFormula compile(String formula, ModelTable table) {
        if (formula == null || table == null) {
            throw new IllegalArgumentException("formula and table are required");
        }
        Matcher matcher = RANDOM_TERM.matcher(formula.replaceAll("\\s+", ""));
        List<RandomEffectTerm> random = new ArrayList<>();
        List<CorrelatedRandomEffectBlock> correlated = new ArrayList<>();
        StringBuffer fixedFormula = new StringBuffer();
        while (matcher.find()) {
            String expression = matcher.group(1);
            boolean explicitlyIndependent = matcher.group(2).equals("||");
            String groupExpression = matcher.group(3);
            for (String groupName : expandedGroups(groupExpression)) {
                List<String> groups = groupLabels(groupName, table);
                addTerms(random, correlated, expression, groupName, groups,
                    explicitlyIndependent, table);
            }
            matcher.appendReplacement(fixedFormula, "");
        }
        matcher.appendTail(fixedFormula);
        if (random.isEmpty() && correlated.isEmpty()) {
            throw new IllegalArgumentException("mixed formula has no random-effect term");
        }
        String cleaned = fixedFormula.toString()
            .replaceAll("\\++", "+")
            .replace("~+", "~")
            .replaceAll("\\+$", "");
        CompiledFormula fixed = Formula.compile(cleaned, table);
        return new CompiledMixedFormula(fixed, random, correlated);
    }

    private static void addTerms(
            List<RandomEffectTerm> destination,
            List<CorrelatedRandomEffectBlock> correlated,
            String expression,
            String groupName,
            List<String> groups,
            boolean independent,
            ModelTable table) {
        String[] tokens = expression.split("\\+", -1);
        boolean intercept = !expression.startsWith("0+");
        List<String> slopes = new ArrayList<>();
        for (String token : tokens) {
            if (token.isEmpty() || token.equals("0")) continue;
            if (token.equals("1")) intercept = true;
            else slopes.add(token);
        }
        if (!independent && intercept && !slopes.isEmpty()) {
            int effects = slopes.size() + 1;
            double[][] design = new double[table.rows()][effects];
            List<String> effectNames = new ArrayList<>(effects);
            effectNames.add("(Intercept)");
            for (int row = 0; row < table.rows(); row++) design[row][0] = 1.0;
            for (int slopeIndex = 0; slopeIndex < slopes.size(); slopeIndex++) {
                String slope = slopes.get(slopeIndex);
                if (!table.isNumeric(slope))
                    throw new IllegalArgumentException(
                        "random slope must be numeric: " + slope);
                effectNames.add(slope);
                double[] values = table.numeric(slope);
                for (int row = 0; row < table.rows(); row++)
                    design[row][slopeIndex + 1] = values[row];
            }
            correlated.add(CorrelatedRandomEffectBlock.of(
                expression + "|" + groupName, groups, effectNames, design));
            return;
        }
        if (intercept)
            destination.add(RandomEffectTerm.randomIntercept(
                "1|" + groupName, groups));
        for (String slope : slopes) {
            if (!table.isNumeric(slope))
                throw new IllegalArgumentException(
                    "random slope must be numeric: " + slope);
            destination.add(RandomEffectTerm.randomSlope(
                "0+" + slope + "|" + groupName, groups,
                table.numeric(slope)));
        }
        if (!intercept && slopes.isEmpty())
            throw new IllegalArgumentException(
                "random-effect expression produces no terms");
    }

    private static List<String> expandedGroups(String expression) {
        String[] nested = expression.split("/", -1);
        if (nested.length == 1) return List.of(expression);
        List<String> result = new ArrayList<>(nested.length);
        String current = nested[0];
        if (current.isEmpty())
            throw new IllegalArgumentException("nested group is empty");
        result.add(current);
        for (int index = 1; index < nested.length; index++) {
            if (nested[index].isEmpty())
                throw new IllegalArgumentException("nested group is empty");
            current += ":" + nested[index];
            result.add(current);
        }
        return result;
    }

    private static List<String> groupLabels(
            String expression, ModelTable table) {
        String[] factors = expression.split(":", -1);
        String[][] values = new String[factors.length][];
        for (int factor = 0; factor < factors.length; factor++) {
            if (!table.isCategorical(factors[factor]))
                throw new IllegalArgumentException(
                    "random-effect group must be categorical: "
                        + factors[factor]);
            values[factor] = table.categorical(factors[factor]);
        }
        String[] combined = new String[table.rows()];
        for (int row = 0; row < table.rows(); row++) {
            StringBuilder label = new StringBuilder();
            for (int factor = 0; factor < values.length; factor++) {
                if (values[factor][row] == null)
                    throw new IllegalArgumentException(
                        "missing random-effect group at row " + row);
                if (factor > 0) label.append(':');
                label.append(values[factor][row]);
            }
            combined[row] = label.toString();
        }
        return Arrays.asList(combined);
    }
}
