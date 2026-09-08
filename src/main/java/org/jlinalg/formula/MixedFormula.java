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

/** Compiler for lme4-style grouped blocks with numeric/categorical slopes and interactions. */
public final class MixedFormula {
    private static final Pattern RANDOM_TERM =
        Pattern.compile("\\(([^|()]+)(\\|\\|?)([^()]+)\\)");

    private MixedFormula() { }

    /**
     * Compiles terms such as {@code y ~ x + (1|family) + (0+x|family)}.
     * The compiled sparse random designs are reused by every subsequent fit.
     */
    public static CompiledMixedFormula compile(String formula, ModelTable table) {
        return compile(formula, table, FormulaOptions.defaults());
    }

    /** Compiles fixed-effect contrasts and residual precision weights as well as random blocks. */
    public static CompiledMixedFormula compile(String formula, ModelTable table, FormulaOptions options) {
        if (formula == null || table == null || options == null) {
            throw new IllegalArgumentException("formula and table are required");
        }
        Matcher matcher = RANDOM_TERM.matcher(formula.replaceAll("\\s+", ""));
        String responseName = formula.split("~", -1)[0].trim();
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
                    explicitlyIndependent, table, responseName, options);
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
        if (cleaned.endsWith("~")) cleaned += "1";
        CompiledFormula fixed = Formula.compile(cleaned, table, options);
        return new CompiledMixedFormula(fixed, random, correlated);
    }

    private static void addTerms(
            List<RandomEffectTerm> destination,
            List<CorrelatedRandomEffectBlock> correlated,
            String expression,
            String groupName,
            List<String> groups,
            boolean independent,
            ModelTable table, String responseName, FormulaOptions options) {
        RandomTerms terms = randomTerms(expression);
        if (!independent) {
            addBlock(destination, correlated, expression + "|" + groupName, groupName, groups,
                randomDesign(terms, table, responseName, options));
        } else {
            // lme4 || splits formula TERMS. A factor's full indicator columns
            // remain a correlated block, while numeric terms remain scalar.
            if (terms.intercept()) destination.add(RandomEffectTerm.randomIntercept("1|" + groupName, groups));
            for (List<String> term : terms.terms()) {
                String name = String.join(":", term);
                addBlock(destination, correlated, "0+" + name + "|" + groupName, groupName, groups,
                    randomDesign(new RandomTerms(false, List.of(term)), table, responseName, options));
            }
        }
    }

    private record RandomTerms(boolean intercept, List<List<String>> terms) { }
    private record Encoded(List<String> names, double[] values) {
        int columns() { return names.size(); }
    }

    private static RandomTerms randomTerms(String expression) {
        boolean intercept = true;
        List<String> variableOrder = new ArrayList<>();
        List<List<String>> expanded = new ArrayList<>();
        for (String token : expression.replace("-1", "+-1").split("\\+")) {
            if (token.isEmpty()) continue;
            if (token.equals("0") || token.equals("-1")) { intercept = false; continue; }
            if (token.equals("1")) { intercept = true; continue; }
            String[] factors = token.split("\\*", -1);
            if (factors.length > 10) throw new IllegalArgumentException("random interaction expands too many terms");
            for (String factor : factors) for (String variable : factor.split(":", -1)) {
                if (variable.isEmpty()) throw new IllegalArgumentException("empty random interaction variable");
                if (!variableOrder.contains(variable)) variableOrder.add(variable);
            }
            for (int mask = 1; mask < (1 << factors.length); mask++) {
                List<String> term = new ArrayList<>();
                for (int i = 0; i < factors.length; i++) if ((mask & (1 << i)) != 0)
                    for (String variable : factors[i].split(":", -1)) if (!term.contains(variable)) term.add(variable);
                term.sort(java.util.Comparator.comparingInt(variableOrder::indexOf));
                if (!expanded.contains(term)) expanded.add(List.copyOf(term));
            }
        }
        // R orders main effects before interactions, retaining encounter order within degree.
        expanded.sort(java.util.Comparator.comparingInt(List::size));
        if (!intercept && expanded.isEmpty()) throw new IllegalArgumentException("random expression produces no columns");
        return new RandomTerms(intercept, List.copyOf(expanded));
    }

    private static Encoded randomDesign(RandomTerms specification, ModelTable table,
            String responseName, FormulaOptions options) {
        List<Encoded> parts = new ArrayList<>();
        if (specification.intercept()) {
            double[] ones = new double[table.rows()]; Arrays.fill(ones, 1);
            parts.add(new Encoded(List.of("(Intercept)"), ones));
        }
        boolean firstFactor = !specification.intercept();
        List<List<String>> previous = new ArrayList<>();
        for (List<String> term : specification.terms()) {
            Encoded interaction = null;
            for (String variable : term) {
                List<String> marginal = new ArrayList<>(term); marginal.remove(variable);
                boolean contrast = marginal.isEmpty() || previous.stream().anyMatch(t -> t.containsAll(marginal));
                boolean categorical = table.isCategorical(variable);
                boolean full = !contrast || (categorical && firstFactor);
                if (categorical) firstFactor = false;
                Encoded encoded = encodeRandom(variable, full, table, responseName, options);
                interaction = interaction == null ? encoded : interactRandom(interaction, encoded, table.rows());
            }
            parts.add(interaction); previous.add(term);
        }
        int columns = parts.stream().mapToInt(Encoded::columns).sum();
        double[] design = new double[table.rows() * columns];
        List<String> names = new ArrayList<>(); int start = 0;
        for (Encoded part : parts) {
            names.addAll(part.names());
            for (int r = 0; r < table.rows(); r++)
                System.arraycopy(part.values(), r * part.columns(), design, r * columns + start, part.columns());
            start += part.columns();
        }
        return new Encoded(List.copyOf(names), design);
    }

    private static Encoded encodeRandom(String name, boolean full, ModelTable table,
            String responseName, FormulaOptions options) {
        if (!full || !table.isCategorical(name)) {
            // Reuse the established numeric/treatment/sum contrast encoder.
            CompiledFormula contrast = Formula.compile(responseName + "~0+" + name, table,
                new FormulaOptions(options.contrastCoding(), null));
            return new Encoded(contrast.coefficientNames(), contrast.design());
        }
        String[] values = table.categorical(name);
        List<String> levels = new ArrayList<>(new java.util.LinkedHashSet<>(Arrays.asList(values)));
        if (levels.contains(null)) throw new IllegalArgumentException("missing categorical random slope: " + name);
        if (levels.size() < 2) throw new IllegalArgumentException("categorical random slope requires at least two levels: " + name);
        double[] design = new double[table.rows() * levels.size()];
        for (int r = 0; r < table.rows(); r++) design[r * levels.size() + levels.indexOf(values[r])] = 1;
        return new Encoded(levels.stream().map(level -> name + level).toList(), design);
    }

    private static Encoded interactRandom(Encoded left, Encoded right, int rows) {
        int columns = left.columns() * right.columns();
        double[] design = new double[rows * columns]; List<String> names = new ArrayList<>();
        // R model.matrix varies the first variable fastest in an interaction.
        for (int j = 0; j < right.columns(); j++) for (int i = 0; i < left.columns(); i++) {
            names.add(left.names().get(i) + ":" + right.names().get(j));
            for (int r = 0; r < rows; r++) design[r * columns + j * left.columns() + i] =
                left.values()[r * left.columns() + i] * right.values()[r * right.columns() + j];
        }
        return new Encoded(List.copyOf(names), design);
    }

    private static void addBlock(List<RandomEffectTerm> scalar, List<CorrelatedRandomEffectBlock> blocks,
            String name, String groupName, List<String> groups, Encoded encoded) {
        if (encoded.columns() == 1) {
            String effect = encoded.names().get(0);
            scalar.add(effect.equals("(Intercept)") ? RandomEffectTerm.randomIntercept("1|" + groupName, groups)
                : RandomEffectTerm.randomSlope("0+" + effect + "|" + groupName, groups, encoded.values()));
        } else {
            double[][] design = new double[groups.size()][encoded.columns()];
            for (int r = 0; r < groups.size(); r++)
                System.arraycopy(encoded.values(), r * encoded.columns(), design[r], 0, encoded.columns());
            blocks.add(CorrelatedRandomEffectBlock.of(name, groups, encoded.names(), design));
        }
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
