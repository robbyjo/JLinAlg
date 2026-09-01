/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/** RAM-style observed-variable path model. */
public final class SemModel {
    enum Kind { REGRESSION, VARIANCE, COVARIANCE }
    record Element(Kind kind, int first, int second, String label,
                   double start, boolean fixed) { }

    private final List<String> variables;
    private final List<Element> elements;
    private final List<String> freeLabels;
    private final Map<String, Integer> freeIndex;
    private final boolean[] varianceParameter;

    private SemModel(Builder builder) {
        this.variables = List.copyOf(builder.variables);
        List<Element> values = new ArrayList<>(builder.elements);
        for (int variable = 0; variable < variables.size(); variable++) {
            final int target = variable;
            boolean present = values.stream().anyMatch(value ->
                value.kind() == Kind.VARIANCE && value.first() == target);
            if (!present) values.add(new Element(Kind.VARIANCE, variable, variable,
                variables.get(variable) + "~~" + variables.get(variable), 1.0, false));
        }
        this.elements = List.copyOf(values);
        List<String> labels = new ArrayList<>();
        Map<String, Boolean> type = new HashMap<>();
        for (Element value : elements) {
            if (value.fixed()) continue;
            boolean variance = value.kind() == Kind.VARIANCE;
            Boolean previous = type.putIfAbsent(value.label(), variance);
            if (previous != null && previous != variance) {
                throw new IllegalArgumentException(
                    "an equality label cannot mix variance and nonvariance parameters");
            }
            if (!labels.contains(value.label())) labels.add(value.label());
        }
        this.freeLabels = List.copyOf(labels);
        Map<String, Integer> indices = new HashMap<>();
        for (int index = 0; index < labels.size(); index++) indices.put(labels.get(index), index);
        this.freeIndex = Map.copyOf(indices);
        this.varianceParameter = new boolean[labels.size()];
        for (Element value : elements) {
            if (!value.fixed() && value.kind() == Kind.VARIANCE)
                varianceParameter[freeIndex.get(value.label())] = true;
        }
    }

    public static Builder builder(String... variables) { return new Builder(variables); }
    public List<String> variables() { return variables; }
    public List<String> freeParameterLabels() { return freeLabels; }
    public int freeParameterCount() { return freeLabels.size(); }
    List<Element> elements() { return elements; }
    int freeIndex(String label) { return freeIndex.get(label); }
    boolean varianceParameter(int index) { return varianceParameter[index]; }

    /** Fluent builder; repeated labels impose equality constraints. */
    public static final class Builder {
        private final List<String> variables;
        private final Map<String, Integer> index;
        private final List<Element> elements = new ArrayList<>();
        private final HashSet<String> positions = new HashSet<>();

        private Builder(String[] variableNames) {
            if (variableNames == null || variableNames.length < 2)
                throw new IllegalArgumentException("SEM requires at least two variables");
            this.variables = List.of(variableNames.clone());
            this.index = new HashMap<>();
            for (int position = 0; position < variables.size(); position++) {
                String value = variables.get(position);
                if (value == null || value.isBlank() || index.put(value, position) != null)
                    throw new IllegalArgumentException("SEM variable names must be unique and nonblank");
            }
        }

        public Builder regression(String outcome, String predictor, double start) {
            return regression(outcome + "~" + predictor, outcome, predictor, start);
        }
        public Builder regression(
                String label, String outcome, String predictor, double start) {
            return add(Kind.REGRESSION, outcome, predictor, label, start, false);
        }
        public Builder fixedRegression(String outcome, String predictor, double value) {
            return add(Kind.REGRESSION, outcome, predictor, "fixed", value, true);
        }
        public Builder variance(String variable, double start) {
            return variance(variable + "~~" + variable, variable, start);
        }
        public Builder variance(String label, String variable, double start) {
            if (!(start > 0.0)) throw new IllegalArgumentException("variance start must be positive");
            return add(Kind.VARIANCE, variable, variable, label, start, false);
        }
        public Builder fixedVariance(String variable, double value) {
            if (!(value > 0.0)) throw new IllegalArgumentException("variance must be positive");
            return add(Kind.VARIANCE, variable, variable, "fixed", value, true);
        }
        public Builder covariance(String first, String second, double start) {
            return covariance(first + "~~" + second, first, second, start);
        }
        public Builder covariance(
                String label, String first, String second, double start) {
            return add(Kind.COVARIANCE, first, second, label, start, false);
        }
        public Builder fixedCovariance(String first, String second, double value) {
            return add(Kind.COVARIANCE, first, second, "fixed", value, true);
        }

        private Builder add(Kind kind, String first, String second,
                String label, double start, boolean fixed) {
            Integer firstIndex = index.get(first);
            Integer secondIndex = index.get(second);
            if (firstIndex == null || secondIndex == null || label == null || label.isBlank()
                    || !Double.isFinite(start))
                throw new IllegalArgumentException("invalid SEM element");
            String position = kind + ":" + firstIndex + ":" + secondIndex;
            if (kind == Kind.COVARIANCE && firstIndex > secondIndex)
                position = kind + ":" + secondIndex + ":" + firstIndex;
            if (!positions.add(position)) throw new IllegalArgumentException("duplicate SEM element");
            elements.add(new Element(kind, firstIndex, secondIndex, label, start, fixed));
            return this;
        }

        public SemModel build() { return new SemModel(this); }
    }
}
