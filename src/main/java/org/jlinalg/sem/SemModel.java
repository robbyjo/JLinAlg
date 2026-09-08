/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/** RAM path model. Data columns follow {@link #variables()}, excluding latent variables. */
public final class SemModel {
    enum Kind { REGRESSION, VARIANCE, COVARIANCE, INTERCEPT }
    record Element(Kind kind, int first, int second, String label,
                   double start, boolean fixed) { }

    private final List<String> variables;
    private final List<String> allVariables;
    private final boolean meanStructure;
    private final List<Element> elements;
    private final List<String> freeLabels;
    private final Map<String, Integer> freeIndex;
    private final boolean[] varianceParameter;

    private SemModel(Builder builder) {
        this.allVariables = List.copyOf(builder.variables);
        this.variables = allVariables.subList(0, builder.observedCount);
        this.meanStructure = builder.meanStructure;
        List<Element> values = new ArrayList<>(builder.elements);
        for (int variable = 0; variable < allVariables.size(); variable++) {
            final int target = variable;
            boolean present = values.stream().anyMatch(value ->
                value.kind() == Kind.VARIANCE && value.first() == target);
            if (!present) values.add(new Element(Kind.VARIANCE, variable, variable,
                allVariables.get(variable) + "~~" + allVariables.get(variable), 1.0, false));
            if (meanStructure && variable < variables.size() && values.stream().noneMatch(value ->
                    value.kind() == Kind.INTERCEPT && value.first() == target))
                values.add(new Element(Kind.INTERCEPT, variable, variable,
                    allVariables.get(variable) + "~1", 0.0, false));
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
    public List<String> latentVariables() { return allVariables.subList(variables.size(), allVariables.size()); }
    public boolean hasMeanStructure() { return meanStructure; }
    List<String> allVariables() { return allVariables; }
    /** Copy this specification, retaining constraints and starting values. */
    public Builder toBuilder() { return new Builder(this); }
    public List<String> freeParameterLabels() { return freeLabels; }
    public int freeParameterCount() { return freeLabels.size(); }
    List<Element> elements() { return elements; }
    int freeIndex(String label) { return freeIndex.get(label); }
    boolean varianceParameter(int index) { return varianceParameter[index]; }
    SemModel freeElement(Kind kind, String first, String second, String label) {
        Builder b=toBuilder();
        Integer i=b.index.get(first),j=b.index.get(second);
        if(i==null||j==null)throw new IllegalArgumentException("unknown modification variable");
        double start=0;
        for(Element e:elements) if(e.kind()==kind && ((e.first()==i && e.second()==j)
                || (kind==Kind.COVARIANCE && e.first()==j && e.second()==i))) {
            if(!e.fixed())throw new IllegalArgumentException("modification must release a fixed or omitted parameter");
            start=e.start();b.elements.remove(e);
        }
        b.positions.clear();
        List<Element> keep=new ArrayList<>(b.elements);b.elements.clear();
        for(Element e:keep)b.add(e.kind(),allVariables.get(e.first()),allVariables.get(e.second()),e.label(),e.start(),e.fixed());
        b.add(kind,first,second,label,start,false);
        return b.build();
    }

    /** Fluent builder; repeated labels impose equality constraints. */
    public static final class Builder {
        private final List<String> variables;
        private final Map<String, Integer> index;
        private final List<Element> elements = new ArrayList<>();
        private final HashSet<String> positions = new HashSet<>();
        private final int observedCount;
        private boolean meanStructure;

        private Builder(String[] variableNames) {
            if (variableNames == null || variableNames.length < 2)
                throw new IllegalArgumentException("SEM requires at least two variables");
            this.variables = new ArrayList<>(List.of(variableNames.clone()));
            this.observedCount = variableNames.length;
            this.index = new HashMap<>();
            for (int position = 0; position < variables.size(); position++) {
                String value = variables.get(position);
                if (value == null || value.isBlank() || index.put(value, position) != null)
                    throw new IllegalArgumentException("SEM variable names must be unique and nonblank");
            }
        }

        private Builder(SemModel model) {
            this(model.variables.toArray(String[]::new));
            for (String latent : model.latentVariables()) latent(latent);
            this.meanStructure = model.meanStructure;
            for (Element e : model.elements) add(e.kind(), model.allVariables.get(e.first()),
                model.allVariables.get(e.second()), e.label(), e.start(), e.fixed());
        }

        public Builder latent(String... names) {
            for (String name : names) {
                if (name == null || name.isBlank() || index.containsKey(name))
                    throw new IllegalArgumentException("latent names must be unique and nonblank");
                index.put(name, variables.size()); variables.add(name);
            }
            return this;
        }
        public Builder loading(String indicator, String latent, double start) {
            return regression(indicator, latent, start);
        }
        public Builder loading(String label, String indicator, String latent, double start) {
            return regression(label, indicator, latent, start);
        }
        public Builder fixedLoading(String indicator, String latent, double value) {
            return fixedRegression(indicator, latent, value);
        }
        /** Adds free observed intercepts where unspecified; latent intercepts default to zero. */
        public Builder meanStructure() { meanStructure = true; return this; }
        public Builder intercept(String variable, double start) {
            return intercept(variable + "~1", variable, start);
        }
        public Builder intercept(String label, String variable, double start) {
            meanStructure = true;
            return add(Kind.INTERCEPT, variable, variable, label, start, false);
        }
        public Builder fixedIntercept(String variable, double value) {
            meanStructure = true;
            return add(Kind.INTERCEPT, variable, variable, "fixed", value, true);
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
            if (!(value >= 0.0)) throw new IllegalArgumentException("fixed variance must be nonnegative");
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
            if(kind==Kind.COVARIANCE && firstIndex.equals(secondIndex))
                throw new IllegalArgumentException("use variance for a diagonal covariance element");
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
