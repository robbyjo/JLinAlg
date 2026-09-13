/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Small expression language: no scripting engine, reflection, or executable input. */
final class EnrichmentExpression {
    private interface Node { Object eval(Map<String, String> row, int count); }
    private final List<String> tokens = new ArrayList<>();
    private final Set<String> columns;
    private int position;
    private final Node root;

    EnrichmentExpression(String expression, Set<String> columns) {
        this.columns = columns;
        for (int i = 0; i < expression.length();) {
            char c = expression.charAt(i);
            if (Character.isWhitespace(c)) { i++; continue; }
            int start = i++;
            if (c == '\'' || c == '"' || c == '`') {
                while (i < expression.length() && expression.charAt(i) != c) i++;
                if (i == expression.length()) throw new IllegalArgumentException("unterminated expression quote");
                i++;
            } else if (Character.isDigit(c) || c == '.') {
                while (i < expression.length() && (Character.isDigit(expression.charAt(i)) || expression.charAt(i) == '.')) i++;
                if (i < expression.length() && (expression.charAt(i) == 'e' || expression.charAt(i) == 'E')) {
                    i++;
                    if (i < expression.length() && (expression.charAt(i) == '+' || expression.charAt(i) == '-')) i++;
                    while (i < expression.length() && Character.isDigit(expression.charAt(i))) i++;
                }
            } else if (Character.isLetter(c) || c == '_') {
                while (i < expression.length() && (Character.isLetterOrDigit(expression.charAt(i)) || "_.".indexOf(expression.charAt(i)) >= 0)) i++;
            } else if (i < expression.length() && Set.of("&&", "||", "<=", ">=", "==", "!=").contains(expression.substring(start, i + 1))) i++;
            tokens.add(expression.substring(start, i));
        }
        root = or();
        if (position != tokens.size()) throw new IllegalArgumentException("unexpected expression token: " + tokens.get(position));
    }
    Boolean test(Map<String, String> row, int count) { return bool(root.eval(row, count)); }
    private Node or() {
        Node node = and();
        while (take("||")) { Node left = node, right = and(); node = (r, n) -> {
            Boolean a = bool(left.eval(r,n));
            if (Boolean.TRUE.equals(a)) return Boolean.TRUE;
            Boolean b = bool(right.eval(r,n));
            if (Boolean.TRUE.equals(a) || Boolean.TRUE.equals(b)) return Boolean.TRUE;
            if (a == null || b == null) return null;
            return Boolean.FALSE;
        }; }
        return node;
    }
    private Node and() {
        Node node = compare();
        while (take("&&")) { Node left = node, right = compare(); node = (r, n) -> {
            Boolean a = bool(left.eval(r,n));
            if (Boolean.FALSE.equals(a)) return Boolean.FALSE;
            Boolean b = bool(right.eval(r,n));
            if (Boolean.FALSE.equals(a) || Boolean.FALSE.equals(b)) return Boolean.FALSE;
            if (a == null || b == null) return null;
            return Boolean.TRUE;
        }; }
        return node;
    }
    private Node compare() {
        Node left = sum();
        if (position == tokens.size() || !Set.of("<", "<=", ">", ">=", "==", "!=").contains(tokens.get(position))) return left;
        String op = tokens.get(position++); Node right = sum();
        return (r, n) -> {
            Object a = left.eval(r,n), b = right.eval(r,n);
            if (a == null || b == null) return null;
            int c;
            if (op.equals("==") || op.equals("!=")) c = a instanceof Number || b instanceof Number
                ? Double.compare(number(a), number(b)) : a.toString().compareTo(b.toString());
            else c = Double.compare(number(a), number(b));
            return switch(op) { case "<" -> c < 0; case "<=" -> c <= 0; case ">" -> c > 0;
                case ">=" -> c >= 0; case "==" -> c == 0; default -> c != 0; };
        };
    }
    private Node sum() {
        Node node = product();
        while (position < tokens.size() && Set.of("+", "-").contains(tokens.get(position))) {
            String op = tokens.get(position++); Node left = node, right = product();
            node = arithmetic(left, right, op);
        }
        return node;
    }
    private Node product() {
        Node node = atom();
        while (position < tokens.size() && Set.of("*", "/").contains(tokens.get(position))) {
            String op = tokens.get(position++); Node left = node, right = atom();
            node = arithmetic(left, right, op);
        }
        return node;
    }
    private static Node arithmetic(Node left, Node right, String op) {
        return (r,n) -> {
            Object a = left.eval(r,n), b = right.eval(r,n);
            if (a == null || b == null) return null;
            double x = number(a), y = number(b);
            return finite(switch(op) { case "+" -> x+y; case "-" -> x-y; case "*" -> x*y; default -> x/y; });
        };
    }
    private Node atom() {
        if (take("!")) { Node child = atom(); return (r,n) -> { Boolean v = bool(child.eval(r,n)); return v == null ? null : !v; }; }
        if (take("-")) { Node child = atom(); return (r,n) -> { Object v=child.eval(r,n); return v == null ? null : -number(v); }; }
        if (take("(")) { Node child = or(); expect(")"); return child; }
        if (position == tokens.size()) throw new IllegalArgumentException("incomplete selection expression");
        String token = tokens.get(position++);
        if (token.startsWith("'") || token.startsWith("\"")) return (r,n) -> token.substring(1, token.length()-1);
        if (token.equals("true") || token.equals("false")) return (r,n) -> Boolean.valueOf(token);
        if (token.equals("bonferroni")) {
            Node alpha = (r,n) -> 0.05, count = (r,n) -> n;
            if (take("(")) { alpha = sum(); if (take(",")) count = sum(); expect(")"); }
            Node a = alpha, m = count;
            return (r,n) -> {
                double av = number(a.eval(r,n)), mv = number(m.eval(r,n));
                if (!(av > 0 && av <= 1) || mv < 1 || mv != Math.rint(mv))
                    throw new IllegalArgumentException("bonferroni requires alpha in (0,1] and a positive integer test count");
                return av / mv;
            };
        }
        if (Set.of("abs", "is_missing", "is_finite").contains(token) && take("(")) {
            Node child = sum(); expect(")");
            return (r,n) -> {
                Object v = child.eval(r,n);
                if (token.equals("is_missing")) return v == null;
                if (token.equals("is_finite")) { if(v==null) return false; try { number(v); return true; } catch(IllegalArgumentException e) { return false; } }
                return v == null ? null : Math.abs(number(v));
            };
        }
        if (Character.isDigit(token.charAt(0)) || token.charAt(0) == '.') {
            double value = number(token); return (r,n) -> value;
        }
        String column = token.startsWith("`") ? token.substring(1,token.length()-1) : token;
        if (!columns.contains(column)) throw new IllegalArgumentException("expression column is absent: " + column);
        return (r,n) -> { String value = r.get(column); return missing(value) ? null : value; };
    }
    static boolean missing(String value) { return value == null || Set.of("", "NA", "NAN", "NULL", ".").contains(value.trim().toUpperCase(java.util.Locale.ROOT)); }
    private boolean take(String token) { if(position < tokens.size() && tokens.get(position).equals(token)) { position++; return true; } return false; }
    private void expect(String token) { if(!take(token)) throw new IllegalArgumentException("expected expression token: " + token); }
    private static Boolean bool(Object value) {
        if(value == null) return null;
        if(value instanceof Boolean b) return b;
        throw new IllegalArgumentException("selection expression must return a boolean");
    }
    private static double number(Object value) {
        if (value == null) throw new IllegalArgumentException("missing numeric expression value");
        try { return finite(value instanceof Number n ? n.doubleValue() : Double.parseDouble(value.toString())); }
        catch (NumberFormatException e) { throw new IllegalArgumentException("nonnumeric expression value: " + value, e); }
    }
    private static double finite(double value) { if(!Double.isFinite(value)) throw new IllegalArgumentException("nonfinite numeric expression value"); return value; }
}
