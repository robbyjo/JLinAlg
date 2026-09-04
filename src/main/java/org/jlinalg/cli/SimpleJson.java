/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Minimal JSON reader for dependency-free command-line service clients. */
final class SimpleJson {
    private final String text;
    private int offset;

    private SimpleJson(String text) {
        this.text = text;
    }

    static Object parse(String text) throws IOException {
        SimpleJson reader = new SimpleJson(text);
        Object result = reader.value();
        reader.space();
        if (reader.offset != text.length()) reader.fail("trailing content");
        return result;
    }

    private Object value() throws IOException {
        space();
        if (offset >= text.length()) return fail("expected a value");
        return switch (text.charAt(offset)) {
            case '{' -> object();
            case '[' -> array();
            case '"' -> string();
            case 't' -> literal("true", Boolean.TRUE);
            case 'f' -> literal("false", Boolean.FALSE);
            case 'n' -> literal("null", null);
            default -> number();
        };
    }

    private Map<String, Object> object() throws IOException {
        offset++;
        Map<String, Object> result = new LinkedHashMap<>();
        space();
        if (take('}')) return result;
        while (true) {
            space();
            if (offset >= text.length() || text.charAt(offset) != '"')
                return fail("expected an object key");
            String key = string();
            space();
            if (!take(':')) return fail("expected ':'");
            result.put(key, value());
            space();
            if (take('}')) return result;
            if (!take(',')) return fail("expected ',' or '}'");
        }
    }

    private List<Object> array() throws IOException {
        offset++;
        List<Object> result = new ArrayList<>();
        space();
        if (take(']')) return result;
        while (true) {
            result.add(value());
            space();
            if (take(']')) return result;
            if (!take(',')) return fail("expected ',' or ']'");
        }
    }

    private String string() throws IOException {
        offset++;
        StringBuilder result = new StringBuilder();
        while (offset < text.length()) {
            char value = text.charAt(offset++);
            if (value == '"') return result.toString();
            if (value != '\\') {
                if (value < 0x20) return fail("control character in string");
                result.append(value);
                continue;
            }
            if (offset >= text.length()) return fail("incomplete escape");
            char escaped = text.charAt(offset++);
            switch (escaped) {
                case '"', '\\', '/' -> result.append(escaped);
                case 'b' -> result.append('\b');
                case 'f' -> result.append('\f');
                case 'n' -> result.append('\n');
                case 'r' -> result.append('\r');
                case 't' -> result.append('\t');
                case 'u' -> result.append(unicode());
                default -> { return fail("invalid escape"); }
            }
        }
        return fail("unterminated string");
    }

    private char unicode() throws IOException {
        if (offset + 4 > text.length()) return fail("incomplete unicode escape");
        try {
            char result = (char) Integer.parseInt(
                text.substring(offset, offset + 4), 16);
            offset += 4;
            return result;
        } catch (NumberFormatException exception) {
            return fail("invalid unicode escape");
        }
    }

    private Number number() throws IOException {
        int start = offset;
        if (take('-')) { /* sign */ }
        while (offset < text.length()
                && Character.isDigit(text.charAt(offset))) offset++;
        if (take('.')) while (offset < text.length()
                && Character.isDigit(text.charAt(offset))) offset++;
        if (offset < text.length()
                && (text.charAt(offset) == 'e' || text.charAt(offset) == 'E')) {
            offset++;
            if (offset < text.length()
                    && (text.charAt(offset) == '+'
                        || text.charAt(offset) == '-')) offset++;
            while (offset < text.length()
                    && Character.isDigit(text.charAt(offset))) offset++;
        }
        String token = text.substring(start, offset);
        try {
            if (!token.contains(".") && !token.contains("e")
                    && !token.contains("E")) return Long.valueOf(token);
            return Double.valueOf(token);
        } catch (NumberFormatException exception) {
            return fail("invalid number");
        }
    }

    private Object literal(String token, Object result) throws IOException {
        if (!text.startsWith(token, offset)) return fail("invalid literal");
        offset += token.length();
        return result;
    }

    private boolean take(char value) {
        if (offset < text.length() && text.charAt(offset) == value) {
            offset++;
            return true;
        }
        return false;
    }

    private void space() {
        while (offset < text.length()
                && Character.isWhitespace(text.charAt(offset))) offset++;
    }

    private <T> T fail(String message) throws IOException {
        throw new IOException("invalid JSON at offset " + offset + ": " + message);
    }
}
