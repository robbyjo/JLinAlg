/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.ServiceLoader;
import org.jlinalg.pipeline.OmicsTransform;
import org.jlinalg.pipeline.OmicsTransformProvider;
import org.jlinalg.pipeline.OmicsTransforms;

/** Parser for deterministic row-wise omics transform pipelines. */
final class TransformParser {
    private TransformParser() { }

    static OmicsTransform parse(
            List<String> specifications, List<Path> pluginPaths) {
        Map<String, OmicsTransformProvider> providers =
            providers(pluginPaths);
        List<OmicsTransform> transforms = new ArrayList<>();
        for (String specification : specifications) {
            int equals = specification.indexOf('=');
            if (equals < 0 || !specification.substring(0, equals).trim()
                    .equalsIgnoreCase("<omics>"))
                throw new IllegalArgumentException(
                    "transform target must be <omics>: " + specification);
            for (String stage : specification.substring(equals + 1)
                    .split("\\|", -1))
                transforms.add(stage(stage.trim(), providers));
        }
        return transforms.isEmpty() ? OmicsTransforms.identity()
            : OmicsTransforms.compose(
                transforms.toArray(OmicsTransform[]::new));
    }

    private static OmicsTransform stage(
            String text, Map<String, OmicsTransformProvider> providers) {
        String lower = text.toLowerCase(Locale.ROOT);
        if (lower.equals("identity") || lower.equals("identity()"))
            return OmicsTransforms.identity();
        if (lower.equals("log1p") || lower.equals("log1p()"))
            return OmicsTransforms.log1p();
        if (lower.equals("zscore") || lower.equals("zscore()")
                || lower.equals("zscore(ddof=1)"))
            return OmicsTransforms.zScore();
        if (lower.equals("int") || lower.equals("int()")
                || lower.equals("int(method=blom)"))
            return OmicsTransforms.rankInverseNormal();
        if (lower.startsWith("mvalue("))
            return OmicsTransforms.mValue(parameter(lower, "epsilon", 1e-6));
        if (lower.startsWith("winsor(")) {
            double symmetric = parameter(lower, "p", Double.NaN);
            double low = Double.isNaN(symmetric)
                ? parameter(lower, "lower", 0.0) : symmetric;
            double high = Double.isNaN(symmetric)
                ? parameter(lower, "upper", 1.0) : 1.0 - symmetric;
            return OmicsTransforms.winsorize(low, high);
        }
        if (lower.startsWith("log("))
            return OmicsTransforms.shiftedLog(
                parameter(lower, "offset", 0.0));
        if (lower.startsWith("expr("))
            throw new IllegalArgumentException(
                "custom expression transforms are not enabled in this build");
        String name = lower.contains("(")
            ? lower.substring(0, lower.indexOf('(')) : lower;
        OmicsTransformProvider provider = providers.get(name);
        if (provider != null)
            return provider.create(parameters(text));
        throw new IllegalArgumentException("unknown transform stage: " + text);
    }

    private static Map<String, OmicsTransformProvider> providers(
            List<Path> paths) {
        List<URL> urls = new ArrayList<>();
        for (Path path : paths) {
            try {
                urls.add(path.toAbsolutePath().normalize()
                    .toUri().toURL());
            } catch (MalformedURLException exception) {
                throw new IllegalArgumentException(
                    "invalid transform plugin path: " + path, exception);
            }
        }
        ClassLoader loader = urls.isEmpty()
            ? TransformParser.class.getClassLoader()
            : new URLClassLoader(urls.toArray(URL[]::new),
                TransformParser.class.getClassLoader());
        Map<String, OmicsTransformProvider> result = new LinkedHashMap<>();
        for (OmicsTransformProvider provider
                : ServiceLoader.load(OmicsTransformProvider.class, loader)) {
            String name = provider.name().trim().toLowerCase(Locale.ROOT);
            if (name.isEmpty() || result.put(name, provider) != null)
                throw new IllegalArgumentException(
                    "duplicate or blank transform provider name: " + name);
        }
        return Map.copyOf(result);
    }

    private static Map<String, String> parameters(String call) {
        int open = call.indexOf('(');
        if (open < 0) return Map.of();
        int close = call.lastIndexOf(')');
        if (close < open)
            throw new IllegalArgumentException(
                "malformed transform call: " + call);
        Map<String, String> result = new LinkedHashMap<>();
        String body = call.substring(open + 1, close).trim();
        if (body.isEmpty()) return Map.of();
        for (String token : body.split(",", -1)) {
            String[] pair = token.trim().split("=", 2);
            if (pair.length != 2 || pair[0].isBlank()
                    || result.put(pair[0].trim(), pair[1].trim()) != null)
                throw new IllegalArgumentException(
                    "custom transform parameters must be unique name=value pairs");
        }
        return Map.copyOf(result);
    }

    private static double parameter(
            String call, String name, double defaultValue) {
        int open = call.indexOf('(');
        int close = call.lastIndexOf(')');
        if (open < 0 || close < open)
            throw new IllegalArgumentException(
                "malformed transform call: " + call);
        String body = call.substring(open + 1, close);
        for (String token : body.split(",", -1)) {
            String[] pair = token.trim().split("=", 2);
            if (pair.length == 2 && pair[0].trim().equals(name)) {
                try {
                    return Double.parseDouble(pair[1].trim());
                } catch (NumberFormatException exception) {
                    throw new IllegalArgumentException(
                        "invalid transform parameter: " + token, exception);
                }
            }
        }
        return defaultValue;
    }
}
