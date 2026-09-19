/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import org.yaml.snakeyaml.LoaderOptions;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.SafeConstructor;

/** Safe, declarative configuration shared by every CLI command. */
final class ProjectConfiguration {
    private ProjectConfiguration() { }
    static Yaml yaml() {
        LoaderOptions options = new LoaderOptions();
        options.setAllowDuplicateKeys(false);
        options.setMaxAliasesForCollections(0);
        options.setCodePointLimit(8_000_000);
        return new Yaml(new SafeConstructor(options));
    }
    record Resolved(String[] arguments, Map<String,Object> provenance) { }
    static Resolved resolve(String[] args, Path project, Path home) throws IOException {
        List<String> cli = new ArrayList<>();
        Path local = home.resolve(".jlinalg/config.yaml"), config = project.resolve("jlinalg.yaml");
        boolean disabled = false, explicitLocal = false, explicitProject = false;
        for (int i=0; i<args.length; i++) {
            switch (args[i]) {
                case "--no-config" -> disabled = true;
                case "--local-config", "--config" -> {
                    String flag = args[i];
                    if (++i == args.length || args[i].startsWith("--"))
                        throw new IllegalArgumentException(flag + " requires a path");
                    Path path = project.resolve(args[i]).toAbsolutePath().normalize();
                    if (flag.equals("--config")) { config=path; explicitProject=true; }
                    else { local=path; explicitLocal=true; }
                }
                default -> cli.add(args[i]);
            }
        }
        if (disabled && (explicitLocal || explicitProject))
            throw new IllegalArgumentException("--no-config cannot be combined with configuration paths");
        String command = !cli.isEmpty() && !cli.get(0).startsWith("-") ? cli.get(0)
            : cli.contains("--enrichment") ? "enrichment" : "association";
        boolean show = command.equals("config");
        if (show) {
            if (cli.size()!=3 || !cli.get(1).equals("--command"))
                throw new IllegalArgumentException("Usage: config --command COMMAND [--config FILE] [--local-config FILE]");
            command = cli.get(2);
        }
        Map<String,Object> merged = new LinkedHashMap<>();
        List<Map<String,String>> sources = new ArrayList<>();
        if (!disabled) {
            load(local, explicitLocal, merged, sources);
            load(config, explicitProject, merged, sources);
        }
        Map<String,Object> options = new LinkedHashMap<>(mapping(merged.get("defaults"), "defaults"));
        Map<String,Object> commands = mapping(merged.get("commands"), "commands");
        options.putAll(mapping(commands.get(command), "commands."+command));
        Set<String> explicit = new HashSet<>();
        if (!show) for (String value:cli) if(value.startsWith("--")) explicit.add(value.substring(2).split("=",2)[0]);
        List<String> forwarded = new ArrayList<>();
        boolean positional = !show && !cli.isEmpty() && !cli.get(0).startsWith("-");
        if (positional) forwarded.add(cli.get(0));
        for (var entry:options.entrySet()) {
            String key=entry.getKey();
            if (!key.matches("[a-z][a-z0-9-]*") || Set.of("config","local-config","no-config").contains(key))
                throw new IllegalArgumentException("Invalid configured option: " + key);
            if (explicit.contains(key) || entry.getValue()==null) continue;
            List<?> values=entry.getValue() instanceof List<?> list ? list : List.of(entry.getValue());
            for(Object value:values) {
                if (isSwitch(command,key)) {
                    if(!(value instanceof Boolean flag))
                        throw new IllegalArgumentException("Switch " + key + " requires a YAML boolean");
                    if(flag) forwarded.add("--"+key);
                } else {
                    if(value instanceof Map<?,?> || value instanceof List<?>)
                        throw new IllegalArgumentException("Option " + key + " requires scalar values or {path: ...}");
                    forwarded.add("--"+key); forwarded.add(String.valueOf(value));
                }
            }
        }
        if (!show) forwarded.addAll(cli.subList(positional?1:0,cli.size()));
        Map<String,Object> provenance = new LinkedHashMap<>();
        provenance.put("schema_version",1); provenance.put("command",command);
        provenance.put("sources",sources); provenance.put("arguments",forwarded);
        if(show) provenance.put("configured_options",options);
        return new Resolved(show ? new String[]{"config"} : forwarded.toArray(String[]::new),provenance);
    }
    // Keep arity command-specific: e.g. --joint is a switch for mr-estimate
    // but a valued boolean for twas/pwas. Unknown options remain parser errors.
    private static boolean isSwitch(String command,String key) {
        if(Set.of("help","overwrite","no-log","resume","version").contains(key))return true;
        return switch(command) {
            case "association" -> Set.of("conditional-gwas-summary","dry-run","explain").contains(key);
            case "meta-analysis","meta-regression","iv-regression","glm-predict" -> key.equals("no-intercept");
            case "arima-regression" -> Set.of("no-intercept","smooth","parameter-uncertainty").contains(key);
            case "beta-regression","penalized-regression" -> Set.of("no-intercept","no-standardize").contains(key);
            case "confounders","batch-adjust" -> Set.of("write-adjusted","center","no-center","scale","add-mean","nonparametric","mean-only").contains(key);
            case "rare-meta" -> Set.of("leave-variant-out","leave-cohort-out","cohort-results").contains(key);
            case "coloc" -> key.equals("no-trim");
            case "mr-estimate" -> key.equals("joint");
            default -> false;
        };
    }
    private static void load(Path file, boolean required, Map<String,Object> merged,
            List<Map<String,String>> sources) throws IOException {
        if (!Files.exists(file)) {
            if(required) throw new IOException("Configuration is absent: "+file);
            return;
        }
        Map<String,Object> data;
        try(var reader=Files.newBufferedReader(file)) { data=mapping(yaml().load(reader),file.toString()); }
        if(!Integer.valueOf(1).equals(data.get("schema_version")))
            throw new IllegalArgumentException("Configuration requires schema_version: 1 in "+file);
        for(String key:data.keySet()) if(!Set.of("schema_version","defaults","commands").contains(key))
            throw new IllegalArgumentException("Unknown configuration section: "+key);
        resolvePaths(data,file.toAbsolutePath().getParent());
        merge(merged,data);
        sources.add(Map.of("path",file.toAbsolutePath().normalize().toString(),"sha256",FollowupSupport.hash(file)));
    }
    static Map<String,Object> mapping(Object value,String name) {
        if(value==null) return new LinkedHashMap<>();
        if(!(value instanceof Map<?,?> map)) throw new IllegalArgumentException(name+" must be a mapping");
        Map<String,Object> result=new LinkedHashMap<>();
        for(var e:map.entrySet()) {
            if(!(e.getKey() instanceof String key)) throw new IllegalArgumentException(name+" requires string keys");
            result.put(key,e.getValue());
        }
        return result;
    }
    private static Object pathValue(Object value,Path base) {
        if(value instanceof Map<?,?> raw) {
            Map<String,Object> map=mapping(raw,"option");
            if(map.containsKey("path")) {
                if(map.size()!=1 || !(map.get("path") instanceof String path))
                    throw new IllegalArgumentException("Path values must be {path: STRING}");
                return base.resolve(path).normalize().toString();
            }
            resolvePaths(map,base); return map;
        }
        if(value instanceof List<?> list) return list.stream().map(v->pathValue(v,base)).toList();
        return value;
    }
    private static void resolvePaths(Map<String,Object> map,Path base) {
        map.replaceAll((key,value)->pathValue(value,base));
    }
    private static void merge(Map<String,Object> target,Map<String,Object> source) {
        for(var e:source.entrySet()) {
            if(e.getValue() instanceof Map<?,?> && target.get(e.getKey()) instanceof Map<?,?>) {
                Map<String,Object> child=mapping(target.get(e.getKey()),e.getKey());
                merge(child,mapping(e.getValue(),e.getKey())); target.put(e.getKey(),child);
            } else target.put(e.getKey(),e.getValue());
        }
    }
}
