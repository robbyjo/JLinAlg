/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.raremetal;

import java.io.*;
import java.nio.file.Path;
import java.util.*;

/** Model declarations for quantitative rare-variant summaries. Missing historical
 * declarations are unknown, never evidence that a model has been calibrated.
 * Version 1 describes Gaussian scores; it does not enable additional trait models.
 */
public final class ScoreModelMetadata {
    private static final List<String> REQUIRED=List.of("SummaryMetadataVersion", "TraitType", "TraitId",
        "TraitUnits", "TraitTransformation", "NullModel", "GenotypeModel", "EffectScale",
        "ScoreCalibration", "CovarianceModel", "CovarianceScale");
    private static final Set<String> RECOGNIZED;
    static {
        Set<String> keys=new HashSet<>(REQUIRED);
        keys.addAll(List.of("AnalyzedSamples", "GenomeBuild", "InverseNormal", "BinaryTrait"));
        RECOGNIZED=Set.copyOf(keys);
    }
    private final Map<String,String> values;
    private final List<String> columns;
    private ScoreModelMetadata(Map<String,String> values,List<String> columns) {
        this.values=Map.copyOf(values);this.columns=List.copyOf(columns);
    }
    static ScoreModelMetadata read(Path path)throws IOException {
        Map<String,String> values=new HashMap<>();List<String> columns=List.of();
        try(BufferedReader reader=RareMetalStudy.open(path)) {
            for(String line;(line=reader.readLine())!=null;) {
                if(line.isBlank())continue;
                String clean=line.replaceFirst("^#+", "");
                if(clean.startsWith("CHROM")) {columns=List.of(clean.trim().split("\\s+"));break;}
                if(!line.startsWith("#"))break;
                int equals=line.indexOf('=');
                if(line.startsWith("##")&&equals>2) {
                    String key=line.substring(2,equals),value=line.substring(equals+1).trim();
                    if(RECOGNIZED.contains(key)) {
                        if(value.isBlank()||values.putIfAbsent(key,value)!=null)
                            throw new IOException("blank or duplicate metadata "+key+" in "+path);
                    }
                }
            }
        }
        ScoreModelMetadata result=new ScoreModelMetadata(values,columns);
        result.validateQuantitative(false);
        return result;
    }
    /** Declared value, or null when the source did not declare it. */
    public String value(String key){return values.get(key);}
    /** Immutable recognized declarations, excluding unknown producer headers. */
    public Map<String,String> declarations(){return values;}
    List<String> columns(){return columns;}
    /** Whether every version-1 model field is explicitly declared. */
    public boolean complete(){return values.keySet().containsAll(REQUIRED);}
    /** Audit label; legacy assumptions must be checked against cohort provenance. */
    public String status(){return complete()?"declared-quantitative":"legacy-assumed-quantitative";}
    /** Reject unsupported declarations even when legacy metadata is permitted. */
    public void validateQuantitative(boolean requireDeclared)throws IOException {
        allow("SummaryMetadataVersion", "1");
        allow("TraitType", "quantitative");
        allow("NullModel", "unrelated-Gaussian", "related-Gaussian-REML");
        allow("GenotypeModel", "diploid-additive");
        allow("EffectScale", "trait-per-alt-allele");
        allow("ScoreCalibration", "asymptotic-normal");
        allow("CovarianceModel", "model-based");
        allow("CovarianceScale", "per-sample", "score");
        allow("InverseNormal", "ON", "OFF");
        allow("BinaryTrait", "False", "false", "0");
        String transform=value("TraitTransformation"),inverse=value("InverseNormal");
        if(transform!=null&&inverse!=null&&!transform.equals(inverse.equals("ON")?"inverse-normal":"none"))
            throw new IOException("TraitTransformation conflicts with InverseNormal");
        if(requireDeclared&&!complete()) {
            List<String> missing=REQUIRED.stream().filter(k->!values.containsKey(k)).toList();
            throw new IOException("strict model metadata requires "+String.join(",",missing));
        }
    }
    private void allow(String key,String... allowed)throws IOException {
        String value=value(key);
        if(value!=null&&!Arrays.asList(allowed).contains(value))
            throw new IOException("unsupported "+key+"="+value+"; rare-meta supports quantitative Gaussian scores with asymptotic-normal calibration only");
    }
    private String transformation(){
        String result=value("TraitTransformation");
        return result!=null?result:value("InverseNormal")==null?null:value("InverseNormal").equals("ON")?"inverse-normal":"none";
    }
    /** Check effect comparability across independent cohorts. Gaussian OLS and
     * Gaussian REML may coexist; their within-cohort covariance models differ. */
    public void requireCompatible(ScoreModelMetadata other)throws IOException {
        for(String key:List.of("TraitType","TraitId","TraitUnits","GenotypeModel","EffectScale",
                "ScoreCalibration","CovarianceModel"))match(key,value(key),other.value(key));
        match("TraitTransformation",transformation(),other.transformation());
    }
    /** Score and covariance files must describe the same cohort and null fit.
     * Missing old headers remain unknown; conflicting declarations always fail. */
    void requireSameFilePair(ScoreModelMetadata other)throws IOException {
        for(String key:RECOGNIZED)match(key,value(key),other.value(key));
        match("TraitTransformation",transformation(),other.transformation());
    }
    ScoreModelMetadata withMissingFrom(ScoreModelMetadata other)throws IOException {
        requireSameFilePair(other);
        Map<String,String> merged=new HashMap<>(other.values);merged.putAll(values);
        return new ScoreModelMetadata(merged,columns);
    }
    private static void match(String key,String left,String right)throws IOException {
        if(left!=null&&right!=null&&!left.equals(right))
            throw new IOException("incompatible metadata "+key+": "+left+" / "+right);
    }
    static String gaussianHeaders(String traitId,String units,String nullModel)throws IOException {
        Map<String,String> fields=new LinkedHashMap<>();
        fields.put("SummaryMetadataVersion","1");fields.put("TraitType","quantitative");
        if(traitId!=null)fields.put("TraitId",traitId);
        fields.put("TraitUnits",units);fields.put("TraitTransformation","none");
        fields.put("NullModel",nullModel);fields.put("GenotypeModel","diploid-additive");
        fields.put("EffectScale","trait-per-alt-allele");fields.put("ScoreCalibration","asymptotic-normal");
        fields.put("CovarianceModel","model-based");fields.put("CovarianceScale","per-sample");
        StringBuilder result=new StringBuilder();
        for(var field:fields.entrySet()) {
            String value=field.getValue();
            if(value==null||value.isBlank()||!value.equals(value.trim())||value.chars().anyMatch(Character::isISOControl))
                throw new IOException("metadata "+field.getKey()+" must be nonblank text without control characters or outer whitespace");
            result.append("##").append(field.getKey()).append('=').append(value).append('\n');
        }
        return result.toString();
    }
}
