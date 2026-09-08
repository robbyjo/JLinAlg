/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

/** Tests R model.matrix coding, covariance parameterization, and original-coordinate predictions. */
public final class MixedCategoricalFormulaTest {
    private static List<String> lines(String name) throws Exception {
        try (var in = MixedCategoricalFormulaTest.class.getResourceAsStream("/r-reference/"+name)) {
            return new String(java.util.Objects.requireNonNull(in,name).readAllBytes(),StandardCharsets.UTF_8).lines().toList();
        }
    }
    private static ModelTable table(String name) throws Exception {
        List<String> lines = lines(name); int n = lines.size()-1;
        String[] header = lines.get(0).split("\t");
        String[][] data = new String[n][];
        for (int r=0;r<n;r++) data[r] = lines.get(r+1).split("\t");
        ModelTable.Builder table = ModelTable.builder(n);
        for (int c=0;c<header.length;c++) {
            if (List.of("f","h","g").contains(header[c])) {
                String[] values = new String[n];
                for (int r=0;r<n;r++) values[r] = data[r][c];
                table.categorical(header[c],values);
            } else {
                double[] values = new double[n];
                for (int r=0;r<n;r++) values[r] = Double.parseDouble(data[r][c]);
                table.numeric(header[c],values);
            }
        }
        return table.build();
    }

    @Test void fortyCategoricalDesignsMatchRIncludingMarginalInteractions() throws Exception {
        ModelTable table = table("mixed-categorical-design-data.tsv");
        Map<String,List<String[]>> cases = new LinkedHashMap<>();
        for (String line : lines("mixed-categorical-designs.tsv").subList(1,lines("mixed-categorical-designs.tsv").size())) {
            String[] fields = line.split("\t");
            cases.computeIfAbsent(fields[0]+"\t"+fields[1],ignored -> new ArrayList<>()).add(fields);
        }
        assertEquals(40,cases.size());
        for (var entry : cases.entrySet()) {
            String[] description = entry.getKey().split("\t");
            var compiled = MixedFormula.compile("y~x+("+description[0]+"|g)",table,
                new FormulaOptions(ContrastCoding.valueOf(description[1]),null));
            var block = compiled.correlatedRandomEffects().get(0);
            assertEquals(entry.getValue().stream().map(f -> f[2]).toList(),block.effectNames(),entry.getKey());
            double[] design = block.effectDesign(); int columns = block.effectCount();
            for (int c=0;c<columns;c++) {
                String[] expected = entry.getValue().get(c)[3].split(",");
                for (int r=0;r<table.rows();r++) assertEquals(Double.parseDouble(expected[r]),design[r*columns+c],1e-12,
                    entry.getKey()+" column="+c+" row="+r);
            }
        }
    }

    @Test void fittedCategoricalCovariancesAndModesMatchLme4() throws Exception {
        ModelTable table = table("mixed-categorical-fit-data.tsv");
        String[] names = {"treatment","indicators","sum","interaction","factorial"};
        String[] expressions = {"1+f","0+f","1+f","0+f:x","1+f*x"};
        for (int index=0;index<names.length;index++) {
            String name = names[index];
            var formula = MixedFormula.compile("y~x+offset(o)+("+expressions[index]+"|g)",table,
                new FormulaOptions(name.equals("sum")?ContrastCoding.SUM:ContrastCoding.TREATMENT,"w"));
            Properties reference = new Properties();
            try (var in = getClass().getResourceAsStream("/r-reference/mixed-categorical-"+name+".properties")) { reference.load(in); }
            long started = System.nanoTime();
            var result = formula.fitSparseUnstructured(RemlOptions.builder().maximumIterations(1000).build(),BackendPolicy.CPU);
            var fit = result.correlatedFit();
            System.out.printf("categorical %s seconds=%.4f ll=%.12f evaluations=%d converged=%s%n",name,
                (System.nanoTime()-started)/1e9,fit.logLikelihood(),result.evaluations(),fit.converged());
            assertTrue(fit.converged(),name);
            assertEquals(Double.parseDouble(reference.getProperty("logLik")),fit.logLikelihood(),2e-6,name);
            near(reference,"beta0",fit.beta()[0],2e-5); near(reference,"beta1",fit.beta()[1],2e-5);
            near(reference,"residual",fit.residualVariance(),2e-5);
            var random = fit.randomEffects().get(0);
            for (int c=0;c<random.covariance().length;c++) near(reference,"cov"+c,random.covariance()[c],1e-3);
            for (int c=0;c<random.effectNames().size();c++) near(reference,"mode"+c,random.modes()[c],5e-4);
            for (int r=0;r<12;r++) near(reference,"fitted"+r,fit.fittedValues()[r],5e-4);
        }
    }

    private static void near(Properties reference,String key,double value,double tolerance) {
        double expected = Double.parseDouble(reference.getProperty(key));
        assertEquals(expected,value,tolerance*Math.max(1,Math.abs(expected)),key);
    }

    @Test void doubleBarsSplitTermsWhileFactorIndicatorColumnsStayCorrelated() throws Exception {
        ModelTable table = table("mixed-categorical-fit-data.tsv");
        var fit = MixedFormula.compile("y~x+(1+f*x||g)",table);
        assertEquals(2,fit.randomEffects().size());
        assertEquals("1|g",fit.randomEffects().get(0).name());
        assertEquals("0+x|g",fit.randomEffects().get(1).name());
        assertEquals(2,fit.correlatedRandomEffects().size());
        assertEquals(List.of("fc","fa","fb"),fit.correlatedRandomEffects().get(0).effectNames());
        assertEquals(List.of("fc:x","fa:x","fb:x"),fit.correlatedRandomEffects().get(1).effectNames());
        var explicit = MixedFormula.compile("y~x+(1|g)+(0+f|g)+(0+x|g)+(0+f:x|g)",table);
        for (int b=0;b<2;b++) assertArrayEquals(explicit.correlatedRandomEffects().get(b).effectDesign(),
            fit.correlatedRandomEffects().get(b).effectDesign());
    }

    @Test void missingCategoricalRandomSlopeIsRejectedBeforeFit() {
        var table = ModelTable.builder(4).numeric("y",1,2,3,4).categorical("g","g","g","h","h")
            .categorical("f","a",null,"b","a").build();
        assertThrows(IllegalArgumentException.class,() -> MixedFormula.compile("y~1+(f|g)",table));
        assertThrows(IllegalArgumentException.class,() -> MixedFormula.compile("y~1+(0+f|g)",table));
    }

    public static void main(String[] args) throws Exception {
        var test = new MixedCategoricalFormulaTest();
        test.fortyCategoricalDesignsMatchRIncludingMarginalInteractions();
        test.doubleBarsSplitTermsWhileFactorIndicatorColumnsStayCorrelated();
        test.missingCategoricalRandomSlopeIsRejectedBeforeFit();
        test.fittedCategoricalCovariancesAndModesMatchLme4();
    }
}
