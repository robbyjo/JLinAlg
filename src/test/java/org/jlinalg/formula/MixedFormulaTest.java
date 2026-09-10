/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class MixedFormulaTest {
    @Test
    void completeCasesAlignFixedRandomWeightOffsetAndGroupingRows() {
        int n = 24;
        double[] y = new double[n], x = new double[n], z = new double[n];
        double[] weights = new double[n], offset = new double[n];
        double[] unused = new double[n];
        String[] group = new String[n];
        for (int row = 0; row < n; row++) {
            x[row] = row / 10.0;
            z[row] = row % 3 - 1;
            weights[row] = 1 + row % 2;
            offset[row] = .1 * Math.sin(row);
            unused[row] = row;
            group[row] = "g" + row / 4;
            y[row] = 2 + .4 * x[row] + .2 * z[row] + offset[row];
        }
        y[1] = Double.NaN;
        x[2] = Double.NaN;
        z[3] = Double.NaN;
        offset[4] = Double.NaN;
        weights[5] = Double.NaN;
        group[6] = null;
        unused[7] = Double.NaN;
        ModelTable table = ModelTable.builder(n)
            .numeric("y", y).numeric("x", x).numeric("z", z)
            .numeric("w", weights).numeric("o", offset)
            .numeric("unused", unused)
            .categorical("g", group).build();
        FormulaOptions formulaOptions =
            new FormulaOptions(ContrastCoding.TREATMENT, "w");

        assertThrows(IllegalArgumentException.class, () ->
            MixedFormula.compile("y~x+offset(o)+(1+z||g)", table,
                formulaOptions));
        CompiledMixedFormula compiled = MixedFormula.compile(
            "y~x+offset(o)+(1+z||g)", table,
            MixedFormulaOptions.builder()
                .formulaOptions(formulaOptions)
                .missingDataPolicy(MissingDataPolicy.OMIT)
                .build());

        int[] expected = new int[n - 6];
        expected[0] = 0;
        for (int row = 7; row < n; row++) expected[row - 6] = row;
        assertArrayEquals(expected, compiled.retainedRows());
        assertEquals(n, compiled.originalRows());
        assertEquals(6, compiled.omittedRows());
        assertEquals(n - 6, compiled.fixed().rows());
        for (var term : compiled.randomEffects()) {
            assertEquals(n - 6, term.observations());
        }
        var fit = compiled.fitSparse(
            RemlOptions.builder().initialVariances(1, 1, 1).build(),
            BackendPolicy.CPU);
        assertEquals(n - 6, fit.fittedValues().length);
        assertEquals(3, fit.varianceComponents().length);
    }

    @Test
    void mapsFormulaGroupToPedigreePrecisionAndKeepsOrdinaryTerms() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("parent1"),
            PedigreeIndividual.founder("parent2"),
            new PedigreeIndividual("member1", "parent1", "parent2"),
            new PedigreeIndividual("member2", "parent1", "parent2")));
        int n = 24;
        double[] y = new double[n];
        String[] individual = new String[n], batch = new String[n];
        for (int row = 0; row < n; row++) {
            individual[row] = row % 2 == 0 ? "member1" : "member2";
            batch[row] = "b" + row / 4;
            y[row] = 3 + (row % 2 == 0 ? -.8 : .8)
                + .3 * Math.sin(row / 4.0) + .05 * Math.cos(row);
        }
        ModelTable table = ModelTable.builder(n).numeric("y", y)
            .categorical("individual", individual)
            .categorical("batch", batch).build();
        CompiledMixedFormula compiled = MixedFormula.compile(
            "y~1+(1|individual)+(1|batch)", table,
            MixedFormulaOptions.builder()
                .pedigree("individual", pedigree)
                .build());

        assertEquals(1, compiled.pedigreeRandomEffects().size());
        assertEquals(1, compiled.randomEffects().size());
        assertEquals(pedigree.individualIds(), compiled.pedigreeRandomEffects()
            .get(0).randomEffect().coefficientNames());
        var fit = compiled.fitSparse(
            RemlOptions.builder().initialVariances(1, 1, 1)
                .maximumIterations(300).build(), BackendPolicy.CPU);
        assertTrue(fit.converged());
        assertEquals(3, fit.varianceComponents().length);
    }

    @Test
    void pedigreeMappingRejectsSlopeAndFitsJointUnstructuredCovariance() {
        List<PedigreeIndividual> members=new java.util.ArrayList<>();
        for(int i=0;i<8;i++)members.add(PedigreeIndividual.founder("i"+i));
        Pedigree pedigree=Pedigree.of(members);
        int n=80;double[] y=new double[n],x=new double[n];
        String[] id=new String[n],g=new String[n];
        for(int row=0;row<n;row++){
            int individual=row%8,group=row/8;
            id[row]="i"+individual;g[row]="g"+group;x[row]=row%5-2;
            y[row]=2+.3*x[row]+.7*Math.sin(individual)
                +.5*Math.cos(group)*(1+.2*x[row])+.08*Math.sin(2.3*row);
        }
        ModelTable table=ModelTable.builder(n).numeric("y",y).numeric("x",x)
            .categorical("id",id).categorical("g",g).build();
        MixedFormulaOptions options = MixedFormulaOptions.builder()
            .pedigree("id", pedigree).build();
        assertThrows(IllegalArgumentException.class, () ->
            MixedFormula.compile("y~x+(1+x|id)", table, options));
        CompiledMixedFormula joint = MixedFormula.compile(
            "y~x+(1|id)+(1+x|g)", table, options);
        var fit=joint.fitSparse(RemlOptions.builder()
            .initialVariances(1,1,1).maximumIterations(500).build(),
            BackendPolicy.CPU);
        assertTrue(fit.converged());
        assertEquals(3,fit.varianceComponents().length);
        assertTrue(fit.varianceComponents()[1]>0);
    }

    @Test
    void offsetDoesNotConsumeTheFollowingRandomTermAndZeroInterceptSlopesStayCorrelated() {
        ModelTable table = ModelTable.builder(4).numeric("y",1,2,3,4).numeric("x",0,1,2,3)
            .numeric("z",1,3,2,4).numeric("o",.1,.2,.3,.4)
            .categorical("g","a","a","b","b").build();
        var fit = MixedFormula.compile("y~x+offset(o)+(0+x+z|g)",table);
        assertEquals(java.util.List.of("x","z"),fit.correlatedRandomEffects().get(0).effectNames());
        assertEquals(0,fit.randomEffects().size());
        assertEquals(.1,fit.fixed().offset()[0]);
    }

    @Test
    void weightedOffsetsRestoreBothConditionalAndMarginalPredictions() {
        int n = 40; double[] y = new double[n], w = new double[n], offset = new double[n];
        String[] group = new String[n];
        for (int r = 0; r < n; r++) {
            group[r] = "g"+(r/4); w[r] = 1+(r%3)*.5; offset[r] = .7*Math.sin(r);
            y[r] = 3+Math.sin(r/4)*2+.2*Math.cos(r*2.7)+offset[r];
        }
        var table = ModelTable.builder(n).numeric("y",y).numeric("w",w).numeric("o",offset)
            .categorical("g",group).build();
        var formula = MixedFormula.compile("y~offset(o)+(1|g)",table,
            new FormulaOptions(ContrastCoding.TREATMENT,"w"));
        var sparse = formula.fitSparse(RemlOptions.defaults(),BackendPolicy.CPU);
        var fit = formula.fit(RemlOptions.defaults(),BackendPolicy.CPU);
        for (int r = 0; r < n; r++) {
            assertEquals(y[r],sparse.fittedValues()[r]+sparse.residuals()[r],1e-10);
            assertEquals(y[r]-offset[r]-fit.beta()[0],fit.reml().residuals()[r],1e-10);
        }
        var interval = formula.profileFixedEffect(0,.95,-1,7,
            RemlOptions.builder().maximumIterations(200).build(),BackendPolicy.CPU);
        assertTrue(interval.lowerFound() && interval.upperFound());
    }
    @Test
    void compilesSparseRandomInterceptWithoutRuntimeParsing() {
        ModelTable table = ModelTable.builder(9)
            .numeric("y", 0, 1, 2, 4, 5, 6, 8, 9, 10)
            .categorical("subject", "a", "a", "a", "b", "b", "b",
                "c", "c", "c")
            .build();
        CompiledMixedFormula model = MixedFormula.compile(
            "y ~ 1 + (1 | subject)", table);

        assertTrue(model.randomEffects().get(0).sparse());
        assertEquals(3, model.randomEffects().get(0).coefficients());
        assertEquals(5.0, model.fit(
            RemlOptions.builder().initialVariances(10, 2).build(),
            BackendPolicy.CPU).reml().fixedEffects()[0], 1e-10);
    }

    @Test
    void expandsIndependentDoubleBarAndNestedGroupingShorthand() {
        ModelTable table = ModelTable.builder(8)
            .numeric("y", 1, 2, 3, 4, 5, 6, 7, 8)
            .numeric("x", 0, 1, 0, 1, 0, 1, 0, 1)
            .categorical("site", "a", "a", "a", "a", "b", "b", "b", "b")
            .categorical("subject", "u", "u", "v", "v", "u", "u", "v", "v")
            .build();

        CompiledMixedFormula independent = MixedFormula.compile(
            "y ~ x + (1 + x || site)", table);
        assertEquals(2, independent.randomEffects().size());
        assertEquals("1|site", independent.randomEffects().get(0).name());
        assertEquals("0+x|site", independent.randomEffects().get(1).name());

        CompiledMixedFormula nested = MixedFormula.compile(
            "y ~ x + (1 | site/subject)", table);
        assertEquals(2, nested.randomEffects().size());
        assertEquals("1|site:subject", nested.randomEffects().get(1).name());
        assertTrue(nested.fitSparse(
            RemlOptions.builder().initialVariances(1, 1, 1).build(),
            BackendPolicy.CPU).randomCoefficientCount() > 0);
    }

    @Test
    void compilesAndFitsCorrelatedInterceptSlopeBlock() {
        ModelTable table = ModelTable.builder(12)
            .numeric("y", 1.0, 2.0, 2.8, 4.1, 2.1, 3.4,
                4.2, 5.6, 0.4, 1.8, 2.4, 3.7)
            .numeric("x", 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3)
            .categorical("subject", "a", "a", "a", "a",
                "b", "b", "b", "b", "c", "c", "c", "c")
            .build();
        CompiledMixedFormula model = MixedFormula.compile(
            "y ~ x + (1 + x | subject)", table);

        assertEquals(1, model.correlatedRandomEffects().size());
        var fit = model.fitCorrelated(
            RemlOptions.builder().maximumIterations(100).build(),
            BackendPolicy.CPU);
        assertEquals(2, fit.beta().length);
        assertEquals(4, fit.randomEffects().get(0).covariance().length);
        assertTrue(Double.isFinite(fit.randomEffects().get(0)
            .correlation(0, 1)));
    }
}
