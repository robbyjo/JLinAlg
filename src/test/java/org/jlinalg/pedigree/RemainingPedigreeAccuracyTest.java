/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.*;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;
class RemainingPedigreeAccuracyTest {
    static List<PedigreeIndividual> entries(){return List.of(PedigreeIndividual.founder("a"),PedigreeIndividual.founder("b"),
        new PedigreeIndividual("c","a","b"),new PedigreeIndividual("d","a","b"),new PedigreeIndividual("e","c","d"),new PedigreeIndividual("f","e",null));}
    @Test void sparseConstructionRejectsCyclesWithoutDenseRelationshipMatrix() {
        var cycle=List.of(new PedigreeIndividual("a","b",null),new PedigreeIndividual("b","a",null));
        assertThrows(IllegalArgumentException.class,()->PedigreeRandomEffectTerm.ofUninbred("animal",List.of("a"),cycle));
        assertThrows(IllegalArgumentException.class,()->PedigreeRandomEffectTerm.ofSparse("animal",List.of("a"),cycle,new double[2]));
    }
    @Test void inbreedingOneParentAndArbitraryOrderingMatchRInverse() throws Exception {
        var p=new Properties();try(var reader=Files.newBufferedReader(Path.of("src/test/resources/r-reference/remaining-mixed/reference.properties"))){p.load(reader);}
        double[] expectedA=numbers(p,"pedigree.A"),expectedQ=numbers(p,"pedigree.precision");
        List<PedigreeIndividual> source=entries();int[] order={4,2,5,1,0,3};List<PedigreeIndividual> shuffled=new ArrayList<>();
        for(int i:order)shuffled.add(source.get(i));var pedigree=Pedigree.of(shuffled);
        var term=PedigreeRandomEffectTerm.ofSparse("animal",List.of("f","e","f"),shuffled,pedigree.inbreedingCoefficients());
        double[] q=new double[36];var precision=term.precision();int[] starts=precision.rowStarts(),columns=precision.columnIndices();double[] values=precision.values();
        for(int i=0;i<6;i++)for(int j=starts[i];j<starts[i+1];j++)q[6*i+columns[j]]=values[j];
        double[] a=pedigree.relationshipMatrix();
        for(int i=0;i<6;i++)for(int j=0;j<6;j++) {
            assertEquals(expectedA[6*order[i]+order[j]],a[6*i+j],1e-14);
            assertEquals(expectedQ[6*order[i]+order[j]],q[6*i+j],2e-14);
            double product=0;for(int k=0;k<6;k++)product+=a[6*i+k]*q[6*k+j];assertEquals(i==j?1:0,product,2e-14);
        }
        assertEquals(6,term.randomEffect().coefficients());assertArrayEquals(new int[]{2,0,2},term.randomEffect().columnIndices());
    }
    @Test void unobservedAncestorsAndRepeatedRecordsRetainDenseSparseLikelihoodEquivalence() {
        var pedigree=Pedigree.of(entries());int n=60;double[] y=new double[n],x=new double[n],basis=new double[n*n];List<String> ids=new ArrayList<>();
        String[] observed={"c","d","e","f"};
        for(int i=0;i<n;i++){ids.add(observed[i%4]);x[i]=1;y[i]=new double[]{0,2,7,12}[i%4]+(i%3);}
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)basis[i*n+j]=pedigree.relationship(ids.get(i),ids.get(j));
        var term=PedigreeRandomEffectTerm.of("animal",ids,pedigree);var options=new GlmmLaplaceOptions(100,100,1e-8,1,1e-8,100,null);
        var sparse=SparseGlmmLaplace.fitWithPrecision(y,x,n,1,GlmFamilies.poisson(),List.of(term.randomEffect()),List.of(term.precision()),null,null,options,BackendPolicy.CPU);
        var dense=GlmmLaplace.fit(y,x,n,1,GlmFamilies.poisson(),List.of(new VarianceComponent("animal",n,basis)),null,null,options,BackendPolicy.CPU);
        assertTrue(sparse.converged());assertTrue(dense.converged());assertEquals(6,sparse.componentCoefficients("animal").length);
        assertEquals(dense.marginalLogLikelihood(),sparse.marginalLogLikelihood(),2e-7);
        assertArrayEquals(dense.beta(),sparse.beta(),2e-6);assertArrayEquals(dense.varianceComponents(),sparse.varianceComponents(),2e-6);
        assertArrayEquals(dense.componentPredictor("animal"),sparse.componentPredictor("animal"),2e-6);
    }
    private static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
}
