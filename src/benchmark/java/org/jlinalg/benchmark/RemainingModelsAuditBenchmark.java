/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.*;
import java.util.*;
import java.util.function.Supplier;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.survival.*;
import org.jlinalg.mediation.*;
import org.jlinalg.nonlinear.*;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.*;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.reml.*;

/** Accuracy-gated, checksum-consuming timings for the remaining-model audit. */
public final class RemainingModelsAuditBenchmark {
 private RemainingModelsAuditBenchmark() { }
 private static final Path ROOT=Path.of("src/test/resources/remaining-model-audit");
 private static double[][] data(String name)throws Exception {
  var lines=Files.readAllLines(ROOT.resolve(name+".csv"));double[][] out=new double[lines.size()][];
  for(int i=0;i<out.length;i++)out[i]=Arrays.stream(lines.get(i).split(",")).mapToDouble(Double::parseDouble).toArray();return out;
 }
 private static double[] ref(String model,String quantity)throws Exception {
  List<Double> values=new ArrayList<>();
  for(String line:Files.readAllLines(ROOT.resolve("reference.csv"))) {String[] a=line.split(",");if(a[0].equals(model)&&a[1].equals(quantity))values.add(Double.parseDouble(a[3]));}
  return values.stream().mapToDouble(Double::doubleValue).toArray();
 }
 private static void time(String name,int repeats,Supplier<double[]> run,double[] reference,double tolerance) {
  double[] check=run.get();double error=0;
  for(int j=0;j<check.length;j++)error=Math.max(error,Math.abs(check[j]-reference[j]));
  if(!Double.isFinite(error)||error>tolerance)throw new AssertionError(name+" reference error="+error);
  for(int i=0;i<10;i++)run.get();
  double checksum=0;long start=System.nanoTime();
  for(int i=0;i<repeats;i++)for(double value:run.get())checksum+=value;
  double ms=(System.nanoTime()-start)/1e6/repeats;
  System.out.printf(Locale.ROOT,"%s,%d,%.9f,%.12g,%.12g%n",name,repeats,ms,error,checksum);
 }
 public static void main(String[] args)throws Exception {
  System.out.println("model,repeats,milliseconds_per_fit,max_reference_error,checksum");
  double[][] c=data("cox");double[] start=new double[c.length],stop=start.clone(),offset=start.clone();boolean[] event=new boolean[c.length];int[] strata=new int[c.length];double[][] cx=new double[c.length][2];
  for(int i=0;i<c.length;i++){start[i]=c[i][0];stop[i]=c[i][1];event[i]=c[i][2]==1;strata[i]=(int)c[i][3];cx[i][0]=c[i][4];cx[i][1]=c[i][5];offset[i]=c[i][6];}
  var survival=new CoxSurvivalData(start,stop,event,strata);
  for(CoxTies ties:CoxTies.values()) {
   String name="cox_"+ties.toString().toLowerCase(Locale.ROOT);
   time(name,50,()->{var f=CoxRegression.fit(survival,cx,offset,CoxOptions.defaults().withTies(ties),BackendPolicy.CPU);if(!f.converged())throw new AssertionError(name+" did not converge");return f.beta();},ref(name,"beta"),1e-7);
  }
  double[][] gd=data("gamma_frailty");double[] gt=new double[gd.length];boolean[] ge=new boolean[gd.length];double[][] gx=new double[gd.length][1];List<String> groups=new ArrayList<>();
  for(int i=0;i<gd.length;i++){gt[i]=gd[i][0];ge[i]=gd[i][1]==1;gx[i][0]=gd[i][2];groups.add("g"+(int)gd[i][3]);}
  double[] gr=new double[5];gr[0]=ref("gamma_frailty","beta")[0];System.arraycopy(ref("gamma_frailty","mode"),0,gr,1,4);
  var gammaSurvival=CoxSurvivalData.rightCensored(gt,ge);
  time("gamma_frailty",20,()->{var f=CoxGammaFrailty.fitAtVariance(gammaSurvival,gx,groups,null,.4,CoxGammaFrailtyOptions.defaults(),BackendPolicy.CPU);if(!f.converged())throw new AssertionError("gamma frailty did not converge");double[] values=new double[5];values[0]=f.beta()[0];System.arraycopy(f.logFrailties(),0,values,1,4);return values;},gr,2e-6);
  double[][] d=data("mediation");double[] y=new double[d.length],t=y.clone(),m=y.clone();double[][] cov=new double[d.length][1];List<String> gs=new ArrayList<>();
  for(int i=0;i<d.length;i++){y[i]=d[i][0];t[i]=d[i][1];m[i]=d[i][2];cov[i][0]=d[i][3];gs.add("g"+(int)d[i][4]);}
  var term=RandomEffectTerm.randomIntercept("g",gs);
  time("mediation",50,()->{var f=MediationAnalysis.fit(y,t,m,cov,OlsOptions.defaults(),BackendPolicy.CPU);return effects(f.aPath(),f.bPath(),f.directEffect(),f.totalEffect(),f.indirectEffect());},ref("mediation","effects"),2e-6);
  time("mediation_mixed",20,()->{var f=MediationAnalysis.fitMixed(y,t,m,cov,List.of(term),RemlOptions.defaults(),.95,BackendPolicy.CPU);if(!f.converged())throw new AssertionError("mixed mediation did not converge");return effects(f.aPath(),f.bPath(),f.directEffect(),f.totalEffect(),f.indirectEffect());},ref("mediation_mixed","effects"),2e-6);
  double[][] n=data("nonlinear");double[] ny=Arrays.stream(n).mapToDouble(row->row[0]).toArray();
  time("nonlinear",50,()->{var f=NonlinearFixedModel.fit(ny,new double[]{1,.1},mean(n),NonlinearModelOptions.defaults(),BackendPolicy.CPU);if(!f.converged())throw new AssertionError("nonlinear did not converge");return f.beta();},ref("nonlinear","beta"),1e-7);
  double[][] p=data("nonlinear_pedigree");double[] py=Arrays.stream(p).mapToDouble(row->row[0]).toArray();List<String> ids=new ArrayList<>();for(double[] row:p)ids.add(""+(char)('a'+(int)row[2]));
  var pedigree=Pedigree.of(List.of(PedigreeIndividual.founder("a"),PedigreeIndividual.founder("b"),new PedigreeIndividual("c","a","b"),new PedigreeIndividual("d","a","b"),PedigreeIndividual.founder("e"),PedigreeIndividual.founder("f")));
  var animal=PedigreeRandomEffectTerm.of("animal",ids,pedigree);
  double[] pb=ref("nonlinear_pedigree","beta"),pv=ref("nonlinear_pedigree","variance");
  time("nonlinear_pedigree",20,()->{var f=NonlinearMixedModel.fitPedigree(py,new double[]{1,.1},mean(p),List.of(animal),List.of(),RemlOptions.builder().varianceEstimation(VarianceEstimation.ML).varianceBounds(1e-8,100).build(),new NonlinearModelOptions(100,1e-7,1e-10,1,20),BackendPolicy.CPU);if(!f.converged())throw new AssertionError("pedigree nonlinear did not converge");var v=f.linearizedModel().varianceComponents();return new double[]{f.parameters()[0],f.parameters()[1],Math.log(v[0]),Math.log(v[1])};},new double[]{pb[0],pb[1],Math.log(pv[0]),Math.log(pv[1])},2e-5);
 }
 private static NonlinearMeanFunction mean(double[][] d) {
  return (b,i)->{double e=Math.exp(b[1]*d[i][1]);return new NonlinearMeanFunction.Evaluation(b[0]*e,new double[]{e,d[i][1]*b[0]*e});};
 }
 private static double[] effects(MediationEffect a,MediationEffect b,MediationEffect direct,MediationEffect total,MediationEffect indirect) {
  return new double[]{a.estimate(),b.estimate(),direct.estimate(),total.estimate(),indirect.estimate(),indirect.standardError()};
 }
}
