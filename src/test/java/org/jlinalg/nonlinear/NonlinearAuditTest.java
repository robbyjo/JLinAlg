package org.jlinalg.nonlinear;
import static org.junit.jupiter.api.Assertions.*;
import static org.jlinalg.nonlinear.RemainingAuditData.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.pedigree.*;
import org.jlinalg.reml.*;
import org.junit.jupiter.api.Test;
public class NonlinearAuditTest {
 @Test void largeResponseOffsetCannotCertifyAnUnfittedSlope() {
  for(double shift:new double[]{0,1e12,1e15}) {
   double[] y=new double[100];double xy=0,xx=0;
   for(int i=0;i<y.length;i++) {
    double x=i%9-4;y[i]=shift+2*x+.25*(i%2==0?1:-1);
    xy+=x*(y[i]-shift);xx+=x*x;
   }
   double optimum=xy/xx,referenceSse=0;
   for(int i=0;i<y.length;i++){double r=y[i]-(shift+optimum*(i%9-4));referenceSse+=r*r;}
   var f=NonlinearFixedModel.fit(y,new double[]{0},(b,i)->new NonlinearMeanFunction.Evaluation(shift+b[0]*(i%9-4),new double[]{i%9-4}),NonlinearModelOptions.defaults(),BackendPolicy.CPU);
   assertTrue(f.objective()<=referenceSse+1e-12,"must reach the analytic solution's evaluated SSE at shift="+shift);
   if(shift==0) {assertTrue(f.converged());assertEquals(optimum,f.beta()[0],1e-14);}
   else assertFalse(f.converged(),"a quantized, nonstationary plateau is not convergence");
   if(shift==1e15) {
    assertEquals(optimum,f.beta()[0],1e-14);assertEquals(6.25,f.objective(),0);
    assertEquals(2,f.iterations(),"do not drift across an equal-SSE plateau");
   }
  }
 }
 @Test void stationarityIsInvariantToResponseUnits() {
  for(double units:new double[]{1e-12,1,1e12}) {
   double[] y=new double[100];double xy=0,xx=0;
   for(int i=0;i<y.length;i++){double x=i%9-4;y[i]=units*(2*x+.25*(i%2==0?1:-1));xy+=x*(y[i]/units);xx+=x*x;}
   var f=NonlinearFixedModel.fit(y,new double[]{0},(b,i)->new NonlinearMeanFunction.Evaluation(units*b[0]*(i%9-4),new double[]{units*(i%9-4)}),NonlinearModelOptions.defaults(),BackendPolicy.CPU);
   assertTrue(f.converged());assertEquals(xy/xx,f.beta()[0],1e-14);
   assertEquals(6.2485207100591716,f.objective()/units/units,1e-12);
  }
 }
 @Test void scaledParametersCannotFalselyConverge() {
  var f=NonlinearFixedModel.fit(new double[]{1,2,3},new double[]{0},(b,i)->new NonlinearMeanFunction.Evaluation(b[0]*1e12,new double[]{1e12}),NonlinearModelOptions.defaults(),BackendPolicy.CPU);
  assertTrue(f.converged());assertEquals(2e-12,f.beta()[0],1e-25);assertEquals(2,f.objective(),1e-12);
 }
 @Test void covarianceUsesReturnedJacobianAndIterationBudgetIsExact() {
  var f=NonlinearFixedModel.fit(new double[]{3,4,5},new double[]{1},(b,i)->new NonlinearMeanFunction.Evaluation(b[0]*b[0],new double[]{2*b[0]}),new NonlinearModelOptions(1,1e-8,1e-10,1,20),BackendPolicy.CPU);
  assertFalse(f.converged());assertEquals(1,f.iterations());assertEquals(2.5,f.beta()[0],1e-14);
  assertEquals(Math.sqrt(f.residualVariance()/(12*f.beta()[0]*f.beta()[0])),f.standardErrors()[0],1e-14);
 }
 @Test void noisyExponentialMatchesNls() throws Exception {
  double[][] d=data("nonlinear");double[] y=Arrays.stream(d).mapToDouble(r->r[0]).toArray();
  var f=NonlinearFixedModel.fit(y,new double[]{1,.1},(b,i)->new NonlinearMeanFunction.Evaluation(b[0]*Math.exp(b[1]*d[i][1]),new double[]{Math.exp(b[1]*d[i][1]),d[i][1]*b[0]*Math.exp(b[1]*d[i][1])}),NonlinearModelOptions.defaults(),BackendPolicy.CPU);
  assertTrue(f.converged());for(int j=0;j<2;j++){assertEquals(ref("nonlinear","beta",j),f.beta()[j],1e-7);assertEquals(ref("nonlinear","se",j),f.standardErrors()[j],1e-8);}
  assertEquals(ref("nonlinear","sse",0),f.objective(),1e-10);
 }
 @Test void pedigreeNonlinearMlMatchesIndependentDenseLikelihood() throws Exception {
  double[][] d=data("nonlinear_pedigree");double[] y=Arrays.stream(d).mapToDouble(r->r[0]).toArray();
  var pedigree=Pedigree.of(List.of(PedigreeIndividual.founder("a"),PedigreeIndividual.founder("b"),new PedigreeIndividual("c","a","b"),new PedigreeIndividual("d","a","b"),PedigreeIndividual.founder("e"),PedigreeIndividual.founder("f")));
  List<String> ids=new ArrayList<>();for(double[] row:d)ids.add(""+(char)('a'+(int)row[2]));
  var term=PedigreeRandomEffectTerm.of("animal",ids,pedigree);
  var f=NonlinearMixedModel.fitPedigree(y,new double[]{1,.1},(b,i)->new NonlinearMeanFunction.Evaluation(b[0]*Math.exp(b[1]*d[i][1]),new double[]{Math.exp(b[1]*d[i][1]),d[i][1]*b[0]*Math.exp(b[1]*d[i][1])}),List.of(term),List.of(),RemlOptions.builder().varianceEstimation(VarianceEstimation.ML).varianceBounds(1e-8,100).build(),new NonlinearModelOptions(100,1e-7,1e-10,1,20),BackendPolicy.CPU);
  System.out.println("pedigree beta="+Arrays.toString(f.parameters())+" variance="+Arrays.toString(f.linearizedModel().varianceComponents())+" ll="+f.linearizedModel().logLikelihood()+" converged="+f.converged()+" iterations="+f.iterations());
  assertTrue(f.converged());for(int j=0;j<2;j++){assertEquals(ref("nonlinear_pedigree","beta",j),f.parameters()[j],2e-5);assertEquals(ref("nonlinear_pedigree","variance",j),f.linearizedModel().varianceComponents()[j],2e-5);}
  assertEquals(ref("nonlinear_pedigree","ll",0),f.linearizedModel().logLikelihood(),2e-6);
  assertArrayEquals(f.parameters(),f.linearizedModel().beta(),1e-6);
 }
 @Test void invalidParametersAndZeroResidualDfFailExplicitly() {
  NonlinearMeanFunction mean=(b,i)->new NonlinearMeanFunction.Evaluation(b[0],new double[]{1});
  assertThrows(IllegalArgumentException.class,()->NonlinearFixedModel.fit(new double[]{1,2},new double[]{Double.NaN},mean,NonlinearModelOptions.defaults(),BackendPolicy.CPU));
  assertThrows(IllegalArgumentException.class,()->NonlinearFixedModel.fit(new double[]{1,2},new double[]{0,0},mean,NonlinearModelOptions.defaults(),BackendPolicy.CPU));
 }
 @Test void mixedStepControlsAreHonoredAndUnfinishedStepsDoNotConverge() {
  double[] y=new double[40];List<String> groups=new ArrayList<>();
  for(int i=0;i<40;i++){y[i]=2+.2*i+Math.sin(i*1.7)+2*Math.sin(i/5);groups.add("g"+(i/5));}
  var term=org.jlinalg.mixed.RandomEffectTerm.randomIntercept("g",groups);
  NonlinearMeanFunction mean=(b,i)->new NonlinearMeanFunction.Evaluation(b[0]+b[1]*i,new double[]{1,i});
  var full=NonlinearMixedModel.fit(y,new double[]{0,0},mean,List.of(term),RemlOptions.defaults(),NonlinearModelOptions.defaults(),BackendPolicy.CPU);
  var shortFit=NonlinearMixedModel.fit(y,new double[]{0,0},mean,List.of(term),RemlOptions.defaults(),new NonlinearModelOptions(1,1e-8,1e-10,.1,20),BackendPolicy.CPU);
  assertTrue(full.converged());assertFalse(shortFit.converged());assertEquals(1,shortFit.iterations());
  for(int j=0;j<2;j++)assertEquals(.1*full.parameters()[j],shortFit.parameters()[j],1e-8);
 }
 public static void main(String[] args)throws Exception{run(NonlinearAuditTest.class);}
}
