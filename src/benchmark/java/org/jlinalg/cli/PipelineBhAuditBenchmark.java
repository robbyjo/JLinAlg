package org.jlinalg.cli;
import java.nio.file.*;
import java.util.*;
/** End-to-end disk BH timings: write, sort, adjust, publish; checksum outside timer. */
public class PipelineBhAuditBenchmark {
 public static void main(String[] args)throws Exception {
  int n=20000;Path root=Path.of("build/pipeline-audit");
  Files.createDirectories(root);
  System.out.println("method,rep,n,milliseconds,checksum");
  for(int rep=-3;rep<7;rep++) {
   Path output=Files.createTempFile(root,"bh-timing-",".tsv");Files.delete(output);
   long start=System.nanoTime();
   try(var bh=new ExternalBh(output,false)) {
    bh.writeHeader(List.of("id","p"));
    for(int i=0;i<n;i++) {
     double p=i%101==0?Double.NaN:((i*7919L)%100003)/100003.0;
     bh.write(List.of(Integer.toString(i),Double.toString(p)),p);
    }
    bh.finish();
   }
   double ms=(System.nanoTime()-start)/1e6,checksum=0;
   List<String> lines=Files.readAllLines(output);
   if(lines.size()!=n+1)throw new AssertionError("lost rows");
   for(int i=1;i<lines.size();i++){
    String[] f=lines.get(i).split("\t",-1);
    if(Integer.parseInt(f[0])!=i-1)throw new AssertionError("row order");
    if(!f[2].isEmpty())checksum+=Double.parseDouble(f[2])*(1+(i-1)%17);
   }
   Files.delete(output);
   if(!Double.isFinite(checksum)||Math.abs(checksum-176709.28000988538)>1e-7)throw new AssertionError("BH checksum differs from frozen R reference: "+checksum);
   System.out.printf(Locale.ROOT,"java-disk-bh,%d,%d,%.9f,%.17g%n",rep,n,ms,checksum);
  }
 }
}
