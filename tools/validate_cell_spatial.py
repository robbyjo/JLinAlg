"""Run every cell/spatial CLI operation and independent raw-input R/NumPy comparisons.

Requires Python numpy and R limma. Writes a fresh timestamped audit directory.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jar", default="build/cli/jlinalg-0.3.6.jar")
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--r-library", required=True)
    parser.add_argument("--out")
    args = parser.parse_args()
    output = ROOT / (args.out or "build/cell-spatial-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S"))
    output.mkdir(parents=True, exist_ok=False)
    example=ROOT/"examples/cell-spatial"
    common=[]
    for f in ("counts", "cells", "samples", "features"):
        common += ["--"+f, str(example/(f+".tsv"))]
    rargs=["--rscript",args.rscript,"--r-library",args.r_library]
    design=["--reference","control","--tested","case"]
    geometry=["--units","micrometer","--coordinate-system","synthetic-v1","--resolution","cell","--radius","10.1"]
    timings={}
    def run(name,command,method,*extra,paired=False,inputs_root=None):
        dest=output/name
        inputs=list(common)
        if inputs_root:
            for f in ("counts","cells","samples","features"): inputs[inputs.index("--"+f)+1]=str(inputs_root/(f+".tsv"))
        if paired: inputs[inputs.index("--samples")+1]=str(example/"paired-samples.tsv")
        cmd=["java","-jar",str((ROOT/args.jar).resolve()),command,"--no-config","--method",method,*inputs,*map(str,extra),"--out",str(dest)]
        start=time.perf_counter()
        p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=240)
        timings[name]=time.perf_counter()-start
        (output/(name+".log")).write_text(p.stdout+p.stderr,encoding="utf-8")
        if p.returncode: raise RuntimeError(name+": "+p.stdout+p.stderr)
        assert (dest/"manifest.yaml").is_file() and Path(str(dest)+".config.yaml").is_file()
        return dest
    run("qc","single-cell","qc")
    bulk=run("pseudobulk","single-cell","pseudobulk")
    state=run("state","single-cell","state",*design,*rargs)
    pair=run("paired-state","single-cell","state",*design,"--paired","true",*rargs,paired=True)
    abundance=run("abundance","single-cell","abundance",*design)
    run("paired-abundance","single-cell","abundance",*design,"--paired","true",paired=True)
    run("pathways","single-cell","pathways",*design,"--gene-sets",example/"gene-sets.tsv")
    rep=run("representation","single-cell","representation",*rargs,"--components",6,"--clusters",3)
    run("graph","spatial","graph",*geometry)
    knn=run("knn","spatial","graph",*geometry,"--graph","knn","--k",3)
    adjacency_geometry=geometry[:-2]+["--graph","adjacency","--edges",knn/"edges.tsv"]
    run("adjacency","spatial","graph",*adjacency_geometry)
    spatial=run("autocorrelation","spatial","autocorrelation",*geometry,"--permutations",99,"--feature-list",example/"selected-features.tsv")
    run("neighborhoods","spatial","neighborhoods",*geometry,"--permutations",99)
    run("compare","spatial","compare",*geometry,*design)
    run("gradient","spatial","gradient",*geometry[:-2],*design,"--distance-column","distance","--feature-list",example/"selected-features.tsv")

    # Independent raw sparse aggregation: compare each sample/population/feature sum.
    cells={r["obs_id"]:r for r in read(example/"cells.tsv")}
    genes=[r["feature_id"] for r in read(example/"features.tsv")]
    sums={}
    for r in read(example/"counts.tsv"):
        c=cells[r["obs_id"]]; key=(c["sample_id"],c["cell_type"],r["feature_id"])
        sums[key]=sums.get(key,0)+int(r["count"])
    meta={r["aggregate_id"]:r for r in read(bulk/"pseudobulk-samples.tsv")}
    for r in read(bulk/"pseudobulk-counts.tsv"):
        for key,m in meta.items():
            assert float(r[key])==sums.get((m["sample_id"],m["cell_type"],r["feature_id"]),0)
    # Independent NumPy least squares for abundance; independent R tails below.
    x=np.array([[float(v) for k,v in r.items() if k!="sample_id"] for r in read(abundance/"design.tsv")])
    for kind in ("B","T"):
        y=np.array([float(r["asin_sqrt"]) for r in read(abundance/"proportions.tsv") if r["population"]==kind])
        beta=np.linalg.lstsq(x,y,rcond=None)[0]; df=len(y)-x.shape[1]
        variance=np.sum((y-x@beta)**2)/df; se=np.sqrt(np.linalg.inv(x.T@x)[1,1]*variance)
        r=next(r for r in read(abundance/"results.tsv") if r["population"]==kind)
        np.testing.assert_allclose([float(r["effect"]),float(r["se"])],[beta[1],se],rtol=1e-9,atol=1e-12)
    # Independent dense W definition of Moran/Geary, all sections/features.
    sparse={(r["obs_id"],r["feature_id"]):float(r["count"]) for r in read(example/"counts.tsv")}
    libraries={oid:sum(sparse.get((oid,g),0) for g in genes) for oid in cells}
    all_edges=read(spatial/"edges.tsv")
    for r in read(spatial/"autocorrelation.tsv"):
        ids=[oid for oid,c in cells.items() if c["section_id"]==r["section_id"]]
        index={oid:i for i,oid in enumerate(ids)}; n=len(ids); w=np.zeros((n,n))
        for e in all_edges:
            if e["source"] in index:
                a,b=index[e["source"]],index[e["target"]]; w[a,b]=w[b,a]=1
        y=np.array([np.log1p(1e4*sparse.get((oid,r["feature_id"]),0)/libraries[oid]) for oid in ids]); z=y-y.mean()
        moran=n*(z@w@z)/(w.sum()*(z@z))
        geary=(n-1)*np.sum(w*(y[:,None]-y[None,:])**2)/(2*w.sum()*(z@z))
        np.testing.assert_allclose([float(r["moran_i"]),float(r["geary_c"])],[moran,geary],rtol=1e-10,atol=1e-12)
    # Refit directly from raw counts and metadata in independent R script, not adapter inputs.
    reference=output/"reference.R"
    reference.write_text('''a <- commandArgs(TRUE)
.libPaths(c(a[1],.libPaths())); library(limma)
root <- a[2]; dest <- a[3]
read <- function(f) read.delim(file.path(root,f),check.names=FALSE,stringsAsFactors=FALSE)
cells <- read("cells.tsv"); raw <- read("counts.tsv"); features <- read("features.tsv")
merged <- merge(raw,cells[,c("obs_id","sample_id","cell_type")],by="obs_id")
sample <- read("samples.tsv")
composition <- table(factor(cells$sample_id,levels=sample$sample_id),cells$cell_type)
x <- model.matrix(~factor(condition,levels=c("control","case"))+factor(batch),sample)
abundance <- do.call(rbind,lapply(colnames(composition),function(kind) {
  y <- asin(sqrt(composition[,kind]/rowSums(composition)))
  fit <- summary(lm(y~x-1))$coefficients[2,]
  data.frame(population=kind,effect=fit[1],se=fit[2],p=fit[4])
}))
write.table(abundance,file.path(dest,"reference-abundance.tsv"),sep="\\t",quote=FALSE,row.names=FALSE)
for(paired in c(FALSE,TRUE)) {
  sample <- read(if(paired) "paired-samples.tsv" else "samples.tsv")
  sample$condition <- factor(sample$condition,levels=c("control","case"))
  x <- if(paired) model.matrix(~condition+factor(donor_id),sample) else model.matrix(~condition+factor(batch),sample)
  for(kind in sort(unique(cells$cell_type))) {
    subset <- merged[merged$cell_type==kind,]
    counts <- unclass(xtabs(count ~ factor(feature_id,levels=features$feature_id)+factor(sample_id,levels=sample$sample_id),subset))
    eligible <- rowSums(counts>=10)>=2
    v <- voom(counts[eligible,],x,lib.size=colSums(counts),normalize.method="none",plot=FALSE)
    fit <- eBayes(lmFit(v,x))
    z <- data.frame(feature_id=rownames(counts)[eligible],effect=fit$coefficients[,2],se=fit$stdev.unscaled[,2]*sqrt(fit$s2.post),p=fit$p.value[,2])
    write.table(z,file.path(dest,paste0("reference-",paired,"-",kind,".tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
  }
}
capture.output(sessionInfo(),file=file.path(dest,"reference-session.txt"))
''',encoding="utf-8")
    subprocess.run([args.rscript,"--vanilla",str(reference),args.r_library,str(example),str(output)],check=True,cwd=ROOT,timeout=120)
    for flag,dest in (("FALSE",state),("TRUE",pair)):
        actual={(r["population"],r["feature_id"]):r for r in read(dest/"results.tsv")}
        for kind in ("B","T"):
            for r in read(output/f"reference-{flag}-{kind}.tsv"):
                np.testing.assert_allclose([float(actual[kind,r["feature_id"]][key]) for key in ("effect","se","p")],[float(r[key]) for key in ("effect","se","p")],rtol=1e-8,atol=1e-11)
    actual={r["population"]:r for r in read(abundance/"results.tsv")}
    for r in read(output/"reference-abundance.tsv"):
        np.testing.assert_allclose([float(actual[r["population"]][key]) for key in ("effect","se","p")],[float(r[key]) for key in ("effect","se","p")],rtol=1e-8,atol=1e-11)
    assert len(read(rep/"embedding.tsv"))==len(cells)
    assert any(r["feature_id"]=="GENE01" and r["population"]=="T" and float(r["effect"])>0 for r in read(state/"results.tsv"))

    # Identifier round trip through both R branches, including values R ordinarily coerces.
    unusual=output/"unusual-identifiers"; unusual.mkdir()
    feature_map={g:g for g in genes}; feature_map.update({genes[0]:"NA",genes[1]:"001",genes[2]:"1",genes[3]:'gene"quoted'})
    sample_map={f"S{i}":f"{i:03}" for i in range(1,9)}
    cell_map={oid:f"{i:05}" for i,oid in enumerate(cells)}
    for f in ("counts","cells","samples","features"):
        rows=read(example/(f+".tsv"))
        for r in rows:
            if "feature_id" in r:r["feature_id"]=feature_map[r["feature_id"]]
            if "obs_id" in r:r["obs_id"]=cell_map[r["obs_id"]]
            if "sample_id" in r:r["sample_id"]=sample_map[r["sample_id"]]
        with (unusual/(f+".tsv")).open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
    renamed=run("identifier-state","single-cell","state",*design,*rargs,inputs_root=unusual)
    renamed_rep=run("identifier-representation","single-cell","representation",*rargs,"--components",6,"--clusters",3,inputs_root=unusual)
    renamed_results={(r["population"],r["feature_id"]):r for r in read(renamed/"results.tsv")}
    for r in read(state/"results.tsv"):
        actual=renamed_results[r["population"],feature_map[r["feature_id"]]]
        np.testing.assert_allclose([float(actual[k]) for k in ("effect","se","p")],[float(r[k]) for k in ("effect","se","p")],rtol=1e-12,atol=1e-13)
    assert [r["obs_id"] for r in read(renamed_rep/"embedding.tsv")]==list(cell_map.values())
    assert {r["feature_id"] for r in read(renamed_rep/"variable-features.tsv")}==set(feature_map.values())
    filtered=run("no-eligible-genes","single-cell","state",*design,*rargs,"--min-gene-samples",9)
    assert all(r["status"]=="low_counts" and r["p"]=="NA" for r in read(filtered/"results.tsv"))
    assert all(r["status"]=="not_tested" and r["output_directory"]!="NA" for r in read(filtered/"population-status.tsv"))
    report={"operations":len(timings),"observations":len(cells),"samples":8,"features":len(genes),"seconds":timings,"validation":"raw aggregation, independent limma refits (independent and paired), NumPy/R abundance, dense-W spatial statistics"}
    (output/"validation.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2));print(output)

if __name__=="__main__":main()
