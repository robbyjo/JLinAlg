# Independent fixtures: stats::phyper/p.adjust, limma and BiasedUrn; no JLinAlg calls.
# Run from the repository root with R 4.6.x, limma 3.68.5 and BiasedUrn 2.0.12.
.libPaths(c("build/r-library", .libPaths()))
Sys.setlocale("LC_COLLATE", "C")
stopifnot(as.character(packageVersion("limma")) == "3.68.5",
          as.character(packageVersion("BiasedUrn")) == "2.0.12")
options(digits=17)
dest <- "src/test/resources/enrichment"
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x, name) write.table(x, file.path(dest,name), sep="\t", row.names=FALSE, quote=FALSE)

cases <- rbind(c(2,1,2,2,2), c(20,80,15,8,3), c(200,800,100,45,2),
               c(30,70,20,5,.2), c(30,70,95,28,3), c(100,900,100,25,.5),
               c(50,50,50,40,10), c(5,995,100,4,8), c(10,90,30,8,1),
               c(20,80,90,12,.1), c(100,100,100,90,.1))
colnames(cases) <- c("red","white","draws","observed","odds")
expected <- apply(cases,1,function(x) BiasedUrn::pWNCHypergeo(x[4]-1,x[1],x[2],x[3],x[5],precision=1e-12,lower.tail=FALSE))
write_tsv(data.frame(cases,p=expected),"wallenius.tsv")
x <- c(0,0,1,1,0,1,0,0,0,1,1,1,0,0,1,0)
write_tsv(data.frame(x=x,pwf=limma::tricubeMovingAverage(x)),"smoother.tsv")

# Explicit normalized probe-to-gene edges: mapping is shared input, statistics are independent.
probes <- data.frame(probe=character(),genes=character(),selected=integer())
for(g in 1:40) for(j in seq_len(1+g%%5)) {
  genes <- sprintf("G%02d",g)
  if(j==1 && g%%3==0) genes <- paste(genes,sprintf("G%02d",1+g%%40),sep=";")
  selected <- as.integer((g%%7==0 && j==1) || (g<=8 && j==2))
  probes <- rbind(probes,data.frame(probe=sprintf("p%03d",nrow(probes)+1),genes=genes,selected=selected))
}
write_tsv(probes,"probes.tsv")
sets <- list(early=sprintf("G%02d",1:10), spread=sprintf("G%02d",seq(7,35,7)),
             zero=c("G11","G12"), late=sprintf("G%02d",30:40), multi=c("G06","G07","G21","G22"))
writeLines(vapply(names(sets),function(id) paste(c(id,id,sets[[id]]),collapse="\t"),character(1)),file.path(dest,"sets.gmt"))

# Execute the published missMethyl functions at a pinned revision. The sole replaced
# function is the array-manifest adapter, because the fixture already supplies gene IDs.
ref <- "build/enrichment-reference"
dir.create(ref,recursive=TRUE,showWarnings=FALSE)
commit <- "c5d0a518aaacd98c9f31c4b6ca48239246363ad4"
for(name in c("gometh.R","gsameth.R","getMappedEntrezIDs.R")) {
  file <- file.path(ref,name)
  if(!file.exists(file)) download.file(paste0("https://raw.githubusercontent.com/bioc/missMethyl/",commit,"/R/",name),file,mode="wb")
  source(file)
}
edges <- do.call(rbind,lapply(seq_len(nrow(probes)),function(i) data.frame(cpg=probes$probe[i],entrezid=strsplit(probes$genes[i],";",fixed=TRUE)[[1]],group="Body")))
.getFlatAnnotation <- function(...) edges
sig <- probes$probe[probes$selected==1]
mapped <- getMappedEntrezIDs(sig,probes$probe)
pwf <- .estimatePWF(mapped$de,as.vector(mapped$equiv))
ans <- gsameth(sig,probes$probe,collection=sets)
weighted <- vapply(sets,function(s) sum(mapped$fract.counts$frac[mapped$sig.eg %in% s]),numeric(1))
odds <- vapply(sets,function(s) mean(pwf[mapped$universe %in% s])/mean(pwf[!mapped$universe %in% s]),numeric(1))
write_tsv(data.frame(term=names(sets),weighted=weighted,odds=odds,p=ans[names(sets),"P.DE"],BH=ans[names(sets),"FDR"]),"gsameth.tsv")
writeLines(c(paste("R",getRversion()),paste("limma",packageVersion("limma")),
             paste("BiasedUrn",packageVersion("BiasedUrn")),paste("missMethyl reference",commit),
             "LC_COLLATE=C; custom normalized edges; equiv.cpg=TRUE; fract.counts=TRUE; prior.prob=TRUE"),file.path(dest,"reference.txt"))
