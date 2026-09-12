# Synthetic independent cohorts with valid, block-diagonal score covariances.
root <- "build/rare-meta-benchmark"
dir.create(root,recursive=TRUE,showWarnings=FALSE)
nvariants <- 20000L
n <- 1000L
header <- c("##ProgramName=RareMetalWorker","##Version=4.14.0","##Samples=1000","##AnalyzedSamples=1000",
 "##Families=1000","##AnalyzedFamilies=1000","##Founders=1000","##AnalyzedFounders=1000",
 "##TraitSummaries\tmin\t25th\tmedian\t75th\tmax\tmean\tvariance",
 "##Trait\t-3\t-1\t0\t1\t3\t0\t1",
 "#CHROM\tPOS\tREF\tALT\tN_INFORMATIVE\tFOUNDER_AF\tALL_AF\tINFORMATIVE_ALT_AC\tCALL_RATE\tHWE_PVALUE\tN_REF\tN_HET\tN_ALT\tU_STAT\tSQRT_V_STAT\tALT_EFFSIZE\tPVALUE")
v <- 19.8
for (cohort in 1:2) {
  file <- file(file.path(root,paste0("cohort",cohort,".score")),"w")
  writeLines(header,file)
  for(i in seq_len(nvariants)) {
    u <- sqrt(v)*(sin(i*.731+cohort)+cos(i*.23))
    writeLines(paste(1,sprintf("%d",i*100L),"A","G",n,.01,.01,20,1,1,980,20,0,
      format(u,digits=17),format(sqrt(v),digits=17),format(u/v,digits=17),
      format(2*pnorm(-abs(u/sqrt(v))),digits=17),sep="\t"),file)
  }
  close(file)
  cv<-file(file.path(root,paste0("cohort",cohort,".cov")),"w")
  writeLines(c("##ProgramName=RareMetalWorker","##Version=4.14.0","#CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES"),cv)
  for(i in seq_len(nvariants))writeLines(paste(1,sprintf("%d",i*100L),sprintf("%d",i*100L),v/n,sep="\t"),cv)
  close(cv)
}
# Binary connections keep file lists LF-only for the native Linux reference.
lf <- function(lines,path) { con<-file(path,"wb");on.exit(close(con));writeLines(lines,con) }
lf(c("cohort\tscores", "A\tcohort1.score.gz", "B\tcohort2.score.gz"),file.path(root,"cohorts.tsv"))
lf(c("cohort1.score.gz","cohort2.score.gz"),file.path(root,"summaryfiles"))
lf(c("cohort1.cov.gz","cohort2.cov.gz"),file.path(root,"covfiles"))
