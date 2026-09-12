# Regenerate independent score-summary fixtures with R 4.6 / SKAT / CompQuadForm.
.libPaths(c("build/r-library", .libPaths()))
library(SKAT)
library(CompQuadForm)
target <- "src/test/resources/raremetal/independent-scores.tsv"
cases <- list(
  diagonal=list(u=c(2,-3,4),v=diag(c(4,9,16))),
  correlated=list(u=c(3,-2,4),v=matrix(c(4,1,.5,1,3,-.2,.5,-.2,2),3)),
  positive=list(u=c(5,4,3),v=matrix(c(4,1,.5,1,3,-.2,.5,-.2,2),3)),
  tail=list(u=c(12,-7,10),v=diag(c(1,2,3)))
)
rho <- c(0,.25,.5,.75,1)
out <- list()
for (name in names(cases)) {
  u<-cases[[name]]$u; v<-cases[[name]]$v; q<-sum(u*u)
  eigenvalues<-eigen(v,symmetric=TRUE)$values
  d<-davies(q,eigenvalues,acc=1e-8,lim=1000000)
  # Imhof provides an independent integration check when Davies is near zero.
  im<-imhof(q,eigen(v,symmetric=TRUE)$values,epsabs=1e-12,epsrel=1e-12,limit=1000000)
  # Independent rank-three spherical integration: radius^2 ~ chi-square(3)
  # is independent of a uniform direction. Integrate survival, not 1-CDF.
  sphere<-integrate(Vectorize(function(z) {
    integrate(function(phi) pchisq(q/(eigenvalues[1]*(1-z*z)*cos(phi)^2+
      eigenvalues[2]*(1-z*z)*sin(phi)^2+eigenvalues[3]*z*z),3,lower.tail=FALSE),
      0,pi/2,rel.tol=1e-10,abs.tol=1e-50,subdivisions=1000)$value*2/pi
  }),0,1,rel.tol=1e-10,abs.tol=1e-50,subdivisions=1000)$value
  qr<-matrix((1-rho)*q+rho*sum(u)^2,nrow=1)
  skato<-suppressWarnings(SKAT:::SKAT_Optimal_Get_Pvalue(qr,chol(v),rho,"optimal.adj"))
  out[[name]]<-data.frame(case=name,score=paste(u,collapse=","),covariance=paste(as.vector(t(v)),collapse=","),
    burden_beta=sum(u)/sum(v),burden_se=1/sqrt(sum(v)),burden_p=2*pnorm(-abs(sum(u)/sqrt(sum(v)))),
    skat_q=q,skat_davies=d$Qq,skat_ifault=d$ifault,skat_imhof=im$Qq,skato=skato$p.value,skat_spherical=sphere)
}
write.table(do.call(rbind,out),target,sep="\t",quote=FALSE,row.names=FALSE)
writeLines(trimws(capture.output(sessionInfo()),which="right"),"src/test/resources/raremetal/R-session.txt")
