# Independent R fixtures: polar integration of a bivariate standard normal.
# This integrates the rejection region directly, not JLinAlg's conditional
# noncentral gamma series. Base R only; no SKAT moment approximation.
options(digits=17)
rho <- c(0,.25,.5,.75,1)
polar <- function(f) {
  knots<-seq(0,2*pi,length.out=65)
  sum(vapply(seq_len(64),function(i)integrate(Vectorize(f),knots[i],knots[i+1],
      rel.tol=1e-8,abs.tol=1e-10,subdivisions=3000)$value,0.))/(2*pi)
}
skato <- function(u,V) {
  E <- eigen(V,symmetric=TRUE); L <- E$vectors %*% diag(sqrt(E$values))
  kernels <- lapply(rho,function(r) (1-r)*diag(2)+r*matrix(1,2,2))
  spectra <- lapply(kernels,function(K)eigen(t(L)%*%K%*%L,symmetric=TRUE)$values)
  tail <- function(q,e) {
    if(q==0)return(1)
    scale <- exp(-q/(2*max(e)))
    scale*polar(function(t) {
      v <- sum(e*c(cos(t),sin(t))^2)
      if(v<=0)return(0)
      exp(-q/(2*v)+q/(2*max(e)))
    })
  }
  p <- mapply(function(K,e)tail(drop(t(u)%*%K%*%u),e),kernels,spectra)
  minimum <- min(p)
  critical <- vapply(spectra,function(e)uniroot(function(q)tail(q,e)/minimum-1,
      c(0,max(e)*qchisq(minimum,2,lower.tail=FALSE)),tol=1e-11)$root,0.)
  adjusted <- minimum*polar(function(t) {
    a <- L%*%c(cos(t),sin(t))
    ratios <- mapply(function(K,q)q/drop(t(a)%*%K%*%a),kernels,critical)
    exp(-min(ratios)/2-log(minimum))
  })
  c(minimum,adjusted)
}
cases <- list(moderate=list(c(2,-1),matrix(c(2,.3,.3,1),2)),
              rare=list(c(10,-3),matrix(c(2,.3,.3,1),2)),
              equal=list(c(7,4),diag(2)),
              correlated=list(c(3,1),matrix(c(2,1.1,1.1,1),2)))
rows <- lapply(names(cases),function(name){
  x<-cases[[name]];p<-skato(x[[1]],x[[2]])
  data.frame(case=name,u=paste(x[[1]],collapse=","),v=paste(as.vector(t(x[[2]])),collapse=","),minimum=p[1],adjusted=p[2])
})
# Rank-three identity covariance: a separate exact convolution of independent
# chi-square(1) burden and chi-square(2) orthogonal residual.
convolution<-function(q,r) {
  if(r==0)return(pchisq(q,3,lower.tail=FALSE))
  if(r==1)return(pchisq(q/3,1,lower.tail=FALSE))
  b<-1+2*r;limit<-sqrt(q/b);scale<-exp(-q/(2*b))
  scale*integrate(function(t)sqrt(2/pi)*exp(-t*t/2-(q-b*t*t)/(2*(1-r))+q/(2*b)),
    0,limit,rel.tol=1e-10,abs.tol=1e-12)$value+2*pnorm(limit,lower.tail=FALSE)
}
for(u in list(c(2,-1,1),c(9,4,-1))) {
  q<-sapply(rho,function(r)(1-r)*sum(u*u)+r*sum(u)^2)
  mp<-min(mapply(convolution,q,rho))
  cv<-sapply(rho,function(r)uniroot(function(q)convolution(q,r)/mp-1,c(0,3*qchisq(mp,3,lower.tail=FALSE)),tol=1e-11)$root)
  upper<-sqrt(cv[length(cv)]/3)
  fun<-function(t){bound<-min((cv[rho<1]-3*rho[rho<1]*t*t)/(1-rho[rho<1]))-t*t
    sqrt(2/pi)*exp(-t*t/2-max(0,bound)/2-log(mp))}
  knots<-seq(0,upper,length.out=65)
  integral<-sum(sapply(1:64,function(i)integrate(Vectorize(fun),knots[i],knots[i+1],rel.tol=1e-8,abs.tol=1e-10)$value))
  adj<-mp*integral+2*pnorm(upper,lower.tail=FALSE)
  rows[[length(rows)+1]]<-data.frame(case=paste0('rank3_',length(rows)),u=paste(u,collapse=','),v=paste(as.vector(diag(3)),collapse=','),minimum=mp,adjusted=adj)
}
write.table(do.call(rbind,rows),'src/test/resources/raremetal/advanced-skato.tsv',sep='\t',row.names=FALSE,quote=FALSE)
# Noncentral rank-one reference includes small tails and large shifts.
nc <- expand.grid(q=c(.5,8,60,120),shift=c(0,.5,5,20))
# Exact rank-one identity avoids pchisq's noncentral upper-tail cancellation.
nc$p <- mapply(function(q,shift)pnorm(sqrt(q)-sqrt(shift),lower.tail=FALSE)+
                  pnorm(-sqrt(q)-sqrt(shift)),nc$q,nc$shift)
write.table(nc,'src/test/resources/raremetal/noncentral-tail.tsv',sep='\t',row.names=FALSE,quote=FALSE)
# VT: covariance of the standardized singleton and two-variant burden.
V<-matrix(c(2,.3,.3,1),2);u<-c(2,-1);W<-rbind(c(1,0),c(1,1))
B<-W%*%V%*%t(W);S<-diag(1/sqrt(diag(B)));C<-S%*%B%*%S
L<-t(chol(C));maxz<-max(abs(S%*%W%*%u))
vt<-polar(function(t){a<-L%*%c(cos(t),sin(t));exp(-maxz^2/(2*max(a^2)))})
writeLines(format(vt,digits=17),'src/test/resources/raremetal/advanced-vt.txt')
# Unequal, noncentral two-term mixture: integrate one normal coordinate and
# evaluate the other with two normal tails. No gamma-series implementation.
mixed<-data.frame(q=c(5,15,40,80))
mixed$p<-sapply(mixed$q,function(q){
  scale<-exp(-q/(2*1.2))
  f<-function(x){remaining<-(q-.8*(x+sqrt(.5))^2)/1.2
    tail<-if(remaining<=0)1 else pnorm(sqrt(remaining)-sqrt(2),lower.tail=FALSE)+pnorm(-sqrt(remaining)-sqrt(2))
    dnorm(x)*tail/scale}
  knots<-sort(unique(c(seq(-12,12,by=.125),pmax(-12,pmin(12,c(-sqrt(q/.8),sqrt(q/.8))-sqrt(.5))))))
  scale*sum(sapply(seq_len(length(knots)-1),function(i)integrate(Vectorize(f),knots[i],knots[i+1],rel.tol=1e-8,abs.tol=1e-8)$value))
})
write.table(mixed,'src/test/resources/raremetal/noncentral-mixture.tsv',sep='\t',row.names=FALSE,quote=FALSE)
