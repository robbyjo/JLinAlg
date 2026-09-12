options(digits=17)
suppressPackageStartupMessages(library(survival))
out <- 'src/test/resources/conditional-score-export'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(260912)
n <- 180
x <- rnorm(n); g <- matrix(rbinom(n*4,2,.35),n,4)
g[,2] <- ifelse(seq_len(n)%%3==0,g[,1],g[,2])
o <- seq(-.2,.2,length.out=n)
binary <- rbinom(n,1,plogis(-.4+.3*x+.8*g[,1]+o))
count <- rpois(n,exp(.3+.2*x+.35*g[,1]+o))
start <- rep(c(0,.2),length.out=n)
stop <- start + ceiling(rexp(n,exp(.15*x+.3*g[,1]+o))*5)/5
event <- rbinom(n,1,.7)
g[c(4,17,61),3] <- NA
d <- data.frame(IID=paste0('S',1:n),x=x,o=o,binary=binary,count=count,start=start,stop=stop,event=event)
d$x[11] <- NA
write.table(d,file.path(out,'phenotype.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
con<-file(file.path(out,'variants.vcf'),'w')
writeLines(c('##fileformat=VCFv4.2','##contig=<ID=1,length=1000000>',
 '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Alternate allele dosage">',
 paste(c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',d$IID),collapse='\t')),con)
for(j in 1:4)writeLines(paste(c(1,100*j,paste0('g',j),'A','G','.','PASS','.','DS',ifelse(is.na(g[,j]),'.',g[,j])),collapse='\t'),con)
close(con)
keep<-!is.na(d$x);d<-d[keep,];g<-g[keep,,drop=FALSE];raw<-g
for(j in 1:4)g[is.na(g[,j]),j]<-mean(g[,j],na.rm=TRUE)
emit<-function(con,key,v)writeLines(paste0(key,'=',paste(format(v,digits=17,scientific=TRUE,trim=TRUE),collapse=',')),con)
for(family in c('binomial','poisson','gaussian'))for(conditioned in c(FALSE,TRUE)) {
 y<-if(family=='binomial')d$binary else d$count
 X<-cbind(1,d$x);targets<-1:4
 if(conditioned){X<-cbind(X,g[,1]);targets<-2:4}
 fit<-glm.fit(X,y,family=get(family)(),offset=d$o,control=glm.control(epsilon=1e-12,maxit=100))
 stopifnot(fit$converged)
 w<-if(family=='binomial')fit$fitted.values*(1-fit$fitted.values) else if(family=='poisson')fit$fitted.values else rep(1,nrow(X))
 dispersion<-if(family=='gaussian')sum((y-fit$fitted.values)^2)/fit$df.residual else 1
 G<-g[,targets,drop=FALSE]
 residual<-G-X%*%solve(crossprod(X,w*X),crossprod(X,w*G))
 U<-as.vector(crossprod(residual,y-fit$fitted.values))/dispersion;V<-crossprod(residual,w*residual)/dispersion
 con<-file(file.path(out,paste0(family,if(conditioned)'-conditional' else '','.properties')),'w')
 emit(con,'u',U);emit(con,'v',as.vector(t(V)));emit(con,'nullBeta',fit$coefficients)
 if(family=='binomial') {
  emit(con,'nCases',sum(y==1));emit(con,'nControls',sum(y==0))
  emit(con,'calledCases',colSums(!is.na(raw[y==1,targets,drop=FALSE])))
  emit(con,'calledControls',colSums(!is.na(raw[y==0,targets,drop=FALSE])))
  emit(con,'acCases',colSums(raw[y==1,targets,drop=FALSE],na.rm=TRUE))
  emit(con,'acControls',colSums(raw[y==0,targets,drop=FALSE],na.rm=TRUE))
 }
 close(con)
}
for(right in c(FALSE,TRUE))for(ties in c('efron','breslow'))for(conditioned in c(FALSE,TRUE)) {
 X<-matrix(d$x,ncol=1);targets<-1:4
 if(conditioned){X<-cbind(X,g[,1]);targets<-2:4}
 response<-if(right)Surv(d$stop,d$event) else Surv(d$start,d$stop,d$event)
 fit<-coxph(response~X+offset(d$o),ties=ties,control=coxph.control(eps=1e-11))
 Z<-cbind(X,g[,targets,drop=FALSE]);p<-ncol(X)
 at<-coxph(response~Z+offset(d$o),ties=ties,
     init=c(coef(fit),rep(0,length(targets))),control=coxph.control(iter.max=0))
 H<-solve(vcov(at));s<-colSums(residuals(at,type='score'));a<-(p+1):ncol(Z);b<-1:p
 U<-s[a]-H[a,b,drop=FALSE]%*%solve(H[b,b,drop=FALSE],s[b])
 V<-H[a,a,drop=FALSE]-H[a,b,drop=FALSE]%*%solve(H[b,b,drop=FALSE],H[b,a,drop=FALSE])
 con<-file(file.path(out,paste0('cox-',if(right)'right-' else '',ties,if(conditioned)'-conditional' else '','.properties')),'w')
 emit(con,'u',as.vector(U));emit(con,'v',as.vector(t(V)));emit(con,'nullBeta',coef(fit));close(con)
}
writeLines(c(R.version.string,paste('survival',packageVersion('survival'))),file.path(out,'versions.txt'))
