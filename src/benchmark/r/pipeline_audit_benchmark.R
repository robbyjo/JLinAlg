.libPaths(c('build/r-library',.libPaths()))
options(digits=17)
out <- 'src/benchmark/resources/pipeline-audit'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
n<-20000;i<-0:(n-1)
p<-ifelse(i%%101==0,NA,((i*7919)%%100003)/100003)
bh<-function(){q<-rep(NA_real_,n);q[!is.na(p)]<-p.adjust(p[!is.na(p)],'BH');q}
rows<-lapply(-3:6,function(rep){start<-proc.time()[3];for(k in 1:20)q<-bh();
 data.frame(method='R-memory-bh',rep=rep,n=n,milliseconds=1000*(proc.time()[3]-start)/20,
 checksum=sum(q*(1+i%%17),na.rm=TRUE))})
write.csv(do.call(rbind,rows),file.path(out,'r-bh-raw.csv'),row.names=FALSE)
n<-400;m<-128;i<-0:(n-1)
y<-sin(.13*i)+.3*cos(.07*i);x<-cbind(1,cos(.07*i))
g<-sapply(0:(m-1),function(j)ifelse((i+j)%%97==0,NA,sin(.023*(i+1)*(j+1))+.2*cos(.031*(i+1)*(j+1))))
scan<-function(){t(vapply(1:m,function(j){z<-g[,j];z[is.na(z)]<-mean(z,na.rm=TRUE)
 fit<-lm.fit(cbind(x,z),y);rss<-sum(fit$residuals^2);v<-chol2inv(qr.R(fit$qr))*rss/fit$df.residual
 beta<-fit$coefficients[3];se<-sqrt(v[3,3]);c(beta=beta,se=se,p=2*pt(-abs(beta/se),fit$df.residual))},numeric(3)))}
rows<-lapply(-5:8,function(rep){start<-proc.time()[3];for(k in 1:25)result<-scan();
 data.frame(method='R-lm.fit',rep=rep,n=n,markers=m,milliseconds=1000*(proc.time()[3]-start)/25,
 checksum=sum((result[,1]+result[,2]+result[,3])*(1:m)))})
write.csv(do.call(rbind,rows),file.path(out,'r-ols-raw.csv'),row.names=FALSE)
write.table(data.frame(marker=0:(m-1),scan()),'src/test/resources/pipeline-audit/ols-scan-reference.tsv',sep='\t',row.names=FALSE,quote=FALSE)
print(do.call(rbind,rows))
tail2 <- function() integrate(function(t)exp(-20/(2*(cos(t)^2+10*sin(t)^2))),
    0,pi/2,rel.tol=1e-12,abs.tol=1e-13)$value*2/pi
rows<-lapply(-5:8,function(rep){start<-proc.time()[3];checksum<-0
 for(k in 1:5000)checksum<-checksum+tail2()
 data.frame(method='R-rank-two',rep=rep,n=5000,milliseconds=1000*(proc.time()[3]-start)/5000,checksum=checksum)})
write.csv(do.call(rbind,rows),file.path(out,'r-kernel-raw.csv'),row.names=FALSE)

# Same positive gamma-mixture identity as Java, including coefficient generation
# inside every timed call. Fixed 600-term truncation has PGF omitted mass <1e-12
# of this q=100 tail; independently assert it rather than timing unconverged work.
tail3 <- function() {
 lambda<-c(1,.5,.25);beta<-min(lambda);terms<-600L;r<-1-beta/lambda
 coeff<-numeric(terms+1L);coeff[1]<-prod(sqrt(beta/lambda))
 g<-vapply(1:terms,function(j)sum(r^j)/2,0)
 for(k in 1:terms)coeff[k+1]<-sum(g[1:k]*rev(coeff[1:k]))/k
 p<-sum(coeff*pchisq(100/beta,df=length(lambda)+2*(0:terms),lower.tail=FALSE))
 z<-1+(1/max(r)-1)/2
 logBound<-sum(log(beta/lambda)/2-log1p(-r*z)/2)-(terms+1)*log(z)
 stopifnot(logBound<log(p)+log(1e-12),abs(p/2.5052827321161237e-23-1)<1e-10)
 p
}
rows<-lapply(-5:8,function(rep){start<-proc.time()[3];checksum<-0
 for(k in 1:25)checksum<-checksum+tail3()
 data.frame(method='R-rank-three',rep=rep,n=25,milliseconds=1000*(proc.time()[3]-start)/25,checksum=checksum)})
write.csv(do.call(rbind,rows),file.path(out,'r-rank-three-raw.csv'),row.names=FALSE)
