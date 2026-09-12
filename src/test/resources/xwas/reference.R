# Independent base-R QR/optim references and reproducible synthetic examples.
# Run from the repository root: Rscript src/test/resources/xwas/reference.R
options(digits=17)
dir.create('examples/xwas', recursive=TRUE, showWarnings=FALSE)
dest <- 'src/test/resources/xwas'
write_tab <- function(x, name) write.table(x, file.path('examples/xwas',name), sep='\t', quote=FALSE, row.names=FALSE)
write_matrix <- function(x, labels, name) {
  colnames(x) <- labels
  write_tab(data.frame(row=labels,x,check.names=FALSE),name)
}
expect <- function(x,name) write.table(data.frame(value=as.vector(x)),file.path(dest,name),sep='\t',quote=FALSE,row.names=FALSE)
set.seed(421)
n <- 1200; nt <- 4; B <- 30; M <- 100000; N <- 10000
ld <- 1+19*((seq_len(n)*17)%%97)/96
wld <- 1+ld/2
lambda <- c(.5,.6,.7,.8)
S <- tcrossprod(lambda)+diag(c(.2,.3,.15,.25))
z <- t(vapply(seq_len(n),function(i) as.vector(t(chol(diag(nt)+N/M*ld[i]*S))%*%rnorm(nt)),numeric(nt)))
tab <- data.frame(variant=paste0('rs',seq_len(n)),ea='A',oa='C',ld_score=ld,weight_ld=wld)
traits <- LETTERS[1:nt]
for(i in seq_len(nt)){tab[[paste0('z_',traits[i])]]<-z[,i];tab[[paste0('n_',traits[i])]]<-N}
write_tab(tab,'ldsc.tsv')
blocks <- floor((seq_len(n)-1)*B/n)
clamp <- function(x,lo,hi) max(lo,min(hi,x))
wls <- function(y,x,w,keep=rep(TRUE,n)) {
  unname(lm.wfit(cbind(x,1)[keep,,drop=FALSE],y[keep],w[keep])$coefficients)
}
reg <- function(i,j,h1=NULL,h2=NULL) {
  y<-z[,i]*z[,j];x<-N*ld/M
  coef<-c((mean(y)-(i==j))/mean(x),as.numeric(i==j))
  for(step in 0:2){
    if(i==j) variance<-2*(coef[2]+clamp(coef[1],0,1)*N*pmax(ld,1)/M)^2 else {
      a<-h1$coef[2]+clamp(h1$coef[1],0,1)*N*pmax(ld,1)/M
      b<-h2$coef[2]+clamp(h2$coef[1],0,1)*N*pmax(ld,1)/M
      c<-coef[2]+clamp(coef[1],-1,1)*N*pmax(ld,1)/M
      variance<-a*b+c*c
    }
    w<-1/(pmax(wld,1)*variance);coef<-wls(y,x,w)
  }
  deletes<-t(vapply(0:(B-1),function(b)wls(y,x,w,blocks!=b),numeric(2)))
  list(coef=coef,deletes=deletes)
}
fits<-matrix(vector('list',nt*nt),nt,nt)
for(i in 1:nt)fits[[i,i]]<-reg(i,i)
for(i in 1:(nt-1))for(j in (i+1):nt)fits[[i,j]]<-reg(i,j,fits[[i,i]],fits[[j,j]])
estimated<-matrix(0,nt,nt);ints<-estimated;deletes<-matrix(0,B,nt*(nt+1)/2);k<-1
for(i in 1:nt)for(j in i:nt){estimated[i,j]<-estimated[j,i]<-fits[[i,j]]$coef[1];ints[i,j]<-ints[j,i]<-fits[[i,j]]$coef[2];deletes[,k]<-fits[[i,j]]$deletes[,1];k<-k+1}
V<-(B-1)*cov(deletes)*(B-1)/B # cov uses divisor B-1; desired (B-1)/B * sum
expect(t(estimated),'ldsc-S.tsv');expect(t(V),'ldsc-V.tsv');expect(t(ints),'ldsc-intercepts.tsv')

# Independent factor covariance fixture: SPD genetic S and correlated known sampling V.
pairs<-unlist(lapply(1:nt,function(i)paste0(traits[i],':',traits[i:nt])))
S[1,2]<-S[2,1]<-S[1,2]+.008
S[2,4]<-S[4,2]<-S[2,4]-.004
m<-length(pairs);Vf<-outer(1:m,1:m,function(i,j).2^abs(i-j))*0.0001
write_matrix(S,traits,'factor-S.tsv');write_matrix(Vf,pairs,'factor-V.tsv')
vec<-function(a)unlist(lapply(1:nrow(a),function(i)a[i,i:ncol(a)]))
mu<-function(par)vec(tcrossprod(par[1:nt])+diag(exp(par[nt+1:nt])))
obj<-function(par){r<-vec(S)-mu(par);as.numeric(crossprod(r,solve(Vf,r)))}
gradient<-function(par){
  jac<-sapply(seq_along(par),function(i){d<-rep(0,length(par));d[i]<-1e-5;(mu(par+d)-mu(par-d))/(2e-5)})
  as.vector(-2*t(jac)%*%solve(Vf,vec(S)-mu(par)))
}
opt<-optim(c(lambda,log(c(.2,.3,.15,.25))),obj,gr=gradient,method='BFGS',control=list(reltol=1e-15,maxit=10000))
stopifnot(opt$convergence==0)
par<-opt$par
J<-sapply(seq_along(par),function(i){d<-rep(0,length(par));d[i]<-1e-5;(mu(par+d)-mu(par-d))/(2e-5)})
pc<-solve(t(J)%*%solve(Vf,J));delta<-diag(c(rep(1,nt),exp(par[nt+1:nt])));pc<-delta%*%pc%*%delta
expect(c(par[1:nt],exp(par[nt+1:nt]),obj(par)),'factor-fit.tsv');expect(t(pc),'factor-covariance.tsv')
corr<-matrix(.15,nt,nt);diag(corr)<-1
write_matrix(corr,traits,'sampling-correlation.tsv')
betas<-c(.05,.04,.09,.07);ses<-c(.01,.015,.02,.012);C<-diag(ses)%*%corr%*%diag(ses)
b<-as.numeric(crossprod(par[1:nt],solve(C,betas))/crossprod(par[1:nt],solve(C,par[1:nt])))
se<-as.numeric(1/sqrt(crossprod(par[1:nt],solve(C,par[1:nt]))));r<-betas-b*par[1:nt];q<-as.numeric(crossprod(r,solve(C,r)))
expect(c(b,se,q),'factor-association.tsv')
gt<-data.frame(variant='rsFollowup',ea='A',oa='C');for(i in 1:nt){gt[[paste0('beta_',traits[i])]]<-betas[i];gt[[paste0('se_',traits[i])]]<-ses[i]};write_tab(gt,'factor-gwas.tsv')

R<-matrix(c(1,.2,-.1,.2,1,.3,-.1,.3,1),3)
sdg<-c(.6,.7,.5);weights<-c(.4,-.2,.1);zs<-c(2,-1,3)
a<-weights*sdg;variance<-as.numeric(crossprod(a,R%*%a));stat<-sum(a*zs)/sqrt(variance)
expect(c(stat,2*pnorm(-abs(stat)),variance),'twas.tsv')
write_matrix(R,c('v1','v2','v3'),'ld.tsv')
write_tab(data.frame(variant=c('v1','v2','v3'),ea='A',oa='C',sd=sdg),'reference.tsv')
write_tab(data.frame(variant=c('v1','v2','v3'),ea=c('A','C','A'),oa=c('C','A','C'),z=c(2,1,3)),'gwas.tsv')
write_tab(data.frame(model=c('gene_tissue1','gene_tissue1','gene_tissue1','gene_tissue2','gene_tissue2'),variant=c('v1','v2','v3','v1','v2'),ea='A',oa='C',weight=c(weights,.1,.3)),'weights.tsv')
W<-cbind(a,c(.1,.3,0)*sdg);G<-t(W)%*%R%*%W;u<-as.vector(t(W)%*%zs)
expect(as.numeric(crossprod(u,solve(G,u))),'twas-joint.tsv')

# Ridge comparator uses augmented QR directly, independently of coordinate descent.
x<-matrix(rnorm(160),80,2);y<-1+2*x[,1]-.7*x[,2]+rnorm(80,sd=.3)
train<-data.frame(sample=paste0('train',1:80),y=y,x1=x[,1],x2=x[,2]);write_tab(train,'train.tsv')
xt<-matrix(rnorm(80),40,2);yt<-1+2*xt[,1]-.7*xt[,2]+rnorm(40,sd=.3)
write_tab(data.frame(sample=paste0('test',1:40),y=yt,x1=xt[,1],x2=xt[,2]),'test.tsv')
sx<-sqrt(colMeans(sweep(x,2,colMeans(x))^2));standard<-scale(x,center=TRUE,scale=sx)
ridge<-qr.solve(rbind(standard,sqrt(80*.1)*diag(2)),c(y-mean(y),0,0))/sx
intercept<-mean(y)-sum(colMeans(x)*ridge);pred<-as.vector(intercept+xt%*%ridge)
cal<-coef(lm(yt~pred));expect(c(intercept,ridge,mean(y),sqrt(mean((yt-pred)^2)),1-sum((yt-pred)^2)/sum((yt-mean(y))^2),cal),'score.tsv')
write_tab(data.frame(variant=c('v1','v2'),ea=c('A','G'),oa=c('C','T'),weight=c(.3,-.1)),'polygenic-weights.tsv')
write_tab(data.frame(variant=c('v1','v2'),ea=c('C','G'),oa=c('A','T')),'target-alleles.tsv')
write_tab(data.frame(sample=paste0('external',1:8),v1=c(0,0,1,1,1,2,2,2),v2=c(0,1,0,1,2,0,1,2),case_status=c(1,1,0,1,0,0,0,1)),'target-dosages.tsv')
cat('Independent xWAS fixtures generated\n')
