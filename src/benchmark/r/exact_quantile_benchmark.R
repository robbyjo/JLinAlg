.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(quantreg))
options(digits=17)
args <- commandArgs(trailingOnly=TRUE)
repeats <- if (length(args)) as.integer(args[1]) else 100L
stopifnot(repeats>0)
root <- 'src/benchmark/resources/exact-quantile-benchmark'
dir.create(root,showWarnings=FALSE,recursive=TRUE)
refs <- data.frame(case=character(),tau=double(),objective=double(),coefficient=integer(),value=double())
cat('case,method,repeats,milliseconds_per_fit,objective_difference,checksum\n')
for (shape in list(c(256,4),c(1024,6),c(4096,8))) {
  n <- shape[1];p <- shape[2];i <- 0:(n-1)
  X <- matrix(1,n,p);X[,2]<-(i-(n-1)/2)/(n/4)
  for (j in 2:(p-1)) X[,j+1] <- sin(i*(.17+j*.13))+cos(i*(.07+j*.03))*.2
  beta <- c(1.2,-.35,sapply(2:(p-1),function(j) .5*sin(j)))
  y <- drop(X%*%beta)+(.35+.12*abs(X[,2]))*(sin(i*1.731)+.4*cos(i*2.37))
  y[i%%97==0] <- y[i%%97==0]+3
  name <- paste0('n',n,'p',p);tau<-.75
  write.table(cbind(y,X),file.path(root,paste0(name,'.csv')),sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)
  objective <- function(b) {e<-y-drop(X%*%b);sum(ifelse(e>=0,tau*e,(tau-1)*e))}
  reference <- rq.fit.br(X,y,tau=tau)
  target <- objective(reference$coefficients)
  refs <- rbind(refs,data.frame(case=name,tau=tau,objective=target,coefficient=0:(p-1),value=reference$coefficients))
  for (method in c('br','fnb')) {
    run <- if (method=='br') function() rq.fit.br(X,y,tau=tau)$coefficients else function() rq.fit.fnb(X,y,tau=tau,eps=1e-9)$coefficients
    error <- abs(objective(run())-target);stopifnot(is.finite(error),error<1e-6)
    for (k in 1:10) run()
    checksum <- 0
    started <- Sys.time()
    for(k in seq_len(repeats)) checksum<-checksum+sum(run())
    elapsed <- as.numeric(difftime(Sys.time(),started,units='secs'))
    cat(sprintf('%s,%s,%d,%.9f,%.12g,%.12g\n',name,method,repeats,elapsed*1000/repeats,error,checksum))
  }
}
write.table(refs,file.path(root,'reference.csv'),sep=',',row.names=FALSE,quote=FALSE)
