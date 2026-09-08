.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(quantreg))
options(digits=17)
root <- 'src/test/resources/exact-quantile'
dir.create(root,showWarnings=FALSE,recursive=TRUE)
refs <- data.frame(case=character(),tau=double(),objective=double(),unique=integer())
coefficients <- data.frame(case=character(),tau=double(),coefficient=integer(),value=double())
save_case <- function(name,y,X,taus=c(.01,.1,.5,.9,.99),unique=TRUE) {
  write.table(cbind(y,X),file.path(root,paste0(name,'.csv')),sep=',',quote=FALSE,row.names=FALSE,col.names=FALSE)
  for (tau in taus) {
    fit <- suppressWarnings(rq.fit.br(X,y,tau=tau))
    e <- y-drop(X%*%fit$coefficients)
    objective <- sum(ifelse(e>=0,tau*e,(tau-1)*e))
    refs <<- rbind(refs,data.frame(case=name,tau=tau,objective=objective,unique=as.integer(unique)))
    coefficients <<- rbind(coefficients,data.frame(case=name,tau=tau,coefficient=seq_len(ncol(X))-1L,value=fit$coefficients))
  }
}
i <- 0:100; x <- (i-50)/20
save_case('regular',2+.3*x-.2*sin(i*.7)+sin(i*2.1)*(1+.1*abs(x)),cbind(1,x,sin(i*.7)))
i <- 0:79; x <- (i%%9)-4
save_case('ties',round(2+.4*x+sin(i*1.7),0),cbind(1,x),unique=FALSE)
save_case('even',c(-3,-2,0,1,4,7,8,10),matrix(1,8,1),taus=c(.25,.5,.75),unique=FALSE)
save_case('median',c(1,2,3),matrix(1,3,1),taus=.5)
X <- cbind(1,c(-2,0,1,3),c(1,0,3,2),c(0,1,2,0))
save_case('saturated',c(3,-2,8,1),X,taus=c(.1,.5,.9))
i <- 0:126; x <- (i-63)/30
y <- .7*x-.4*cos(i*.53)+sin(i*1.9)
y[c(5,32,100)] <- c(15,-20,25)
save_case('no_intercept',y,cbind(x,cos(i*.53)))
i <- 0:44; x <- (i-22)/10
save_case('perfect',1.25-.75*x,cbind(1,x),taus=c(.1,.5,.9))
write.table(refs,file.path(root,'reference.csv'),sep=',',quote=FALSE,row.names=FALSE)
write.table(coefficients,file.path(root,'coefficients.csv'),sep=',',quote=FALSE,row.names=FALSE)
cat('Generated',nrow(refs),'exact quantile fits using quantreg',as.character(packageVersion('quantreg')),'\n')
