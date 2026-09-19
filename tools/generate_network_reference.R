# Independent reference fixtures for native Java methods. Run after generate_followup_examples.py.
# Optional first argument specifies the library containing glmnet.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) .libPaths(c(args[1],.libPaths()))
root <- "src/test/resources/network-reference"
dir.create(root,recursive=TRUE,showWarnings=FALSE)
a <- as.matrix(read.delim("examples/followup/discovery.tsv",row.names=1))[,1:6]
b <- as.matrix(read.delim("examples/followup/replication.tsv",row.names=1))[,1:6]
std <- function(x) sweep(sweep(x,2,colMeans(x)),2,sqrt(colMeans(sweep(x,2,colMeans(x))^2)),"/")
x <- std(a)
coef <- matrix(0,ncol(x),ncol(x),dimnames=list(colnames(x),colnames(x)))
for(j in seq_len(ncol(x))) {
  fit <- glmnet::glmnet(x[,-j,drop=FALSE],x[,j],lambda=.15,alpha=1,intercept=FALSE,
                        standardize=FALSE,control=list(thresh=1e-18,maxit=1000000))
  beta <- as.numeric(coef(fit))[-1]
  # Independently solve the active-set KKT equations, then check inactive inequalities.
  xx <- x[,-j,drop=FALSE]; active <- which(beta != 0)
  beta[active] <- solve(crossprod(xx[,active,drop=FALSE])/nrow(x),
    crossprod(xx[,active,drop=FALSE],x[,j])/nrow(x) - .15*sign(beta[active]))
  gradient <- as.numeric(crossprod(xx,x[,j]-xx%*%beta)/nrow(x))
  stopifnot(max(abs(gradient[active]-.15*sign(beta[active])))<1e-10,
            all(abs(gradient[-active])<=.15+1e-10))
  coef[j,-j] <- beta
}
write.table(data.frame(feature=rownames(coef),coef),file.path(root,"lasso.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
ca <- cor(a); cb <- cor(b); pairs <- which(upper.tri(ca),arr.ind=TRUE)
z <- (atanh(ca[pairs])-atanh(cb[pairs]))/sqrt(1/(nrow(a)-3)+1/(nrow(b)-3))
p <- 2*pnorm(-abs(z))
write.table(data.frame(source=colnames(a)[pairs[,1]],target=colnames(a)[pairs[,2]],r_a=ca[pairs],r_b=cb[pairs],z=z,p=p,bh=p.adjust(p,"BH")),file.path(root,"differential.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(sample_id=rownames(a),a),file.path(root,"a.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(sample_id=rownames(b),b),file.path(root,"b.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(root,"sessionInfo.txt"))
