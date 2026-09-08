.libPaths(c('build/r-library',.libPaths()))
options(digits=17)
out <- 'src/test/resources/pipeline-audit'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
y <- c(1,2,2,4,5,7,7,9)
g <- rep(0:3,each=2)
fit <- summary(lm(y~g))$coefficients[2,]
x <- 0:19
strong <- x + rep(c(1,-1),10)*1e-7
strongfit <- summary(lm(strong~x))$coefficients[2,]
values <- c(ols_beta=fit[1],ols_se=fit[2],ols_t=fit[3],ols_p=fit[4],
            strong_se=strongfit[2], dosage_var=var(c(1-1e-8,1,1+1e-8)),
            single_kernel=pchisq(3.5/2,1,lower.tail=FALSE))
# Independent polar-coordinate integral: R^2 ~ chi-square(2), angle uniform.
tail2 <- function(q,a=1,b=10) integrate(function(t)
    exp(-q/(2*(a*cos(t)^2+b*sin(t)^2))),0,pi/2,
    rel.tol=1e-12,abs.tol=1e-13)$value*2/pi
values <- c(values,rank2_tail=tail2(20),rank2_critical=uniroot(function(q) tail2(q)-.05,c(0,200),tol=1e-11)$root)
values <- c(values,rank1_q50=pchisq(50,1,lower.tail=FALSE),rank1_q200=pchisq(200,1,lower.tail=FALSE),rank1_q700=pchisq(700,1,lower.tail=FALSE))
write.table(data.frame(metric=names(values),value=unname(values)),file.path(out,'reference.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
p <- c(0,.01,.04,.01,.2,1,NA,.9)
# Failed rows are not part of the tested family; make that denominator explicit.
q <- rep(NA_real_,length(p));q[!is.na(p)] <- p.adjust(p[!is.na(p)],'BH')
write.table(data.frame(p=p,q=q),file.path(out,'bh-reference.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NaN')
print(values)
