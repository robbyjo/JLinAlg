.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(lmerTest))
options(digits=17)
out <- 'src/benchmark/resources/r-reference'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
g <- rep(0:1499,each=8); j <- rep(0:7,1500); x <- (j-3.5)/2
y <- 2+.8*x+1.3*sin(g*.71)+.7*cos(g*.113)+(.4*sin(g*.33)+.2*cos(g*.47))*x+.3*sin((seq_along(g))*1.731)
d <- data.frame(y=y,x=x,g=sprintf('g%04d',g),h='h0',w=1,o=0)
write.table(d,file.path(out,'mixed-scalability.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
ctrl <- lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000,rhoend=1e-9))
fit <- lmerTest::lmer(y~x+(1+x|g),d,control=ctrl)
times <- replicate(3,system.time(lmerTest::lmer(y~x+(1+x|g),d,control=ctrl))['elapsed'])
s <- coef(summary(fit))
values <- c(logLik=as.numeric(logLik(fit)),beta0=fixef(fit)[1],beta1=fixef(fit)[2],
  df0=s[1,3],df1=s[2,3],fitSeconds=median(times))
names(values) <- c('logLik','beta0','beta1','df0','df1','fitSeconds')
writeLines(c(paste0('# R ',getRversion(),'; lme4 ',packageVersion('lme4'),'; lmerTest ',packageVersion('lmerTest')),
  paste(names(values),format(values,digits=17,scientific=TRUE),sep='=')),file.path(out,'mixed-scalability.properties'))
print(values)
# pbkrtest KR is deliberately not run at n=12000: its observation inverse is
# dense. Java's memory-bounded KR is checked against balanced-model DF=1499.
