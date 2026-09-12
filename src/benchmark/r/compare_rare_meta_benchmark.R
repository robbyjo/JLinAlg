args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop("Supply JLinAlg single.tsv and RAREMETAL singlevar.results")
j <- read.delim(args[1],check.names=FALSE)
r <- read.delim(args[2],comment.char="#",header=FALSE,check.names=FALSE)
id <- paste(r[[1]],r[[2]],r[[3]],r[[4]],sep=":")
stopifnot(nrow(j)==20000L,nrow(r)==20000L,!anyDuplicated(id),!anyDuplicated(j$feature_id))
j <- j[match(id,j$feature_id),]
stopifnot(!anyNA(j$feature_id))
delta <- c(beta_absolute=max(abs(j$beta-r[[8]])),
           se_absolute=max(abs(j$se-r[[9]])),
           p_absolute=max(abs(j$p_value-r[[11]])),
           p_relative=max(abs(j$p_value-r[[11]])/r[[11]]))
print(delta,digits=12)
# Native output has six significant digits; use corresponding rounding tolerance.
stopifnot(delta[1]<1e-6,delta[2]<1e-6,delta[3]<1e-6,delta[4]<1e-5)
