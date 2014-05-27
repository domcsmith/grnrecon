#########################################################################
# Finding the correct null in the MI context
# Dominic Smith 2014/05/27
# Theory based on Bradley Efron's paper <<>>
#########################################################################

library("quantreg")
library("splines")
library("stats")

setwd("~/academic/projects/grnrecon/cmi/code/")
data.df <-read.table("../data/microarray/output/GSE2350_mi.out", sep="\t", nrows=10000)

y <- hist(data.df[,3],
          breaks=50,
          plot=FALSE)$counts

x <- hist(data.df[,3],
          breaks=50,
          plot=FALSE)$mid

spl <- model.matrix(y ~ bs(x, df=15))


hist(data.df[,3],
     breaks=50,
     plot=TRUE)

for (tau in 3/4) {
    fit <- rq(y ~ bs(x,df=15), tau=tau)
    dens.fit <- spl %*% fit$coef
    lines(x, dens.fit)
}

cbind(x, dens.fit)
