#########################################################################
# Finding the correct null in the MI context
# Dominic Smith 2014/05/27
# Theory based on Bradley Efron's paper <<>>
#########################################################################

library("quantreg")
library("splines")
library("stats")

setwd("~/academic/projects/grnrecon/cmi/code/")
data.df <-read.table("../data/microarray/output/GSE2350_mi.out", sep="\t", nrows=100000)

# transform data to z scale
z <- (data.df[,3]-mean(data.df[,3]))/sd(data.df[,3])

y <- hist(z,
          breaks=500,
          plot=FALSE,
          freq=FALSE)$density

x <- hist(z,
          breaks=500,
          plot=FALSE,
          freq=FALSE)$mid

spl <- model.matrix(y ~ bs(x, df=25))


hist(z,
     breaks=50,
     plot=TRUE,
     freq=FALSE,
     col="gray",
     ylim=c(0,1.5))

for (tau in 1/2) {
    fit <- rq(y ~ bs(x,df=25), tau=tau)
    dens.fit <- spl %*% fit$coef
    lines(x, dens.fit, col="red")
}

# find centre of main peak as estimator for mean of empirical null (delta0)
delta0 <- x[which(abs(predict(fit))==max(abs(predict(fit))))]
quad.data <- subset(cbind(x,log(dens.fit)), x > delta0-1.5 & x < delta0+1.5)
quadfit <- lm(quad.data[,2]~quad.data[,1]+I(quad.data[,1]^2))
sigma0 <- (-2*quadfit[[1]][3])^(-1/2)

# find half width as estimator for sd of empirical null (sigma0)
#curve(quadfit[[1]][1]+quadfit[[1]][2]*x+quadfit[[1]][3]*x^2, from = -2, to=2, type="l")
#points(x, log(dens.fit), ylim=c(-8,10))


z0 <- rnorm(1000, delta0, sigma0)
lines(density(z0))
