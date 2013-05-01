pkgname <- "sn"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('sn')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("T.Owen")
### * T.Owen

flush(stderr()); flush(stdout())

### Name: T.Owen
### Title: Owen's function
### Aliases: T.Owen
### Keywords: math

### ** Examples

owen <- T.Owen(1:10, 2)



cleanEx()
nameEx("ais")
### * ais

flush(stderr()); flush(stdout())

### Name: ais
### Title: Australian Institute of Sport data
### Aliases: ais
### Keywords: datasets

### ** Examples

data(ais, package="sn")
attach(ais)
pairs(ais[,c(3:4,10:13)], main = "AIS data")
plot(Wt~sport)



cleanEx()
nameEx("cp.to.dp")
### * cp.to.dp

flush(stderr()); flush(stdout())

### Name: cp.to.dp
### Title: Conversion between equivalent parametrizations
### Aliases: cp.to.dp dp.to.cp
### Keywords: distribution

### ** Examples

cp <- dp.to.cp(c(30,30,2,4))
dp <- cp.to.dp(cp)



cleanEx()
nameEx("dmsn")
### * dmsn

flush(stderr()); flush(stdout())

### Name: dmsn
### Title: Multivariate skew-normal distribution
### Aliases: dmsn pmsn rmsn
### Keywords: distribution multivariate

### ** Examples

x <- seq(-3,3,length=15)
xi <- c(0.5, -1)
Omega <- diag(2)
Omega[2,1] <- Omega[1,2] <- 0.5
alpha <- c(2,-6)
pdf <- dmsn(cbind(x,2*x-1), xi, Omega, alpha)
rnd <- rmsn(10,  xi, Omega, alpha)
p1 <- pmsn(c(2,1), xi, Omega, alpha)
p2 <- pmsn(c(2,1), xi, Omega, alpha, abseps=1e-12, maxpts=10000)



cleanEx()
nameEx("dmst")
### * dmst

flush(stderr()); flush(stdout())

### Name: dmst
### Title: Multivariate skew-t distribution
### Aliases: dmst pmst rmst
### Keywords: distribution multivariate

### ** Examples

x <- seq(-4,4,length=15)
xi <- c(0.5, -1)
Omega <- diag(2)
Omega[2,1] <- Omega[1,2] <- 0.5
alpha <- c(2,2)
pdf <- dmst(cbind(x,2*x-1), xi, Omega, alpha, df=5)
rnd <- rmst(10,  xi, Omega, alpha, 6)
p1 <- pmst(c(2,1), xi, Omega, alpha, df=5)
p2 <- pmst(c(2,1), xi, Omega, alpha, df=5, abseps=1e-12, maxpts=10000)



cleanEx()
nameEx("dsn")
### * dsn

flush(stderr()); flush(stdout())

### Name: dsn
### Title: Skew-Normal Distribution
### Aliases: dsn psn qsn rsn
### Keywords: distribution

### ** Examples

pdf <- dsn(seq(-3,3,by=0.1), shape=3)
cdf <- psn(seq(-3,3,by=0.1), shape=3)
qu <- qsn(seq(0.1,0.9,by=0.1), shape=-2)
rn <- rsn(100, 5, 2, 5)



cleanEx()
nameEx("dsn2.plot")
### * dsn2.plot

flush(stderr()); flush(stdout())

### Name: dsn2.plot
### Title: Plot of Bivariate Skew-normal Density Function
### Aliases: dsn2.plot
### Keywords: distribution

### ** Examples

x <- y <- seq(-5, 5, length=35)
dsn2.plot(x, y, c(-1,2), diag(c(1,2.5)), c(2,-3))



cleanEx()
nameEx("dst")
### * dst

flush(stderr()); flush(stdout())

### Name: dst
### Title: Skew-t Distribution
### Aliases: dst pst qst rst
### Keywords: distribution

### ** Examples

pdf <- dst(seq(-4,4,by=0.1), shape=3, df=5)
rnd <- rst(100, 5, 2, -5, 8)
q <- qst(c(0.25,0.5,0.75), shape=3, df=5)
pst(q, shape=3, df=5)  # must give back c(0.25,0.5,0.75)




cleanEx()
nameEx("dst2.plot")
### * dst2.plot

flush(stderr()); flush(stdout())

### Name: dst2.plot
### Title: Plot of bivariate skew-t density function
### Aliases: dst2.plot
### Keywords: distribution

### ** Examples

x <- y <- seq(-5, 5, length=35)
dst2.plot(x, y, c(-1,2), diag(c(1,2.5)), c(2,-3), df=5)



cleanEx()
nameEx("frontier")
### * frontier

flush(stderr()); flush(stdout())

### Name: frontier
### Title: Simulated sample from a skew-normal distribution
### Aliases: frontier
### Keywords: datasets

### ** Examples

data(frontier, package="sn")
a <- sn.2logL.profile(y=frontier)
a <- sn.2logL.profile(y=frontier, param.range=c(0.8,1.6,10,30),
        use.cp=FALSE, npts=11)



cleanEx()
nameEx("gamma1.to.lambda")
### * gamma1.to.lambda

flush(stderr()); flush(stdout())

### Name: gamma1.to.lambda
### Title: Converts skewness to shape parameter of skew-normal distribution
### Aliases: gamma1.to.lambda
### Keywords: distribution

### ** Examples

gamma1.to.lambda(seq(-0.95, 0.95, length=11))



cleanEx()
nameEx("msn.affine")
### * msn.affine

flush(stderr()); flush(stdout())

### Name: msn.affine
### Title: Affine transformation of a multivariate skew-normal or skew-t
###   variable
### Aliases: msn.affine mst.affine
### Keywords: multivariate distribution

### ** Examples

dp<- list(xi=c(1,1,2), Omega=toeplitz(1/1:3), alpha=c(3,-1,2))
A <- matrix(c(1,-1,1,3,0,-2), 2, 3, byrow=TRUE) 
dp1 <- msn.affine(dp, 1:2, A)
#
dp$df <- 5
dp2<-  mst.affine(dp,,A[1,,drop=FALSE])
dp3<-  mst.affine(dp,,A[1,,drop=FALSE], drop=FALSE)
if(zapsmall(dp2$scale^2 - dp3$Omega)) print("something wrong here!")



cleanEx()
nameEx("msn.cond.plot")
### * msn.cond.plot

flush(stderr()); flush(stdout())

### Name: msn.cond.plot
### Title: Plot of the density of a conditional skew-normal variate
### Aliases: msn.cond.plot
### Keywords: multivariate distribution

### ** Examples

Omega <- diag(3)+0.5*outer(rep(1,3),rep(1,3))
a<- msn.cond.plot(rep(0,3), Omega, 1:3, 3, -0.75)



cleanEx()
nameEx("msn.conditional")
### * msn.conditional

flush(stderr()); flush(stdout())

### Name: msn.conditional
### Title: Cumulants and distribution of a skew-normal variate after
###   conditioning
### Aliases: msn.conditional
### Keywords: multivariate distribution

### ** Examples

Omega <- diag(3)+0.5*outer(rep(1,3),rep(1,3))
a<- msn.conditional(rep(0,3), Omega, 1:3, 3, -0.75)



cleanEx()
nameEx("msn.fit")
### * msn.fit

flush(stderr()); flush(stdout())

### Name: msn.fit
### Title: Fitting multivariate skew-normal distributions
### Aliases: msn.fit
### Keywords: distribution regression

### ** Examples

data(ais, package="sn")
attach(ais)
# a simple-sample case
b <- msn.fit(y=cbind(Ht,Wt))
#
# a regression case:
a <- msn.fit(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-6))
#
# refine the previous outcome
a1 <- msn.fit(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-9), start=a$dp)



cleanEx()
nameEx("msn.marginal")
### * msn.marginal

flush(stderr()); flush(stdout())

### Name: msn.marginal
### Title: Marginal components of a multivariate skew-normal distribution
### Aliases: msn.marginal
### Keywords: multivariate distribution

### ** Examples

xi <- c(10,0,-30)
Omega <- 5*diag(3)+outer(1:3,1:3)
alpha <- c(1,-3,5)
msn.marginal(xi,Omega,alpha,c(3,1))
msn.marginal(dp=list(xi=xi,Omega=Omega,alpha=alpha), comp=3)



cleanEx()
nameEx("msn.mle")
### * msn.mle

flush(stderr()); flush(stdout())

### Name: msn.mle
### Title: Maximum likelihood estimation for a multivariate skew-normal
###   distribution
### Aliases: msn.mle
### Keywords: distribution regression

### ** Examples

data(ais, package="sn")
attach(ais)
# a simple-sample case
a <- msn.mle(y=cbind(Ht,Wt))
#
# a regression case:
b  <- msn.mle(X=cbind(1,Ht,Wt), y=ssf)
b1 <- msn.mle(X=cbind(1,Ht,Wt), y=ssf, algorithm="Nelder-Mead")
b2 <- msn.mle(X=cbind(1,Ht,Wt), y=ssf, start=b1$dp)



cleanEx()
nameEx("msn.quantities")
### * msn.quantities

flush(stderr()); flush(stdout())

### Name: msn.quantities
### Title: Quantities related to the multivariate skew-normal distribution.
### Aliases: msn.quantities
### Keywords: multivariate distribution

### ** Examples

Omega <- 5*diag(3)+outer(1:3,1:3)
msn.quantities(c(0,0,1), Omega, c(-2,2,3))



cleanEx()
nameEx("mst.fit")
### * mst.fit

flush(stderr()); flush(stdout())

### Name: mst.fit
### Title: Fitting multivariate skew-t distributions
### Aliases: mst.fit
### Keywords: distribution regression

### ** Examples

data(ais, package="sn")
attach(ais)
# a simple-sample case
b <- mst.fit(y=cbind(Ht,Wt))
#
# a regression case:
a <- mst.fit(X=cbind(1,Ht,Wt), y=bmi)
#
# refine the previous outcome
a1 <- mst.fit(X=cbind(1,Ht,Wt), y=bmi, start=a$dp)



cleanEx()
nameEx("mst.mle")
### * mst.mle

flush(stderr()); flush(stdout())

### Name: mst.mle
### Title: Maximum likelihood estimation for a (multivariate) skew-t
###   distribution
### Aliases: mst.mle st.mle
### Keywords: distribution regression

### ** Examples

data(ais, package="sn")
attach(ais)
X.mat <- model.matrix(~lbm+sex)
b <- sn.mle(X.mat, bmi)
# 
b <- mst.mle(y=cbind(Ht,Wt))
#
# a multivariate regression case:
a <- mst.mle(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-6))
#
# refine the previous outcome
a1 <- mst.mle(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-9), start=a$dp)



cleanEx()
nameEx("sample.centralmoments")
### * sample.centralmoments

flush(stderr()); flush(stdout())

### Name: sample.centralmoments
### Title: Sample centralmoments
### Aliases: sample.centralmoments
### Keywords: univar

### ** Examples

data(ais, package='sn')
mom <- sample.centralmoments(ais[,"bmi"])
st.cumulants.inversion(cum=c(mom[1:3], mom[4]-3*mom[2]^2))
# parameters of the fitted ST distribution



cleanEx()
nameEx("sn.2logL.profile")
### * sn.2logL.profile

flush(stderr()); flush(stdout())

### Name: sn.2logL.profile
### Title: Twice profile relative negative loglikelihood for skew-normal
###   models
### Aliases: sn.2logL.profile
### Keywords: distribution

### ** Examples

data(ais, package="sn")
attach(ais)
a <- sn.2logL.profile(y=bmi)
## Not run: 
##D a <- sn.2logL.profile(y=bmi, use.cp=FALSE, param.range=c(3,6,1,5))
##D a <- sn.2logL.profile(X=cbind(1,lbm), y=bmi, param.range=c(0.5,0.95), npts=31)
##D #
##D data(frontier, package="sn")
##D a <- sn.2logL.profile(y=frontier, param.range=c(0.8,2, 2,30),
##D         use.cp=FALSE, npts=16)
##D 	
## End(Not run)



cleanEx()
nameEx("sn.Einfo")
### * sn.Einfo

flush(stderr()); flush(stdout())

### Name: sn.Einfo
### Title: Expected Fisher information for SN distribution parameters
### Aliases: sn.Einfo
### Keywords: distribution

### ** Examples

info <- sn.Einfo(dp=c(0,1,5), n=3)
#
data(ais, package="sn")
M <- model.matrix(~ais$"Ht")
mle <- sn.mle(X=M, y=ais$"Wt", plot.it=FALSE)
info <- sn.Einfo(cp=mle$cp, x=M)



cleanEx()
nameEx("sn.cumulants")
### * sn.cumulants

flush(stderr()); flush(stdout())

### Name: sn.cumulants
### Title: Cumulants of the skew-normal distribution
### Aliases: sn.cumulants
### Keywords: distribution

### ** Examples

sn.cumulants(shape=c(0,2.5,5,10), n=5)
sn.cumulants(dp=c(10,3,-8), n=6)



cleanEx()
nameEx("sn.em")
### * sn.em

flush(stderr()); flush(stdout())

### Name: sn.em
### Title: Fitting Skew-normal variables using the EM algorithm
### Aliases: sn.em
### Keywords: regression distribution

### ** Examples

data(ais, package="sn")
attach(ais)
#
a<-sn.em(y=bmi)
#
a<-sn.em(X=cbind(1,lbm,lbm^2),y=bmi)
#
M<-model.matrix(~lbm+I(ais$sex))
b<-sn.em(M,bmi)
#
fit <- sn.em(y=bmi, fixed=c(NA, 2, 3), l.eps=0.001)



cleanEx()
nameEx("sn.mle")
### * sn.mle

flush(stderr()); flush(stdout())

### Name: sn.mle
### Title: Maximum likelihood estimation for skew-normal models
### Aliases: sn.mle
### Keywords: regression distribution

### ** Examples

data(ais, package="sn")
attach(ais)
a<-sn.mle(y=bmi)
#
a<-sn.mle(X=cbind(1,lbm),y=bmi)
#
b<-sn.mle(X=model.matrix(~lbm+sex), y=bmi)



cleanEx()
nameEx("sn.mle.grouped")
### * sn.mle.grouped

flush(stderr()); flush(stdout())

### Name: sn.mle.grouped
### Title: Maximum likelihood estimation of SN and ST distribution for
###   grouped data
### Aliases: sn.mle.grouped st.mle.grouped
### Keywords: distribution

### ** Examples

data(ais, package="sn")
attach(ais)
breaks<- c(130,160, seq(170, 190, by=2.5), 200, 230)
f <- cut(Ht[sex=="female"], breaks = breaks)
freq <- tabulate(f, length(levels(f)))
b1 <- sn.mle.grouped(breaks, freq)
b2 <- st.mle.grouped(breaks, freq, start=c(b1$end,log(5)))
print(b2$dp)
#
us.income <- c(0,seq(from=0.2, to=1.8, by=0.1), 2.0, 2.5, 5.0, Inf)
mid <- (us.income[-1]+us.income[-length(us.income)])/2
mid[length(mid)] <- 6.5
cum.freq<- c(1.78, 3.25, 5.56, 8.16, 11.12, 14.21, 17.54, 20.78, 24.00,
             27.52, 30.77, 34.21, 37.56, 40.70, 44.41, 47.85, 51.22, 
             57.60, 72.12, 96.40, 100) / 100
freq<- round(diff(c(0,cum.freq*34660)))
a <- st.mle.grouped(breaks=log(us.income), freq, trace=TRUE,
        start=c(1.2, log(0.9), -2.1, log(20)))
print(a$dp)



cleanEx()
nameEx("sn.mmle")
### * sn.mmle

flush(stderr()); flush(stdout())

### Name: sn.mmle
### Title: Modified maximum likelihood estimation for skew-normal ans
###   skew-t models
### Aliases: sn.mmle st.mmle
### Keywords: regression distribution

### ** Examples

data(ais, package="sn")
attach(ais)
a <-  sn.mmle(y=bmi)
#
M <- model.matrix(~lbm+sex)
b <- sn.mmle(M,bmi)



cleanEx()
nameEx("st.2logL.profile")
### * st.2logL.profile

flush(stderr()); flush(stdout())

### Name: st.2logL.profile
### Title: Twice profile relative negative loglikelihood for skew-t models
### Aliases: st.2logL.profile
### Keywords: distribution

### ** Examples

data(ais, package="sn")
attach(ais)
a <- st.2logL.profile(y=bmi, xlab="alpha", ylab="log(df)")
## Not run: 
##D a <- st.2logL.profile(y=bmi, fixed.comp=4, fixed.values=log(c(1,25)), npts=26)
##D a <- st.2logL.profile(X=cbind(1,lbm), y=bmi, fixed.comp=5,  
##D          fixed.values=log(c(5,25)), xlab="log(df)", npts=26)
##D a <- st.2logL.profile(X=cbind(1,Ht), y=Wt, fixed.comp=c(4,5),
##D          fixed.values=cbind(c(-1,5), log(c(2,25))),
##D          xlab="alpha", ylab="log(df)", npts=12)
##D 	 
## End(Not run)



cleanEx()
nameEx("st.cumulants")
### * st.cumulants

flush(stderr()); flush(stdout())

### Name: st.cumulants
### Title: Cumulants of the skew-t distribution
### Aliases: st.cumulants st.cumulants.inversion
### Keywords: distribution

### ** Examples

st.cumulants(shape=c(0,3,9), df=5)
cum <- st.cumulants(dp=c(10, 2, -8, 5.2))
st.cumulants.inversion(cum)
#
data(ais, package='sn')
mom <- sample.centralmoments(ais[,"bmi"])
st.cumulants.inversion(cum=c(mom[1:3],mom[4]-3*mom[2]^2))
# parameters of the ST distribution fitted by method of moments



cleanEx()
nameEx("zeta")
### * zeta

flush(stderr()); flush(stdout())

### Name: zeta
### Title: Function 'log(2*pnorm(x))' and its derivatives
### Aliases: zeta
### Keywords: math

### ** Examples

y <- zeta(2,seq(-20,20,by=0.5))
#
for(k in 0:5) curve(zeta(k,x), from=-1.5, to=5, col = k+2, add = k > 0)
legend(3.5, -0.5, legend=as.character(0:5), col=2:7, lty=1)



### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
