# S-plus library for the Skew-Normal (SN) distribution 
# Author: A.Azzalini <azzalini@stat.unipd.it> 
# Home-page: http://www.stat.unipd.it/~azzalini/SN
# major updates: 29/8/1997, 10/12/1997, 1/10/1998, 12/10/1998, 1/4/99
# This is version 0.22 of the library (19-Sept-2000) for R 1.0.1.

dsn <- function(x, location=0, scale=1, shape=0)
  2*dnorm((x-location)/scale)*pnorm((shape*(x-location)/scale))/scale

psn <- function(x, location=0, scale=1, shape=0) 
  pnorm((x-location)/scale) - 2*T.Owen((x-location)/scale, shape)
 
rsn <- function(n=1, location=0, scale=1, shape=0){
  u1 <- rnorm(n)
  u2 <- rnorm(n)
  id <- (u2 > shape*u1)
  u1[id] <- (-u1[id])
  return(location+scale*u1)
  }

qsn <- function(p, location=0, scale=1, shape=0, tol=1.e-8){
  na <- is.na(p) | (p<0) | (p>1)
  zero <- (p==0)
  one  <- (p==1)
  p <- replace(p,(na|zero|one),0.5)
  cum <- as.vector(sn.cumulants(shape,4))
  g1<-cum[3]/cum[2]^(3/2)
  g2<-cum[4]/cum[2]^2
  x <- qnorm(p)
  x <- x+(x^2-1)*g1/6+x*(x^2-3)*g2/24-x*(2*x^2-5)*g1^2/36
  x <- cum[1]+sqrt(cum[2])*x # a Cornish-Fisher start
  max.err <- 1
  while(max.err>tol){
    x1 <- x-(psn(x,0,1,shape)-p)/dsn(x,0,1,shape)
    max.err <- max(abs(x1-x)/(1+abs(x)))
    x <- x1
    }
  x <- replace(x,na,NA)
  x <- replace(x,zero,-Inf)
  x <- replace(x,one,Inf)
  return(location+scale*x)
}


sn.cumulants <- function(shape=0, n=4)
 {
   cumulants.half.norm <- function(n=4){
     n <- max(n,2)
     n <- as.integer(2*ceiling(n/2))
     half.n  <-  as.integer(n/2)
     m <- 0:(half.n-1)
     a <- sqrt(2/pi)/(gamma(m+1)*2^m*(2*m+1))
     signs <- rep(c(1,-1),half.n)[1:half.n]
     a <- as.vector(rbind(signs*a,rep(0,half.n)))
     coeff <- rep(a[1],n)
     for (k in 2:n) {
        ind <- 1:(k-1)
        coeff[k] <- a[k]-sum(ind*coeff[ind]*a[rev(ind)]/k)
        }
     kappa <- coeff*gamma((1:n)+1)
     kappa[2] <- 1+kappa[2]
     return(kappa)
    }
  # delta<-delta.of(shape)
  delta <- shape/sqrt(1+shape^2)
  kv <- cumulants.half.norm(n)
  if(length(kv)>n) kv<-kv[-(n+1)]
  kv[2] <- kv[2]-1
  kappa <- outer(delta,1:n,"^")*matrix(rep(kv,length(shape)),ncol=n,byrow=T)
  kappa[,2] <- kappa[,2]+1
  kappa
}

# lambda.of <- function(delta) delta/sqrt(1-delta^2)

# delta.of <- function(lambda) {
#   inf <- (abs(lambda)==Inf)
#   delta <-lambda/sqrt(1+lambda^2)
#   delta[inf] <- sign(lambda[inf])
#   delta
#}

T.Owen <- function(h, a, jmax=50, cut.point=6)
{
 T.int <-function(h,a,jmax,cut.point)
   {
     fui<- function(h,i) (h^(2*i))/((2^i)*gamma(i+1)) 
     seriesL <- seriesH <- NULL
     i  <- 0:jmax
     low<- (h<=cut.point)
     hL <- h[low]
     hH <- h[!low]
     L  <- length(hL)
     if (L>0) {
       b    <- outer(hL,i,fui)
       cumb <- apply(b,1,cumsum)
       b1   <- exp(-0.5*hL^2)*t(cumb)
       matr <- matrix(1,jmax+1,L)-t(b1)
       jk   <- rep(c(1,-1),jmax)[1:(jmax+1)]/(2*i+1)
       matr <- t(matr*jk) %*%  a^(2*i+1)
       seriesL  <- (atan(a)-as.vector(matr))/(2*pi)
     }
     if (length(hH) >0) 
       seriesH <- atan(a)*exp(-0.5*(hH^2)*a/atan(a))*
                  (1+0.00868*(hH^4)*a^4)/(2*pi)
     series <- c(seriesL,seriesH)
     id <- c((1:length(h))[low],(1:length(h))[!low]) 
     series[id] <- series  # re-sets in original order
     series
  }
  if(!is.vector(a) | length(a)>1) stop("a must be a vector of length 1")
  if(!is.vector(h)) stop("h must be a vector")
  aa <- abs(a)    
  ah <- abs(h)
  if(aa==Inf) return(0.5*pnorm(-ah))
  if(aa==0)   return(rep(0,length(h)))
  na  <- is.na(h)
  inf <- (ah==Inf)
  ah  <- replace(ah,(na|inf),0)
  if(aa<=1)
    owen <- T.int(ah,aa,jmax,cut.point)
  else
    owen<-0.5*pnorm(ah)+pnorm(aa*ah)*(0.5-pnorm(ah))- 
               T.int(aa*ah,(1/aa),jmax,cut.point)
  owen <- replace(owen,na,NA)
  owen <- replace(owen,inf,0)
  return(owen*sign(a))
}


pnorm2 <- function(x,y,rho){
  if(length(c(x,y,rho))>3) stop("non-scalar arguments")
  if(x==0 & y==0) return(0.25+asin(rho)/(2*pi))
  p <- 0.5*(pnorm(x)+pnorm(y))
  if(x==0) p <- p-0.25*sign(y) else p <- p-T.Owen(x,(y-rho*x)/(x*sqrt(1-rho^2)))
  if(y==0) p <- p-0.25*sign(x) else p <- p-T.Owen(y,(x-rho*y)/(y*sqrt(1-rho^2)))
  if((x*y<0) | ((x*y==0) & (x+y)<0)) p <- p-0.5
  return(p)
  }  

cp.to.dp <- function(param){
  # converts centred parameters cp=(mu,sigma,gamma1)
  # to direct parameters dp=(xi,omega,lambda)
  # Note:  mu can be m-dimensional, the other must be scalars
  b <- sqrt(2/pi)
  m <- length(param)-2
  gamma1 <- param[m+2]
  if(abs(gamma1)>0.9952719) stop("abs(gamma1)>0.9952719 ")
  A <- sign(gamma1)*(abs(2*gamma1/(4-pi)))^(1/3)
  delta <- A/(b*sqrt(1+A^2))
  lambda <- delta/sqrt(1-delta^2)
  E.Z <- b*delta
  sd.Z <- sqrt(1-E.Z^2)
  location    <- param[1:m]
  location[1] <- param[1]-param[m+1]*E.Z/sd.Z
  scale <- param[m+1]/sd.Z
  dp <- c(location,scale,lambda)
  names(dp)[(m+1):(m+2)]<-c("scale","shape")
  dp
  }

dp.to.cp <- function(param){
# converts 'direct' dp=(xi,omega,lambda) to 'centred' cp=(mu,sigma,gamma1)
  m<-length(param)-2
  omega<-param[m+1]
  lambda<-param[m+2]
  mu.Z <- lambda*sqrt(2/(pi*(1+lambda^2)))
  s.Z <- sqrt(1-mu.Z^2)
  gamma1<- 0.5*(4-pi)*(mu.Z/s.Z)^3
  sigma <- omega*s.Z
  mu    <- param[1:m]
  mu[1] <- param[1]+sigma*mu.Z/s.Z
  cp <- c(mu,sigma,gamma1)
  names(cp)[(m+1):(m+2)]<-c("s.d.","skewness")
  cp
}

zeta <- function(k,x){# k integer \in (0,4)
  k <- as.integer(k)
  na <- is.na(x)
  x <- replace(x,na,0)
  if(any(abs(x)==Inf)) stop("Inf not allowed")
  # funzionerebbe per k=0 e 1, ma non per k>1
  ok <- (-35<x)
  if(k==0)  
    {ax <- (-x[!ok])
    ay <- (-0.918938533204673)-0.5*ax^2-log(ax)
    y  <- rep(NA,length(x))
    y[ok] <- log(2*pnorm(x[ok]))
    y[!ok]<- ay
    }
  else {if(k==1) {y  <- (-x)*(1+1/x^2)
          y[ok]<-dnorm(x[ok])/pnorm(x[ok]) }
    else { if(k==2)  y<-(-zeta(1,x)*(x+zeta(1,x)))
      else{ if(k==3)  y<-(-zeta(2,x)*(x+zeta(1,x))-zeta(1,x)*(1+zeta(2,x)))
        else{ if(k==4)  
           y<-(-zeta(3,x)*(x+2*zeta(1,x))-2*zeta(2,x)*(1+zeta(2,x)))
        else {warning("k>4"); y<-rep(NA,length(x)) }}}}}
  replace(y,na,NA)
}


sn.em <-function(X, y, fixed, p.eps=1e-4, l.eps=1.e-2, trace=F, data=F){
#
#  1/10/1998 (elaborando dal em.lm.sn del 2-12-97)
#
#  EM per caso con uno/due/tre parametri ignoti, parametrizzando in modo 
#  "diretta" con (xi, omega, lambda); internamente usa peraltro 'delta'.
#  Le componenti ignote sono i termini NA di fixed,  ma per semplicita' 
#  assumiamo che un NA implica che le componenti alla sua sx sono NA
#  (e quindi il primo elemento di fixed e` sempre NA).
#
  n<-length(y)
  if(missing(X)) X<-matrix(rep(1,n),n,1)
  nc<-ncol(X)
  if(missing(fixed)) fixed <- rep(NA,3)
  if(all(!is.na(fixed))) stop("all parameter are fixed") 
  if(is.na(fixed[3])) iter<-(1-log(l.eps,10)) else iter<-1 
  qrX<-qr(X)
  beta<-qr.coef(qrX,y)
  xi  <- m <-qr.fitted(qrX,y)
  omega  <- fixed[2] 
  lambda <- fixed[3]
  # delta  <- delta.of(lambda)
  delta <- lambda/sqrt(1+lambda^2)
  s<-sqrt(sum((y-xi)^2)/n)
  if(is.na(fixed[3])) {
    gamma1 <- sum((y-m)^3)/(n*s^3)
    a  <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
    delta<-sqrt(pi/2)*a/sqrt(1+a^2)
    if(abs(delta)>=1) delta<-sign(delta)/(1+1/n)
    # lambda<-lambda.of(delta)
    lambda<-delta/sqrt(1-delta^2)
    }
  mean.Z <- sqrt(2/pi)*delta
  sd.Z <- sqrt(1-mean.Z^2)
  if(is.na(fixed[2])) omega  <- s/sd.Z
  if(is.na(fixed[1])) xi     <- m-s*mean.Z/sd.Z
  old.par   <- c(beta,omega,lambda)
  diverge   <- 1
  incr.logL <- Inf
  logL      <- -Inf
  while(diverge>p.eps | incr.logL>l.eps){
    # E-step
    v  <- (y-xi)/omega
    p  <- zeta(1,lambda*v)
    u1 <- omega*(delta*v+p*sqrt(1-delta^2))
    u2<-omega^2*((delta*v)^2+(1-delta^2)+p*v*delta*sqrt(1-delta^2))
    # M-step
    for(i in 1:iter){
      beta<-qr.coef(qrX,y-delta*u1)
      xi <- qr.fitted(qrX,y-delta*u1)
      d  <- y-xi
      Q  <- sum(d^2-2*delta*d*u1+u2)
      if(is.na(fixed[2])) omega <-sqrt(Q/(2*n*(1-delta^2)))
      r  <- 2*sum(d*u1)/Q
      if(is.na(fixed[3])) delta<-(sqrt((2*r)^2+1)-1)/(2*r)
      }
    # convergence?    # lambda<-lambda.of(delta)
    lambda<-delta/sqrt(1-delta^2)
    param<-c(beta,omega,lambda)
    names(param)[(nc+1):(nc+2)]<-c("scale","shape")
    diverge<-sum(abs(param-old.par)/(1+abs(old.par)))/(nc+2)
    old.par<-param
    a<-sum(log(dsn(y,xi,omega,lambda)))
    incr.logL<- a-logL
    logL<-a
    if(trace) print(c(param,logL),digits=5)
  }
  cp <- dp.to.cp(param)
  result<-list(dp=param, cp=cp, logL=logL)
  if(data) result$data <- list(X=X, y=y, residuals=d/omega)
  result
}

#-------------------

gamma1.to.lambda<- function(gamma1){
  max.gamma1 <- 0.5*(4-pi)*(2/(pi-2))^1.5
  na <- (abs(gamma1)>max.gamma1)
  if(any(na)) warning("NAs generated") 
  gamma1<-replace(gamma1,na,NA)
  a    <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  delta<- sqrt(pi/2)*a/sqrt(1+a^2)
  lambda<-delta/sqrt(1-delta^2)
  as.vector(lambda)
}

# Examples:
#  a<-sn.2logL.profile(y=otis)
#  a<-sn.2logL.profile(y=otis, use.cp=F)
#  a<-sn.2logL.profile(X=cbind(1,lbm), y=bmi, npts=50)
#  a<-sn.2logL.profile(y=frontier,param.range=c(0.8,1.6,10,30),use.cp=F,npts=11)

sn.2logL.profile<-function(X=matrix(rep(1,n)), y, 
      param.range=c(sqrt(var(y))*c(2/3, 3/2), -0.95, 0.95),
      use.cp=T, npts= 31 %/% d, plotit=T, ...)
{# plot 1D or 2D profile deviance (=-2logL) using either parameters
   # if(plotit & !exists(.Device)) stop("Device not active")
   n<-length(y)
   d<- round(length(param.range)/2)
   if((d!=1)&(d!=2)) stop(" length(param.range) must be either 2 or 4")
   if(d==1){
      param1 <- seq(param.range[1],param.range[2],length=npts)
      llik <- param2 <- rep(NA,npts)}
   else{ 
      param1 <- seq(param.range[1],param.range[2],length=npts)
      param2 <- seq(param.range[3],param.range[4],length=npts)
      llik   <- matrix(NA,npts,npts)}
   if(use.cp){
      if(d==1){
        gamma1<-param1
        sigma <-param2 
        xlab <- "gamma1"
        ylab <- ""}
      else {
        sigma <-param1
        gamma1<-param2
        xlab <- "sigma"
        ylab <- "gamma1"
        }
      if(max(abs(gamma1))>0.9952719) stop("abs(gamma1)>0.9952719")
      lambda <- gamma1.to.lambda(gamma1)
      sc<-sqrt(1-(2/pi)*lambda^2/(1+lambda^2))      
      }
   else{ # use dp 
      if(d==1) {
        lambda<-param1
        omega<-param2
        xlab <- "lambda"
        ylab <- ""}
      else {
         omega<-param1 
	 sc <- rep(1,npts)
         lambda<-param2
         xlab <- "omega"
         ylab <- "lambda"
         }
      }
   cat(c("Running until",npts,":"))
   for(i in 1:npts){
     cat(" ");cat(i)
     if(d==1) {
       a<-sn.em(X, y, fixed=c(NA,NA,lambda[i]), ...)       
       llik[i]<-a$logL
       # print(c(i,lambda[i],a$logL))
       }
     else{
     for(j in 1:npts){     
       a<-sn.em(X, y, fixed=c(NA,param1[i]/sc[j],lambda[j]), ...)
       llik[i,j] <- a$logL
       # print( c(i,j, param1[i]/sc[j], lambda[j], a$logL))
     }}
   }
  cat("\n")
  #if(plot)
  f<-2*(llik-max(llik))	
  if(plotit){
    if(d==1) plot(param1, f, type="l", 
            xlab=xlab, ylab="profile relative 2(logL)")
    else contour(param1, param2, f, labcex=0.75, 
            xlab=xlab, ylab=ylab,
            levels=-c(0.57, 1.37, 2.77, 4.6, 5.99, 9.2))
            # qchisq(c(0.25,0.5,0.75,0.90,0.95,0.99),2)
    title(main=paste("dataset:", deparse(substitute(y)),
        "\nProfile relative 2(logLikelihood)", sep= " "))	
  }
  invisible( list(param1=param1, param2=param2,
      param.names=c(xlab,ylab), two.logL=f, maximum=2*max(llik)))
}

  
sn.mle <- function(X, y, cp, plotit=T, trace=F, method="L-BFGS-B",
               control=list(iter.max=100, abs.tol=1e-5)) 
{
  xlab<-deparse(substitute(y))
  if(!is.null(dim(y))) {
    if(min(dim(y))==1) y<-as.vector(y)
    else stop("y must be a vector")
    }
  n<-length(y)
  if(missing(X)) X<-as.matrix(rep(1,n))
  m<-ncol(X)
  if(missing(cp)) {
    qrX <- qr(X)
    s <- sqrt(sum(qr.resid(qrX, y)^2)/n)
    gamma1 <- sum(qr.resid(qrX, y)^3)/(n*s^3)
    if(abs(gamma1)>0.99527) gamma1<- sign(gamma1)*0.95
    cp <- c(qr.coef(qrX,y), s, gamma1)
    }
  else{ 
    if(length(cp)!= (m+2)) stop("ncol(X)+2 != length(cp)")}
  opt<- optim(cp, fn=sn.dev, gr=sn.dev.gh, method=method,
         lower=c(-rep(Inf,m),1e-10,-0.99527), 
         upper=c(rep(Inf,m),Inf,0.99527), 
         control=control, X=X, y=y, trace=trace, hessian=F)
  cp <- opt$par
  cat(paste("Message from optimization routine:", opt$message,"\n"))
  logL <- (-opt$value)/2
  info <- attr(sn.dev.gh(cp, X, y, trace=F, hessian=T),"hessian")/2
  se <- sqrt(diag(solve(info)))
  # if(plotit & !exists(".Device")) warning("Device not active")
  # if(exists(".Device") & plotit) {
  if(plotit) {
    dp0<-cp.to.dp(cp)
    if(all(X==rep(1,n))) 
      y0<-y        
    else {
      y0<- as.vector(y - X %*% dp0[1:m])
      dp0<-c(0,dp0[m+1],dp0[m+2])
      xlab<-"residuals"
      }
    if(exists("sm.density",mode="function"))
      {
      a<-sm.density(x=y0,h=hnorm(y0)/1.5, xlab=xlab, lty=3)
      x<-a$eval.points 
      }
    else 
      {
      x<-seq(min(pretty(y0,10)),max(pretty(y0,10)),length=100)
      magic <- n^(1/3)+sqrt(n)
      hist(y0, prob=T, nclass=magic, xlim=c(min(x),max(x)), xlab=xlab, main=xlab)      
      }
    if(n<101) points(y0,rep(0,n),pch=1)
    # title(deparse(substitute(y)))
    lines(x,dsn(x,dp0[1],dp0[2],dp0[3]))
  }  
  list(call=match.call(), cp=cp, logL=logL, se=se, info=info, optim=opt)
}


sn.dev <- function(cp, X, y, trace=F)
{ # -2*logL for centred parameters  
  m <- ncol(X)
  dp <- as.vector(cp.to.dp(cp))
  location <- X %*% as.matrix(dp[1:m])
  scale <- dp[m+1]
  # AVOID: logL <- sum(log(dsn(y,location,dp[m+1],dp[m+2])))
  z <- (y-location)/scale
  nlogL <- (length(y)*log(2.506628274631*scale) + 0.5*sum(z^2)
            - sum(zeta(0,dp[m+2]*z)))
  if(trace) {cat("sn.dev: (cp,dev)="); print(c(cp,2*nlogL))}
  return(2*nlogL) 
}

sn.dev.gh <- function(cp, X, y, trace=F, hessian=F)
{
  # computes gradient and hessian of dev=-2*logL for centred parameters 
  # (and observed information matrix);
  m  <- ncol(X)
  n  <- nrow(X)
  np <- m+2
  score <- rep(NA,np)
  info  <- matrix(NA,np,np)
  beta <- cp[1:m]
  sigma <- cp[m+1]
  gamma1 <- cp[m+2]
  lambda <- gamma1.to.lambda(gamma1)
  mu <- as.vector(X %*% as.matrix(beta))
  d  <- y-mu
  r  <- d/sigma
  E.Z<- lambda*sqrt(2/(pi*(1+lambda^2)))
  s.Z<- sqrt(1-E.Z^2)
  z  <- E.Z+s.Z*r
  p1 <- as.vector(zeta(1,lambda*z))
  p2 <- as.vector(zeta(2,lambda*z))
  omega<- sigma/s.Z
  w    <- lambda*p1-E.Z
  DE.Z <- sqrt(2/pi)/(1+lambda^2)^1.5
  Ds.Z <- (-E.Z/s.Z)*DE.Z
  Dz   <- DE.Z + r*Ds.Z
  DDE.Z<- (-3)*E.Z/(1+lambda^2)^2
  DDs.Z<- -((DE.Z*s.Z-E.Z*Ds.Z)*DE.Z/s.Z^2+E.Z*DDE.Z/s.Z)
  DDz  <- DDE.Z + r*DDs.Z
  score[1:m] <- omega^(-2)*t(X) %*% as.matrix(y-mu-omega*w) 
  score[m+1] <- (-n)/sigma+s.Z*sum(d*(z-p1*lambda))/sigma^2
  score[m+2] <- n*Ds.Z/s.Z-sum(z*Dz)+sum(p1*(z+lambda*Dz))
  Dg.Dl <-1.5*(4-pi)*E.Z^2*(DE.Z*s.Z-E.Z*Ds.Z)/s.Z^4
  score[m+2] <- score[m+2]/Dg.Dl  # convert deriv wrt lamda to gamma1 
  gradient <- (-2)*score
  if(hessian){
     info[1:m,1:m] <- omega^(-2) * t(X) %*% ((1-lambda^2*p2)*X)
     info[1:m,m+1] <- info[m+1,1:m] <- 
            s.Z* t(X) %*% as.matrix((z-lambda*p1)+d*(1-lambda^2*p2)*
            s.Z/sigma)/sigma^2
     info[m+1,m+1] <- (-n)/sigma^2+2*s.Z*sum(d*(z-lambda*p1))/sigma^3 +
            s.Z^2*sum(d*(1-lambda^2*p2)*d)/sigma^4
     tmp <- as.matrix(Ds.Z*w+s.Z*(p1+lambda^2*p2*Dz-DE.Z) )
     info[1:m,m+2] <- info[m+2,1:m] <- 
            t(X)%*%(-2*Ds.Z*d/omega+Ds.Z*w+s.Z*(p1+lambda*p2*(z+lambda*Dz)
            -DE.Z))/sigma 
     info[m+1,m+2] <- info[m+2,m+1] <- 
            -sum(d*(Ds.Z*(z-lambda*p1)+s.Z*(Dz-p1-p2*lambda*(z+lambda*Dz))
             ))/sigma^2
     info[m+2,m+2] <- n*(-DDs.Z*s.Z+Ds.Z^2)/s.Z^2+sum(Dz^2+z*DDz)-
            sum(p2*(z+lambda*Dz)^2)- sum(p1*(2*Dz+lambda*DDz))
     info[np,] <- info[np,]/Dg.Dl # convert info wrt lamda to gamma1 
     info[,np] <- info[,np]/Dg.Dl
     }
  attr(gradient,"hessian") <- 2*info
  if(trace) {cat("sn.dev.gh: gradient="); print(-2*score)}
  return(gradient)
}

# version 0.2, Oct 1998

dmsn <- function(x, xi=rep(0,k), Omega, alpha)
{# Density of Multivariate SN rv with parameters (xi, Omega, alpha) 
 # evaluated at x, which is either a k-vector or n x k matrix
  scale <- sqrt(diag(Omega))
  if(is.vector(x)) {n <-1;         k <- length(x)} 
        else       {n <-dim(x)[1]; k <- dim(x)[2]}
  X     <- t(matrix(x,nrow=n,ncol=k))-xi
  z     <- X/scale
  Q     <- apply((solve(Omega)%*% X)* X,2,sum) #diagonal of (x Omega^(-1) x^T)
  # Det   <- as.numeric(det.Hermitian(as.Matrix(Omega),logarithm=F)$modulus)
  d <- diag(qr(Omega)[[1]])
  Det <- prod(abs(d))
  pdf   <- 2*exp(-0.5*Q)*pnorm(t(z)%*%as.matrix(alpha))/sqrt((2*pi)^k*Det)
  pdf
}


rmsn <- function(n=1, xi=rep(0,k), Omega, alpha){
# generates SN_k(xi,Omega,alpha) variates using transformation method
  k <- ncol(Omega)
  Z <- msn.quantities(xi,Omega,alpha)
  y <- matrix(rnorm(n*k),n,k) %*% chol(Z$Psi)
  # each row of y is N_k(0,Psi)
  abs.y0 <- abs(rnorm(n))  
  abs.y0<-matrix(rep(abs.y0,k),ncol=k)
  delta <- Z$delta
  z <- delta * t(abs.y0) +  sqrt(1-delta^2) * t(y)
  y <- t(xi+Z$omega*z)
  return(y)
  }


plot.dsn2 <- function(x, y, xi, Omega, alpha, ...)
{# plot bivariate density SN_2(xi,Omega,alpha) computed at (x,y) grid
  if(any(dim(Omega)!=c(2,2))) stop("dim(Omega) != c(2,2)")
  nx <- length(x)
  ny <- length(y)
  xoy <- cbind(rep(x,ny), as.vector(matrix(y,nx,ny,byrow=T)))
  X <- matrix(xoy, nx*ny, 2, byrow=F)
  pdf<-dmsn(X, xi, Omega, alpha)
  pdf<-matrix(pdf, nx, ny)
  contour(x, y, pdf, ...)
  invisible(list(x=x,y=x,density=pdf,xi=xi,Omega=Omega,alpha=alpha))
}

msn.quantities <- function(xi,Omega,alpha){
# 21-12/1997; computes various quantities related to SN_k(xi,Omega,alpha)
  Diag <- function(x) diag(x,nrow=length(x),ncol=length(x))
  k <- length(alpha)
  if(length(xi)!=k | any(dim(Omega)!=c(k,k))) 
       stop("dimensions of arguments do not match")
  omega <- sqrt(diag(Omega))
  O.cor <- Diag(1/omega) %*% Omega %*% Diag(1/omega)
  tmp <- as.vector(sqrt(1 + t(as.matrix(alpha))%*%O.cor%*%alpha)) 
  delta<- as.vector(O.cor %*%alpha)/tmp
  lambda<- delta/sqrt(1-delta^2)
  D <-diag(sqrt(1+lambda^2))
  Psi <- D %*% (O.cor-outer(delta,delta)) %*% D
  Psi <- (Psi+t(Psi))/2
  O.inv <- solve(Omega)
  oi<- sqrt(diag(O.inv))
  O.pcor <- Diag(1/oi)%*% (-O.inv) %*% Diag(1/oi)
  diag(O.pcor) <- rep(1,k)
  muZ <- delta*sqrt(2/pi)
  muY <- xi+omega*muZ
  Sigma <- Diag(omega) %*% (O.cor-outer(muZ,muZ)) %*% Diag(omega) 
  Sigma <- (Sigma+t(Sigma))/2
  cv <- muZ/sqrt(1-muZ^2)
  gamma1 <- 0.5*(4-pi)*cv^3
  list(xi=xi, Omega=Omega, alpha=alpha, omega=omega,  mean=muY, variance=Sigma,
       Omega.conc=O.inv, Omega.cor=O.cor, Omega.pcor=O.pcor,
       lambda=lambda, Psi=Psi,  delta=delta, skewness=gamma1)
}

msn.conditional <- function(xi, Omega, alpha, fixed.comp, fixed.values)
{
# conditional Multivariate SN  (6/11/1997).
# Given a rv Y ~ SN_k(xi,Omega,alpha), this function computes cumulants of
# conditrional distribution, given that the fixed.com take on fixed.values;
# then it finds MSN with matching cumlants.  
  Diag <- function(x) diag(x,nrow=length(x),ncol=length(x))
  msqrt <- function(A) Diag(sqrt(diag(as.matrix(A))))
  imsqrt<- function(A) Diag(1/sqrt(diag(as.matrix(A))))
  k<-length(alpha) 
  h<-length(fixed.comp)
  if(any(dim(Omega)!=c(k,k)) | length(xi)!=k | h!=length(fixed.values))
        stop("dimensions of parameters do not match")
  fc <- fixed.comp
  O  <- as.matrix(Omega)
  O11<- O[fc,fc, drop=F]
  O12<- O[fc,-fc, drop=F]
  O21<- O[-fc,fc, drop=F]
  O22<- O[-fc,-fc, drop=F]
  o22<- sqrt(diag(O22))
  inv.O11 <- solve(O11)
  xi1 <- xi[fc, drop=F]
  xi2 <- xi[-fc, drop=F]
  alpha1 <- matrix(alpha[fc])
  alpha2 <- matrix(alpha[-fc]) 
  O22.1 <- O22 - O21 %*% inv.O11 %*% O12
  O22.b <- imsqrt(O22) %*% O22.1 %*% imsqrt(O22)
  xi.c  <- xi2 + as.vector(O21 %*% inv.O11 %*% matrix(fixed.values-xi1))
  a     <- sqrt(1+as.vector(t(alpha2) %*% O22.b %*% alpha2))
  alpha.b <- (alpha1 + msqrt(O11) %*% inv.O11 %*% O12 %*% (alpha2/o22))/a  
  d2    <- as.vector(O22.b %*% alpha2)/a
  x0    <- sum(alpha.b * (fixed.values-xi1)/sqrt(diag(O11)))
  E.c   <- xi.c + zeta(1,x0)*o22*d2
  var.c <- O22.1+zeta(2,x0)*outer(o22*d2,o22*d2)
  gamma1<- zeta(3,x0)*d2^3/diag(O22.b+zeta(2,x0)*d2^2)^1.5
  cum   <- list(as.vector(E.c),var.c,gamma1)
  # cumulants are computed; now choose SN distn to fit them
  a     <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  E.z   <- a/sqrt(1+a^2)
  delta <- E.z*sqrt(pi/2) 
  omega <- sqrt(diag(var.c)/(1-E.z^2))
  O.new <- var.c+outer(omega*E.z,omega*E.z) 
  xi.new<- E.c-omega*E.z
  B   <- Diag(1/omega)
  m   <- as.vector(solve(B %*% O.new %*% B) %*% as.matrix(delta))
  a   <- m/sqrt(1-sum(delta*m))
  # cum2<- msn.cumulants(xi.new,O.new,a)
  list(cumulants=cum, fit=list(xi=xi.new, Omega=O.new, alpha=a, delta=delta))
}


msn.marginal <- function(xi,Omega,alpha,comp)
{# calcola parametri della marginale associata a comp di un SN_k 
  Diag <- function(x) diag(x, nrow=length(x), ncol=length(x))
  xi <- as.vector(xi)
  comp <- as.integer(comp)
  alpha <- as.vector(alpha)
  k <- length(alpha)
  if(length(comp)<k){
    if(any(comp>k | comp<1)) stop("comp makes no sense")
    scale<- sqrt(diag(Omega))
    O   <- Diag(1/scale) %*% Omega %*% Diag(1/scale)
    O11 <- O[comp,comp, drop=F]
    O12 <- O[comp,-comp, drop=F]
    O21 <- O[-comp,comp, drop=F]
    O22 <- O[-comp,-comp, drop=F]
    alpha1<- as.matrix(alpha[comp, drop=F])
    alpha2<- as.matrix(alpha[-comp, drop=F])
    O22.1 <- O22 - O21 %*% solve(O11) %*% O12
    a.sum <- as.vector(t(alpha2) %*% O22.1 %*% alpha2)
    a.new <- as.vector(alpha1+solve(O11) %*% O12 %*% alpha2)/sqrt(1+a.sum)
    O.new <- Diag(scale[comp]) %*% O11 %*% Diag(scale[comp])
    result<- list(xi=xi[comp], Omega=O.new, alpha=a.new)
  }
  else {
   if(any(sort(comp)!=(1:k))) stop("comp makes no sense")
   result <-
    list(xi=xi[comp], Omega=as.matrix(Omega[comp,comp, drop=F]), alpha=alpha[comp]) 
   }
  result
}

plot.msn.cond <- function(xi, Omega, alpha, fixed.comp, fixed.values, n=35)
{# fa il grafico di Y_2|Y_1; assumiamo che dim(Y_2)= 2 
  msn.pdf2.aux <- function(x,y,xi,Omega,alpha,fc,fv)
   {
     nx <- length(x)
     ny <- length(y)
     FV <- matrix(rep(fv,nx*ny), nx*ny, length(fv), byrow=T)
     X <- matrix(NA, nx*ny, length(alpha))
     X[,fc] <- FV
     xoy <- cbind(rep(x,ny), as.vector(matrix(y,nx,ny,byrow=T)))
     X[,-fc] <- matrix(xoy, nx*ny, 2, byrow=F)
     pdf<-dmsn(X,xi,Omega,alpha)
     matrix(pdf,nx,ny)
   }
  dsn2 <- function(x,y,d1,d2,omega)
    {
     u <- (x*(d1-omega*d2)+y*(d2-omega*d1))/
          sqrt((1-omega^2-d1^2-d2^2+2*omega*d1*d2)*(1-omega^2))
     pdfn2 <- exp(-0.5*(x^2-2*omega*x*y+y^2)/(1-omega^2))/
              (2*pi*sqrt(1-omega^2))
     2*pdfn2*pnorm(u)
    }
  fc <- fixed.comp
  fv <- fixed.values
  cond <- msn.conditional(xi,Omega,alpha,fc,fv)
  xi.c <- cond$fit$xi
  O.c  <- cond$fit$Omega
  a.c  <- cond$fit$alpha
  if(any(dim(O.c)!=c(2,2))) stop("length(alpha)-length(fixed.com)!=2")
  scale1<-sqrt(as.vector(O.c[1,1]))
  scale2<-sqrt(as.vector(O.c[2,2]))
  delta <- cond$fit$delta
  omega <-as.vector(O.c[1,2])/(scale1*scale2)
  x<-seq(xi.c[1]-3*scale1, xi.c[1]+3*scale1, length=n)
  y<-seq(xi.c[2]-3*scale2, xi.c[2]+3*scale2, length=n)
  plot(x,y,type="n", main="Conditional multivariate SN pdf")
  z1<-(x-xi.c[1])/scale1
  z2<-(y-xi.c[2])/scale2
  pdf.fit<-outer(z1,z2,dsn2,d1=delta[1],d2=delta[2],omega=omega)/
                (scale1*scale2)
  cond$pdf<-list(x=x,y=y,f.fitted=pdf.fit)
  levels<-pretty(pdf.fit,5)
  contour(x,y,pdf.fit,levels=levels,add=T,col=2)
  # fino a qui per il calcolo della densità approx;
  # ora otteniamo quella esatta
  numer <- msn.pdf2.aux(x, y, xi, Omega, alpha, fc, fv)
  marg  <- msn.marginal(xi, Omega, alpha, fc)
  denom <- dmsn(fv, marg$xi, marg$Omega, marg$alpha)
  pdf.exact<- numer/as.vector(denom)
  contour(x, y, pdf.exact, add=T, levels=levels, col=1, lty=4, labcex=0)
  legend(x[1],y[length(y)],c("approx","exact"), lty=c(1,4),col=c(2,1))
  cond$pdf$f.exact<-pdf.exact
  cond$rel.error<-summary((pdf.fit-pdf.exact)/pdf.exact)
  cond$abs.error<-summary(abs(pdf.fit-pdf.exact))
  invisible(cond)
}


msn.moment.fit <- function(y){
# 31-12-1997: simple fit of MSN distribution usign moments
#  cat("MLE con metodo momenti\n")
  Diag <- function(x) diag(x,nrow=length(x),ncol=length(x))
  y     <- as.matrix(y)
  k     <- ncol(y)
  m.y   <- apply(y,2,mean)
  var.y <- var(y)
  y0    <- (t(y)-m.y)/sqrt(diag(var.y))
  gamma1<- apply(y0^3,1,mean)
  out   <- (abs(gamma1)>0.99527)
  gamma1[out] <- sign(gamma1[out])*0.995
  a     <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  delta <- sqrt(pi/2)*a/sqrt(1+a^2)
  m.z   <- delta*sqrt(2/pi) 
  omega <- sqrt(diag(var.y)/(1-m.z^2))
  Omega <- var.y+outer(omega*m.z,omega*m.z) 
  xi    <- m.y-omega*m.z
  O.cor <- Diag(1/omega) %*% Omega %*% Diag(1/omega)
  O.cor <- (t(O.cor)+O.cor)/2
  O.inv <- solve(O.cor)
  tmp   <- as.vector(1-t(delta) %*% O.inv %*% delta)
  if(tmp<=0) {tmp <- 0.0001; admissible <- F} else admissible<-T
  alpha <-as.vector(O.inv%*%delta)/sqrt(tmp)
  list(xi=xi, Omega=Omega, alpha=alpha, Omega.cor=O.cor, omega=omega, 
       delta=delta, skewness=gamma1, admissible=admissible) 
}

msn.fit <- function(X, y, freq, plotit=T, trace=F, ... )
{
  y.name <- deparse(substitute(y))
  y.names<- dimnames(y)[[2]] 
  y <- as.matrix(y)
  colnames(y)<-y.names
  k <- ncol(y)
  if(missing(freq)) freq<-rep(1,nrow(y))  
  n <- sum(freq)
  if(missing(X)) {X<-rep(1,nrow(y)); missing.X<-TRUE}
  X <- as.matrix(X)
  m <- ncol(X)
  if(length(dimnames(y)[[2]])==0) 
      dimnames(y) <- list(NULL, outer("V",as.character(1:k),paste,sep=""))
  qrX <- qr(X)
  mle<- msn.mle(X=X, y=y, freq=freq, trace=trace, ...)
  mle$call <- match.call()
  # print(mle$nlminb$message)
  beta  <- mle$dp$beta
  Omega <- mle$dp$Omega
  alpha <- mle$dp$alpha
  omega <- sqrt(diag(Omega))
  xi    <- X %*% beta
  if(plotit & all(freq==rep(1,length(y)))) {
    if(missing.X) { 
      y0  <-y 
      xi0 <- apply(xi,2,mean)} 
    else  {
      y0  <- y-xi 
      xi0 <- rep(0,k)
      }
    if(k>1) {
       pairs(y0, labels=y.names,
        panel=function(x, y, Y, y.names, xi, Omega, alpha){
          for(i in 1:length(alpha)){
            # if(y.names[i]==deparse(substitute(x))) Ix<-i
            # if(y.names[i]==deparse(substitute(y))) Iy<-i
            if(all(Y[,i]==x)) Ix<-i
            if(all(Y[,i]==y)) Iy<-i
            }
          points(x,y)
          marg<-msn.marginal(xi,Omega,alpha,c(Ix,Iy))
          xi.marg<-marg$xi
          Omega.marg<-marg$Omega
          alpha.marg<-marg$alpha     
          x1 <- seq(min(x), max(x), length=30)
          x2 <- seq(min(y), max(y), length=30)
          plot.dsn2(x1, x2, xi.marg, Omega.marg, alpha.marg, add=T, col=2)
        },  # end "panel" function
      Y=y0, y.names=y.names, xi=xi0, Omega=Omega, alpha=alpha)}
    else{ # plot for case k=1
      y0<-as.vector(y0)
      x<-seq(min(pretty(y0,10)),max(pretty(y0,10)),length=100)
      if(missing.X){
         dp0<-c(xi0,omega,alpha)
         xlab<-y.name}
      else {
         dp0<-c(0,omega,alpha)
         xlab <- "residuals"}
      hist(y0, prob=T, nclass=n^(1/3)+sqrt(n), xlab=xlab, ylab="density")
      lines(x, dsn(x,dp0[1],dp0[2],dp0[3]))
      if(length(y)<101) points(y0, rep(0,n), pch=1)
      title(y.name)
      }
    cat("Press <Enter> to continue..."); readline()
    par(mfrow=c(1,2))
    pp <- qchisq((1:n)/(n+1),k)
    Xb <- qr.fitted(qrX,y)
    res<- qr.resid(qrX,y)
    rad.n  <- apply((y-Xb) * ((y-Xb) %*% solve(var(res))), 1, sum)
    rad.sn <- apply((y-xi) * ((y-xi) %*% solve(Omega)), 1, sum)
    plot(pp, sort(rad.n), pch=1, ylim=c(0,max(rad.n,rad.sn)),
        xlab="Percentiles of chi-square distribution", 
        ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for normal distribution", sub=y.name)
    plot(pp, sort(rad.sn), pch=1, ylim=c(0,max(rad.n,rad.sn)),
        xlab="Percentiles of chi-square distribution", 
        ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for skew-normal distribution", sub=y.name)
    cat("Press <Enter> to continue..."); readline()
    plot((1:n)/(n+1), sort(pchisq(rad.n,k)), xlab="", ylab="")
    abline(0,1,lty=3)
    title(main="PP-plot for normal distribution", sub=y.name)
    plot((1:n)/(n+1), sort(pchisq(rad.sn,k)), xlab="", ylab="")
    abline(0,1,lty=3)
    title(main="PP-plot for skew-normal distribution", sub=y.name)
    par(mfrow=c(1,1))
    } # end ploting
  dev.norm<- msn.dev(c(qr.coef(qrX,y),rep(0,k)), X, y, freq)
  test <- dev.norm + 2*mle$logL
  p.value <-  1-pchisq(test,k)
  if(trace) {
    cat("LRT for normality (test-function, p-value): ")
    print(c(test,p.value))
    }
  mle$test.normality <- list(LRT=test, p.value=p.value)
  invisible(mle)
}


msn.mle <-function(X, y, freq, start, trace=F, method="L-BFGS-B",
                control=list(iter.max=150, x.tol=1e-8) )
{
  y <- as.matrix(y)
  if(missing(X)) X <- rep(1,nrow(y))
  if(missing(freq)) freq <- rep(1,nrow(y))
  X <- as.matrix(X) 
  k <- ncol(y)  
  n <- sum(freq)
  m <- ncol(X)
  qrX <- qr(X)
  y.names<-dimnames(y)[[2]] 
  x.names<-dimnames(X)[[2]]
  if(missing(start)) {
      beta  <- qr.coef(qrX, y)
      res   <- qr.resid(qrX, y)
      a     <- msn.moment.fit(res)
      Omega <- a$Omega
      omega <- a$omega
      alpha <- a$alpha
      if(!a$admissible) alpha<-alpha/(1+max(abs(alpha)))
      beta[1,] <- beta[1,]-omega*a$delta*sqrt(2/pi)  
     }
  else{
    beta  <- start$beta
    Omega <- start$Omega
    alpha <- start$alpha
    omega <- sqrt(diag(Omega)) 
    }
  al.om <-alpha/omega
  if(trace){ 
    cat("Initial parameters:\n")
    print(cbind(t(beta),al.om,Omega))
    }
  param<- c(beta,al.om)
  dev <- msn.dev(param,X,y,freq) 
  opt <- optim(param, fn=msn.dev, gr=msn.dev.grad, method=method,
        control=control, X=X, y=y, freq=freq, trace=trace)      
  # if(trace) 
  cat(paste("Message from optimization routine:", opt$message,"\n"))
  logL <- (-opt$value)/2
  beta <- matrix(opt$par[1:(m*k)],m,k)
  dimnames(beta)[2] <- list(y.names)
  dimnames(beta)[1] <- list(x.names)
  al.om  <- opt$par[(m*k+1):(m*k+k)]
  xi    <- X %*% beta
  Omega <- t(y-xi) %*% (freq*(y-xi))/n
  omega <- sqrt(diag(Omega))
  alpha <- al.om*omega
  param <- cbind(omega,alpha)
  dimnames(Omega) <- list(y.names,y.names)
  dimnames(param)[1] <- list(y.names)
  info <- num.deriv(opt$par, FUN="msn.dev.grad", X=X, y=y, freq=freq)/2
  se      <- sqrt(diag(solve(info)))
  se.beta <- matrix(se[1:(m*k)],m,k)
  se.alpha<- se[(m*k+1):(m*k+k)]*omega
  dimnames(se.beta)[2]<-list(y.names)
  dimnames(se.beta)[1]<-list(x.names)
  se  <- list(beta=se.beta, alpha=se.alpha, info=info)
  dp  <- list(beta=beta, Omega=Omega, alpha=alpha)
  list(call=match.call(), dp=dp, logL=logL, se=se, optim=opt)
}
 
msn.dev<-function(param, X, y, freq, trace=F){
  k <- ncol(y)
  # if(missing(freq)) freq<-rep(1,nrow(y))
  n <- sum(freq)
  m <- ncol(X)
  beta<-matrix(param[1:(m*k)],m,k)
  al.om<-param[(m*k+1):(m*k+k)]
  y0 <- y-X %*% beta
  Omega <- (t(y0) %*% (y0*freq))/n  
  d <- diag(qr(2*pi*Omega)[[1]])
  logDet <- sum(log(abs(d)))
  dev <- n*logDet-2*sum(zeta(0,y0 %*% al.om)*freq)+n*k
  if(trace) { 
    cat("\nmsn.dev:",dev,"\n","parameters:"); 
    print(rbind(beta,al.om))
    }
  dev
}

msn.dev.grad <- function(param, X, y, freq, trace=F){
  k <- ncol(y)
  # if(missing(freq)) freq<-rep(1,nrow(y))
  n <- sum(freq)
  m <- ncol(X)
  beta<-matrix(param[1:(m*k)],m,k)
  al.om<-param[(m*k+1):(m*k+k)]
  y0 <- y-X %*% beta
  Omega <- (t(y0) %*% (freq*y0))/n
  p1 <- zeta(1,as.vector(y0 %*% al.om))
  Dbeta <- t(X)%*% (y0*freq) %*%solve(Omega) - 
            outer(as.vector(t(X*freq)%*%p1), al.om)
  Dal.om <- as.vector(t(y0*freq) %*% p1)
  if(trace){
    cat("gradient:\n")
    print(rbind(Dbeta,Dal.om))}
  -2*c(Dbeta,Dal.om)
}

num.deriv <- function(coefficients, FUN, ...)
{# da rm.fit: derivate seconde numeriche, se FUN da` il gradiente
   FUN <- get(FUN, inherit = T)
   values <- FUN(coefficients, ...)
   p <- length(values)
   H <- matrix(0, p, p)
   h <- rep(0, p)
   delta <- cbind((abs(coefficients) + 1e-10) * 1e-5, rep(1e-06, p))
   delta <- apply(delta, 1, max)
   for(i in 1:p) {
   	h[i] <- delta[i]
   	new.values <- FUN(coefficients + h, ...)
   	H[, i] <- (new.values - values)/delta[i]
   	h[i] <- 0
   }
   (H+t(H))/2
}

.First.lib <- function(library, pkg)
{
    if(version$major == 0 |(version$major == 1 && version$minor < 0.1))
        stop("This package requires R 1.0.1 or later")
    assign(".sm.home", file.path(library, pkg),
           pos=match("package:sn", search()))
    cat("Library `sn'; Copyright (C) 1998-2000 A.Azzalini\n")
    cat("type `help(SN,package=sn)' for summary information\n")
    invisible()
}

