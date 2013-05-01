\name{dmsn}
\alias{dmsn}
\alias{pmsn}
\alias{rmsn}
\title{
Multivariate skew-normal distribution
}
\description{
Probability density function, distribution function and random number 
generation for the multivariate skew-normal (MSN) distribution.
}

\usage{
dmsn(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, log=FALSE)
pmsn(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, ...)
rmsn(n=1, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL)
}
\arguments{
\item{x}{
for \code{dmsn}, this is either a vector of length \code{d}, 
where \code{d=length(alpha)}, or a matrix with \code{d} columns, 
giving the coordinates of the point(s) where the density must 
be evaluated;
for \code{pmsn}, only a vector of length \code{d} is allowed.
}
\item{xi}{
a numeric vector of length \code{d}
representing the location parameter of the distribution; 
see Background.
In a call to \code{dmsn}, \code{xi} can be a matrix; 
in this case, its dimensions must agree with those of \code{x}.
}
\item{Omega}{
a symmetric positive-definite matrix of dimension \code{(d,d)};
see Background.
}
\item{alpha}{
a numeric vector which regulates the shape of the density; see Background.
}
\item{dp}{
a list with three elements named \code{xi}, \code{Omega} and \code{alpha}
containing quantities as described above. If \code{dp} is specified, 
individual parameters must not be specified. 
}
\item{n}{
a numeric value which represents the number of random vectors
to be drawn.
}
\item{log}{ 
logical (default value: \code{FALSE}); if TRUE, log-densities are returned.
}
\item{...}{
additional parameters passed to \code{pmnorm}
}}
\value{
A vector of density values (\code{dmsn}), or a single probability 
(\code{pmsn}) or a matrix of random  points (\code{rmsn}).
}
\details{Typical usages are
\preformatted{%
dmsn(x, xi=rep(0,length(alpha)), Omega, alpha, log=FALSE)
dmsn(x, dp=, log=FALSE)
pmsn(x, xi=rep(0,length(alpha)), Omega, alpha, ...)
pmsn(x, dp=)
rmsn(n=1, xi=rep(0,length(alpha)), Omega, alpha)
rmsn(n=1, dp=)
}
  The positive-definiteness of \code{Omega} is not tested for
  efficiency reasons. Function \code{pmsn} requires \code{pmnorm}
  from package \code{mnormt};
  the accuracy of its computation can be controlled via use of \code{...}
}
\section{Background}{
The multivariate skew-normal distribution is discussed by
Azzalini and Dalla Valle (1996); the \code{(Omega,alpha)}
parametrization adopted here is the one of Azzalini and Capitanio (1999).
 Notice that the location vector \code{xi}
does not represent the mean vector of the distribution  and similarly
\code{Omega} is not \emph{the} covariance matrix of the distribution,
although it is \emph{a} covariance matrix. 
}

\references{
Azzalini, A. and Dalla Valle, A. (1996).
The multivariate skew-normal distribution.
\emph{Biometrika}
\bold{83}, 715--726.


Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{dsn}},  \code{\link{dmst}}, \code{\link[mnormt]{dmnorm}}
}
\examples{
x <- seq(-3,3,length=15)
xi <- c(0.5, -1)
Omega <- diag(2)
Omega[2,1] <- Omega[1,2] <- 0.5
alpha <- c(2,-6)
pdf <- dmsn(cbind(x,2*x-1), xi, Omega, alpha)
rnd <- rmsn(10,  xi, Omega, alpha)
p1 <- pmsn(c(2,1), xi, Omega, alpha)
p2 <- pmsn(c(2,1), xi, Omega, alpha, abseps=1e-12, maxpts=10000)
}
\keyword{distribution}
\keyword{multivariate}
