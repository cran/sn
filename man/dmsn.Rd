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
dmsn(x, xi=rep(0,d), Omega, alpha)
pmsn(x, xi=rep(0,d), Omega, alpha, ...)
rmsn(n=1, xi=rep(0,d), Omega, alpha)
}
\arguments{
\item{x}{
for \code{dmsn}, this is either a vector of length \code{d} or a matrix 
with \code{d} columns (where \code{d=length(alpha)}), giving 
the coordinates of the point(s) where the density must be avaluated;
for \code{pmsn}, only a vector of length \code{d} is allowed.
}
\item{Omega}{
a positive-definite covariance matrix of dimension \code{(d,d)}.
}
\item{alpha}{
a numeric vector which regulates the shape of the density.
}
\item{xi}{
a numeric vector of lenght \code{d}, or a matrix with \code{d} columns,
representing the location parameter of the distribution.
If \code{xi} is a matrix, its dimensions must agree with those of \code{x}.
}
\item{n}{
a numeric value which represents the number of random vectors
to be drawn.
}
\item{...}{
additional parameters passed to \code{pmvnorm}
}}
\value{
A vector of density values (\code{dmsn}), or a single probability 
(\code{pmsn}) or a matrix of random  points (\code{rmsn}).
}
\section{Background}{
The multivariate skew-normal distribution is discussed by
Azzalini and Dalla Valle (1996); the \code{(Omega,alpha)} parametrization
adopted here is the one of Azzalini and Capitanio (1999).
}
\details{
  The positive-definiteness of \code{Omega} is not tested for
  efficiency reasons. Function
  \code{pmsn} requires  \code{pmvnorm} from \code{library(mvtnorm)};
  the accuracy of its computation can be controlled via use of \code{...}
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
\code{\link{dsn}}, \code{\link{msn.fit}}, \code{\link{dmst}},
\code{\link[mvtnorm]{pmvnorm}}
}
\examples{
x <- seq(-3,3,length=15)
xi <- c(0.5, -1)
Omega <- diag(2)
Omega[2,1] <- Omega[1,2] <- 0.5
alpha <- c(2,2)
pdf <- dmsn(cbind(x,2*x-1), xi, Omega, alpha)
rnd <- rmsn(10,  xi, Omega, alpha)
library(mvtnorm)                        # only once in the session
cdf <- pmsn(c(2,1), xi, Omega, alpha)
}
\keyword{distribution}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
