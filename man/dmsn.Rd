\name{dmsn}
\alias{dmsn}
\alias{rmsn}
\title{
Multivariate skew-normal distribution
}
\description{
Probability density function and random number generation
for the multivariate skew-normal (MSN) distribution.
}
\usage{
dmsn(x, xi=rep(0,k), Omega, alpha)
rmsn(n=1, xi=rep(0,k), Omega, alpha)
}
\arguments{
\item{x}{
either a vector of length \code{k} or a matrix with \code{k} columns, where
\code{k} is \code{length(alpha)}, giving the coordinates of the point(s)
where the density must be avaluated.
}
\item{Omega}{
a covariance matrix of dimension \code{(k,k)}.
}
\item{alpha}{
a numeric vector which regulates the shape of the density.
}
\item{xi}{
a numeric vector of lenght \code{k}, or a matrix with \code{k} columns,
representing the location parameter of the distribution.
If \code{xi} is a matrix, its dimensions must agree with those of \code{x}.
}
\item{n}{
a numeric value which represents the number of random vectors
to be drawn.
}}
\value{
A vector of density values (\code{dmsn}), or a matrix of random 
points (\code{rmsn}).
}
\section{BACKGROUND}{
The multivariate skew-normal distribution is discussed by
Azzalini and Dalla Valle (1996); the \code{(Omega,alpha)} parametrization
adopted here is the one of Azzalini and Capitanio (1999).
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
\code{\link{dsn}}, \code{\link{msn.fit}}, \code{\link{msn.quantities}} }
}
\examples{
x <- seq(-3,3,length=30)
pdf <- dmsn(cbind(x,0), c(0,0), diag(2), c(2,3))
#
rnd <- rmsn(50,  c(0,0), diag(2), c(2,3))
}
\keyword{distribution}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
