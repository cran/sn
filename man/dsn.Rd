\name{dsn}
\alias{dsn}
\alias{psn}
\alias{qsn}
\alias{rsn}
\title{
Skew-Normal Distribution
}
\description{
Density function, distribution function, quantiles and random number
generation for the skew-normal (SN) distribution.
}
\synopsis{
dsn(x, location = 0, scale = 1, shape = 0, dp = NULL, log = FALSE)
psn(x, location = 0, scale = 1, shape = 0, dp = NULL, engine, ...)
qsn(p, location = 0, scale = 1, shape = 0, dp = NULL, tol = 1e-8, engine, ...) 
rsn(n = 1, location = 0, scale = 1, shape = 0, dp = NULL)
}
\usage{
dsn(x, location=0, scale=1, shape=0, log=FALSE)
dsn(x, dp=, log=FALSE)
psn(x, location=0, scale=1, shape=0, engine, ...)
psn(x, dp=, engine, ...)
qsn(p, location=0, scale=1, shape=0, tol=1e-8, ...)
qsn(x, dp=, ...)
rsn(n=1, location=0, scale=1, shape=0)
rsn(x, dp=)
}
\arguments{
\item{x}{
  vector of quantiles. Missing values (\code{NA}s) and \code{Inf}'s
  are allowed.
}
\item{p}{
  vector of probabilities. Missing values (\code{NA}s) are allowed.
}
\item{location}{
  vector of location parameters.
}
\item{scale}{
  vector of (positive) scale parameters.
}
\item{shape}{
  vector of shape parameters.
  With \code{psn} and \code{qsn}, it must be of length 1.
}
\item{dp}{
a vector of length 3, whose elements represent location, scale (positive) 
and shape, respectively. If \code{dp} is specified, the individual
parameters cannot be set.
}
\item{n}{
  sample size.
}
\item{tol}{
  a scalar value which regulates the accuracy of the result of \code{qsn}.
}
\item{log}{
  logical flag used in \code{dsn} (default \code{FALSE}).
  When \code{TRUE}, the logarithm of the density values is returned.
}
\item{engine}{
  character string to select the computing engine, which is either
  \code{"T.Owen"} or \code{"biv.nt.prob"} (the latter from package 
  \code{mnormt}).
  If the parameter is missing, a default selection rule is applied.
}
\item{...}{
  additional parameters passed to \code{T.Owen}
}}
\value{
density (\code{dsn}), probability (\code{psn}),
quantile (\code{qsn}) or  random sample (\code{rsn})
from the skew-normal distribution with given \code{location}, \code{scale}
and \code{shape} parameters.
}
\section{Background}{
The family of skew-normal distributions is an extension of the normal
family, via the introdution of a \code{shape} parameter which regulates
skewness; when \code{shape=0}, the skew-normal distribution reduces to the
normal one.  The density of the SN distribution in the "standard" case
having \code{location=0} and \code{scale=1}  is
\code{2*dnorm(x)*pnorm(shape*x)}.
A multivariate version of the distribution exists.
See the reference below for additional information.
}
\section{Details}{
\code{psn} and \code{qsn} make use either of function \code{\link{T.Owen}}
or \code{\link[mnormt:dmt]{biv.nt.prob}}
}
\references{
Azzalini, A. (1985).
A class of distributions which includes the normal ones.
\emph{Scand. J. Statist.}
\bold{12}, 171-178.
}
\seealso{
 \code{\link{dmsn}}, \code{\link{dst}},  \code{\link{T.Owen}},
 \code{\link[mnormt:dmt]{biv.nt.prob}}
}
\examples{
pdf <- dsn(seq(-3,3,by=0.1), shape=3)
cdf <- psn(seq(-3,3,by=0.1), shape=3)
qu <- qsn(seq(0.1,0.9,by=0.1), shape=-2)
rn <- rsn(100, 5, 2, 5)
}
\keyword{distribution}
