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
\usage{
dsn(x, location=0, scale=1, shape=0)
psn(q, location=0, scale=1, shape=0)
qsn(p, location=0, scale=1, shape=0, tol=1e-8)
rsn(n=1, location=0, scale=1, shape=0)
}
\arguments{
\item{x}{
vector of quantiles. Missing values (\code{NA}s) are allowed.
}
\item{q}{
vector of quantiles. Missing values (\code{NA}s) are allowed.
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
vector of shape parameters. With \code{psn} and `qsn", it must be of
length 1.
}
\item{n}{
sample size.
}
\item{tol}{
a scal value which regulates the accuracy of the result.
}}
\value{
density (\code{dsn}), probability (\code{psn}),
quantile (\code{qsn}) or  random sample (\code{rsn})
from the skew-normal distribution with given \code{location}, \code{scale}
and \code{shape} parameters.
}
\section{BACKGROUND}{
The family of skew-normal distributions is an extension of the normal
family, via the introdution of a \code{shape} parameter which regulates
skewness; when \code{shape=0}, the skew-normal distribution reduces to the
normal one.  The density of the SN distribution when \code{location=0}
and \code{scale=1} is \code{2*dnorm(x)*pnorm(shape*x)}.
A multivariate version of the distribution exists.
See the references below for additional information.
}
\details{
\code{psn} make use of function \code{T.Owen}
}
\references{
Azzalini, A. (1985).
A class of distributions which includes the normal ones.
\emph{Scand. J. Statist.}
\bold{12}, 171-178.


Azzalini, A. and Dalla Valle, A. (1996).
The multivariate skew-normal distribution.
\emph{Biometrika}
\bold{83}, 715--726.
}
\seealso{
\code{\link{T.Owen}}, \code{\link{dmsn}}, \code{\link{sn.mle}}
}
\examples{
pdf <- dsn(seq(-3,3,by=0.1), shape=3)
cdf <- psn(seq(-3,3,by=0.1), shape=3)
qu <- qsn(seq(0.1,0.9,by=0.1), shape=-2)
rn <- rsn(100, 5, 2, 5)
}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
