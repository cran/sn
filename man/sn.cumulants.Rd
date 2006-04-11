\name{sn.cumulants}
\alias{sn.cumulants}
\title{
Cumulants of the skew-normal distribution
}
\description{
Cumulants of the skew-normal distribution.
}
\synopsis{
sn.cumulants(location = 0, scale = 1, shape = 0, dp = NULL, n = 4) 
}
\usage{
sn.cumulants(location = 0, scale = 1, shape = 0,  n = 4)
sn.cumulants(dp=, n = 4)
}
\arguments{
\item{location}{
 location parameter (vector)
 }
\item{scale}{
scale parameter (vector)
}
\item{shape}{
shape parameter (vector)
}
\item{dp}{
 a vector of three elements, whose elements are \code{(location, scale, shape)}
 respectively. If \code{dp} is specified, then the individual
 parameters must not be.
}
\item{n}{
  a scalar integer of the maximal order or cumulants required
}}
\value{
the cumulants up to order \code{n} of the skew-normal distribution
with \code{location=0}, \code{scale=1} and \code{shape} as selected.
}
\details{
The moment generating function (hence the cumulant generating function)
of the distribution is given in the refence below.
The computations method used is proved analytically up to \code{n=3} but it is
seen to behave correctly up to the order which was checked (\code{n=8}).
}
\references{
Azzalini, A. (1985).
A class of distributions which includes the normal ones.
\emph{Scand. J. Statist.}
\bold{12}, 171-178.
}
\seealso{
\code{\link{dsn}},\code{\link{zeta}}
}
\examples{
sn.cumulants(shape=c(0,2.5,5,10), n=5)
sn.cumulants(dp=c(10,3,-8), n=6)
}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
