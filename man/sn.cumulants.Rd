\name{sn.cumulants}
\alias{sn.cumulants}
\title{
cumulants of the skew-normal distribution
}
\description{
cumulants of the skew-normal distribution
}
\usage{
sn.cumulants(shape=0, n=4)
}
\arguments{
\item{shape}{
a vector of shape parameter
}
\item{n}{
a scalar integer
}}
\value{
the cumulants up to order \code{n} of the skew-normal distribution
with \code{location=0}, \code{scale=1} and \code{shape} as selected
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
cum <- sn.cumulants(c(0,1,2,5,10),4)
}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
