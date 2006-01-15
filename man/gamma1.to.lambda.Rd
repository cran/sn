\name{gamma1.to.lambda}
\alias{gamma1.to.lambda}
\title{
Converts skewness to shape parameter of skew-normal distribution
}
\description{
For a given value of the index of skewness (standardized third
cumulant), the function finds the corresponding shape parameter
of a skew-normal distribution.
}
\usage{
gamma1.to.lambda(gamma1)
}
\arguments{
\item{gamma1}{
a numeric vector of indices of skewness.
}}
\value{
A numeric vector of the corresponding shape parameters.
}
\details{
Feasible values for input must have \code{abs(gamma1)<0.5*(4-pi)*(2/(pi-2))^1.5},
which is about 0.99527.
If some values of \code{gamma1} are not in the feasible region, a warning
message is issued, and \code{NA}s are returned.


See the reference below for the expression of the index of skewnnes
of a skew-normal distribution.
}
\references{
Azzalini, A. (1985).
A class of distributions which includes the normal ones.
\emph{Scand. J. Statist.}
\bold{12}, 171-178.
}
\seealso{
\code{\link{dsn}}
}
\examples{
}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
