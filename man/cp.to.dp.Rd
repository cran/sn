\name{cp.to.dp}
\alias{cp.to.dp}
\alias{dp.to.cp}
\title{
Conversion between equivalent parametrizations  
}
\description{
  Convert direct parameters (DP) to centred parameters (CP) of the
  one-dimensional skew-normal distribution and \emph{vice versa}
}
\usage{
cp.to.dp(param)
dp.to.cp(param)
}
\arguments{
\item{param}{
a vector of length at least three.
If \code{lenght(param)} is \code{m+2}, then the first \code{m} components refer
to the regression coefficients (or the location parameter, in
case \code{m} is 1), and the remaining two components refer to scale and
shape, respectively; their role is preserved across parametrizations.
}}
\value{
a vector of the same length of \code{param}, representing \code{param} in the
alternative parametrization; \code{cp.to.dp} converts centred to direct
parameters, \code{dp.to.cp} converts direct to centred parameters.
}
\details{For a description of the two parametrizations, 
see the reference below.
}
\references{
Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{sn.mle}}, \code{\link{sn.em}}
}
\examples{
cp <- dp.to.cp(c(30,30,2,4))
dp <- cp.to.dp(cp)
}
\keyword{distribution}
