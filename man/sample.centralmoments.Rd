\name{sample.centralmoments}
\alias{sample.centralmoments}
\title{Sample centralmoments}
\description{
  Computes sample central moments up to a given order and the
  first moment from the origin}
}
\usage{
sample.centralmoments(x, w = rep(1, length(x)), order=4)
}
\arguments{
  \item{x}{a vector of sample values}
  \item{w}{an optional vector of  weights}
  \item{order}{the maximal order of the central moments to
    be computed; it must be a positive integer (default value 4)
}}
\value{
  A vector containing the first sample central moments, 
  in position  \code{[2:order]}, and  the first moment from the
  origin,   in the first position of the returned vector
}
\details{
  \code{NA}'s are allowed but removed. Averaging of appropriate
  quantities is actually performed by \code{weighted.mean}
  }
\author{Adelchi Azzalini}
\note{
  The second component of the returned vector (if \code{order>1})
  gives the sample variance; notice that it differs from the value
  returned by \code{var(x)}, since this gives the corrected sample
  variance.
  
  Used in conjunction with \code{st.cumulants.inversion},  this 
  function allows to fit a skew-t distribution by the methods
  of moments; see the example below.
}
\seealso{\link{st.cumulants.inversion}, \link{weighted.mean}}
\examples{
data(ais, package='sn')
mom <- sample.centralmoments(ais[,"bmi"])
st.cumulants.inversion(cum=c(mom[1:3], mom[4]-3*mom[2]^2))
# parameters of the fitted ST distribution
}
\keyword{univar}
