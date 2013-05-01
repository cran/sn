\name{msn.quantities}
\alias{msn.quantities}
\title{
Quantities related to the multivariate skew-normal distribution.
}
\description{
Computes mean vector, variance matrix and other relevant quantities
of a given multivariate skew-normal distribution.
}
\usage{msn.quantities(xi=rep(0,length(alpha)), Omega, alpha, dp = NULL)}

\arguments{
\item{xi}{
numeric vector giving the location parameter, of length \code{d}, say.
Missing values are not allowed.
}
\item{Omega}{
a covariance matrix of size \code{d} by \code{d}.
Missing values are not allowed.
}
\item{alpha}{
numeric vector of shape parameter of length \code{d}.
Missing values are not allowed.
}
\item{dp}{
  a list with three components named \code{xi}, \code{Omega},
  \code{alpha}, containing quantities as described above.
  If \code{dp} is set, then the component parameters must not be.
}}
\value{
A list containing the following components:

\item{xi}{
the input parameter \code{xi}.
}
\item{Omega}{
the input parameter  \code{Omega}.
}
\item{alpha}{
the input parameter  \code{alpha}.
}
\item{omega}{
vector of scale parameters.
}
\item{mean}{
 the mean value of the distribution (vector)
}
\item{variance}{
variance-covariance matrix of the distribution.
}
\item{Omega.conv}{
concentration matrix associated to \code{Omega}, i.e. its inverse.
}
\item{Omega.cor}{
correlation matrix associated to \code{Omega}.
}
\item{Omega.pcor}{
partial correlations matrix associated to \code{Omega}.
}
\item{lambda}{
shape parameters of the marginal distributions
}
\item{Psi}{
correlation matrix of the equivalent \code{(lambda,Psi)} parametrization.
}
\item{delta}{
the parameter \code{delta} which determines the shape of the marginal 
distributions; this is related to \code{lambda}
}
\item{skewness}{
numeric vector with marginal indices of skewness (the standardised
third cumulant).
}}
\details{Typical usages are
\preformatted{%
msn.quantities(xi=rep(0,length(alpha)), Omega, alpha)
msn.quantities(dp= )
}
The meaning of the parameters is explained in the references below;
see especially Azzalini and Capitanio (1999).
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
\code{\link{dmsn}}
}
\examples{
Omega <- 5*diag(3)+outer(1:3,1:3)
msn.quantities(c(0,0,1), Omega, c(-2,2,3))
}
\keyword{multivariate}
\keyword{distribution}

