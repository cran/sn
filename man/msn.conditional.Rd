\name{msn.conditional}
\alias{msn.conditional}
\title{
Cumulants and distribution of a skew-normal variate after conditioning
}
\description{
Finds cumulants up to 3rd order of a multivariate skew-normal
distribution conditionally on the values taken on by some of
its components, and finds a multivariate skew-normal distribution
with the same cumulants.
}
\usage{
msn.conditional(xi, Omega, alpha, fixed.comp, fixed.values)
}
\arguments{
\item{xi}{
a numeric vector of length \code{k}, say, giving the location parameter.
}
\item{Omega}{
a covariance matrix of dimension \code{(k,k)}.
}
\item{alpha}{
a numeric vector of length \code{k}, which regulates the shape of the density.
}
\item{fixed.comp}{
a vector containing a subset of \code{1:k} which selects the components
whose values are to be fixed; it must be of length \code{k-2}.
}
\item{fixed.values}{
a numeric vector of values taken on by the components \code{fixed.comp};
it must be of the same length of \code{fixed.comp}.
}}
\value{
a list containing the following elements:

\item{cumulants}{
a list containing mean vector, variance matrix, and indices of 
skewness of the conditional distribution.
}
\item{fit}{
a list containing the parameters of the fitted skew-normal
distribution in the \code{(xi,Omega,alpha)} parametrization, plus
the vector \code{delta}.
}}
\details{
See the reference below for details and background.
}
\references{
Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{plot.msn.cond}}, \code{\link{msn.marginal}}
}
\examples{
Omega <- diag(3)+0.5*outer(rep(1,3),rep(1,3))
a<- msn.conditional(rep(0,3), Omega, 1:3, 3, -0.75)
}
\keyword{multivariate}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
