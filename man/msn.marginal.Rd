\name{msn.marginal}
\alias{msn.marginal}
\title{
Marginal compontents of a multivariate skew-normal distribution
}
\description{
Computes the marginal distribution of a subset of components of a
multivariate skew-normal distribution.
}
\usage{
msn.marginal(xi, Omega, alpha, comp)
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
\item{comp}{
a vector containing a subset of \code{1:k} which selects the components
whose values are to be fixed. A permutation of \code{1:k} is allowed, and
the components of comp' do not need to be sorted.
}}
\value{
a list containing the \code{xi, Omega, alpha} parameters for the marginal
distribution.
}
\details{
See Azzalini and Capitanio (1999) for details.
}
\references{
Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{dmsn}}, \code{\link{msn.conditional}}
}
\examples{
xi <- c(10,0,-30)
Omega <- 5*diag(3)+outer(1:3,1:3)
alpha <- c(1,-3,5)
marg31 <- msn.marginal(xi,Omega,alpha,c(3,1))
}
\keyword{multivariate}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
