\name{msn.marginal}
\alias{msn.marginal}
\title{
Marginal components of a multivariate skew-normal distribution
}
\description{
Computes the marginal distribution of a subset of components of a
multivariate skew-normal distribution.
}
\synopsis{
msn.marginal(xi=NULL, Omega=NULL, alpha=NULL, comp=1:d, dp=NULL)
}
\usage{
msn.marginal(xi, Omega, alpha, comp)
msn.marginal(dp=, comp)
}
\arguments{
\item{xi}{
  a numeric vector of length \code{d}, say, giving the location parameter.
}
\item{Omega}{
  a covariance matrix of dimension \code{(d,d)}.
}
\item{alpha}{
  a numeric vector of length \code{d}, which regulates the shape of the density.
}
\item{comp}{
a vector containing a subset of \code{1:d} selecting the components of the
marginal distribution. A permutation of \code{1:d} is allowed, and
the components of \code{comp} do not need to be sorted.
}
\item{dp}{
  a list containing the components \code{xi}, \code{Omega}, \code{alpha}
}}
\value{
A list containing  components \code{xi, Omega, alpha} with the 
parameters of the marginal distribution. If \code{length(comp)} is  equal to
\code{r}, say, then the new components are of size \code{r, (r,r), r}, 
respectively.
}
\background{
See Azzalini and Capitanio (1999) for background information.
}
\references{
Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B} \bold{61}, 579--602.
}
\seealso{
\code{\link{dmsn}}, \code{\link{msn.conditional}}, \code{\link{msn.affine}}
}
\examples{
xi <- c(10,0,-30)
Omega <- 5*diag(3)+outer(1:3,1:3)
alpha <- c(1,-3,5)
marg31 <- msn.marginal(xi,Omega,alpha,c(3,1))
}
\keyword{multivariate}
\keyword{distribution}

