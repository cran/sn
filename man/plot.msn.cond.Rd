\name{plot.msn.cond}
\alias{plot.msn.cond}
\title{
Plot of the density of a conditional skew-normal variate
}
\description{
Plot of the exact and of the approximate density function of
a multivariate skew-normal variate conditionally on the values
taken on by some components.
}
\usage{
plot.msn.cond(xi, Omega, alpha, fixed.comp, fixed.values, n=35)
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
}
\item{n}{
an integer value which determines the grid size of the density
computations and plot.
}}
\value{
a list containing the following elements:

\item{cumulants,}{
two lists as returned by \code{msn.conditional}
}
\item{pdf}{
a list containing the coordinates \code{x} and \code{y} of the points where the 
densities have been evaluated, and the matrices \code{f.exact} and \code{f.fitted}
of the exact and fitted conditional desntities.
}
\item{rel.error,}{
summary statistics of relative and absolute error of the approximation
}}
\section{Side Effects}{
A contour plot of the exact and approximate densities is produced
on a graphical device.
}
\details{
}
\references{
Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{msn.conditional}}, \code{\link{dmsn}}
}
\examples{
Omega <- diag(3)+0.5*outer(rep(1,3),rep(1,3))
a<- plot.msn.cond(rep(0,3), Omega, 1:3, 3, -0.75)
}
\keyword{multivariate}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
