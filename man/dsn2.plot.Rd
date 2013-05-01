\name{dsn2.plot}
\alias{dsn2.plot}
\title{
Plot of Bivariate Skew-normal Density Function
}
\description{
Produces a contour plot of the density function of a bivariate 
skew-normal variate.
}

\usage{dsn2.plot(x, y, xi, Omega, alpha, dp = NULL, ...) }

\arguments{
\item{x}{
vector of values of the first component.
}
\item{y}{
vector of values of the second component.
}
\item{xi}{
a vector of length 2 containing the location parameter.
}
\item{Omega}{
a 2 by 2 matrix containing a covariance matrix.
}
\item{alpha}{
a vector of length 2 containing the shape parameter.
}
\item{dp}{
  a list with components named \code{xi, Omega, alpha}, containing
  quantities as described above. If this parameter is set, the
  individual parameters must not be.
  }
\item{...}{
additional parameters to be passed to \code{contour}.
}}

\value{
A list containing the original input parameters plus a matrix
containing the density function evaluated at the grid formed 
by the \code{x} and \code{y} values. 
}

\details{Typical usages are
\preformatted{%
dsn2.plot(x, y, xi, Omega, alpha, ...)
dsn2.plot(x, y, dp=, ...)
}
  The density function is evalutate at the grid of points whose
  coordinates are given by vectors \code{x} and \code{y}.
  The actual computation is done by the  function \code{dmsn}.
  A contour level plot is produced on the graphical window. 
}
\section{Background}{
The multivariate skew-normal distribution is discussed by
Azzalini and Dalla Valle (1996); the \code{(Omega,alpha)} parametrization
adopted here is the one of Azzalini and Capitanio (1999).
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
x <- y <- seq(-5, 5, length=35)
dsn2.plot(x, y, c(-1,2), diag(c(1,2.5)), c(2,-3))
}
\keyword{distribution}

