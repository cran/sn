\name{msn.moment.fit}
\alias{msn.moment.fit}
\title{
Fitting multivariate skew-normal distributions with method of moments
}
\description{
Fits a multivariate skew-normal (MSN) distribution to data 
using the method of moments.
}
\usage{
msn.moment.fit(y)
}
\arguments{
\item{y}{
a matrix or a vector.  If \code{y} is a matrix, its rows refer to
observations, and its columns to components of the multivariate
distribution. If \code{y} is a vector, it is converted to a one-column
matrix, and a scalar skew-normal distribution is fitted.
}}
\value{
A list containing the following components:

\item{xi}{
a vector with the location parameter. 
}
\item{Omega}{
a variance matrix  representing the association parameter.
}
\item{alpha}{
a vector of shape  parameters.
}
\item{omega}{
vector of scale parameters corresponding to \code{Omega}.
}
\item{delta}{
the parameter \code{delta} which determines the shape of the marginal 
distributions.
}
\item{skewness}{
numeric vector with marginal indices of skewness (the standardised
third cumulant).
}
\item{admissible}{
a logical value indicating if the estimated parameters lie in 
the admissible region.
}}

\details{
This function is used by \code{msn.mle} to obtain preliminary
estimates, if starting values are not provided. 
After removing the regression component, estimated by ordinary least squares, 
\code{msn.moment.fit} is used to obtain preliminary estimate of 
the other parameters from the least squares residuals.

Although the function accepts a vector \code{y} as input, the use of
\code{sn.mle} is recommended in the scalar case.
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
\code{\link{msn.mle}}, \code{\link{sn.mle}}
}
\keyword{distribution}
