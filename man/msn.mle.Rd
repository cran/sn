\name{msn.mle}
\alias{msn.mle}
\title{
Maximum likelihood estimation for a multivariate skew-normal distribution
}
\description{
Fits a multivariate skew-normal (MSN) distribution to data, or fits a
linear regression model with multivariate skew-normal errors,
using maximum likelihood estimation.
}
\usage{
msn.mle(X, y, freq, start, trace=FALSE,  method="BFGS", control=list(iter.max=150, x.tol=1e-8) )
}
\arguments{
\item{y}{
a matrix or a vector.  If \code{y} is a matrix, rows refer to
observations, and columns to components of the multivariate
distribution. If \code{y} is a vector, it is converted to a one-column
matrix, and a scalar skew-normal distribution is fitted.
}
\item{X}{
a matrix of covariate values.
If missing, a one-column matrix of 1's is created; otherwise,
it must have the same number of rows of \code{y}.
}
\item{freq}{
a vector of weights.
If missing, a one-column matrix of 1's is created; otherwise
it must have the same number of rows of \code{y}.
}
\item{start}{
a list contaning the components \code{beta},\code{Omega}, \code{alpha},
of the type described below. The \code{dp} component of the returned
list from a previous call has the required format.
}
\item{trace}{
logical value which controls printing of the algorithm convergence.
If \code{trace=TRUE}, details are printed. Default value is \code{F}.
}
\item{method}{
this parameter is just passed to the optimizer \code{optim}; 
see the documentation of this function for its usage. Default value is
\code{"BFGS"}.
}
\item{control}{
this parameter is passed to the optimizer \code{optim}; 
see the documentation of this function for its usage.
}}
\value{
A list containing the following components:

\item{call}{
a string containing the calling statement.
}
\item{dp}{
a list containing the direct parameters \code{beta}, \code{Omega}, \code{alpha}.
Here, \code{beta} is a matrix of regression coefficients with
\code{dim(beta)=c(nrow(X),ncol(y))}, \code{Omega} is a covariance matrix of
order \code{ncol(y)}, \code{alpha} is a vector of shape parameters of length
\code{ncol(y)}.
}
\item{se}{
a list containing the components \code{beta}, \code{alpha}, \code{info}.
Here, \code{beta} and \code{alpha} are the standard errors for the
corresponding point estimates;
\code{info} is the observed information matrix for the working parameter,
as explained below.
}
\item{optim}{
the list returned by the optimizer \code{optim}; see the documentation
of this function for explanation of its components.
}}
\details{
The parameter \code{freq} is intended for use with grouped data,
setting the values of \code{y} equal to the central values of the
cells; in this case the resulting estimate is an approximation
to the exact maximum likelihood estimate. If \code{freq} is not
set, exact maximum likelihood estimation is performed.


The working parameter used in the maximization stage is
\code{c(beta,alpha/omega)}, since a profile `deviance' \code{-2*loglikelihood}
for this parameter is actually used;
see Azzalini and Capitanio (1999) for details.
The optimizer \code{optim} is called, supplying the gradient of
the profile deviance.


The function can take a vector \code{y} as input; however the use of
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
\code{\link{dmsn}},\code{\link{sn.mle}},\code{\link{msn.fit}}, \code{\link{optim}}}
}
\examples{
data(ais,package=sn)
attach(ais)
# a simple-sample case
b <- msn.mle(y=cbind(Ht,Wt))
#
# a regression case:
a <- msn.mle(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-6))
#
# refine the previous outcome
a1 <- msn.mle(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-9), start=a$dp)
}
\keyword{distribution}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
