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
msn.mle(X, y, freq, start, trace=FALSE,  
        algorithm=c("nlminb", "Nelder-Mead", "BFGS", "CG",  "SANN"),
        control=list())
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
a list containing the components \code{beta},\code{Omega}, \code{alpha},
of the type described below. The \code{dp} component of the returned
list from a previous call has the required format.
}
\item{trace}{
logical value which controls printing of the algorithm convergence.
If \code{trace=TRUE}, details are printed. Default value is \code{FALSE}.
}
\item{algorithm}{
a character string which selects the numerical optimization procedure 
used to maximize the log-likelihood function. If this string is set 
equal to \code{"nlminb"}, then this function is called; in all other cases,
\code{optim} is called, with \code{method} set equal to the given string.
Default value is \code{"nlminb"}.
}
\item{control}{
this parameter is passed to the  optimizer selected via \code{algorithm};
see the documentation of \code{nlminb} or \code{optim} for its usage.
}}
\value{
A list containing the following components:

\item{call}{
a string containing the calling statement.
}
\item{dp}{
  a list containing terms named \code{beta}, \code{Omega}, \code{alpha},
  where \code{beta} is a matrix of regression coefficients with
  \code{dim(beta)=c(nrow(X),ncol(y))}, \code{Omega} is a covariance
  matrix   of order \code{ncol(y)}, \code{alpha} is a vector of shape
  parameters of length \code{ncol(y)}.
}
\item{se}{
  a list containing the components \code{beta}, \code{alpha}, \code{info}.
  Here, \code{beta} and \code{alpha} are the standard errors for the
  corresponding point estimates;
  \code{info} is the observed information matrix for the working parameter,
  as explained below.
}
\item{algorithm}{
  the list returned by the chose optimizer, either \code{nlminb}
  or  \code{optim}, plus an item with the \code{name} of the selected
  algorithm; see the documentation of either \code{nlminb}
  or  \code{optim}  for explanation  of the other components.
}}
\details{
The parameter \code{freq} is intended for use with grouped data,
setting the values of \code{y} equal to the central values of the
cells; in this case the resulting estimate is an approximation
to the exact maximum likelihood estimate. If \code{freq} is not
set, exact maximum likelihood estimation is performed.

The working parameter used in the maximization stage is
\code{c(beta,alpha/omega)}, since a profile `deviance'
\code{-2*loglikelihood} for this parameter is actually used;
see Azzalini and Capitanio (1999, section 6.1) for details.
The selected optimizer (\code{nlminb} or \code{optim}) is called,
supplying the gradient of the profile deviance. In case the
optimizer is  \code{optim}),  the gradient may or may not
be used, depending on which specific method has been selected.

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
\code{\link{dmsn}},\code{\link{msn.fit}}, \code{\link{nlminb}}, 
\code{\link{optim}}  }
}
\examples{
data(ais, package="sn")
attach(ais)
# a simple-sample case
a <- msn.mle(y=cbind(Ht,Wt))
#
# a regression case:
b  <- msn.mle(X=cbind(1,Ht,Wt), y=ssf)
b1 <- msn.mle(X=cbind(1,Ht,Wt), y=ssf, algorithm="Nelder-Mead")
b2 <- msn.mle(X=cbind(1,Ht,Wt), y=ssf, start=b1$dp)
}
\keyword{distribution}
\keyword{regression}
