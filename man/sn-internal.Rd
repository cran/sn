\name{sn-internal}
\alias{msn.dev}
\alias{msn.dev.grad}
\alias{sn.dev}
\alias{sn.dev.gh}
\alias{msn.moment.fit}
\alias{num.deriv}
\alias{mst.dev}
\alias{mst.dev.grad}
\alias{st.dev.fixed}
\alias{sn.logL.grouped}
\alias{st.logL.grouped}
\alias{sn.SFscore}
\title{Internal sn functions}
\description{
 Internal  functions of package \emph{sn}
}
\usage{
msn.dev(param, X, y, freq, trace=FALSE)
msn.dev.grad(param, X, y, freq, trace=FALSE)
msn.moment.fit(y)
mst.dev(param, X, y, freq, fixed.df=NA, trace=FALSE)
mst.dev.grad(param, X, y, freq, fixed.df=NA, trace=FALSE)
num.deriv(coefficients, FUN, ...)
st.dev.fixed(free.param, X, y, freq, trace=FALSE, fixed.comp=NA, fixed.values=NA)
sn.dev(cp, X, y, trace=FALSE)
sn.dev.gh(cp, X, y, trace=FALSE, hessian=FALSE)
sn.logL.grouped(param, breaks, freq, trace=FALSE)
st.logL.grouped(param, breaks, freq, trace=FALSE)
sn.SFscore(shape, data, trace=FALSE)
}
\arguments{
\item{param,cp, coefficients, shape}{
a numeric vector of parameter values.
}

\item{X}{
a matrix of explanatory variables; must have \code{col(X)} equal to
\code{length(y)}. Missing values (\code{NA}) are not allowed.
If \code{X} is missing, a one-column matrix of 1's is created.
}
\item{y}{
a numeric vector or matrix. Missing values (\code{NA}s) are not allowed.
}
   
\item{freq}{
a vector of frequencies.
}
\item{trace}{
logical value which controls printing of the algorithm convergence.
If \code{trace=TRUE}, details are printed. Default value is \code{FALSE}.
}
\item{free.param}{
a vector of suitably re-parametrized parameters, not to be kept fixed during 
iteration.
}
\item{fixed.comp}{
a vector containing the subset of the parameters for which the
profile log-likelihood function is required; it can be of length 1 or 2.
}
\item{fixed.values}{
a numeric vector of values or a matrix with two columns, giving the
range spanned by the selected parameters. 
}
\item{fixed.df}{
a scalar value contaning the degrees of freedom (df), if these must
be taked as fixed, or \code{NA} (deafult value) if df is a parameter
to be estimated.
}
\item{breaks}{
a vector contaning the cut points of the groups, given
in ascending order. The last value can be \code{Inf}, the
first one can be \code{-Inf}
}
}
\value{
  A numeric value or a vector or a list.
}

\details{
  These functions are not intended to be called directly by the user.
}
\keyword{internal}

