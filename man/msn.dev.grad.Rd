\name{msn.dev.grad}
\alias{msn.dev.grad}
\title{
Gradient of the deviance for multivariate skew-normal distributions
}
\description{
Computes the gradient of twice the negative profile log-likelihood 
for multivariate regression models with errors having a multivariate
skew-normal distribution.
}
\usage{
msn.dev.grad(param, X, y, freq, trace=F)
}
\arguments{
\item{param}{
a numeric vector of parameter values.
}
\item{X}{
a matrix of explanatory variables.
Missing values (\code{NA}s) are not allowed.
}
\item{y}{
a matrix contaning the response variable.
Missing values (\code{NA}s) are not allowed.
}
\item{freq}{
a vector of weights, with length equal to \code{nrow(y)}.
}
\item{trace}{
logical value.
If \code{trace=T}, details are printed. Default value is \code{F}.
}}
\value{
a vector contaning the gradient
}
\details{
This is the gradient function of \code{msn.dev}, 
which is the objective function of \code{msn.mle};
refer to its documentation for additional details.
}
\seealso{
\code{\link{msn.mle}},\code{\link{msn.dev}}}
}
\keyword{multivariate}
\keyword{regression}


