\name{sn.dev.gh}
\alias{sn.dev.gh}
\title{
Derivatives of deviance for skew-normal distributions
}
\description{
Computes first and second derivatives of twice the negative log-likelihood 
(the "deviance")
for regression models with errors having a skew-normal distribution.
}
\usage{
sn.dev.gh(cp, X, y, trace=F, hessian=F)
}
\arguments{
\item{cp}{
a vector of initial values for the centred parameters,
with \code{length(cp)=ncol(X)+2}
}
\item{X}{
a matrix of explanatory variables. Missing values (\code{NA}s) are not allowed.
}
\item{y}{
a vector contaning the observed variable.
Missing values (\code{NA}s) are not allowed.
}
\item{trace}{
logical value.
If \code{trace=T}, details are printed. Default value is \code{F}.
}
\item{hessian}{
logical value.
If \code{hessian=T}, the Hessian is computed. Default value is \code{F}.
}
\value{
a vector containing the gradient; 
if \code{hessian=T}, this is given an attribute containing the Hessian matrix
}
\details{
This function produces of derivatives of the objective function of \code{sn.mle}, 
whose documentation gives additional details.
}
\seealso{
\code{\link{sn.mle}}, \code{\link{sn.dev}},\code{\link{msn.mle}}
}
\keyword{regression}

