\name{sn.dev}
\alias{sn.dev}
\title{
Deviance for skew-normal distributions
}
\description{
Computes twice the negative log-likelihood (the "deviance")
for regression models with errors having a skew-normal distribution.
}
\usage{
sn.dev(cp, X, y, trace=F)
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
}}
\value{
a numeric value
}
\details{
This is the objective function of \code{sn.mle}, whose documentation
gives additional details.
}
\seealso{
\code{\link{sn.mle}},\code{\link{msn.mle}}
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
