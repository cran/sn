\name{msn.dev}
\alias{msn.dev}
\title{
Deviance for multivariate skew-normal distributions
}
\description{
Computes twice the negative profile log-likelihood (essentially,
a form of the deviance)
for multivariate regression models with errors having a multivariate
skew-normal distribution.
}
\usage{
msn.dev(param, X, y, freq, trace=FALSE)
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
If \code{trace=TRUE}, details are printed. Default value is \code{F}.
}}
\value{
A numeric value.
}

\details{
This is the objective function of \code{msn.mle}, whose documentation
gives additional details.
}
\seealso{
\code{\link{msn.mle}},\code{\link{msn.fit}},\code{\link{sn.mle}}
}
\keyword{multivariate}
\keyword{regression}
