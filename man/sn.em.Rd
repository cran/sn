\name{sn.em}
\alias{sn.em}
\title{
Fitting Skew-normal variables using the EM algorithm
}
\description{
Fits a skew-normal (SN) distribution to data, or fits a linear regression
model with skew-normal errors, using the EM algorithm to locate the MLE
estimate. The estimation procedure can be global or it can fix some
components of the parameters vector.
}
\usage{
sn.em(X, y, fixed, p.eps=0.0001, l.eps=0.01, trace=FALSE, data=FALSE)
}
\arguments{
\item{y}{
a vector contaning the observed variable. This is the response
variable in case of linear regression.
}
\item{X}{
a matrix of explanatory variables. If \code{X} is missing, then a one-column
matrix of all 1's is created. If \code{X} is supplied, and an intercept term
is required, then it must include a column of 1's.
}
\item{fixed}{
a vector of length 3, indicating which components of the
parameter vector must be regarded as fixed. In \code{fixed=c(NA,NA,NA)},
which is the default setting, a global maximization is performed.
If the 3rd component is given a value, then maximization is performed
keeping that value fixed for the shape parameter. If the 3rd and 2nd
parameters are fixed, then the scale and the shape parameter are
kept fixed. No other patterns of the fixed values are allowed.
}
\item{p.eps}{
numerical value which regulates the parameter convergence tolerance.
}
\item{l.eps}{
numerical value which regulates the log-likelihood convergence tolerance.
}
\item{trace}{
logical value which controls printing of the algorithm convergence.
If \code{trace=TRUE}, details are printed. Default value is \code{F}.
}
\item{data}{
logical value. If \code{data=TRUE}, the returned list includes the original
data. Default value is \code{data=FALSE}.
}}
\value{
a list with the following components:

\item{dp}{
a vector of the direct parameters, as explained in the references below.
}
\item{cp}{
a vector of the centred parameters, as explained in the references below.
}
\item{logL}{
the log-likelihood at congergence.
}
\item{data}{
optionally (if \code{data=TRUE}), a list containing \code{X} and \code{y,} as supplied
on input, and a vector of \code{residuals}, which should have an approximate
SN distribution with \code{location=0} and \code{scale=1}, in the direct
parametrization.
}}
\details{
The function works using the direct parametrization; on convergence,
the output is then given in both parametrizations.


This function is based on the EM algorithm; it is generally quite slow,
but it appears to be very robust.
See \code{sn.mle} for an alternative method, which also returns standard
errors.
}
\section{Background}{
Background information on the SN distribution is given by Azzalini (1985).
See  Azzalini and Capitanio (1999) for a more detailed discussion of
the direct and centred parametrizations.
}
\references{
Azzalini, A. (1985).
A class of distributions which includes the normal ones.
\emph{Scand. J. Statist.}
\bold{12}, 171-178.


Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{dsn}}, \code{\link{sn.mle}}, \code{\link{cp.to.dp}}
}
\examples{
data(ais,package=sn)
attach(ais)
#
a<-sn.em(y=bmi)
#
a<-sn.em(X=cbind(1,lbm,lbm^2),y=bmi)
#
M<-model.matrix(~lbm+I(ais$sex))
b<-sn.em(M,bmi)
#
fit <- sn.em(y=bmi, fixed=c(NA, 2, 3), l.eps=0.001)
}
\keyword{regression}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
