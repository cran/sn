\name{zeta}
\alias{zeta}
\title{
Function `log(2*P(x))' and its derivatives
}
\description{
The function \code{log(2 P(x))} and its derivatives up to order 4,
if P(x) denotes the standard normal distribution function.
}
\usage{
zeta(k, x)
}
\arguments{
\item{k}{
an integer scalar.
}
\item{x}{
a vector. Missing values (\code{NA}s) are allowed, but \code{Inf}s are not.
}}
\value{
a vector; this is the \code{k}-th order derivative evaluated at \code{x}
}
\details{
\code{k} denotes the order of the derivative. If \code{k} is 0, the function
is evaluated, using \code{pnorm(x)} for \code{x>-35},  an asymptotic expansion
otherwise. For \code{k} between 1 and 4, the derivative of that order is
evaluated. For \code{k>4}, a warning message is issued, and a vector of
\code{NA}s is returned.


This function is used by \code{sn.dev} and \code{msn.dev}, among others.
}
\seealso{
\code{\link{sn.mle}}, \code{\link{msn.mle}}
}
\examples{
y <- zeta(2,seq(-20,20,by=0.5))

for(i in 0:4) 
      curve(zeta(i,x),-3, 4, ylim = c(-1,2), col = i+1, add = i > 0)
}
\keyword{math}
% Converted by Sd2Rd version 0.3-3.
