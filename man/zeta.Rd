\name{zeta}
\alias{zeta}
\title{
Function `log(2*pnorm(x))' and its derivatives
}
\description{
  The function \code{log(2*(pnorm(x))} and its derivatives up to order 4.
}
\usage{
zeta(k, x)
}
\arguments{
\item{k}{
  an integer scalar between 0 and 4.
}
\item{x}{
  a vector. Missing values (\code{NA}s)  and \code{Inf}s are allowed
}}
\value{
  a vector giving the \code{k}-th order derivative evaluated at \code{x}
}
\details{
  For \code{k} between 0 and 4, the derivative of  order \code{k}
  of \code{log(2*pnorm(x))} is  evaluated, where the derivative of
  order \code{k=0} refers to  the function itself.
  If \code{k} is not integer, it is converted to integer and a warning
  message is generated.
  If \code{k<0} or \code{k>4},  \code{NULL} is returned.
  This function is used by \code{sn.dev} and \code{msn.dev}, among others.
}
\seealso{
\code{\link{sn.mle}}, \code{\link{msn.mle}}
}
\examples{
y <- zeta(2,seq(-20,20,by=0.5))
#
for(k in 0:4) curve(zeta(k,x), from=-1, to=5, col = k+2, add = k > 0)
legend(3, -0.5, legend=as.character(0:4), col=2:6, lty=1)
}
\keyword{math}
