\name{zeta}
\alias{zeta}
\title{
Function `log(2*pnorm(x))' and its derivatives
}
\description{
  The function \code{log(2*(pnorm(x))} and its derivatives, 
  including inverse Mills ratio.
}
\usage{
zeta(k, x)
}
\arguments{
\item{k}{
  an integer scalar between 0 and 5.
}
\item{x}{
  a vector. Missing values (\code{NA}s)  and \code{Inf}s are allowed
}}
\value{
  a vector giving the \code{k}-th order derivative evaluated at \code{x}
}
\details{
  For \code{k} between 0 and 5, the derivative of  order \code{k}
  of \code{log(2*pnorm(x))} is  evaluated;   the derivative of
  order \code{k=0} refers to  the function itself.
  If \code{k} is not integer, it is converted to integer and a warning
  message is generated.
  If \code{k<0} or \code{k>5},  \code{NULL} is returned.

  The computation for \code{k>1} is reduced to the case \code{k=1}, making use
  of expressions given by Azzalini and Capitanio (1999).  For numerical
  stability, the evaluation of \code{zeta(1,x)} when \code{x < -50} makes use
  of the asymptotic expansion (26.2.13) in Abramowitz and Stegun (1964).

  \code{zeta(1,-x)} equals \code{dnorm(x)/pnorm(-x)} (in principle, apart from
  the asymptotic expansion mentioned above), called the  
  \emph{inverse Mills ratio}.
}
\references{
  Abramowitz, M. and Stegun, I. A., editors (1964).
  \emph{Handbook of Mathematical Functions}. 
  Dover Publications.


Azzalini, A. and Capitanio, A. (1999).
  Statistical applications of the multivariate skew-normal distribution.
  Technical report available at \url{http://azzalini.stat.unipd.it/SN}.
  An abriged version is published in   \emph{J.Roy.Statist.Soc. B}
  \bold{61}, 579--602.

}
\examples{
y <- zeta(2,seq(-20,20,by=0.5))
#
for(k in 0:5) curve(zeta(k,x), from=-1.5, to=5, col = k+2, add = k > 0)
legend(3.5, -0.5, legend=as.character(0:5), col=2:7, lty=1)
}
\keyword{math}
