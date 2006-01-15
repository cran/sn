\name{sn.mmle}
\alias{sn.mmle}
\title{
Modified maximum likelihood estimation for skew-normal models
}
\description{
  Fits a one-dimensional skew-normal (SN) distribution to data, or fits
  a  linear regression model with skew-normal errors, using a
  modification  of maximum likelihood estimation.
}
\usage{
sn.mmle(X, y,  plot.it=TRUE, trace=FALSE, ...)
}
\arguments{
\item{y}{
  a vector contaning the observed variable. This is the response
  variable in case of linear regression.
  Missing values (\code{NA}s) are not allowed.
}
\item{X}{
  a matrix of explanatory variables.
  If \code{X} is missing, then a one-column matrix of all 1's is created.
  If \code{X} has only one column, then it is assumed to be made of 1's.
  Missing values (\code{NA}s) are not allowed.
}
\item{plot.it}{
  logical value, If \code{plot.it=TRUE} (default),
  a plot of the nonparametric estimate of variable \code{y} (or the residuals,
  in the case of regression), and the parametric fit is superimposed.
  See below for details.
}
\item{trace}{
  logical value which controls printing of the algorithm convergence.
  If \code{trace=TRUE}, details are printed. Default value is \code{FALSE}.
}
\item{...}{
  Additional arguments passed to \code{sn.mle}
}
}
\value{
a list containing the following components:

\item{call}{
a string containing the calling statement
}
\item{dp}{
a vector of length \code{ncol(X)+2} with estimates of the direct parameters
}
\item{se}{
a vector of standard errors for the parameters
}
\item{Einfo}{
the expected Fisher information matrix for the parameters
}
}
}
\details{
If plotting operates, the function \code{sm.density} of the library \code{sm}
is searched.  If \code{sm.density} is not found, an histogram is plotted.
}
\section{Background}{
  Maximum likelihood estimation for SN models can lead to estimates
  of the shape parameters on the frontier (that is \code{Inf} for the
  DP parameters); see  Azzalini and Capitanio (1999) for a discussion
  of this aspect.
  To avoid this situation, an alternative estimation criterion is the
  method of Sartori-Firth, which involves first regular maximum estimation
  and subsequent re-estimation of the shape parameter using a modified
  score function; see the references below for a full discussion.
  The effect of this modification is "negligible" for large sample size,
  but it avoids estimates of the frontier of the parameter space.
}
\references{

Azzalini, A. and Capitanio, A. (1999).
  Statistical applications of the multivariate skew-normal distribution.
  \emph{J.Roy.Statist.Soc. B} \bold{61}, 579--602.

Firth, D.  (1993).
  Bias reduction of maximum likelihood estimates.
  \emph{Biometrika}  \bold{80}, 27--38. (Corr: 95V82 p.667).

Sartori, N. (2003).
  Bias prevention of maximum  likelihood estimates: skew normal and skew t
  distributions. Working paper 2003-1,
  Department of Statistical Sciences, University of Padua. Submitted.
}
\seealso{
\code{\link{sn.mle}}, \code{\link{sn.Einfo}}
}
\examples{
data(ais, package="sn")
attach(ais)
a <-  sn.mmle(y=bmi)
#
M <- model.matrix(~lbm+sex)
b <- sn.mmle(M,bmi)
}
\keyword{regression}
\keyword{distribution}

