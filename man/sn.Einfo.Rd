\name{sn.Einfo}
\alias{sn.Einfo}
\title{
Expected Fisher information for SN distribution parameters
}
\description{
  Computes expected Fisher information for parameters of simple
  sample having one-dimensional skew-normal (SN) distribution
  or regression model having SN errors, in the DP and CP parametrizations.
}
\usage{
sn.Einfo(dp=NULL, cp=NULL, n=1, x=NULL)
}
\arguments{
\item{dp, cp}{vector of direct or centred parameters, respectively;
one of the two vectors must be supplied, but not both. See below for
more details.
}
\item{n}{sample size; if this parameter is supplied, then \code{x} must
not be.
}
\item{x}{design matrix of the regression model; if this parameter is
supplied, then \code{n} must not be.
}}
\value{
a list containing the following components:
\item{dp, cp}{
DP and CP parameters; one of the two vectors is the one supplied on
input, the other one  matches the previous one in the alternative
parametrization
}

\item{info.dp, info.cp}{
matrices of Fisher expected information in the two parametrizations
}

\item{se.dp, se.cp}{
vectors of standard errors in the two parametrizations
}

\item{D}{
matrix of derivatives of DP parameters with respect to CP parameters
}
}
\section{Details}{
When \code{x} is not set, then a simple random sample is assumed and a
 matrix \code{x} with a single column of all 1's is constructed; in this
case, the supplied  vector \code{dp} or \code{cp} must have length 3.
If \code{x} is set, then the supplied vector of parameters must have
length \code{ncol(x)+2}.
}
\section{Background}{
  For the description of the DP and CP parametrizations and
  for the expression of the exptected Fisher information of the DP
  parametrizations, see Azzalini (1985). Further discussion is in
  Azzalini and Capitanio (1999).
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
\code{\link{dsn}}, \code{\link{cp.to.dp}}, \code{\link{dp.to.cp}}
}
\examples{
info <- sn.Einfo(dp=c(0,1,5), n=3)
#
data(ais, package="sn")
M <- model.matrix(~ais$"Ht")
mle <- sn.mle(X=M, y=ais$"Wt", plot.it=FALSE)
info <- sn.Einfo(cp=mle$cp, x=M)
}
\keyword{distribution}
