\name{sn.mle}
\alias{sn.mle}
\title{
Maximum likelihood estimation for skew-normal models
}
\description{
Fits a skew-normal (SN) distribution to data, or fits a linear regression
model with skew-normal errors, using maximum likelihood estimation.
}
\usage{
sn.mle(X, y, cp, plot.it=TRUE, trace=FALSE, method="L-BFGS-B",
          control=list(iter.max=100, abs.tol=1e-5))
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
  If \code{X} is supplied, then it must include a column of 1's.
  Missing values (\code{NA}s) are not allowed.
}
\item{cp}{
a vector of initial values for the centred parameters,
with \code{length(cp)=ncol(X)+2}
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
\item{method}{
this parameter is just passed to the optimizer \code{optim}; see the 
documentation of this function for its usage. Default value is 
\code{"L-BFGS-B"}.}
\item{control}{
this parameter is just passed to the optimizer \code{optim}; 
see the documentation of this function for its usage.
}}
\value{
a list containing the following components:

\item{call}{
a string containing the calling statement
}
\item{cp}{
a vector of length \code{ncol(X)+2} with the centred parameters
}
\item{logL}{
the log-likelihood at convergence
}
\item{se}{
a vector of standard errors for the \code{cp} component
}
\item{info}{
the observed information matrix for the \code{cp} component
}
\item{optim}{
the list returned by the optimizer \code{optim}; see the documentation
of this function for explanation of its components.
}}
\section{Side Effects}{
If \code{plot.it=TRUE} and a graphical device is active, a plot is produced,
as described above.
}
\details{
The optimizer \code{optim} is used, supplying the gradient of the log-likelihood.
Convergence is generally fast and reliable, but inspection of
the returned \code{message} from \code{optim} is always appropriate.
In suspect cases, re-run the function changing the starting \code{cp}
vector.

If plotting operates, the function \code{sm.density} of the library \code{sm}
is searched; this library is associated with the book by Bowman and
Azzalini (1997).  If \code{sm.density} is not found, an histogram is plotted.

% To fit a skew-notmal distribution to grouped data by exact
% maximum likelihood estimation, use \code{sn.mle.grouped}.

}
\section{Background}{
Background information on the SN distribution is given by Azzalini (1985).
See  Azzalini and Capitanio (1999) for a more detailed discussion of
the centred parametrization.
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


Bowman, A.W. and Azzalini, A. (1997).
\emph{Applied Smoothing Techniques for Data Analysis:}
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{dsn}}, \code{\link{sn.em}}, \code{\link{msn.mle}}, \code{\link{optim}}
% \code{\link{sn.mle.grouped}}
}
\examples{
data(ais,package=sn)
attach(ais)
a<-sn.mle(y=bmi)
#
a<-sn.mle(X=cbind(1,lbm),y=bmi)
#
b<-sn.mle(X=model.matrix(~lbm+sex), y=bmi)
}
\keyword{regression}
\keyword{distribution}
