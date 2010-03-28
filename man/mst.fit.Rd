\name{mst.fit}
\alias{mst.fit}
\title{
Fitting multivariate skew-t distributions
}
\description{
Fits a multivariate skew-t (MST) distribution to data, or fits a
linear regression model with multivariate skew-t errors,
using maximum likelihood estimation. The outcome is then displayed
in graphical form.
}
\usage{
mst.fit(X, y, freq, start, fixed.df=NA,  plot.it=TRUE, trace=FALSE, ...)
}
\arguments{
\item{y}{
a matrix or a vector.  If \code{y} is a matrix, its rows refer to
observations, and its columns to components of the multivariate
distribution. If \code{y} is a vector, it is converted to a one-column
matrix, and a scalar skew-t distribution is fitted.
}
\item{X}{
a matrix of covariate values.
If missing, a one-column matrix of 1's is created; otherwise,
it must have the same number of rows of \code{y}.
}
\item{freq}{
a vector of weights.
If missing, a vector of 1's is created; otherwise
it must have the same number of rows of \code{y}.
}
\item{fixed.df}{
a scalar value containing the degrees of freedom (df), if these must
be taken as fixed, or \code{NA} (default value) if df is a parameter
to be estimated.
}
\item{start}{
a list containing the components \code{beta},\code{Omega}, \code{alpha},
\code{df} of the type described below. The \code{dp} component of the returned
list from a previous call has the required format.
}
\item{plot.it}{
logical value which controls the graphical output (default=TRUE);
see below for description.
}
\item{trace}{
logical value which controls printing of the algorithm convergence.
If \code{trace=TRUE}, details are printed. Default value is \code{FALSE}.
}
\item{...}{
  additional parameters passed to \code{msn.mle}; in practice, the
  \code{start}, the \code{algorithm} and the \code{control} parameters
  can be passed.
}}
\value{
A list containing the following components:

\item{call}{
a string containing the calling statement.
}
\item{dp}{
a list containing the direct parameters \code{beta}, \code{Omega}, \code{alpha},
\code{df}. Here, \code{beta} is a matrix of regression coefficients with
\code{dim(beta)=c(nrow(X),ncol(y))}, \code{Omega} is a covariance matrix of
order \code{ncol(y)}, \code{alpha} is a vector of shape parameters of length
\code{ncol(y)}, \code{df} is a positive scalar.
}
\item{logL}{
log-likelihood evaluated at \code{dp}.
}
\item{se}{
a list containing the components \code{beta}, \code{alpha}, \code{info}.
Here, \code{beta} and \code{alpha} are the standard errors for the
corresponding point estimates;
\code{info} is the observed information matrix for the working parameter,
as explained below.
}
\item{algorithm}{
  see the documentation of \code{mst.mle} for its explanation
}
\item{test.normality}{
a list with elements \code{test} and \code{p.value}, which are the value
of the likelihood ratio test statistic for normality (i.e. test that
 all components of the shape parameter are 0 and \code{df=Inf}), 
and the corresponding p-value.
}
\item{mahalanobis}{
a list of with elements \code{distance}, \code{prob} and \code{df}, which are 
the Mahalanobis distances of the residuals from the origin, with respect
to the metric associated to the matrix \code{Omega}, and the values
\code{prob} of the associated probabilities computed from the Snedecor's F 
distribution with  degrees of freedom given by the \code{df} vector of length 
two, whose first component equals \code{ncol(y)} and the second component is 
equal to the \code{df} parameter of fitted value ST distribution unless this 
value has been selected by the used via \code{fixed.df}.
}}
\section{Side Effects}{
Graphical output is produced if \code{(plot.it & missing(freq))==TRUE}. 
Three plots are produced, and the programs pauses between each two of them, 
waiting for the <Enter> key to be pressed.

The first plot uses the variable \code{y} if \code{X} is missing, otherwise
it uses the residuals from the regression.
The form of this plot depends on the value of \code{d=ncol(y)};
if \code{d=1}, an histogram is plotted with the fitted distribution
superimposed. If \code{d>1}, a matrix of scatter-plots is produced, with
superimposed the corresponding bivariate densities of the fitted
distribution.

The second plot has two panels, each representing a QQ-plot of
Mahalanobis distances. The first of these refers to the fitting of a
multivariate normal distribution, a standard statistical procedure; 
the second panel gives the corresponding QQ-plot of suitable Mahalanobis 
distances for the multivariate skew-normal fit.

The third plot is similar to the previous one, except that PP-plots
are produced.
}
\details{
For computing the maximum likelihood estimates, \code{mst.fit}
invokes \code{mst.mle}, while \code{mst.fit} displays the results
in graphical form.
See the documentation of \code{mst.mle} for details of the numerical
procedure for maximum likelihood estimation.
}
\note{
  This function may be removed in future versions of the package, and
  (some of) its functionality transferred somewhere else
  }
\section{Background}{
The family of multivariate skew-t distributions is an extension of the 
multivariate Student's t family, via the introduction of a \code{shape} 
parameter which regulates skewness; when \code{shape=0}, the skew-t 
distribution reduces to the regular symmetric \emph{t}-distribution. 
When \code{df=Inf} the distribution reduces to the multivariate skew-normal 
one; see \code{dmsn}. See the reference below for additional information.
}
\references{
Azzalini, A. and Capitanio, A. (2003).
  Distributions generated by perturbation of symmetry 
  with emphasis on a multivariate skew \emph{t} distribution.
 \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389.

}
\seealso{
\code{\link{mst.mle}}, \code{\link{msn.fit}}, \code{\link{dmst}}, \code{\link{dmsn}}
}
\examples{
data(ais, package="sn")
attach(ais)
# a simple-sample case
b <- mst.fit(y=cbind(Ht,Wt))
#
# a regression case:
a <- mst.fit(X=cbind(1,Ht,Wt), y=bmi)
#
# refine the previous outcome
a1 <- mst.fit(X=cbind(1,Ht,Wt), y=bmi, start=a$dp)
}
\keyword{distribution}
\keyword{regression}
