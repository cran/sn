\name{msn.fit}
\alias{msn.fit}
\title{
Fitting multivariate skew-normal distributions
}
\description{
Fits a multivariate skew-normal (MSN) distribution to data, or fits a
linear regression model with multivariate skew-normal errors,
using maximum likelihood estimation. The outcome is then displayed
in graphical form.
}
\usage{
msn.fit(X, y, freq, plot.it=TRUE, trace=FALSE,  \dots )
}
\arguments{
\item{y}{
a matrix or a vector.  If \code{y} is a matrix, its rows refer to
observations, and its columns to components of the multivariate
distribution. If \code{y} is a vector, it is converted to a one-column
matrix, and a scalar skew-normal distribution is fitted.
}
\item{X}{
a matrix of covariate values.
If missing, a one-column matrix of 1's is created; otherwise,
it must have the same number of rows of \code{y}.
}
\item{freq}{
a vector of weights.
If missing, a one-column matrix of 1's is created; otherwise
it must have the same number of rows of \code{y}.
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
additional paramters passed to \code{msn.mle}; in practice, the \code{start}
parameter can be passed.
}}
\value{
A list containing the following components:

\item{call}{
a string containing the calling statement.
}
\item{dp}{
a list containing the direct parameters \code{beta}, \code{Omega}, \code{alpha}.
Here, \code{beta} is a matrix of regression coefficients with
\code{dim(beta)=c(nrow(X),ncol(y))}, \code{Omega} is a covariance matrix of
order \code{ncol(y)}, \code{alpha} is a vector of shape parameters of length
\code{ncol(y)}.
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
\item{optim}{
the list returned by the optimizer \code{optim}; see the documentation
of this function for explanation of its components.
}
\item{test.normality}{
a list of with elements \code{test} and \code{p.value}, which are the value
of the likelihood ratio test statistic for normality (i.e. test that
all components of the shape parameter are 0), and the corresponging
p-value.
}}
\section{Side Effects}{
Graphical output is produced if \code{(plot.it \& missing(freq))=TRUE} and
a suitable device is active. Three plots are produced, and the programs
pauses between each two of them, waiting for the <Enter> key to be pressed.

The first plot uses the variable \code{y} if \code{X} is missing, otherwise
it uses the residuals from the regression.
The form of this plot depends on the value of \code{k=ncol(y)};
if \code{k=1}, an histogram is plotted with the fitted distribution
superimposed. If \code{k>1}, a matrix of scatterplots is produced, with
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
For computing the maximum likelihood estimates, \code{msn.fit}
invokes \code{msn.mle} which does the actual computational work;
then, \code{msn.fit} displays the results in graphical form.
The documentation of \code{msn.mle} gives details of the numerical
procedure for maximum likelihood estimation.


Although the function accepts a vector \code{y} as input, the use of
\code{sn.mle} is recommended in the scalar case.


}
\section{Background}{
The multivariate skew-normal distribution is discussed by
Azzalini and Dalla Valle (1996); the \code{(Omega,alpha)} parametrization
adopted here is the one of Azzalini and Capitanio (1999).
}
\references{
Azzalini, A. and Dalla Valle, A. (1996).
The multivariate skew-normal distribution.
\emph{Biometrika}
\bold{83}, 715--726.


Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{msn.mle}}, \code{\link{mst.fit}}
}
\examples{
data(ais,package=sn)
attach(ais)
# a simple-sample case
b <- msn.fit(y=cbind(Ht,Wt))
#
# a regression case:
a <- msn.fit(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-6))
#
# refine the previous outcome
a1 <- msn.fit(X=cbind(1,Ht,Wt), y=bmi, control=list(x.tol=1e-9), start=a$dp)
}
\keyword{distribution}
\keyword{regression}
