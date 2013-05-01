\name{sn.2logL.profile}
\alias{sn.2logL.profile}
\title{
Twice profile relative negative loglikelihood for skew-normal models
}
\description{
Computation and plot of 1-dimensional and 2-dimensional profile relative
(-2)*loglikelihood for skew-normal regression models.
}
\usage{
sn.2logL.profile(X=matrix(rep(1, n)), y,
    param.range=c(sqrt(var(y)) * c(2/3, 3/2), -0.95, 0.95),
    use.cp=TRUE, npts=51 \%/\% d, plot.it=TRUE, ...)
}
\arguments{
\item{y}{
a numeric vector. Missing values (\code{NA}'s) are not allowed.
}
\item{X}{
a matrix of explanatory variables. It must have \code{col(X)} equal to
\code{length(y)}, and it must include a column of 1's if an intercept
term is required.  Missing values (\code{NA}'s) are not allowed.
If \code{X} is missing, a one-column matrix of 1's is created.
}
\item{param.range}{
a numeric vector of length either 2 or 4. If the length is 2,
the dimensional value \code{d} is set to 1, and
a 1-dimensional profile is computed and plotted, for the shape
or skewness parameter (depending on the parametrization adopted;
see below); in this case the two value represent the minimum
and maximum value for the range of the parameter.
If the length of \code{param.range} is 4, the first two values
determine the range of the scale parameter, the last two give
the range of the shape (or skewness) parameter; in this case, \code{d=2}.
}
\item{use.cp}{
logical value which selects the parametrization adopted.
If \code{use.cp=TRUE} (default value), the centred parametrization is used,
otherwise the direct parametrization is adopted.
}
\item{npts}{
number of points (in the scalar case) or grid size (in the two-dimensional 
case).
}
\item{plot.it}{
logical value which determines if plotting takes place; default is \code{TRUE}.
}
\item{...}{
any additional parameter is passed to \code{sn.em}.
}}
\value{
A list containing the following components

\item{param1}{
vectors of the parameters values where the function has been
evaluated. If \code{d=2}, the second vector contains \code{NA}s.
}
\item{param.names}{
a character vector of two elements with the names of the \code{param1}
and \code{param2}.
}
\item{2logL}{
a vector or a matrix which represents the profile (-2)*loglikelihood;
this is in the "relative" version, i.e. setting the maximum value to be 0.
}
\item{maximum}{
a numeric value with the maximum which has been subtracted to
obtain the "relative" version of \code{2logL}.
}}
\section{Side Effects}{
If \code{plot.it=TRUE}, a plot of the profile twice relative negative 
loglikeliood (called the `deviance') is produced on a graphical device. 
When \code{length(fixed.comp)=1}, a plot of the deviance is produced
as a function of the chosen parameter component.
When \code{length(fixed.comp)=2}, a contour plot of the deviance is produced
with contour lines corresponding to confidence regions of approximate
probability levels \code{c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99)}. 
}
\details{
Likelihood maximization is performed by \code{sn.em}.

See the reference below for explanation of the two possible parametrizations.
}
\references{
Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\seealso{
\code{\link{sn.em}}, \code{\link{sn.mle}}
}
\examples{
data(ais, package="sn")
attach(ais)
a <- sn.2logL.profile(y=bmi)
\dontrun{
a <- sn.2logL.profile(y=bmi, use.cp=FALSE, param.range=c(3,6,1,5))
a <- sn.2logL.profile(X=cbind(1,lbm), y=bmi, param.range=c(0.5,0.95), npts=31)
#
data(frontier, package="sn")
a <- sn.2logL.profile(y=frontier, param.range=c(0.8,2, 2,30),
        use.cp=FALSE, npts=16)
	}
}
\keyword{distribution}
 
