%  file sn/man/summary.SECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{summary.SECdistr}
\alias{summary.SECdistr}
\alias{summary.SECdistrUv}
\alias{summary.SECdistrMv}
\alias{summary,SECdistrUv-method}
\alias{summary,SECdistrMv-method}

\title{Summary of a \acronym{SEC} distribution object}

\description{Produce a summary of an object of class either
\code{"SECdistrUv"} or \code{"SECdistrMv"}, which refer to a univariate or a
multivariate \acronym{SEC} distribution, respectively. Both types of
objects are produced by \code{makeSECditr}. }

\usage{
\S4method{summary}{SECdistrUv}(object, cp.type = "auto", probs)

\S4method{summary}{SECdistrMv}(object, cp.type = "auto")
}

\arguments{
  \item{object}{an object of class \code{"SECdistrUv"} or \code{"SECdistrMv"}.}

  \item{cp.type}{a character string to select the required variance of
    \acronym{CP} parameterization; possible values are \code{"proper"},
    \code{"pseudo"}, \code{"auto"} (default). For a description of these
    codes, see \code{\link{dp2cp}}.}
  
 \item{probs}{in the univariate case, a vector of probabilities for which 
    the corresponding quantiles are required.}
}

%\details{%%  ~~ If necessary, more details than the description above ~~}

\value{A list with the following components:
  \item{family}{name of the family within the \acronym{SEC} class, character}
  \item{dp}{\acronym{DP} parameters, list or vector}
  \item{name}{the name of the distribution, character string}
  \item{compNames}{in the multivariate case, the names of the components}
  \item{cp}{\acronym{CP} parameters, list or vector}
  \item{cp.type}{the name of the selected variant of the \acronym{CP} set}
  \item{aux}{a list with auxiliary ingredients (mode, coefficients of 
       skewness and kurtosis, in the parametric and non-parametric variants, 
       and more).}
}
 
\author{Adelchi Azzalini}

%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{makeSECdistr}} }

\examples{
f3 <- makeSECdistr(dp=c(3,2,5), family="SC")
summary(f3)
print(s <- summary(f3, probs=(1:9)/10))
print(slot(s, "aux")$mode)
#
dp3 <- list(xi=1:3, Omega=toeplitz(1/(1:3)), alpha=c(-3, 8, 5), nu=6)
st3 <- makeSECdistr(dp=dp3, family="ST", compNames=c("U", "V", "W"))
s <- summary(st3)
print(slot(s, "aux")$delta.star)
print(slot(s, "aux")$mardia)
#
dp2 <- list(xi=rep(0,2), Omega=matrix(c(2,2,2,4),2,2), alpha=c(3,-5), tau=-1)
esn2 <- makeSECdistr(dp=dp2, family="ESN", name="ESN-2d")
summary(esn2)
}
\keyword{multivariate}
\keyword{distribution} 