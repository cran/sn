%  file sn/man/dp2cp.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{dp2cp}
\alias{dp2cp}
\alias{cp2dp}
\title{Conversion between parametrizations of a skew-elliptical distribution}

\description{
  Convert direct parameters (\acronym{DP}) to centred parameters 
 (\acronym{CP}) of a skew-elliptical distribution and \emph{vice versa}.}
  
\usage{
dp2cp(dp, family, obj = NULL, cp.type = "proper", upto = NULL) 
cp2dp(cp, family)
}
 
\arguments{

\item{dp}{a vector (in the univariate case) or a list (in the multivariate
  case) as described in \code{\link{makeSECdistr}}; see \sQuote{Details} 
  for an extented form of usage.}

\item{cp}{a vector or a list, in agreement with \code{dp} as for type and
  dimension.}

\item{family}{a characther string, as described in \code{\link{makeSECdistr}}.}

\item{obj}{optionally, an S4 object of class \code{SECdistrUv} or
  \code{SECdistrMv}, as produced by \code{\link{makeSECdistr}} 
  (default value: \code{NULL}). 
  If this argument is not \code{NULL}, then \code{family} and \code{dp} 
  must not be set.}

\item{cp.type}{character string, which has effect only if \code{family="ST"} 
   or \code{"SC"}, otherwise a warning message is generated.  Possible values 
   are \kbd{"proper", "pseudo", "auto"}, which correspond to the \acronym{CP} 
   parameter set, their `pseudo-\acronym{CP}' version and an automatic
   selection based on \code{nu>4}, where \code{nu} represents the degrees of
   freedom of the \acronym{ST} distribution.}

\item{upto}{numeric value (in \code{1:length(dp)}, default=\code{NULL}) to 
  select how many \acronym{CP} components are computed. 
  Default value \code{upto=NULL} is equivalent to \code{length(dp)}.}
 } 

\value{for \code{dp2cp}, a matching vector (in the univariate case) or a list 
   (in the multivariate case) of \code{cp} parameters; for  \code{cp2dp}, 
   a similar object of \code{dp} parameters.}			
			
\section{Details and Background}{For a description of the \acronym{DP}
parameters, see Section \sQuote{Details} of \code{\link{makeSECdistr}}.  The
\acronym{CP} form of parameterization is cumulant-based. For a univariate
distribution, the \acronym{CP} components are the mean value (first cumulant), 
the standard deviation (square root of the 2nd cumulant), the coefficient of 
skewness (3rd standardized cumulant) and,  for the \acronym{ST}, 
the coefficient of excess kurtosis (4th standardized cumulant). 
For a multivariate distribution, there exists an extension based on the \
same logic; its components represent the
vector mean value, the variance matrix, the vector of marginal coefficients of
skewness and, only for the \acronym{ST}, the Mardia's coefficient of excess
kurtosis. The pseudo-\acronym{CP} variant provides an `approximate form' of
\acronym{CP} when not all required cumulants exist; however, this parameter set
is not uniquely invertible to \acronym{DP}. The names of pseudo-\acronym{CP}
components printed in summary output are composed by adding a \code{~} 
after the usual component name; for example, the first one is denoted
\code{mean~}.

Background information is provided by Azzalini and Capitanio (2014).
Specifically, their Section 3.1.4 presents \acronym{CP} in the univariate 
\acronym{SN} case, Section 4.3.4 \acronym{CP} for the \acronym{ST} case and 
the `pseudo-\acronym{CP}' version. Section 5.2.3 presents the multivariate
extension for the \acronym{SN} distribution, Section 6.2.5 for the 
multivariate \acronym{ST} case. 
For a more detailed discussion, see Arellano-Valle and Azzalini (2013).

It is possible to call the functions with \code{dp} or \code{cp} having more
components than those expected for a given family as described above and in
\code{\link{makeSECdistr}}. In the univariate case, this means that \code{dp}
or \code{cp} can be vectors of longer length than indicated earlier. This
occurrence is interpreted in the sense that the additional components after
the first one are regarded as regression coefficients of a \code{selm} model,
and they are transferred unchanged to the matching components of the
transformed parameter set; the motivation is given in Section 3.1.4 of
Azzalini and Capitanio (2014). In the multivariate case, \code{dp[[1]]} and
\code{cp[[1]]} can be matrices instead of vectors; the rows beyond the first
one are transferred unchanged to \code{cp[[1]]} and \code{dp[[1]]},
respectively. }


\references{
Arellano-Valle, R. B. and Azzalini, A. (2013, available on-line 12 June 2011). 
The centred parameterization and related quantities of the skew-\emph{t} 
distribution. \emph{J. Multiv. Analysis} \bold{113}, 73-90. 

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
\emph{The Skew-Normal and Related  Families}. 
Cambridge University Press, IMS Monographs series.
}

\seealso{
  \code{\link{makeSECdistr}}, \code{\link{summary.SECdistr}}, 
  \code{\link{sn.cumulants}}, the \sQuote{Note} at \code{\link{summary.selm}} 
}

\examples{
# univariate case
cp <- dp2cp(c(1, 2222, 3333, 2, 3), "SN")
dp <- cp2dp(cp, "SN")
# notice that 2nd and 3rd component remain unchanged
#
# multivariate case
dp3 <- list(xi=1:3, Omega=toeplitz(1/(1:3)), alpha=c(-3, 8, 5), nu=6)
cp3 <- dp2cp(dp3, "ST")
dp3.back <- cp2dp(cp3, "ST")
}

\keyword{distribution}
