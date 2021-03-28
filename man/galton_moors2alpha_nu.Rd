%  file sn/man/fournum.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2020 Adelchi Azzalini
%---------------------
\name{galton_moors2alpha_nu}
\alias{galton_moors2alpha_nu}
\concept{quantile}
\concept{asymmetry}
\concept{skewness}
\concept{kurtosis}
% 
\title{Mapping of the (Galton-Bowley, Moors)  measures to the 
      (\eqn{alpha}, \eqn{nu}) parameters of a ST distribution}
      
\description{Given a pair of (Galton-Bowley, Moors) measures of skewness and 
kurtosis for a given sample, a pair of values  (\eqn{alpha}, \eqn{nu}) are
found such that a skew-\emph{t} (\acronym{ST}) distribution with these slant 
and tail-weight  parameter has its (Galton-Bowley, Moors) measures equal to 
the input values.
This function is mainly intended for internal package usage.
}
\usage{
galton_moors2alpha_nu(galton, moors, quick = TRUE, move.in = TRUE, verbose = 0, 
  abstol = 1e-04)
}

\arguments{
  \item{galton}{a numeric value, representing a Galton-Bowley measure}
  \item{moors}{a numeric value, representing a Moors measure}
  \item{quick}{a logical value; if \code{TRUE}, a quick mapping is performed}
  \item{move.in}{if the input values \code{(galton, moors)} are outside the 
     feasible \acronym{ST} region, a suitable point within the feasible area 
     is returned}
  \item{verbose}{a numeric value which regulates the amount of printed detail}
  \item{abstol}{the tolerance value of the mapping, only relevant is
      \code{quick=FALSE}}
}
\details{For background information about the  Galton-Bowley's and 
the Moors measures, see the documentation of \code{\link[sn]{fournum}}.
The working of the mapping by described in Azzalini and Salehi (2020).
}
\value{a named vector of length two, with one or more descriptive attributes}

\references{
Azzalini, A. and Salehi, M. (2020).
Some computational aspects of maximum likelihood estimation of the 
skew-\emph{t} distribution.
In: \emph{Computational and Methodological Statistics and Biostatistics},
edited by  Andriëtte Bekker, Ding-Geng Chen and Johannes T. Ferreira.
Springer. DOI: 10.1007/978-3-030-42196-0
% https://www.springer.com/gp/book/9783030421953#aboutBook
}

\author{Adelchi Azzalini}

\note{This function is mainly intended for internal package usage.
Specifically it is used by \code{\link{st.prelimFit}}.}


\seealso{\code{\link[sn]{fournum}}, \code{\link{st.prelimFit}} }
\examples{
galton_moors2alpha_nu(0.5, 3, quick=FALSE)  # input in the feasible area
galton_moors2alpha_nu(0.5, 3)    # very similar output, much more quickly
galton_moors2alpha_nu(0.5, 0.5)  # input outside the feasible area
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{robust} % use one of  RShowDoc("KEYWORDS")
 
