%  file sn/man/SECdistrMv-class.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{SECdistrMv-class}
\Rdversion{1.1}
\docType{class}
\alias{SECdistrMv-class}
\alias{show,SECdistrMv-method}

\title{Class \code{"SECdistrMv"}}

\description{Multivariate skew-elliptically contoured distributions}

\section{Objects from the Class}{
  Objects can be created by a call to function \code{\link{makeSECdistr}},
  when its argument \code{dp} is a list, or by a suitable transformation of 
  some object of this class.}

\section{Slots}{
  \describe{
    \item{\code{family}:}{a character string which identifies the parametric
       family; currently, possible values are: \kbd{"SN"}, \kbd{"ESN"}, 
       \kbd{"ST"}, \kbd{"SC"}.}
    \item{\code{dp}:}{a list of parameters; its length depends on 
       the selected  \code{family}.}
    \item{\code{name}:}{a character string with the name of the multivariate
        variable; it can be an empty string.}
    \item{\code{compNames}:}{a vector of character strings with the names of 
       the component variables.}
    }
  }
    
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SECdistrMv-class")}: \dots }
    \item{plot}{\code{signature(x = "SECdistrMv-class")}: \dots }
    \item{summary}{\code{signature(object = "SECdistrMv-class")}: \dots }
    }
}

\author{Adelchi Azzalini}

\note{See \code{\link{makeSECdistr}} for a detailed description of 
  \code{family} and \code{dp}. If an object of this class is manipulated by
  \code{affineTransSECdistr} or \code{marginalSECdistr}, the returned object 
  is of the same class, unless the transformation leads to a univariate 
  distribution.}

\seealso{
  \code{\linkS4class{SECdistrUv}},  
  \code{\link{plot,SECdistrMv-method}}, 
  \code{\link{summary,SECdistrMv-method}},
  
  \code{\link{affineTransSECdistr}},  \code{\link{marginalSECdistr}}
}
\examples{
  dp0 <- list(xi=1:2, Omega=diag(3:4), alpha=c(3, -5))
  f10 <- makeSECdistr(dp=dp0, family="SN", name="SN-2D", compNames=c("x", "y"))
  show(f10)
  plot(f10)
  summary(f10)
}
\keyword{classes}
