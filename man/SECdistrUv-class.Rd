%  file sn/man/SECdistrUv-class.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{SECdistrUv-class}
\Rdversion{1.1}
\docType{class}
\alias{SECdistrUv-class}
\alias{show,SECdistrUv-method}

\title{Class \code{"SECdistrUv"}}

\description{Univariate skew-elliptically contoured distributions}

\section{Objects from the class}{
Objects can be created by a call to function \code{\link{makeSECdistr}} 
when its argument \code{dp} is a vector.}

\section{Slots}{
  \describe{
    \item{\code{family}:}{a character string which selects the parametric
       family; currently, possible values are: \kbd{"SN"}, \kbd{"ESN"}, 
       \kbd{"ST"}, \kbd{"SC"}. }
    \item{\code{dp}:}{a numeric vector of parameters; its length depends 
       on the selected \code{family}.}
    \item{\code{name}:}{a character string with name of the distribution.}
  }
}

\section{Methods}{
  \describe{   
    \item{show}{\code{signature(object = "SECdistrUv")}: \dots}
    \item{plot}{\code{signature(x = "SECdistrUv")}: \dots }
    \item{summary}{\code{signature(object = "SECdistrUv")}: \dots}
    }
}

\author{Adelchi Azzalini}

\note{See \code{\link{makeSECdistr}} for a detailed description of 
\code{family} and \code{dp}.}

\seealso{
  \code{\linkS4class{SECdistrMv}},
  \code{\link{plot,SECdistrUv-method}}, 
  \code{\link{summary,SECdistrUv-method}}  
}


\examples{
f2 <- makeSECdistr(dp=c(3, 5, -4, 6.5), family="ST", name="My first ST")
show(f2)
plot(f2)
plot(f2, probs=c(1,5,9)/10)
plot(f2, range=c(-30,10), probs=NULL, col=2, main=NULL)
summary(f2)
}
\keyword{classes}
