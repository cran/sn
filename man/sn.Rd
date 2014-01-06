%  file sn/man/sn.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{SN}
\docType{package}
\encoding{UTF-8}
\alias{SN}
\alias{sn-package}
\concept{skew-elliptical distribution}
\concept{skew-normal distribution}

\title{Package \pkg{sn}: overview}
\description{
  The \pkg{sn} package  provides facilities to  define and manipulate 
  probability distributions of the skew-normal (\acronym{SN}) family and 
  some related ones, notably the skew-\eqn{t} (\acronym{ST}) family, 
  and to apply connected statistical methods for data fitting and 
  diagnostics, in the univariate and the multivariate case.
}
 
\section{A substantial upgrade}{The first version of the package has been
written in 1997 (on CRAN since 1998); subsequent versions have evolved 
gradually up to version 0.4-18 in May 2013.  The present \sQuote{version 1} 
of the package is a substantial re-writing of the earlier \sQuote{version 0}.  
Differences between \sQuote{version 0} and \sQuote{version 1} concern the core
computational and graphical part as well as the user interface. 
The S4 protocol for classes and methods has been adopted.

Broadly speaking, the available tools can be divided in two groups: the
probability section and the statistics section.  For a quick start, one 
could look at their key functions, \code{\link{makeSECdistr}} and
\code{\link{selm}}, respectively, and from here explore the rest.
In the probability section, one finds also functions \code{\link{dsn}},
\code{\link{dst}}, \code{\link{dmsn}} and others alike; these functions 
existed also in \sQuote{version 0} and their working is still very much 
the same (not necessarily so their code).

The upgrade to \sQuote{version 1} appears more or less at the time when the
companion book by Azzalini and Capitanio (2014) is published. Although the two
projects are formally separate, they adopt the same notation, terminology 
and logical frame. This matching and the numerous references in the software 
documentation to specific sections of the book for background information 
should facilitate familiarizing with these tools.}

\section{Backward Compatibility}{%
 
There is a partial backward compatibility of \sQuote{version 1} versus
\sQuote{version 0}. Some functions of the older version would work as before
with virtually no change; a wider set arguments is now allowed. Functions
\code{\link{dsn}}, \code{\link{dst}}, \code{\link{dmsn}} and alike fall in
this category: typically, the names of the arguments have been altered, but
they work as before if called with unnamed arguments; similar cases are
\code{\link{msn.mle}}, \code{\link{sn.cumulants}} and \code{\link{T.Owen}}.
Notice, however, that \code{\link{msn.mle}} and other fitting functions have
effectively been subsumed into the more general fitting function
\code{\link{selm}}.

A second group of functions will work with little or even minimal changes.
Specific examples are functions \code{sn.mle} and \code{sn.mle} which have
become \code{\link{sn.mple}} and \code{\link{st.mple}}, with some additional
arguments (again, one can achieve the same result via \code{\link{selm}}) and
\code{dp.to.cp}, which has been replaced by the more general function
\code{\link{dp2cp}}.

Finally, some functions are not there any longer, with no similarly-working
functions in the new version. The more prominent set of cases is represented
by the functions for computing profile log-likelihoods. There is a long-term
plan to re-instate similar facilities, possibly in a more flexible form, but
not in the near future. }


\section{Requirements}{
  \R version 2.15-3 or higher, plus  packages \pkg{mnormt}, 
  \pkg{numDeriv}, \pkg{stats4} and \pkg{methods} in addition to `standard' 
  packages (\pkg{graphics},  etc.)
}

\section{Version}{
The command \code{citation("sn")} indicates, among other information,
the running version of the package.
The most recent version of the package can be obtained from
the WWW page: \url{http://azzalini.stat.unipd.it/SN}
which also contains other related material.
}

\section{Author}{Adelchi Azzalini.
% Dipart. Scienze Statistiche, Università di Padova, Italia.
Please send comments, error reports, etc. to the author whose web page
is \url{http://azzalini.stat.unipd.it/}.
}

\section{Licence}{
This package and its documentation are usable under the terms of 
the \dQuote{GNU General Public License} version 3 or version 2,
as you prefer; a copy of them is available from 
\url{http://www.R-project.org/Licenses/}.
  
While the software is freely usable, it would be appreciated
if a reference is inserted in publications or other work
which makes use of it; for this purpose, see the command
\code{citation("sn")}. 
}

\references{
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.
}

\keyword{multivariate}
\keyword{distribution}
\keyword{univar}
\keyword{regression}
