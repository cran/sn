%  file sn/man/sn-package.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{sn-package}
\docType{package}
\encoding{UTF-8}
\alias{sn-package}
\alias{SN}
\concept{skew-elliptical distribution}
\concept{skew-normal distribution}
\concept{skew-t distribution}
\concept{symmetric distribution}
\concept{symmetry-modulated distribution}

\title{Package \pkg{sn}: development and brief overview}
\description{
  The \pkg{sn} package  provides facilities to  define and manipulate 
  probability distributions of the skew-normal (\acronym{SN}) family and 
  some related ones, notably the skew-\eqn{t} (\acronym{ST}) family, 
  and to apply connected statistical methods for data fitting and 
  diagnostics, in the univariate and the multivariate case.
}
 
\section{Development and basic facts}{%
The first version of the package was written in 1997 (on CRAN since 1998);
subsequent versions have evolved  gradually up to version 0.4-18 in May 2013.  
In January 2014, version 1.0-0 has been uploaded to CRAN; this represented
a substantial re-writing of the earlier \sQuote{version 0.x}.  
Differences between the \sQuote{version 0} and the \sQuote{version 1} series 
are radical; they concern the core computational and graphical part as well 
as the user interface. 
Since version 1.0-0, the S4 protocol for classes and methods has been adopted.

Broadly speaking, the available tools can be divided in two groups: the
probability section and the statistics section.  For a quick start, one 
could look at their key functions, \code{\link{makeSECdistr}} and
\code{\link{selm}}, respectively, and from here explore the rest.
In the probability section, one finds also functions \code{\link{dsn}},
\code{\link{dst}}, \code{\link{dmsn}} and others alike; these functions 
existed also in \sQuote{version 0} and their working is still very much 
the same (not necessarily so their code). Additional information
on the current package structure and the set of commands is available in a
more extended \code{\link[=pkg-overview]{overview}} of the package.

The first instance of the \sQuote{version 1}  series (that is, 1.0-0) has 
appeared at the same time when the companion book by Azzalini and Capitanio
(2014) was published. 
Although the two projects are formally separate, they adopt the same notation, 
terminology and logical frame. This matching and the numerous references 
in the software documentation to specific sections of the book for background 
information should facilitate familiarizing with these tools.

Information on additional and on more recent change to the package is 
provided in \code{NEWS} file, accessible from the package documentation
index page.

A word of explanation is appropriate about the numerous references to 
Azzalini and Capitanio (2014) in the documentation of the package. 
The reason why the documentation often refers to the monograph rather  
than to the original research papers is because the book provides 
a relatively informal summary of material which has been 
elaborated in a number of technical papers, sometimes very technical 
or with information on the point of interest mixed with other material. 
In other words, the motivation behind this policy is readibility, 
not indulgence in self-citation. When one or a few original sources appeared 
to deliver the required information in a compact and accessible form, 
they have been cited directly.  In any case,  the quoted sections of the
book include bibliographic notes which refer back to the original sources.
}


\section{Backward Compatibility of \sQuote{version 1.x-y}}{%
There is a partial backward compatibility of \sQuote{version 1.x-y} versus
\sQuote{version 0-4.18}. 
Some functions of the older version would work as before
with virtually no change; a wider set arguments is now allowed. Functions
\code{\link{dsn}}, \code{\link{dst}}, \code{\link{dmsn}} and alike fall in
this category: in some cases, the names of the arguments have been altered, 
but they work as before if called with unnamed arguments; similar cases are
\code{\link{msn.mle}}, \code{\link{sn.cumulants}} and \code{\link{T.Owen}}.
Notice, however, that \code{\link{msn.mle}} and other fitting functions have
effectively been subsumed into the more comprehensive fitting function
\code{\link{selm}}.

A second group of functions will work with little or even minimal changes.
Specific examples are functions \code{sn.mle} and \code{st.mle} which have
become \code{\link{sn.mple}} and \code{\link{st.mple}}, with some additional
arguments (again, one can achieve the same result via \code{\link{selm}}). 
Another example is constitude by the group of functions \code{dp.to.cp},
\code{cp.to.dp} and  \code{st.cumulants.inversion}, which have been replaced 
by the more general functions \code{\link{dp2cp}} and \code{\link{cp2dp}}; 
one only needs to pay attention to conversion from 3rd and 4th order 
cumulants  to their standardized form in connection with the replacement of 
\code{st.cumulants.inversion}.

Finally, some functions are not there any longer, with no similarly-working
functions in the new version.  This is the case of \code{sn.mle.grouped}
and \code{st.mle.grouped} for maximum likelihood estimation from
grouped data, that is, data recorded as intervals and corresponding 
frequencies.}


\section{Requirements}{
  \R version 2.15-3 or higher, plus  packages \pkg{mnormt}, 
  \pkg{numDeriv}, \pkg{stats4}   in addition to standard
  packages (\pkg{methods}, \pkg{graphics},  etc.)
}

\section{Version}{
The command \code{citation("sn")} indicates, among other information,
the running version of the package.
The most recent version of the package can be obtained from
the web page: \url{http://azzalini.stat.unipd.it/SN}
which also provides related material.

From the above-indicated web page, one can also obtain the package 
\sQuote{sn0} which is essentially the last \sQuote{version 0} 
(that is, 0.4-18) with suitable renaming of  certain  ingredients.
This allows to have both the current and the old package installed 
at the same time.
}

\section{Author}{Adelchi Azzalini.
% Dipart. Scienze Statistiche, Università di Padova, Italia.
Please send comments, error reports \emph{et cetera} to the author, 
whose web page is \url{http://azzalini.stat.unipd.it/}.
}

\section{Licence}{
This package and its documentation are usable under the terms of 
the \dQuote{GNU General Public License} version 3 or version 2,
as you prefer; a copy of them is available from 
\url{https://www.R-project.org/Licenses/}.
  
While the software is freely usable, it would be appreciated
if a reference is inserted in publications or other work
which makes use of it. For the appropriate way of referencing it, 
see the command \code{citation("sn")}. 
}

\seealso{\code{\link[=pkg-overview]{package-overview}} }

\references{
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.
}

\keyword{multivariate}
\keyword{distribution}
\keyword{univar}
\keyword{regression}
