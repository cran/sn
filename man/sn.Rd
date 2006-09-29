\name{SN}
\alias{SN}
\alias{sn-package}
\title{Package `sn': summary information}
\description{
This package provides functions related to the skew-normal (SN) and
the skew-t (ST) probability distributions, both for the univariate 
and for the the multivariate case, including regression models.
}
\section{Functions}{
% The functions of the scalar case section are:
% \code{dsn}, \code{psn}, \code{qsn}, \code{rsn}, \code{T.Owen}, 
% \code{cp.to.dp}, \code{dp.to.cp}, \code{zeta},
% \code{gamma1.to.lambda}, \code{sn.cumulants}, \code{sn.em}, 
% \code{sn.2logL.profile}, \code{sn.mle}, \code{sn.dev}, 
% \code{sn.dev.gh}.
% 
% 
% The functions of the multivariate section are: \code{dmsn},
% \code{rmsn}, \code{plot.dsn2}, \code{msn.quantities},
% \code{msn.conditional}, \code{msn.marginal}, \code{plot.msn.cond},
% \code{msn.fit}, \code{msn.mle}, \code{msn.dev}, \code{msn.dev.grad},
% \code{msn.moment.fit}, \code{num.deriv}.

The package includes several functions, a number of which are intended
as services to other functions, not really to be called directly by the
users. The following ones are those more relevant for practical use. 

\emph{SN distribution:}
\code{\link{dsn}}, \code{\link{psn}}, \code{\link{qsn}},
\code{\link{rsn}}, \code{\link{cp.to.dp}}, \code{\link{dp.to.cp}},
\code{\link{sn.mle}}, \code{\link{sn.em}}, 
%\code{\link{sn.mle.grouped}} 
for the univariate case, and 
\code{\link{dmsn}}, \code{\link{pmsn}}, \code{\link{rmsn}}, 
\code{\link{msn.fit}}, \code{\link{msn.mle}}, \code{\link{dsn2.plot}} 
for the multivariate case. 


\emph{ST distribution:}
\code{\link{dst}}, \code{\link{pst}}, \code{\link{rst}}, \code{\link{st.mle}}, 
%\code{\link{st.mle.grouped}} 
for the univariate case,
and \code{\link{dmst}},\code{\link{pmst}}, \code{\link{rmst}}, 
\code{\link{mst.fit}}, \code{\link{mst.mle}}, \code{\link{dst2.plot}} 
for the multivariate case. 

It is suggested that a user starts by reading the  documentation of 
(some of) these functions.
}
\section{Requirements}{
  R 2.2.0.
  Two functions make use of package \code{mnormt}: \code{pmsn}, \code{pmst}.
}
\section{Version}{
The version level of the package is given by the  command 
\code{print(.sn.version)}.
The most recent version of the package can be obtained from
the WWW page: \url{http://azzalini.stat.unipd.it/SN}
which also contains other related material
}
\section{Author}{
Adelchi Azzalini, Dipart. Scienze Statistiche, Università di Padova, Italia.
Please send comments, error reports, etc. to the author whose WWW page
is \url{http://azzalini.stat.unipd.it/}.
}
\section{Licence}{
This package and its documentation are usable under the terms of the 
"GNU General Public License", a copy of which is distributed with the 
package. While the software is freely usable, it would be appreciated
if a reference is inserted in publications or other work
which makes use of it; for this purpose, see the command
\code{citation("sn")}. 
}
\section{Acknowledgements}{The package has evolved through several 
versions, developed over some years.

For versions up to 0.20, the following people and institutions 
have contributed. Many thanks go to Antonella Capitanio for testing 
the procedures, and to Brian Ripley, Kurt Hornik and  Martin Maechler
for useful  advice on R.
The function \code{num.deriv2} is based on a similar 
function written  by Monica Chiogna. The first version of this software
and part of the associated theoretical work has been developed 
while the author was at the Nuffield College, Oxford, under the Jemolo 
Fellowship scheme; the generous support of the  college is gratefully 
acknowledged. 
Additional support for the development of the theoretical research work has 
been  provided by the "Consiglio Nazionale delle Ricerche" of Italy, grant 
no.97.01331.CT10.

Versions 0.21 to 0.30 have been supported by "MIUR", Italy, under grant 
scheme PRIN 2000. Again, thanks to Antonella Capitanio for additional
testing of those versions of the package.
}
\keyword{multivariate}
\keyword{distribution}
\keyword{regression}
