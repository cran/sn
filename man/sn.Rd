\name{SN}
\alias{SN}
\title{
The library `sn': summary information
}
\description{
This library provides functions related to the skew-normal probability
distribution, both for the univariate and the the multivariate case.
}
\section{Functions}{
The functions of the scalar case section are:
\code{dsn}, \code{psn}, \code{qsn}, \code{rsn}, \code{T.Owen}, \code{cp.to.dp}, 
\code{dp.to.cp}, \code{zeta},
\code{gamma1.to.lambda}, \code{sn.cumulants}, \code{sn.em}, \code{sn.2logL.profile},
\code{sn.mle}, \code{sn.dev}, \code{sn.dev.gh}.


The functions of the multivariate section are:
\code{dmsn}, \code{rmsn}, \code{plot.dsn2}, \code{msn.quantities}, 
\code{msn.conditional},
\code{msn.marginal}, \code{plot.msn.cond}, \code{msn.fit}, \code{msn.mle}, 
\code{msn.dev}, \code{msn.dev.grad}, \code{msn.moment.fit}, \code{num.deriv}.

Some of these functions are there to be used mostly by other functions
of the library.  Those more relevant for practical use by a user are:
\code{\link{dsn}}, \code{\link{psn}}, \code{\link{qsn}},
\code{\link{rsn}}, \code{\link{cp.to.dp}}, \code{\link{dp.to.cp}},
\code{\link{sn.mle}},\code{\link{sn.em}} in the univariate section,
and \code{\link{dmsn}}, \code{\link{rmsn}}, \code{\link{msn.fit}},
\code{\link{msn.mle}}, \code{\link{plot.dsn2}} in the multivariate
section. It is suggested that a user starts by reading the documentation
of (some of) these.
}
\section{Requirements}{
R 1.0.1
}
\section{Version}{
You are using version 0.22.1 (2001-05-17).
The most recent version of the library can be obtained from
the WWW page: \url{http://azzalini.stat.unipd.it/SN}
which also contains other related material
}
\section{Manual}{
There is no manual  besides on-line documentation and its derived versions. 
The papers listed below provide background information. 
}
\section{Author}{
Adelchi Azzalini, Dipart. Scienze Statistiche, Università di Padova, Italia.
Please send comments, error reports, etc. to \email{azzalini@stat.unipd.it}
}
\section{Acknowledgements}{
Many thanks go to Antonella Capitanio for testing the procedures,
and to Brian Ripley and Martin Maechler for useful advice on R. 
The function \code{num.deriv} is based on a similar function written 
by Monica Chiogna. 
This software and part of the associated theoretical work has been developed 
while the author was at the Nuffield College, Oxford, under the Jemolo 
Fellowship scheme; the generous support of the  college is gratefully 
acknowledged. 
Additional support for the development of the theoretical research work has 
been  provided by the "Consiglio Nazionale delle Ricerche" of Italy, grant 
no.97.01331.CT10.
}
\references{
Azzalini, A. (1985).
A class of distributions which includes the normal ones.
\emph{Scand. J. Statist.}
\bold{12}, 171--178.


Azzalini, A. (1986).
Further results on a class of distributions which includes 
the normal ones.
\emph{Statistica}
\bold{46},199--208.


Azzalini, A. and Dalla Valle, A. (1996).
The multivariate skew-normal distribution.
\emph{Biometrika}
\bold{83}, 715--726.


Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579--602.
}
\section{LICENCE}{
This library and its documentation are usable under the terms of the 
"GNU General Public License", a copy of which is distributed with the 
package.
}
\keyword{multivariate}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
