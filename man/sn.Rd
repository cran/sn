\name{SN}
\alias{SN}
\title{
The library `sn': summary information
}
\description{
This library provides functions related to the skew-normal probability
distribution, both for the univariate and the the multivariate case.
}
\section{FUNCTIONS}{
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
}
\section{REQUIREMENTS}{
R 1.0.1
}
\section{VERSION}{
You are using version 0.22 (2000/09/15).
The most recent version of the library can be obtained from
the WWW page: \url{http://www.stat.unipd.it/~azzalini/SN}
which also contains other related material
}
\section{MANUAL}{
There is no manual for this version (and possibly there will never 
be one).
All documentation is on-line, apart for the papers listed below, 
which provide background information. Some functions are not documented.
}
\section{AUTHOR}{
Adelchi Azzalini, Dept Statistical Sciences, University of Padua, Italy.
Please send comments, error reports, etc. to \email{azzalini@stat.unipd.it}
}
\section{ACKNOWLEDGEMENTS}{
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
\bold{12}, 171-178.


Azzalini, A. (1986).
Further results on a class of distributions which includes 
the normal ones.
\emph{Statistica}
\bold{46},199-208.


Azzalini, A. and Dalla Valle, A. (1996).
The multivariate skew-normal distribution.
\emph{Biometrika}
\bold{83}, 715--726.


Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew-normal distribution.
\emph{J.Roy.Statist.Soc. B}
\bold{61}, 579-602.
}
\section{LICENCE}{
This library and its documentation are usable under the terms of the 
"GNU General Public License", a copy of which is distributed with the 
package.
}
\keyword{multivariate}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
