\name{pnorm2}
\alias{pnorm2}
\title{
Bivariate normal integral
}
\description{
Computes the cumulative distribution function of a bivariate
normal distribution with standardized marginals.
}
\usage{
pnorm2(x, y, rho)
}
\arguments{
\item{x}{
a numerical value representing the first coordinate
}
\item{y}{
a numerical value representing the second coordinate
}
\item{rho}{
a numerical value representing the correlation;
it must be between -1 and 1.
}}
\value{
a numerical value with the required probability
}
\details{
The input parameters must be all scalars.


This function is based on function \code{T.Owen}.
}
\seealso{
\code{\link{T.Owen}}
}
\examples{
p <- pnorm2(1.2, 0.5, 0.67)
}
\keyword{distribution}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
