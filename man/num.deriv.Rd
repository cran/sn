\name{num.deriv}
\alias{num.deriv}
\title{
Numerical derivatives
}
\description{
Computes numerical derivatives of a vector-valued function
}
\usage{
num.deriv(coefficients, FUN, ...)
}
\arguments{
\item{coefficients}{
a vector containing the coordinates of the point where the
derivatives must be evaluated
}
\item{FUN}{
a vector-valued function accepting a vector argument of the same length
}
\item{...}{
additional arguments to \code{FUN}
}}
\details{
The function is intended for numerical evaluation of the Hessian
of a function whose gradient is given by \code{FUN}.
Used by library \code{SN} to obtain the second derivative of
\code{msn.dev.grad}
}
\seealso{
\code{\link{msn.dev.grad}}
}
\keyword{math}
