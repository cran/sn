\name{ais}
\alias{ais}
\title{Australian Institute of Sport data}
\usage{data(ais)}
\description{
  Data on 102 male and 100 female athletes collected at the Australian
  Institute of Sport, courtesy of Richard Telford and Ross Cunningham.
}

\format{
  A data frame with 202 observations on 13 variables.
  \tabular{rll}{
    [, 1] \tab sex  \tab sex \cr
    [, 2] \tab sport\tab sport \cr
    [, 3] \tab rcc  \tab red cell count \cr
    [, 4] \tab wcc  \tab white cell count \cr
    [, 5] \tab Hc   \tab Hematocrit \cr
    [, 6] \tab Hg   \tab Hemoglobin \cr
    [, 7] \tab Fe   \tab plasma ferritin concentration \cr
    [, 8] \tab bmi  \tab body mass index, weight/(height)² \cr
    [, 9] \tab Bfat \tab body fat percentage \cr
    [,10] \tab ssf  \tab sum of skin folds \cr
    [,11] \tab lbm  \tab lean body mass \cr
    [,12] \tab Ht   \tab height (cm) \cr
    [,13] \tab Wt   \tab weight (Kg) \cr
  }
}
\source{
  Cook and Weisberg (1994),
  \emph{An Introduction to Regression Graphics}.                   
  John Wiley & Sons, New York.
}
\examples{
data(ais)
attach(ais)
pairs(ais[,c(3:4,10:13)], main = "AIS data")
plot(Wt~sport)
}
\keyword{datasets}
