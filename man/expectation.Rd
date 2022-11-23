\name{expectation}
\alias{expectation}
\docType{data}
\title{
Inflation Expectation Dataset
}
\description{
Dataset
}
\usage{data("expectation")}
\format{
  A data frame with 47 observations on the following variables.
  \describe{
    \item{\code{date}}{date from 2009Q1 to 2020Q4}
    \item{\code{nq_inf_exp}}{Natural log of Next quarter inflation expectation (index)}
    \item{\code{n12m_inf_exp}}{Natural log of Next 12-month inflation expectation (index)}
    \item{\code{hawkish}}{Hawkishness}
    \item{\code{dovish}}{Dovishness}
    \item{\code{nethawkish}}{nethawkishness}
    \item{\code{headline_inf}}{actual headline inflation (percent)}
    \item{\code{core_inf}}{actual core inflation (percent)}
    \item{\code{food_inf}}{actual food inflation (percent)}
  }
}
\details{
The sentiment indices were extracted (text mined) from CBN monetary policy committee ommuniques with the aid of the monetary policy dictionary adopted in Mate et al. (2021). 
}
\source{
....
}
\references{
Mate A, Sebok M, Barczikay T (2021) The effect of central bank communication on sovereign bond yields: The case of Hungary. PLoS ONE 16(2): e0245515. https://doi.org/10.1371/journal.pone.0245515
}
\examples{
data(expectation)
}
\keyword{datasets}
