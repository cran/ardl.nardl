\name{fuel_price}
\alias{fuel_price}
\docType{data}
\title{
Time series data from 2000M01 through 2021M03
}
\description{
Time series data on fuel pump price, oil price and exchange rate
}
\usage{data("fuel_price")}
\format{
  A data frame with 255 observations on the following variables.
  \describe{
    \item{\code{fpp}}{log of Fuel Pump Price Per Litre - Average (PMS)}
    \item{\code{wti}}{log of West Texas Intemediate}
    \item{\code{bdc}}{log of Bureau De-change exchange rate (N/$)}
    \item{\code{date}}{time series from 2000M01 through 2021M03}
  }
}
\source{
fpp is obtained from Central Bank of Nigeria 2021Q1 Statistical Bulletin

wti obtained from Wold bank commodity price (pinksheet)

bdc sourced from http://statistics.cbn.gov.ng/cbn-onlinestats/

}
\examples{
data(fuel_price)
}
\keyword{datasets}
