\name{syg_data}
\alias{syg_data}
\docType{data}
\title{
Time series data on output and unemployment in Canada, Japan and United States from 1981M01 to 2022M07
}
\description{
Data on unemployment and industrial production for Canada, Japan  and United States
}
\usage{data("syg_data")}
\format{
  A data frame with 0 observations on the following variables.
  \describe{
    \item{\code{ca_u}}{Canada (CA) unemployment rate }
    \item{\code{jp_u}}{Japan (JP) unemployment rate }
    \item{\code{us_u}}{United States (US) unemployment rate }
    \item{\code{us_ip}}{Natural log of US industrial production (s.a)}
    \item{\code{jp_ip}}{Natural log of JP industrial production (s.a)}
    \item{\code{ca_ip}}{Natural log of CA industrial production (s.a)}
    \item{\code{time}}{1981M01 - 2022M07}
  }
}
\details{
_U is Unemployment rate (monthly), Total, All persons
_IP is index of industrial production, s.a.

}
\source{
https://stats.oecd.org/index.aspx?lang=en 
}
\examples{
data(syg_data)
}
\keyword{datasets}
