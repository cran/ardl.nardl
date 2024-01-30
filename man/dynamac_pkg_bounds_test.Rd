\name{dynamac_pkg_bounds_test}
\alias{dynamac_pkg_bounds_test}
\title{
A function in 'dynamac' package to Perform Pesaran, Shin, and Smith (2001) cointegration test

}
\description{
This is an Abridged version of Jordan and Philips (2020). Please, see the dynamac package for the full documentation and due reference be given to the authors of dynamac  - Jordan and Philips (2020).
}
\usage{
dynamac_pkg_bounds_test(obs, fstat, tstat = NULL, case, k)
}
\arguments{
  \item{obs}{
  Number of observations
}
  \item{fstat}{
  F-statistics 
}
  \item{tstat}{
  't-statistic of the lagged dependent variable'
}
  \item{case}{
  case I to IV
}
  \item{k}{
  'number of regressors appearing in levels in the estimated model, not including the lagged dependent variable' 
}
}
\value{
\item{Fstat }{The bounds F statistic, lower and upper bounds}
\item{tstat }{The bounds t statistic, lower and upper bounds}
}
\references{
Jordan S, Philips A (2020). _dynamac: Dynamic Simulation and Testing for Single-Equation
      ARDL Models_. R package version 0.1.11 <https://CRAN.R-project.org/package=dynamac>
}
\seealso{
\code{\link[dynamac]{pssbounds}}
}
% \keyword{ ~kwd1 }
