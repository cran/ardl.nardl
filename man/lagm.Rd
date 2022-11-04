\name{lagm}
\alias{lagm}
\title{
Lag a matrix
}
\description{
Take the lag of a variable using the specified number of lag(s).
}
\usage{
lagm(x, k)
}
\arguments{
  \item{x}{
a matrix object
}
  \item{k}{
number of lags
}
}
\value{
\item{x }{dimentional matrix with k}
}
\examples{
x <- cbind(x1 = c(23, 45, 56, 78, 98, 21, 17),
            x2 = c(3, 4, 6, 7, 2, 9, 17))
lagm(x, k = 3)
}
