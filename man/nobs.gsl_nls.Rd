% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nls_methods.R
\name{nobs.gsl_nls}
\alias{nobs.gsl_nls}
\title{Extract the number of observations}
\usage{
\method{nobs}{gsl_nls}(object, ...)
}
\arguments{
\item{object}{An object inheriting from class \code{"gsl_nls"}.}

\item{...}{At present no optional arguments are used.}
}
\value{
Integer number of observations similar to \code{\link[stats]{nobs}}
}
\description{
Returns the number of \emph{observations} from a \code{"gsl_nls"} object.
}
\examples{
## data
set.seed(1)
n <- 25
xy <- data.frame(
  x = (1:n) / n,
  y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
)
## model
obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))

nobs(obj)
}
\seealso{
\code{\link[stats]{nobs}}
}
