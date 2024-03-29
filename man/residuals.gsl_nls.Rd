% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nls_methods.R
\name{residuals.gsl_nls}
\alias{residuals.gsl_nls}
\title{Extract model residuals}
\usage{
\method{residuals}{gsl_nls}(object, type = c("response", "pearson"), ...)
}
\arguments{
\item{object}{An object inheriting from class \code{"gsl_nls"}.}

\item{type}{character; if \code{"response"} the raw residuals are returned, if \code{"pearson"}
the Pearson are returned, i.e. the raw residuals divided by their standard error.}

\item{...}{At present no optional arguments are used.}
}
\value{
Numeric vector of model residuals similar to \code{\link[stats]{residuals}}.
}
\description{
Returns the model residuals from a fitted \code{"gsl_nls"} object.
\code{resid} can also be used as an alias.
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

residuals(obj)
}
\seealso{
\code{\link[stats]{residuals}}
}
