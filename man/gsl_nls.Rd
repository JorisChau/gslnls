\name{gsl_nls}
\alias{gsl_nls}
\alias{gsl_nls.formula}
\alias{gsl_nls.function}
\title{GSL Nonlinear Least Squares fitting}
\usage{
gsl_nls(fn, ...)

\method{gsl_nls}{formula}(
  fn,
  data = parent.frame(),
  start,
  algorithm = c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"),
  control = gsl_nls_control(),
  jac = NULL,
  fvv = NULL,
  trace = FALSE,
  subset,
  weights,
  na.action,
  model = FALSE,
  ...
)

\method{gsl_nls}{function}(
  fn,
  y,
  start,
  algorithm = c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"),
  control = gsl_nls_control(),
  jac = NULL,
  fvv = NULL,
  trace = FALSE,
  weights,
  ...
)
}
\arguments{
\item{fn}{a nonlinear model defined either as a two-sided \link{formula} including variables and parameters,
or as a \link{function} returning a numeric vector, with first argument the vector of parameters to be estimated.
See the individual method descriptions below.}

\item{data}{an optional data frame in which to evaluate the variables in \code{fn} if
defined as a \link{formula}. Can also be a list or an environment, but not a matrix.}

\item{y}{numeric response vector if \code{fn} is defined as a \link{function}, equal in
length to the vector returned by evaluation of the function \code{fn}.}

\item{start}{a named list or named numeric vector of starting estimates. \code{start} is
only allowed to be missing if \code{fn} is a \code{\link{selfStart}} model. If \code{fn}
is a \code{formula}, a naive guess for \code{start} is tried, but this should not be relied on.}

\item{algorithm}{character string specifying the algorithm to use. The following choices are supported:
\itemize{
\item \code{"lm"} Levenberg-Marquadt algorithm (default)
\item \code{"lmaccel"} Levenberg-Marquadt algorithm with geodesic acceleration.
Can be faster than \code{"lm"} but less stable. Stability is controlled by the
\code{avmax} parameter in \code{control}, setting \code{avmax} to zero is analogous
to not using geodesic acceleration.
\item \code{"dogleg"} Powell's dogleg algorithm
\item \code{"ddogleg"} Double dogleg algorithm, an improvement over \code{"dogleg"}
by including information about the Gauss-Newton step while the iteration is still
far from the minimum.
\item \code{"subspace2D"} 2D generalization of the dogleg algorithm. This method
searches a larger subspace for a solution, it can converge more quickly than \code{"dogleg"}
on some problems.
}}

\item{control}{an optional list of control parameters to tune the least squares iterations.
See \code{\link{gsl_nls_control}} for the available control parameters and their default values.}

\item{jac}{either \code{NULL} (default) or a \link{function} returning the \code{n} by \code{p} dimensional Jacobian matrix of
the nonlinear model \code{fn}, where \code{n} is the number of observations and \code{p} the
number of parameters. If a function, the first argument must be the vector of parameters of length \code{p}.
If \code{NULL}, the Jacobian is computed internally using a finite difference approximations.
Can also be \code{TRUE}, in which case \code{jac} is derived symbolically with \code{\link[stats]{deriv}},
this only works if \code{fn} is defined as a (non-selfstarting) formula.}

\item{fvv}{either \code{NULL} (default) or a \link{function} returning an \code{n} dimensional vector containing
the second directional derivatives of the nonlinear model \code{fn}, with \code{n} the number of observations.
This argument is only used if geodesic acceleration is enabled (\code{algorithm = "lmaccel"}).
If a function, the first argument must be the vector of parameters of length \code{p} and the second argument must be the velocity vector
also of length \code{p}. If \code{NULL}, the second directional derivative vector is computed internal
using a finite difference approximation. Can also be \code{TRUE}, in which case \code{fvv} is derived
symbolically with \code{\link[stats]{deriv}}, this only works if \code{fn} is defined as a (non-selfstarting) formula.}

\item{trace}{logical value indicating if a trace of the iteration progress should be printed.
Default is \code{FALSE}. If \code{TRUE}, the residual (weighted) sum-of-squares and the current parameter estimates are printed after each iteration.}

\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.
This argument is only used if \code{fn} is defined as a \link{formula}.}

\item{weights}{an optional numeric vector of (fixed) weights. When present, the objective function is
weighted least squares.}

\item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. The
default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.
The 'factory-fresh' default is \code{\link{na.omit}}. Value \code{\link{na.exclude}} can be useful.
This argument is only used if \code{fn} is defined as a \link{formula}.}

\item{model}{a logical value. If \code{TRUE}, the model frame is returned as part of the object. Defaults to \code{FALSE}.
This argument is only used if \code{fn} is defined as a \link{formula}.}

\item{...}{additional arguments passed to the calls of \code{fn}, \code{jac} and \code{fvv} if
defined as functions.}

}
\value{
If \code{fn} is a \code{formula} returns a list object of class \code{nls}.
If \code{fn} is a \code{function} returns a list object of class \code{gsl_nls}.
See the individual method descriptions for the structures of the returned lists and the generic functions
applicable to objects of both classes.
}
\description{
Determine the nonlinear (weighted) least-squares estimates of the parameters of a
nonlinear model using the GNU Scientific Library (GSL).
}
\section{Methods (by class)}{
\itemize{
\item \code{formula}: If \code{fn} is a \code{formula}, the returned list object is of classes \code{gsl_nls} and \code{nls}.
Therefore, all generic functions applicable to objects of class \code{nls}, such as \code{anova}, \code{coef}, \code{confint},
\code{deviance}, \code{df.residual}, \code{fitted}, \code{formula}, \code{logLik}, \code{nobs}, \code{predict}, \code{print}, \code{profile},
\code{residuals}, \code{summary}, \code{vcov} and \code{weights} are also applicable to the returned list object.
In addition, a method \code{confintd} is available for inference of derived parameters.

\item \code{function}: If \code{fn} is a \code{function}, the first argument must be the vector of parameters and
the function should return a numeric vector containing the nonlinear model evaluations at
the provided parameter and predictor or covariate vectors. In addition, the argument \code{y}
needs to contain the numeric vector of observed responses, equal in length to the numeric
vector returned by \code{fn}. The returned list object is (only) of class \code{gsl_nls}.
Although the returned object is not of class \code{nls}, the following generic functions remain
applicable for an object of class \code{gsl_nls}: \code{anova}, \code{coef}, \code{confint}, \code{deviance},
\code{df.residual}, \code{fitted}, \code{formula}, \code{logLik}, \code{nobs}, \code{predict}, \code{print},
\code{residuals}, \code{summary}, \code{vcov} and \code{weights}. In addition, a method \code{confintd}
is available for inference of derived parameters.
}}

\examples{
# Example 1: exponential model
# (https://www.gnu.org/software/gsl/doc/html/nls.html#exponential-fitting-example)

## data
set.seed(1)
n <- 50
x <- (seq_len(n) - 1) * 3 / (n - 1)
f <- function(A, lam, b, x) A * exp(-lam * x) + b
y <- f(A = 5, lam = 1.5, b = 1, x) + rnorm(n, sd = 0.25)

## model fit
ex1_fit <- gsl_nls(
  fn = y ~ A * exp(-lam * x) + b,                        ## model formula
  data = data.frame(x = x, y = y),                       ## model fit data
  start = c(A = 0, lam = 0, b = 0)                       ## starting values
)
summary(ex1_fit)                                         ## model summary
predict(ex1_fit, interval = "prediction")                ## prediction intervals

## analytic Jacobian 1
gsl_nls(
  fn = y ~ A * exp(-lam * x) + b,                        ## model formula
  data = data.frame(x = x, y = y),                       ## model fit data
  start = c(A = 0, lam = 0, b = 0),                      ## starting values
  jac = function(par) with(as.list(par),                 ## jacobian
    cbind(A = exp(-lam * x), lam = -A * x * exp(-lam * x), b = 1)
  )
)

## analytic Jacobian 2
gsl_nls(
  fn = y ~ A * exp(-lam * x) + b,                        ## model formula
  data = data.frame(x = x, y = y),                       ## model fit data
  start = c(A = 0, lam = 0, b = 0),                      ## starting values
  jac = TRUE                                             ## automatic derivation
)

## self-starting model
gsl_nls(
  fn =  y ~ SSasymp(x, Asym, R0, lrc),                   ## model formula
  data = data.frame(x = x, y = y)                        ## model fit data
)

# Example 2: Gaussian function
# (https://www.gnu.org/software/gsl/doc/html/nls.html#geodesic-acceleration-example-2)

## data
set.seed(1)
n <- 300
x <- seq_len(n) / n
f <- function(a, b, c, x) a * exp(-(x - b)^2 / (2 * c^2))
y <- f(a = 5, b = 0.4, c = 0.15, x) * rnorm(n, mean = 1, sd = 0.1)

## Levenberg-Marquadt (default)
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  trace = TRUE                                          ## verbose output
)

## Levenberg-Marquadt w/ geodesic acceleration 1
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  algorithm = "lmaccel",                                ## algorithm
  trace = TRUE                                          ## verbose output
)

## Levenberg-Marquadt w/ geodesic acceleration 2
## second directional derivative
fvv <- function(par, v, x) {
  with(as.list(par), {
    zi <- (x - b) / c
    ei <- exp(-zi^2 / 2)
    2 * v[["a"]] * v[["b"]] * zi / c * ei + 2 * v[["a"]] * v[["c"]] * zi^2 / c * ei -
      v[["b"]]^2 * a / c^2 * (1 - zi^2) * ei -
      2 * v[["b"]] * v[["c"]] * a / c^2 * zi * (2 - zi^2) * ei -
      v[["c"]]^2 * a / c^2 * zi^2 * (3 - zi^2) * ei
  })
}

## analytic fvv 1
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  algorithm = "lmaccel",                                ## algorithm
  trace = TRUE,                                         ## verbose output
  fvv = fvv,                                            ## analytic fvv
  x = x                                                 ## argument passed to fvv
)

## analytic fvv 2
gsl_nls(
  fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
  data = data.frame(x = x, y = y),                      ## model fit data
  start = c(a = 1, b = 0, c = 1),                       ## starting values
  algorithm = "lmaccel",                                ## algorithm
  trace = TRUE,                                         ## verbose output
  fvv = TRUE                                            ## automatic derivation
)

# Example 3: Branin function
# (https://www.gnu.org/software/gsl/doc/html/nls.html#comparing-trs-methods-example)

## Branin model function
branin <- function(x) {
  a <- c(-5.1 / (4 * pi^2), 5 / pi, -6, 10, 1 / (8 * pi))
  f1 <- x[2] + a[1] * x[1]^2 + a[2] * x[1] + a[3]
  f2 <- sqrt(a[4] * (1 + (1 - a[5]) * cos(x[1])))
  c(f1, f2)
}

## Dogleg minimization w/ model as function
gsl_nls(
  fn = branin,                   ## model function
  y = c(0, 0),                   ## response vector
  start = c(x1 = 6, x2 = 14.5),  ## starting values
  algorithm = "dogleg"           ## algorithm
)
}
\references{
M. Galassi et al., \emph{GNU Scientific Library Reference Manual (3rd Ed.)}, ISBN 0954612078.
}
\seealso{
\code{\link[stats]{nls}}

\url{https://www.gnu.org/software/gsl/doc/html/nls.html}
}
