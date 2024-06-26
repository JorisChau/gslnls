\name{gsl_nls_large}
\alias{gsl_nls_large}
\alias{gsl_nls_large.formula}
\alias{gsl_nls_large.function}
\title{GSL Large-scale Nonlinear Least Squares fitting}
\usage{
gsl_nls_large(fn, ...)

\method{gsl_nls_large}{formula}(
  fn,
  data = parent.frame(),
  start,
  algorithm = c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D", "cgst"),
  control = gsl_nls_control(),
  jac,
  fvv,
  trace = FALSE,
  subset,
  weights,
  na.action,
  model = FALSE,
  ...
)

\method{gsl_nls_large}{function}(
  fn,
  y,
  start,
  algorithm = c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D", "cgst"),
  control = gsl_nls_control(),
  jac,
  fvv,
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
\item \code{"lm"} Levenberg-Marquardt algorithm (default).
\item \code{"lmaccel"} Levenberg-Marquardt algorithm with geodesic acceleration.
Can be faster than \code{"lm"} but less stable. Stability is controlled by the
\code{avmax} parameter in \code{control}, setting \code{avmax} to zero is analogous
to not using geodesic acceleration.
\item \code{"dogleg"} Powell's dogleg algorithm.
\item \code{"ddogleg"} Double dogleg algorithm, an improvement over \code{"dogleg"}
by including information about the Gauss-Newton step while the iteration is still
far from the minimum.
\item \code{"subspace2D"} 2D generalization of the dogleg algorithm. This method
searches a larger subspace for a solution, it can converge more quickly than \code{"dogleg"}
on some problems.
\item \code{"cgst"} Steihaug-Toint Conjugate Gradient algorithm, a generalization of the dogleg algorithm
that avoids solving for the Gauss-Newton step directly, instead using an iterative conjugate gradient algorithm.
The method performs well at points where the Jacobian is singular, and is also suitable for large-scale problems
where factoring the Jacobian matrix is prohibitively expensive.
}}

\item{control}{an optional list of control parameters to tune the least squares iterations and multistart algorithm.
See \code{\link{gsl_nls_control}} for the available control parameters and their default values.}

\item{jac}{a \link{function} returning the \code{n} by \code{p} dimensional Jacobian matrix of
the nonlinear model \code{fn}, where \code{n} is the number of observations and \code{p} the
number of parameters. The first argument must be the vector of parameters of length \code{p}.
Can also be \code{TRUE}, in which case \code{jac} is derived symbolically with \code{\link[stats]{deriv}},
this only works if \code{fn} is defined as a (non-selfstarting) formula. If \code{fn} is a \code{\link{selfStart}} model,
the Jacobian specified in the \code{"gradient"} attribute of the self-start model is used instead.}

\item{fvv}{a \link{function} returning an \code{n} dimensional vector containing
the second directional derivatives of the nonlinear model \code{fn}, with \code{n} the number of observations.
This argument is only used if geodesic acceleration is enabled (\code{algorithm = "lmaccel"}).
The first argument must be the vector of parameters of length \code{p} and the second argument must be the velocity vector
also of length \code{p}. Can also be \code{TRUE}, in which case \code{fvv} is derived
symbolically with \code{\link[stats]{deriv}}, this only works if \code{fn} is defined as a (non-selfstarting) formula.
If the model \link{function} in \code{fn} also returns a \code{"hessian"} attribute (similar to the \code{"gradient"} attribute
in a \code{selfStart} model), this Hessian matrix is used to evaluate the second directional derivatives instead.}

\item{trace}{logical value indicating if a trace of the iteration progress should be printed.
Default is \code{FALSE}. If \code{TRUE}, the residual (weighted) sum-of-squares,
the squared (Euclidean) norm of the current parameter estimates and the condition number of the Jacobian
are printed after each iteration.}

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
Determine the nonlinear least-squares estimates of the parameters of a large
nonlinear model system using the \code{gsl_multilarge_nlinear} module present in
the GNU Scientific Library (GSL).
}
\section{Methods (by class)}{
\itemize{
\item \code{gsl_nls_large(formula)}: If \code{fn} is a \code{formula}, the returned list object is of classes \code{gsl_nls} and \code{nls}.
Therefore, all generic functions applicable to objects of class \code{nls}, such as \code{anova}, \code{coef}, \code{confint},
\code{deviance}, \code{df.residual}, \code{fitted}, \code{formula}, \code{logLik}, \code{nobs}, \code{predict}, \code{print}, \code{profile},
\code{residuals}, \code{summary}, \code{vcov} and \code{weights} are also applicable to the returned list object.
In addition, a method \code{confintd} is available for inference of derived parameters.

\item \code{gsl_nls_large(function)}: If \code{fn} is a \code{function}, the first argument must be the vector of parameters and
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
# Large NLS example
# (https://www.gnu.org/software/gsl/doc/html/nls.html#large-nonlinear-least-squares-example)

## number of parameters
p <- 250

## model function
f <- function(theta) {
  c(sqrt(1e-5) * (theta - 1), sum(theta^2) - 0.25)
}

## jacobian function
jac <- function(theta) {
  rbind(diag(sqrt(1e-5), nrow = length(theta)), 2 * t(theta))
}

## dense Levenberg-Marquardt
\donttest{
gsl_nls_large(
  fn = f,                       ## model
  y = rep(0, p + 1),            ## (dummy) responses
  start = 1:p,                  ## start values
  algorithm = "lm",             ## algorithm
  jac = jac,                    ## jacobian
  control = list(maxiter = 250)
)
}

## dense Steihaug-Toint conjugate gradient
\donttest{
gsl_nls_large(
  fn = f,                       ## model
  y = rep(0, p + 1),            ## (dummy) responses
  start = 1:p,                  ## start values
  jac = jac,                    ## jacobian
  algorithm = "cgst"            ## algorithm
)
}

## sparse Jacobian function
jacsp <- function(theta) {
  rbind(Matrix::Diagonal(x = sqrt(1e-5), n = length(theta)), 2 * t(theta))
}

## sparse Levenberg-Marquardt
gsl_nls_large(
  fn = f,                       ## model
  y = rep(0, p + 1),            ## (dummy) responses
  start = 1:p,                  ## start values
  algorithm = "lm",             ## algorithm
  jac = jacsp,                  ## sparse jacobian
  control = list(maxiter = 250)
)

## sparse Steihaug-Toint conjugate gradient
gsl_nls_large(
  fn = f,                       ## model
  y = rep(0, p + 1),            ## (dummy) responses
  start = 1:p,                  ## start values
  jac = jacsp,                  ## sparse jacobian
  algorithm = "cgst"            ## algorithm
)

}
\references{
M. Galassi et al., \emph{GNU Scientific Library Reference Manual (3rd Ed.)}, ISBN 0954612078.
}
\seealso{
\code{\link{gsl_nls}}

\url{https://www.gnu.org/software/gsl/doc/html/nls.html}
}
