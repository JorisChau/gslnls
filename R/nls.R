#' GSL Nonlinear Least Squares fitting
#'
#' @description
#' Determine the nonlinear least-squares estimates of the parameters of a
#' nonlinear model using the \code{gsl_multifit_nlinear} module present in
#' the GNU Scientific Library (GSL).
#'
#' @section Multi-start algorithm:
#' If \code{start} is a list or matrix of parameter ranges, or contains any missing values, a modified version of the multi-start algorithm described in
#' Hickernell and Yuan (1997) is applied. Note that the \code{start} parameter ranges are only used to bound the domain for the
#' starting values, i.e. the resulting parameter estimates are not constrained to lie within these bounds, use \code{lower} and/or \code{upper} for
#' this purpose instead. Quasi-random starting values are sampled in the unit hypercube from a Sobol sequence if \code{p < 41} or from a Halton sequence (up to \code{p = 1229}) otherwise.
#' The initial starting values are scaled to the specified parameter ranges using an inverse (scaled) logistic function favoring starting values near the center of the
#' (scaled) domain. The trust region method as specified by \code{algorithm} used for the inexpensive and expensive local search (see Algorithm 2.1 of Hickernell
#' and Yuan (1997)) are the same, only differing in the number of search iterations \code{mstart_p} versus \code{mstart_maxiter}, where
#' \code{mstart_p} is typically much smaller than \code{mstart_maxiter}. When a new stationary point is detected, the scaling step from the unit hypercube to
#' the starting value domain is updated using the diagonal of the estimated trust method's scaling matrix \code{D}, which improves optimization performance
#' especially when the parameters live on very different scales. The multi-start algorithm terminates when NSP (number of stationary points)
#' is larger than or equal to \code{mstart_minsp} and NWSP (number of worse stationary points) is larger than \code{mstart_r + sqrt(mstart_r) * NSP},
#' or when the maximum number of major iterations \code{mstart_maxstart} is reached. After termination of the multi-start algorithm, a full
#' single-start optimization is executed starting from the best multi-start solution.
#'
#' @section Missing starting values:
#' If \code{start} contains missing (or infinite) values, the multi-start algorithm is executed without fixed parameter ranges for the missing parameters.
#' The ranges for the missing parameters are initialized to the unit interval and dynamically increased or decreased in each major iteration
#' of the multi-start algorithm. The decision to increase or decrease a parameter range is driven by the minimum and maximum parameter values
#' obtained by the first \code{mstart_q} inexpensive local searches ordered by their squared loss, which typically provide a decent indication of the
#' order of magnitude of the parameter range in which to search for the optimal solution. Note that this procedure is not expected to always
#' return a global minimum of the nonlinear least-squares objective. Especially when the objective function contains many local optima,
#' the algorithm may be unable to select parameter ranges that include the global minimizing solution. In this case, it may help to increase
#' the values of \code{mstart_n}, \code{mstart_r} or \code{mstart_minsp} to avoid early termination of the algorithm at the cost of
#' increased computational effort.
#'
#' @param fn a nonlinear model defined either as a two-sided \link{formula} including variables and parameters,
#' or as a \link{function} returning a numeric vector, with first argument the vector of parameters to be estimated.
#' See the individual method descriptions below.
#' @param data an optional data frame in which to evaluate the variables in \code{fn} if
#' defined as a \link{formula}. Can also be a list or an environment, but not a matrix.
#' @param y numeric response vector if \code{fn} is defined as a \link{function}, equal in
#' length to the vector returned by evaluation of the function \code{fn}.
#' @param start a vector, list or matrix of initial parameter values or parameter ranges. \code{start} is only allowed to be missing
#' if \code{fn} is a \code{\link{selfStart}} model. The following choices are supported:
#' \itemize{
#'  \item a named list or named vector of numeric starting values. If \code{start} has no missing values, a standard single-start optimization
#'  is performed. If \code{start} contains missing values for one or more parameters, a multi-start algorithm (see \sQuote{Details}) with
#'  dynamic starting ranges for the undefined parameters and fixed starting values for the remaining parameters is executed.
#'  If \code{start} is a named list or vector containing \emph{only} missing values, the multi-start algorithm considers dynamically changing starting
#'  ranges for all parameters. Note that there is no guarantee that the optimizing solution is a global minimum of the least-squares objective.
#'  \item a named list with starting parameter ranges in the form of length-2 numeric vectors. Can also be a (\code{2} by \code{p}) named matrix with as columns
#' the numeric starting ranges for the parameters. If \code{start} contains no missing values, a multi-start algorithm with fixed
#' starting ranges for the parameters is executed. Otherwise, if \code{start} contains infinities or missing values (e.g. \code{c(0, Inf)} or \code{c(NA, NA)}),
#' the multi-start algorithm considers dynamically changing starting ranges for the parameters with infinite and/or missing ranges.
#' }
#' @param algorithm character string specifying the algorithm to use. The following choices are supported:
#' \itemize{
#' \item \code{"lm"} Levenberg-Marquardt algorithm (default).
#' \item \code{"lmaccel"} Levenberg-Marquardt algorithm with geodesic acceleration.
#' Stability is controlled by the \code{avmax} parameter in \code{control}, setting \code{avmax}
#' to zero is analogous to not using geodesic acceleration.
#' \item \code{"dogleg"} Powell's dogleg algorithm.
#' \item \code{"ddogleg"} Double dogleg algorithm, an improvement over \code{"dogleg"}
#' by including information about the Gauss-Newton step while the iteration is still
#' far from the minimum.
#' \item \code{"subspace2D"} 2D generalization of the dogleg algorithm. This method
#' searches a larger subspace for a solution, it can converge more quickly than \code{"dogleg"}
#' on some problems.
#' }
#' @param loss character string specifying the loss function to optimize. The following choices are supported:
#' \itemize{
#' \item \code{"default"} default squared loss function.
#' \item \code{"huber"} Huber loss function.
#' \item \code{"barron"} Barron's smooth family of loss functions.
#' \item \code{"biweight"} Tukey's biweight/bisquare loss function.
#' \item \code{"welsh"} Welsh/Leclerc loss function.
#' \item \code{"optimal"} Optimal loss function (Maronna et al. (2006), Section 5.9.1).
#' \item \code{"hampel"} Hampel loss function.
#' \item \code{"ggw"} Generalized Gauss-Weight loss function.
#' \item \code{"lqq"} Linear Quadratic Quadratic loss function.
#' }
#' If a character string, the default tuning parameters as specified by \code{\link{gsl_nls_loss}} are used.
#' Instead, a list as returned by \code{\link{gsl_nls_loss}} with non-default tuning parameters is also accepted.
#' For all choices other than \code{rho = "default"}, iterative reweighted least squares (IRLS) is used to
#' solve the MM-estimation problem.
#' @param control an optional list of control parameters to tune the least squares iterations and multistart algorithm.
#' See \code{\link{gsl_nls_control}} for the available control parameters and their default values.
#' @param lower	a named list or named numeric vector of parameter lower bounds, or an unnamed numeric
#' scalar to be replicated for all parameters. If missing (default), the parameters are unconstrained from below.
#' @param upper a named list or named numeric vector of parameter upper bounds, or an unnamed numeric
#' scalar to be replicated for all parameters. If missing (default), the parameters are unconstrained from above.
#' @param jac either \code{NULL} (default) or a \link{function} returning the \code{n} by \code{p} dimensional Jacobian matrix of
#' the nonlinear model \code{fn}, where \code{n} is the number of observations and \code{p} the
#' number of parameters. If a function, the first argument must be the vector of parameters of length \code{p}.
#' If \code{NULL}, the Jacobian is computed internally using a finite difference approximations.
#' Can also be \code{TRUE}, in which case \code{jac} is derived symbolically with \code{\link[stats]{deriv}},
#' this only works if \code{fn} is defined as a (non-selfstarting) formula. If \code{fn} is a \code{\link{selfStart}} model,
#' the Jacobian specified in the \code{"gradient"} attribute of the self-start model is used instead.
#' @param fvv either \code{NULL} (default) or a \link{function} returning an \code{n} dimensional vector containing
#' the second directional derivatives of the nonlinear model \code{fn}, with \code{n} the number of observations.
#' This argument is only used if geodesic acceleration is enabled (\code{algorithm = "lmaccel"}).
#' If a function, the first argument must be the vector of parameters of length \code{p} and the second argument must be the velocity vector
#' also of length \code{p}. If \code{NULL}, the second directional derivative vector is computed internal
#' using a finite difference approximation. Can also be \code{TRUE}, in which case \code{fvv} is derived
#' symbolically with \code{\link[stats]{deriv}}, this only works if \code{fn} is defined as a (non-selfstarting) formula.
#' If the model \link{function} in \code{fn} also returns a \code{"hessian"} attribute (similar to the \code{"gradient"} attribute
#' in a \code{selfStart} model), this Hessian matrix is used to evaluate the second directional derivatives instead.
#' @param trace logical value indicating if a trace of the iteration progress should be printed.
#' Default is \code{FALSE}. If \code{TRUE}, the residual (weighted) sum-of-squares and the current parameter estimates
#' are printed after each iteration.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' This argument is only used if \code{fn} is defined as a \link{formula}.
#' @param weights an optional numeric vector of (fixed) weights of length \code{n} or an \code{n}-by-\code{n}
#' symmetric positive definite weight matrix. If \code{weights} is a vector or a diagonal matrix, the objective function is weighted least squares.
#' If \code{weights} is a general matrix, the objective function is generalized least squares.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s. The
#' default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.
#' The 'factory-fresh' default is \code{\link{na.omit}}. Value \code{\link{na.exclude}} can be useful.
#' This argument is only used if \code{fn} is defined as a \link{formula}.
#' @param model a logical value. If \code{TRUE}, the model frame is returned as part of the object. Defaults to \code{FALSE}.
#' This argument is only used if \code{fn} is defined as a \link{formula}.
#' @param ... additional arguments passed to the calls of \code{fn}, \code{jac} and \code{fvv} if
#' defined as functions.
#' @return
#' If \code{fn} is a \code{formula} returns a list object of class \code{nls}.
#' If \code{fn} is a \code{function} returns a list object of class \code{gsl_nls}.
#' See the individual method descriptions for the structures of the returned lists and the generic functions
#' applicable to objects of both classes.
#' @useDynLib gslnls, .registration = TRUE
#' @importFrom stats nls numericDeriv deriv as.formula coef deviance df.residual fitted vcov formula getInitial model.weights
#' @importFrom stats pf pt qt setNames sigma nobs hatvalues printCoefmat symnum cooks.distance naprint
#' @seealso \code{\link[stats]{nls}}
#' @seealso \url{https://www.gnu.org/software/gsl/doc/html/nls.html}
#' @seealso \url{https://CRAN.R-project.org/package=robustbase/vignettes/psi_functions.pdf}
#' @references M. Galassi et al., \emph{GNU Scientific Library Reference Manual (3rd Ed.)}, ISBN 0954612078.
#' @references Hickernell, F.J. and Yuan, Y. (1997) \emph{“A simple multistart algorithm for global optimization”}, OR Transactions, Vol. 1 (2).
#' @examples
#' # Example 1: exponential model
#' # (https://www.gnu.org/software/gsl/doc/html/nls.html#exponential-fitting-example)
#'
#' ## data
#' set.seed(1)
#' n <- 25
#' x <- (seq_len(n) - 1) * 3 / (n - 1)
#' f <- function(A, lam, b, x) A * exp(-lam * x) + b
#' y <- f(A = 5, lam = 1.5, b = 1, x) + rnorm(n, sd = 0.25)
#'
#' ## model fit
#' ex1_fit <- gsl_nls(
#'   fn = y ~ A * exp(-lam * x) + b,                        ## model formula
#'   data = data.frame(x = x, y = y),                       ## model fit data
#'   start = c(A = 0, lam = 0, b = 0)                       ## starting values
#' )
#' summary(ex1_fit)                                         ## model summary
#' predict(ex1_fit, interval = "prediction")                ## prediction intervals
#'
#' ## multi-start
#' gsl_nls(
#'   fn = y ~ A * exp(-lam * x) + b,                             ## model formula
#'   data = data.frame(x = x, y = y),                            ## model fit data
#'   start = list(A = c(0, 100), lam = c(0, 10), b = c(-10, 10)) ## fixed starting ranges
#' )
#' ## missing starting values
#' gsl_nls(
#'   fn = y ~ A * exp(-lam * x) + b,                        ## model formula
#'   data = data.frame(x = x, y = y),                       ## model fit data
#'   start = c(A = NA, lam = NA, b = NA)                    ## dynamic starting ranges
#' )
#'
#' ## robust regression
#' gsl_nls(
#'   fn = y ~ A * exp(-lam * x) + b,                      ## model formula
#'   data = data.frame(x = x, y = y),                     ## model fit data
#'   start = c(A = 0, lam = 0, b = 0),                    ## starting values
#'   loss = "barron"                                      ## L1-L2 loss
#' )
#'
#' ## analytic Jacobian 1
#' gsl_nls(
#'   fn = y ~ A * exp(-lam * x) + b,                        ## model formula
#'   data = data.frame(x = x, y = y),                       ## model fit data
#'   start = c(A = 0, lam = 0, b = 0),                      ## starting values
#'   jac = function(par) with(as.list(par),                 ## jacobian
#'     cbind(A = exp(-lam * x), lam = -A * x * exp(-lam * x), b = 1)
#'   )
#' )
#'
#' ## analytic Jacobian 2
#' gsl_nls(
#'   fn = y ~ A * exp(-lam * x) + b,                        ## model formula
#'   data = data.frame(x = x, y = y),                       ## model fit data
#'   start = c(A = 0, lam = 0, b = 0),                      ## starting values
#'   jac = TRUE                                             ## automatic derivation
#' )
#'
#' ## self-starting model
#' gsl_nls(
#'   fn =  y ~ SSasymp(x, Asym, R0, lrc),                   ## model formula
#'   data = data.frame(x = x, y = y)                        ## model fit data
#' )
#'
#' # Example 2: Gaussian function
#' # (https://www.gnu.org/software/gsl/doc/html/nls.html#geodesic-acceleration-example-2)
#'
#' ## data
#' set.seed(1)
#' n <- 100
#' x <- seq_len(n) / n
#' f <- function(a, b, c, x) a * exp(-(x - b)^2 / (2 * c^2))
#' y <- f(a = 5, b = 0.4, c = 0.15, x) * rnorm(n, mean = 1, sd = 0.1)
#'
#' ## Levenberg-Marquardt (default)
#' gsl_nls(
#'   fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
#'   data = data.frame(x = x, y = y),                      ## model fit data
#'   start = c(a = 1, b = 0, c = 1),                       ## starting values
#'   trace = TRUE                                          ## verbose output
#' )
#'
#' ## Levenberg-Marquardt w/ geodesic acceleration 1
#' gsl_nls(
#'   fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
#'   data = data.frame(x = x, y = y),                      ## model fit data
#'   start = c(a = 1, b = 0, c = 1),                       ## starting values
#'   algorithm = "lmaccel",                                ## algorithm
#'   trace = TRUE                                          ## verbose output
#' )
#'
#' ## Levenberg-Marquardt w/ geodesic acceleration 2
#' ## second directional derivative
#' fvv <- function(par, v, x) {
#'   with(as.list(par), {
#'     zi <- (x - b) / c
#'     ei <- exp(-zi^2 / 2)
#'     2 * v[["a"]] * v[["b"]] * zi / c * ei + 2 * v[["a"]] * v[["c"]] * zi^2 / c * ei -
#'       v[["b"]]^2 * a / c^2 * (1 - zi^2) * ei -
#'       2 * v[["b"]] * v[["c"]] * a / c^2 * zi * (2 - zi^2) * ei -
#'       v[["c"]]^2 * a / c^2 * zi^2 * (3 - zi^2) * ei
#'   })
#' }
#'
#' ## analytic fvv 1
#' gsl_nls(
#'   fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
#'   data = data.frame(x = x, y = y),                      ## model fit data
#'   start = c(a = 1, b = 0, c = 1),                       ## starting values
#'   algorithm = "lmaccel",                                ## algorithm
#'   trace = TRUE,                                         ## verbose output
#'   fvv = fvv,                                            ## analytic fvv
#'   x = x                                                 ## argument passed to fvv
#' )
#'
#' ## analytic fvv 2
#' gsl_nls(
#'   fn = y ~ a * exp(-(x - b)^2 / (2 * c^2)),             ## model formula
#'   data = data.frame(x = x, y = y),                      ## model fit data
#'   start = c(a = 1, b = 0, c = 1),                       ## starting values
#'   algorithm = "lmaccel",                                ## algorithm
#'   trace = TRUE,                                         ## verbose output
#'   fvv = TRUE                                            ## automatic derivation
#')
#'
#' # Example 3: Branin function
#' # (https://www.gnu.org/software/gsl/doc/html/nls.html#comparing-trs-methods-example)
#'
#' ## Branin model function
#' branin <- function(x) {
#'   a <- c(-5.1 / (4 * pi^2), 5 / pi, -6, 10, 1 / (8 * pi))
#'   f1 <- x[2] + a[1] * x[1]^2 + a[2] * x[1] + a[3]
#'   f2 <- sqrt(a[4] * (1 + (1 - a[5]) * cos(x[1])))
#'   c(f1, f2)
#' }
#'
#' ## Dogleg minimization w/ model as function
#' gsl_nls(
#'   fn = branin,                   ## model function
#'   y = c(0, 0),                   ## response vector
#'   start = c(x1 = 6, x2 = 14.5),  ## starting values
#'   algorithm = "dogleg"           ## algorithm
#' )
#'
#' # Available example problems
#' nls_test_list()
#'
#' ## BOD regression
#' ## (https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml)
#' (boxbod <- nls_test_problem(name = "BoxBOD"))
#' with(boxbod,
#'      gsl_nls(
#'        fn = fn,
#'        data = data,
#'        start = list(b1 = NA, b2 = NA)
#'      )
#' )
#'
#' ## Rosenbrock function
#' (rosenbrock <- nls_test_problem(name = "Rosenbrock"))
#' with(rosenbrock,
#'      gsl_nls(
#'        fn = fn,
#'        y = y,
#'        start = c(x1 = NA, x2 = NA),
#'        jac = jac
#'      )
#' )
#'
#' @export
gsl_nls <- function (fn, ...) {
  UseMethod("gsl_nls")
}

#' @describeIn gsl_nls
#' If \code{fn} is a \code{formula}, the returned list object is of classes \code{gsl_nls} and \code{nls}.
#' Therefore, all generic functions applicable to objects of class \code{nls}, such as \code{anova}, \code{coef}, \code{confint},
#' \code{deviance}, \code{df.residual}, \code{fitted}, \code{formula}, \code{logLik}, \code{nobs}, \code{predict}, \code{print}, \code{profile},
#' \code{residuals}, \code{summary}, \code{vcov}, \code{hatvalues}, \code{cooks.distance} and \code{weights} are also applicable to the returned list object.
#' In addition, a method \code{confintd} is available for inference of derived parameters.
#' @export
gsl_nls.formula <- function(fn, data = parent.frame(), start,
                            algorithm = c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"),
                            loss = c("default", "huber", "barron", "bisquare", "welsh", "optimal", "hampel", "ggw", "lqq"),
                            control = gsl_nls_control(), lower, upper, jac = NULL, fvv = NULL,
                            trace = FALSE, subset, weights, na.action, model = FALSE, ...) {

  ## adapted from src/library/stats/nls.R
  formula <- as.formula(fn)
  algorithm <- match.arg(algorithm, c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"))

  if(!is.list(data) && !is.environment(data))
    stop("'data' must be a list or an environment")

  mf <- match.call()		# for creating the model frame
  varNames <- all.vars(formula) # parameter and variable names from formula
  ## adjust a one-sided model formula by using 0 as the response
  if (length(formula) == 2L) {
    formula[[3L]] <- formula[[2L]]
    formula[[2L]] <- 0
  }
  ## for prediction we will need to know those which are in RHS
  form2 <- formula; form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)
  mWeights <- missing(weights)

  ## get names of the parameters from the starting values or selfStart model
  pnames <-
    if (missing(start)) {
      if(!is.null(attr(data, "parameters"))) {
        names(attr(data, "parameters"))
      } else { ## try selfStart - like object
        cll <- formula[[length(formula)]]
        if(is.symbol(cll)) { ## replace  y ~ S   by   y ~ S + 0 :
          ## formula[[length(formula)]] <-
          cll <- substitute(S + 0, list(S = cll))
        }
        fn <- as.character(cll[[1L]])
        if(is.null(func <- tryCatch(get(fn), error=function(e)NULL)))
          func <- get(fn, envir=parent.frame()) ## trying "above"
        if(!is.null(pn <- attr(func, "pnames")))
          as.character(as.list(match.call(func, call = cll))[-1L][pn])
      }
    } else if(is.matrix(start)) {
      colnames(start)
    } else {
      names(start)
    }

  if(!is.null(environment(formula))) {
    env <- environment(formula)
  } else {
    env <- parent.frame()
  }

  ## Heuristics for determining which names in formula represent actual
  ## variables :

  ## If it is a parameter it is not a variable
  if(length(pnames))
    varNames <- varNames[is.na(match(varNames, pnames))]

  ## This aux.function needs to be as complicated because
  ## exists(var, data) does not work (with lists or dataframes):
  lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)),
                                   error = function(e) -1L)
  if(length(varNames)) {
    n <- vapply(varNames, lenVar, 0)
    if(any(not.there <- n == -1L)) {
      nnn <- names(n[not.there])
      if(missing(start)) {
        stop("No starting values specified for parameters. Specify 'start' or use a selfStart model")
      } else {
        # has 'start' but forgot some
        stop(gettextf("parameters without starting value in 'data': %s",
                      paste(nnn, collapse=", ")), domain = NA)
      }
    }
  } else { ## length(varNames) == 0
    stop("No parameters to fit and/or no data variables present")
  }

  ## Transform multi-start starting values to matrix
  if(!missing(start)) {
    if(is.list(start) && any(lengths(start) > 1)) {
      if(!all(is.element(lengths(start), c(1L, 2L)))) {
        stop("List elements of 'start' must be of length 1 or 2 to specify (multi-start) parameter values or ranges")
      }
      if(any(len1 <- lengths(start) < 2)) {
        start[len1] <- lapply(start[len1], rep, times = 2L)
      }
      start <- as.matrix(data.frame(start))
    } else if(is.matrix(start) && nrow(start) != 2L) {
      stop("Matrix 'start' must have exactly 2 rows to specify (multi-start) parameter ranges")
    }
    ## Handle infinite/missing start values
    if(is.vector(start)) {
      if(any((napars <- sapply(start, Negate(is.finite))))) {
        start <- vapply(names(start), function(nm) if(napars[[nm]]) c(-0.1, 0.75) else rep(start[[nm]], 2L), numeric(2))
        .has_start <- sapply(!napars, rep, times = 2L)
      } else {
        .has_start <- !napars  ## not used (single-start)
      }
    } else if(is.matrix(start)) {
      if(any(napars <- apply(start, 2L, Negate(is.finite)))) {
        start[1L, napars[1L, ]]  <- -0.1
        start[2L, napars[2L, ]] <- 0.75
      }
      .has_start <- !napars
    } else {
      stop("Unrecognized format for starting values and/or ranges defined in 'start'")
    }
    ## Check start ranges
    if(is.matrix(start)) {
      if(length(pnames) > 1229L)
        stop("GSL quasi-random Halton sequences in multi-start are only available up to 1229 parameters")
      if(any(start[1L, ] > start[2L, ]))
        stop("Multi-start parameter lower bounds cannot be larger than upper bounds")
      else if(!any(start[1L, ] < start[2L, ]))
        start <- start[1L, ]  ## single-start
    }
  }

  ## If its length is a multiple of the response or LHS of the formula,
  ## then it is probably a variable.
  ## This may fail (e.g. when LHS contains parameters):
  respLength <- length(eval(formula[[2L]], data, env))
  if(length(n) > 0L) {
    varIndex <- n %% respLength == 0
    if(is.list(data) && diff(range(n[names(n) %in% names(data)])) > 0) {
      ## 'data' is a list that can not be coerced to a data.frame
      ## (not using varNames, varIndex at all - inconsistency FIXME?)
      mf <- data
      if(!missing(subset))
        warning("argument 'subset' will be ignored")
      if(!missing(na.action))
        warning("argument 'na.action' will be ignored")
      if(missing(start))
        start <- getInitial(formula, data=mf, control=control, trace=trace)
      startEnv <- new.env(hash = FALSE, parent = environment(formula)) # small
      if(is.matrix(start)) {
        for(i in colnames(start)) {
          startEnv[[i]] <- (start[1L, i] + start[2L, i]) / 2
        }
      } else {
        for (i in names(start))
          startEnv[[i]] <- start[[i]]
      }
      rhs <- eval(formula[[3L]], data, startEnv)
      n <- nrow(rhs)
      ## mimic what model.frame.default does
      wts <- if (!mWeights) eval(substitute(weights), data, environment(formula)) else NULL
    }
    else {
      vNms <- varNames[varIndex]
      if(any(nEQ <- vNms != make.names(vNms))) vNms[nEQ] <- paste0("`", vNms[nEQ], "`")
      mf$formula <-  # replace by one-sided linear model formula
        as.formula(paste("~", paste(vNms, collapse = "+")),
                   env = environment(formula))
      mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <-
        mf$fn <- mf$jac <- mf$fvv <- mf$lower <- mf$upper <- mf$loss <- NULL
      ## need stats:: for non-standard evaluation
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights) model.weights(mf) else NULL
    }
    if(!is.null(wts)) {
      stopifnot(
        "'weights' should be numeric" =  is.numeric(wts),
        "missing weights not allowed" = !any(is.na(wts))
      )
      if(is.matrix(wts)) {
        ## error if non-positive definite
        .swts <- t(chol(wts))
      } else {
        stopifnot("non-positive weights not allowed" = all(wts > 0))
        .swts <- sqrt(wts)
      }
    } else {
      .swts <- NULL
    }
  }
  else {
    stop("no data variables present")
  }

  ## set up iteration
  if(missing(start)) {
    start <- getInitial(formula, data=mf, control=control, trace=trace)
    .has_start <- structure(rep(TRUE, times = length(start)), names = names(start)) ## not used (single-start)
  }
  for(var in varNames[!varIndex])
    mf[[var]] <- eval(as.name(var), data, env)
  varNamesRHS <- varNamesRHS[ varNamesRHS %in% varNames[varIndex] ]

  ## parameter constraints
  if(!missing(lower) || !missing(upper)) {
    if(missing(lower)) lower <- structure(rep(-Inf, length(pnames)), names = pnames)
    if(missing(upper)) upper <- structure(rep(Inf, length(pnames)), names = pnames)
    if(is.list(lower)) {
      if(any(lengths(lower) != 1L))
        stop("List elements of 'lower' must have exactly length 1 to specify parameter lower bounds")
      lower <- unlist(lower)
    }
    if(is.list(upper)) {
      if(any(lengths(upper) != 1L))
        stop("List elements of 'upper' must have exactly length 1 to specify parameter upper bounds")
      upper <- unlist(upper)
    }
    if(is.null(names(lower)) && identical(length(lower), 1L) && length(pnames) > 1L)
      lower <- structure(rep(lower, length(pnames)), names = pnames)
    if(is.null(names(upper)) && identical(length(upper), 1L) && length(pnames) > 1L)
      upper <- structure(rep(upper, length(pnames)), names = pnames)
    if(is.null(names(lower)) || any(is.na(match(names(lower), pnames))))
      stop("Failed to match parameter names between 'start' and 'lower'")
    if(is.null(names(upper)) || any(is.na(match(names(upper), pnames))))
      stop("Failed to match parameter names between 'start' and 'upper'")
    if(is.matrix(start))
      startnames <- colnames(start)
    else
      startnames <- names(start)
    .lupars <- matrix(rep(c(-Inf, Inf), times = length(startnames)), nrow = 2L, ncol = length(startnames), dimnames = list(NULL, startnames))
    .lupars[1, match(names(lower), pnames)] <- unname(lower)
    .lupars[2, match(names(upper), pnames)] <- unname(upper)
    if(any(.lupars[1L, ] > .lupars[2L, ]))
      stop("Parameter lower bounds cannot be larger than upper bounds")
    if(is.matrix(start)) {
      if(any(!.has_start[1, ])) {
        start1 <- start[1, !.has_start[1, ]]
        start[1, !.has_start[1, ]] <- pmax(start1, .lupars[1, !.has_start[1, ]])
        start[2, !.has_start[2, ]] <- start[2, !.has_start[2, ]] + (start[1, !.has_start[1, ]] - start1)
      }
      if(any(!.has_start[2, ])) {
        start[2, !.has_start[2, ]] <- pmin(start[2, !.has_start[2, ]], .lupars[2, !.has_start[2, ]])
      }
      if(any(start[1, ] < .lupars[1, ] | start[2, ] > .lupars[2, ]))
        stop("Starting parameter ranges must be contained within 'lower' and 'upper' bounds")
    } else if(!is.matrix(start) && any(vapply(startnames, function(nm) start[[nm]] < .lupars[1, nm] || start[[nm]] > .lupars[2, nm], logical(1))))
      stop("Starting parameters must be contained within 'lower' and/or 'upper' bounds")
    if(all(is.infinite(.lupars))) {
      .lupars <- NULL
    }
  } else {
    .lupars <- NULL
  }

  ## function call
  .fn <- function(par, .data = mf) eval(formula[[3L]], envir = c(as.list(par), .data))
  if(is.matrix(start)) {
    .fcall <- tryCatch(.fn((start[1L, ] + start[2L, ]) / 2), error = function(err) err)
  } else {
    .fcall <- tryCatch(.fn(start), error = function(err) err)
  }
  if(inherits(.fcall, "error"))
    stop(sprintf("failed to evaluate 'fn' at starting values: %s", .fcall$message))
  if(!is.numeric(.fcall))
    stop("'fn' failed to return a numeric vector at starting values")
  if(any(is.na(.fcall)))
    stop("missing values returned by 'fn' at starting values")

  .lhs <- eval(formula[[2L]], envir = mf)

  ## n >= p
  if(length(.lhs) < length(pnames)) {
    stop("negative residual degrees of freedom, cannot fit a model with less observations than parameters")
  }

  ## jac call
  if(!is.function(jac) && !is.null(attr(.fcall, "gradient"))) {
    jac <- function(par, .data = mf) attr(.fn(par, .data), "gradient")
  } else if(isTRUE(jac)) {
    jac <- NULL
    .exprjac <- tryCatch(stats::deriv(formula[[3L]], namevec = pnames), error = function(err) err)
    if(inherits(.exprjac, "error")) {
      warning(sprintf("failed to symbolically derive 'jac': %s", .exprjac$message))
    } else if(is.expression(.exprjac)){
      jac <- function(par, .data = mf) {
        grad <- eval(.exprjac, envir = c(as.list(par), .data))
        attr(grad, "gradient")
      }
    }
  }

  if(is.function(jac)) {
    .jac <- function(par) jac(par, ...)
    if(is.matrix(start)) {
      .dfcall <- tryCatch(.jac((start[1L, ] + start[2L, ]) / 2), error = function(err) err)
    } else {
      .dfcall <- tryCatch(.jac(start), error = function(err) err)
    }
    if(inherits(.dfcall, "error"))
      stop(sprintf("failed to evaluate 'jac' at starting values: %s", .dfcall$message))
    if(!is.numeric(.dfcall) || !is.matrix(.dfcall) || !identical(dim(.dfcall), c(length(.lhs), length(pnames))))
      stop("'jac' failed to return a numeric matrix of expected dimensions at starting values")
    if(any(is.na(.dfcall)))
      stop("missing values returned by 'jac' at starting values")
  } else {
    .jac <- NULL
  }

  ## fvv call
  .fvv <- NULL
  if(identical(algorithm, "lmaccel")) {
    if(!is.function(fvv) && !is.null(attr(.fcall, "hessian"))) {
      fvv <- function(par, v) {
        hess <- attr(.fn(par), "hessian")
        c(matrix(hess, nrow = nrow(hess), ncol = ncol(hess) * ncol(hess)) %*% c(outer(v, v)))
      }
    } else if(isTRUE(fvv)) {
      fvv <- NULL
      .exprfvv <- tryCatch(stats::deriv(formula[[3L]], namevec = pnames, hessian = TRUE), error = function(err) err)
      if(inherits(.exprfvv, "error")) {
        warning(sprintf("failed to symbolically derive 'fvv': %s", .exprfvv$message))
      } else if(is.expression(.exprfvv)){
        fvv <- function(par, v) {
          grad <- eval(.exprfvv, envir = c(as.list(par), mf))
          hess <- attr(grad, "hessian")
          c(matrix(hess, nrow = nrow(hess), ncol = ncol(hess) * ncol(hess)) %*% c(outer(v, v)))
        }
      }
    }

    if(is.function(fvv)) {
      .fvv <- function(par, v) fvv(par, v, ...)
      if(is.matrix(start)) {
        .fvvcall <- tryCatch(.fvv((start[1L, ] + start[2L, ]) / 2, structure(rep(1, ncol(start)), names = colnames(start))),
                             error = function(err) err)
      } else {
        .fvvcall <- tryCatch(.fvv(start, structure(rep(1, length(start)), names = names(start))), error = function(err) err)
      }
      if(inherits(.fvvcall, "error"))
        stop(sprintf("failed to evaluate 'fvv' at starting values: %s", .fvvcall$message))
      if(!is.numeric(.fvvcall) || !identical(length(.fvvcall), length(.lhs)))
        stop("'fvv' failed to return a numeric vector equal in length to 'y' at starting values")
      if(any(is.na(.fvvcall)))
        stop("missing values returned by 'fvv' at starting values")
    }
  }

  ## loss function
  if(is.character(loss)) {
    .loss_config <- gsl_nls_loss(rho = loss)
  } else {
    .loss_config <- do.call(gsl_nls_loss, args = as.list(loss))
  }
  .loss_config$rho <- match(.loss_config$rho, c("default", "huber", "barron", "bisquare", "welsh", "optimal", "hampel", "ggw", "lqq"), nomatch = 1L) - 1L

  ## control arguments
  trace <- isTRUE(trace)
  .ctrl <- do.call(gsl_nls_control, args = if(!missing(control)) as.list(control) else list())
  .ctrl$scale <- match.arg(.ctrl$scale, c("more", "levenberg", "marquardt"))
  .ctrl$solver <- match.arg(.ctrl$solver, c("qr", "cholesky", "svd"))
  .ctrl$fdtype <- match.arg(.ctrl$fdtype, c("forward", "center"))
  stopifnot(
    is.numeric(.ctrl$maxiter), length(.ctrl$maxiter) == 1, .ctrl$maxiter >= 1,
    is.numeric(.ctrl$factor_up), length(.ctrl$factor_up) == 1, .ctrl$factor_up > 0,
    is.numeric(.ctrl$factor_down), length(.ctrl$factor_down) == 1, .ctrl$factor_down > 0,
    is.numeric(.ctrl$avmax), length(.ctrl$avmax) == 1, .ctrl$avmax > 0,
    is.numeric(.ctrl$h_df), length(.ctrl$h_df) == 1, .ctrl$h_df > 0,
    is.numeric(.ctrl$h_fvv), length(.ctrl$h_fvv) == 1, .ctrl$h_fvv > 0,
    is.numeric(.ctrl$xtol), length(.ctrl$xtol) == 1, .ctrl$xtol > 0,
    is.numeric(.ctrl$ftol), length(.ctrl$ftol) == 1, .ctrl$ftol > 0,
    is.numeric(.ctrl$gtol), length(.ctrl$gtol) == 1, .ctrl$gtol > 0,
    is.numeric(.ctrl$mstart_n), length(.ctrl$mstart_n) == 1, .ctrl$mstart_n >= 1,
    is.numeric(.ctrl$mstart_p), length(.ctrl$mstart_p) == 1, .ctrl$mstart_p >= 1,
    is.numeric(.ctrl$mstart_q), length(.ctrl$mstart_q) == 1, .ctrl$mstart_q >= 1,
    is.numeric(.ctrl$mstart_r), length(.ctrl$mstart_r) == 1, .ctrl$mstart_r > 1,
    is.numeric(.ctrl$mstart_s), length(.ctrl$mstart_s) == 1, .ctrl$mstart_s >= 1,
    is.numeric(.ctrl$mstart_tol), length(.ctrl$mstart_tol) == 1, .ctrl$mstart_tol > 0,
    is.numeric(.ctrl$mstart_maxiter), length(.ctrl$mstart_maxiter) == 1, .ctrl$mstart_maxiter >= 1,
    is.numeric(.ctrl$mstart_maxstart), length(.ctrl$mstart_maxstart) == 1, .ctrl$mstart_maxstart >= 1,
    is.numeric(.ctrl$mstart_minsp), length(.ctrl$mstart_minsp) == 1, .ctrl$mstart_minsp >= 1,
    is.numeric(.ctrl$irls_maxiter), length(.ctrl$irls_maxiter) == 1, .ctrl$irls_maxiter >= 1,
    is.numeric(.ctrl$irls_xtol), length(.ctrl$irls_xtol) == 1, .ctrl$irls_xtol > 0
  )
  .ctrl_int <- c(
    as.integer(.ctrl$maxiter),
    isTRUE(trace),
    match(algorithm, c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D")) - 1L,
    match(.ctrl$scale, c("more", "levenberg", "marquardt")) - 1L,
    match(.ctrl$solver, c("qr", "cholesky", "svd")) - 1L,
    match(.ctrl$fdtype, c("forward", "center")) - 1L,
    as.integer(.ctrl$mstart_n),
    as.integer(.ctrl$mstart_p),
    as.integer(.ctrl$mstart_q),
    as.integer(.ctrl$mstart_s),
    as.integer(.ctrl$mstart_maxiter),
    as.integer(.ctrl$mstart_maxstart),
    as.integer(.ctrl$mstart_minsp),
    as.integer(!is.list(start)),
    as.integer(.ctrl$irls_maxiter)
  )
  .ctrl_dbl <- unlist(.ctrl[c("factor_up", "factor_down", "avmax", "h_df", "h_fvv", "xtol", "ftol", "gtol", "mstart_r", "mstart_tol", "irls_xtol")])
  if(!all(.has_start)) {
    .ctrl_dbl[["mstart_r"]] <- 10 * .ctrl_dbl[["mstart_r"]]
  }

  ## optimize
  cFit <- .Call(
    C_nls,
    .fn, .lhs, .jac, .fvv, environment(), start, .swts, .lupars, .ctrl_int, .ctrl_dbl, .has_start, .loss_config,
    PACKAGE = "gslnls"
  )

  ## convert to nls object
  m <- nlsModel(formula, mf, cFit, .swts, jac)

  convInfo <- list(
    isConv = as.logical(!cFit$conv),
    finIter = cFit$niter,
    finTol = cFit$ssrtol,
    nEval = cFit$neval,
    trsName = paste("multifit", cFit$algorithm, sep = "/"),
    stopCode = cFit$conv,
    stopMessage = cFit$status
  )

  nls.out <- list(m = m, data = substitute(data), convInfo = convInfo, call = match.call())

  nls.out$call$algorithm <- algorithm
  nls.out$call$control <- nls.control() ## needed for profiler
  nls.out$call$trace <- trace
  nls.out$call$formula <- nls.out$m$formula() ## needed for external generics
  if(trace) {
    nls.out$partrace <- cFit$partrace[seq_len(cFit$niter + 1L), , drop = FALSE]
    nls.out$devtrace <- cFit$ssrtrace[seq_len(cFit$niter + 1L)]
  }
  nls.out$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  nls.out$control <- .ctrl
  if(model)
    nls.out$model <- mf
  if(!mWeights)
    nls.out$weights <- wts
  if(!is.null(.lupars)) {
    nls.out$lower <- .lupars[1L, ]
    nls.out$upper <- .lupars[2L, ]
  }
  if(!is.null(cFit$irls)) {
    nls.out$irls <- cFit$irls
    nls.out$irls$irls_conv <- as.logical(!cFit$irls$irls_conv)
  }

  class(nls.out) <- c("gsl_nls", "nls")

  return(nls.out)

}

#' @describeIn gsl_nls
#' If \code{fn} is a \code{function}, the first argument must be the vector of parameters and
#' the function should return a numeric vector containing the nonlinear model evaluations at
#' the provided parameter and predictor or covariate vectors. In addition, the argument \code{y}
#' needs to contain the numeric vector of observed responses, equal in length to the numeric
#' vector returned by \code{fn}. The returned list object is (only) of class \code{gsl_nls}.
#' Although the returned object is not of class \code{nls}, the following generic functions remain
#' applicable for an object of class \code{gsl_nls}: \code{anova}, \code{coef}, \code{confint}, \code{deviance},
#' \code{df.residual}, \code{fitted}, \code{formula}, \code{logLik}, \code{nobs}, \code{predict}, \code{print},
#' \code{residuals}, \code{summary}, \code{vcov}, \code{hatvalues}, \code{cooks.distance} and \code{weights}.
#' In addition, a method \code{confintd} is available for inference of derived parameters.
#' @export
gsl_nls.function <- function(fn, y, start,
                             algorithm = c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"),
                             loss = c("default", "huber", "barron", "bisquare", "welsh", "optimal", "hampel", "ggw", "lqq"),
                             control = gsl_nls_control(), lower, upper, jac = NULL, fvv = NULL,
                             trace = FALSE, weights, ...) {

  ## algorithm
  algorithm <- match.arg(algorithm, c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"))

  ## starting values
  .start <- start   ## do not modify start argument
  if(is.list(.start) && any(lengths(.start) > 1)) {
    if(!all(is.element(lengths(.start), c(1L, 2L)))) {
      stop("List elements of 'start' must be of length 1 or 2 to specify (multi-start) parameter values or ranges")
    }
    if(any(len1 <- lengths(.start) < 2)) {
      .start[len1] <- lapply(.start[len1], rep, times = 2L)
    }
    .start <- as.matrix(data.frame(.start))
  } else if(is.matrix(.start) && nrow(.start) != 2L) {
    stop("Matrix 'start' must have exactly 2 rows to specify (multi-start) parameter ranges")
  }
  ## Handle infinite/missing start values
  if(is.vector(.start)) {
    if(any((napars <- sapply(.start, Negate(is.finite))))) {
      .start <- vapply(names(.start), function(nm) if(napars[[nm]]) c(-0.1, 0.75) else rep(.start[[nm]], 2L), numeric(2))
      .has_start <- sapply(!napars, rep, times = 2L)
    } else {
      .has_start <- !napars  ## not used (single-start)
    }
  } else if(is.matrix(.start)) {
    if(any(napars <- apply(.start, 2L, Negate(is.finite)))) {
      .start[1L, napars[1L, ]]  <- -0.1
      .start[2L, napars[2L, ]] <- 0.75
    }
    .has_start <- !napars
  } else {
    stop("Unrecognized format for starting values and/or ranges defined in 'start'")
  }
  ## Check start ranges
  if(is.matrix(.start)) {
    if(ncol(.start) > 1229L)
      stop("GSL quasi-random Halton sequences in multi-start are only available up to 1229 parameters")
    if(any(.start[1L, ] > .start[2L, ]))
      stop("Multi-start parameter lower bounds cannot be larger than upper bounds")
    else if(!any(.start[1L, ] < .start[2L, ]))
      start <- .start <- start[1L, ] ## single-start
  }

  ## function call
  if(!is.numeric(y))
    stop("'y' should be a numeric response vector")

  ## n >= p
  p <- ifelse(is.matrix(.start), ncol(.start), length(.start))
  if(length(y) < p) {
    stop("negative residual degrees of freedom, cannot fit a model with less observations than parameters")
  }

  ## parameter constraints
  if(!missing(lower) || !missing(upper)) {
    if(is.matrix(.start))
      pnames <- colnames(.start)
    else
      pnames <- names(.start)
    if(missing(lower)) lower <- structure(rep(-Inf, length(pnames)), names = pnames)
    if(missing(upper)) upper <- structure(rep(Inf, length(pnames)), names = pnames)
    if(is.list(lower)) {
      if(any(lengths(lower) != 1L))
        stop("List elements of 'lower' must have exactly length 1 to specify parameter lower bounds")
      lower <- unlist(lower)
    }
    if(is.list(upper)) {
      if(any(lengths(upper) != 1L))
        stop("List elements of 'upper' must have exactly length 1 to specify parameter upper bounds")
      upper <- unlist(upper)
    }
    if(is.null(names(lower)) && identical(length(lower), 1L) && length(pnames) > 1L)
      lower <- structure(rep(lower, length(pnames)), names = pnames)
    if(is.null(names(upper)) && identical(length(upper), 1L) && length(pnames) > 1L)
      upper <- structure(rep(upper, length(pnames)), names = pnames)
    if(is.null(names(lower)) || any(is.na(match(names(lower), pnames))))
      stop("Failed to match parameter names between 'start' and 'lower'")
    if(is.null(names(upper)) || any(is.na(match(names(upper), pnames))))
      stop("Failed to match parameter names between 'start' and 'upper'")
    .lupars <- matrix(rep(c(-Inf, Inf), times = length(pnames)), nrow = 2L, ncol = length(pnames), dimnames = list(NULL, pnames))
    .lupars[1, match(names(lower), pnames)] <- unname(lower)
    .lupars[2, match(names(upper), pnames)] <- unname(upper)
    if(any(.lupars[1L, ] > .lupars[2L, ]))
      stop("Parameter lower bounds cannot be larger than upper bounds")
    if(is.matrix(.start)) {
      if(any(!.has_start[1, ])) {
        start1 <- .start[1, !.has_start[1, ]]
        .start[1, !.has_start[1, ]] <- pmax(start1, .lupars[1, !.has_start[1, ]])
        .start[2, !.has_start[2, ]] <- .start[2, !.has_start[2, ]] + (.start[1, !.has_start[1, ]] - start1)
      }
      if(any(!.has_start[2, ])) {
        .start[2, !.has_start[2, ]] <- pmin(.start[2, !.has_start[2, ]], .lupars[2, !.has_start[2, ]])
      }
      if(any(.start[1, ] < .lupars[1, ] | .start[2, ] > .lupars[2, ]))
        stop("Starting parameter ranges must be contained within 'lower' and 'upper' bounds")
    } else if(!is.matrix(.start) && any(vapply(pnames, function(nm) .start[[nm]] < .lupars[1, nm] || .start[[nm]] > .lupars[2, nm], logical(1))))
      stop("Starting parameters must be contained within 'lower' and/or 'upper' bounds")
    if(all(is.infinite(.lupars))) {
      .lupars <- NULL
    }
  } else {
    .lupars <- NULL
  }

  ## function call
  .fn <- function(par) fn(par, ...)
  if(is.matrix(.start) && is.list(start)) {
    .fcall <- tryCatch(.fn(as.list((.start[1L, ] + .start[2L, ]) / 2)), error = function(err) err)
  } else if(is.matrix(.start)) {
    .fcall <- tryCatch(.fn((.start[1L, ] + .start[2L, ]) / 2), error = function(err) err)
  } else {
    .fcall <- tryCatch(.fn(.start), error = function(err) err)
  }
  if(inherits(.fcall, "error"))
    stop(sprintf("failed to evaluate 'fn' at starting values: %s", .fcall$message))
  if(!is.numeric(.fcall) || !identical(length(.fcall), length(y)))
    stop("'fn' failed to return a numeric vector equal in length to 'y' at starting values")
  if(any(is.na(.fcall)))
    stop("missing values returned by 'fn' at starting values")

  ## jac call
  if(!is.function(jac) && !is.null(attr(.fcall, "gradient"))) {
    jac <- function(par, ...) attr(fn(par, ...), "gradient")
  }

  if(is.function(jac)) {
    .jac <- function(par) jac(par, ...)
    if(is.matrix(.start) && is.list(start)) {
      .dfcall <- tryCatch(.jac(as.list((.start[1L, ] + .start[2L, ]) / 2)), error = function(err) err)
    } else if(is.matrix(.start)) {
      .dfcall <- tryCatch(.jac((.start[1L, ] + .start[2L, ]) / 2), error = function(err) err)
    } else {
      .dfcall <- tryCatch(.jac(.start), error = function(err) err)
    }
    if(inherits(.dfcall, "error"))
      stop(sprintf("failed to evaluate 'jac' at starting values: %s", .dfcall$message))
    if(!is.numeric(.dfcall) || !is.matrix(.dfcall) || !identical(dim(.dfcall), c(length(y), p)))
      stop("'jac' failed to return a numeric matrix of expected dimensions at starting values")
    if(any(is.na(.dfcall)))
      stop("missing values returned by 'jac' at starting values")
  } else {
    .jac <- NULL
  }

  ## fvv call
  .fvv <- NULL
  if(identical(algorithm, "lmaccel")) {
    if(!is.function(fvv) && !is.null(attr(.fcall, "hessian"))) {
      fvv <- function(par, v, ...) {
        hess <- attr(fn(par, ...), "hessian")
        c(matrix(hess, nrow = nrow(hess), ncol = ncol(hess) * ncol(hess)) %*% c(outer(v, v)))
      }
    }
    if(is.function(fvv)) {
      .fvv <- function(par, v) fvv(par, v, ...)
      if(is.matrix(.start) && is.list(start)) {
        .fvvcall <- tryCatch(.fvv(as.list((.start[1L, ] + .start[2L, ]) / 2), structure(rep(1, p), names = colnames(.start))), error = function(err) err)
      } else if(is.matrix(.start)) {
        .fvvcall <- tryCatch(.fvv((.start[1L, ] + .start[2L, ]) / 2, structure(rep(1, p), names = colnames(.start))), error = function(err) err)
      } else {
        .fvvcall <- tryCatch(.fvv(.start, structure(rep(1, p), names = names(.start))), error = function(err) err)
      }
      if(inherits(.fvvcall, "error"))
        stop(sprintf("failed to evaluate 'fvv' at starting values: %s", .fvvcall$message))
      if(!is.numeric(.fvvcall) || !identical(length(.fvvcall), length(y)))
        stop("'fvv' failed to return a numeric vector equal in length to 'y' at starting values")
      if(any(is.na(.fvvcall)))
        stop("missing values returned by 'fvv' at starting values")
    }
  }

  ## loss function
  if(is.character(loss)) {
    .loss_config <- gsl_nls_loss(rho = loss)
  } else {
    .loss_config <- do.call(gsl_nls_loss, args = as.list(loss))
  }
  .loss_config$rho <- match(.loss_config$rho, c("default", "huber", "barron", "bisquare", "welsh", "optimal", "hampel", "ggw", "lqq"), nomatch = 1L) - 1L

  ## control arguments
  trace <- isTRUE(trace)
  .ctrl <- do.call(gsl_nls_control, args = if(!missing(control)) as.list(control) else list())
  .ctrl$scale <- match.arg(.ctrl$scale, c("more", "levenberg", "marquardt"))
  .ctrl$solver <- match.arg(.ctrl$solver, c("qr", "cholesky", "svd"))
  .ctrl$fdtype <- match.arg(.ctrl$fdtype, c("forward", "center"))
  stopifnot(
    is.numeric(.ctrl$maxiter), length(.ctrl$maxiter) == 1, .ctrl$maxiter >= 1,
    is.numeric(.ctrl$factor_up), length(.ctrl$factor_up) == 1, .ctrl$factor_up > 0,
    is.numeric(.ctrl$factor_down), length(.ctrl$factor_down) == 1, .ctrl$factor_down > 0,
    is.numeric(.ctrl$avmax), length(.ctrl$avmax) == 1, .ctrl$avmax > 0,
    is.numeric(.ctrl$h_df), length(.ctrl$h_df) == 1, .ctrl$h_df > 0,
    is.numeric(.ctrl$h_fvv), length(.ctrl$h_fvv) == 1, .ctrl$h_fvv > 0,
    is.numeric(.ctrl$xtol), length(.ctrl$xtol) == 1, .ctrl$xtol > 0,
    is.numeric(.ctrl$ftol), length(.ctrl$ftol) == 1, .ctrl$ftol > 0,
    is.numeric(.ctrl$gtol), length(.ctrl$gtol) == 1, .ctrl$gtol > 0,
    is.numeric(.ctrl$mstart_n), length(.ctrl$mstart_n) == 1, .ctrl$mstart_n >= 1,
    is.numeric(.ctrl$mstart_p), length(.ctrl$mstart_p) == 1, .ctrl$mstart_p >= 1,
    is.numeric(.ctrl$mstart_q), length(.ctrl$mstart_q) == 1, .ctrl$mstart_q >= 1,
    is.numeric(.ctrl$mstart_r), length(.ctrl$mstart_r) == 1, .ctrl$mstart_r > 1,
    is.numeric(.ctrl$mstart_s), length(.ctrl$mstart_s) == 1, .ctrl$mstart_s >= 1,
    is.numeric(.ctrl$mstart_tol), length(.ctrl$mstart_tol) == 1, .ctrl$mstart_tol > 0,
    is.numeric(.ctrl$mstart_maxiter), length(.ctrl$mstart_maxiter) == 1, .ctrl$mstart_maxiter >= 1,
    is.numeric(.ctrl$mstart_maxstart), length(.ctrl$mstart_maxstart) == 1, .ctrl$mstart_maxstart >= 1,
    is.numeric(.ctrl$mstart_minsp), length(.ctrl$mstart_minsp) == 1, .ctrl$mstart_minsp >= 1,
    is.numeric(.ctrl$irls_maxiter), length(.ctrl$irls_maxiter) == 1, .ctrl$irls_maxiter >= 1,
    is.numeric(.ctrl$irls_xtol), length(.ctrl$irls_xtol) == 1, .ctrl$irls_xtol > 0
  )
  .ctrl_int <- c(
    as.integer(.ctrl$maxiter),
    isTRUE(trace),
    match(algorithm, c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D")) - 1L,
    match(.ctrl$scale, c("more", "levenberg", "marquardt")) - 1L,
    match(.ctrl$solver, c("qr", "cholesky", "svd")) - 1L,
    match(.ctrl$fdtype, c("forward", "center")) - 1L,
    as.integer(.ctrl$mstart_n),
    as.integer(.ctrl$mstart_p),
    as.integer(.ctrl$mstart_q),
    as.integer(.ctrl$mstart_s),
    as.integer(.ctrl$mstart_maxiter),
    as.integer(.ctrl$mstart_maxstart),
    as.integer(.ctrl$mstart_minsp),
    as.integer(!is.list(start)),
    as.integer(.ctrl$irls_maxiter)
  )
  .ctrl_dbl <- unlist(.ctrl[c("factor_up", "factor_down", "avmax", "h_df", "h_fvv", "xtol", "ftol", "gtol", "mstart_r", "mstart_tol", "irls_xtol")])
  if(!all(.has_start)) {
    .ctrl_dbl[["mstart_r"]] <- 10 * .ctrl_dbl[["mstart_r"]]
  }

  ## weights
  if(!missing(weights)) {
    stopifnot(
      "'weights' should be numeric" =  is.numeric(weights),
      "missing weights not allowed" =  !any(is.na(weights))
    )
    if(is.matrix(weights)) {
      stopifnot(
        "dimensions of 'weights' should be equal to length of 'y'" =  all(dim(weights) == length(y))
      )
      ## error if non-positive definite
      .swts <- t(chol(weights))
    } else {
      stopifnot(
        "'weights' should be equal in length to 'y'" = identical(length(weights), length(y)),
        "non-positive weights not allowed" =  all(weights > 0)
      )
      .swts <- sqrt(weights)
    }
  } else {
    .swts <- NULL
  }

  ## optimize
  cFit <- .Call(
    C_nls,
    .fn, y, .jac, .fvv, environment(), .start, .swts, .lupars, .ctrl_int, .ctrl_dbl, .has_start, .loss_config,
    PACKAGE = "gslnls"
  )

  m <- gslModel(
    fn, y, cFit,
    if(!all(.has_start) && is.list(start)) as.list(.start[1L, ]) else if(!all(.has_start) || is.matrix(start)) .start[1L, ] else .start,
    .swts, jac, ...
  )

  ## mimick nls object
  convInfo <- list(
    isConv = as.logical(!cFit$conv),
    finIter = cFit$niter,
    finTol = cFit$ssrtol,
    nEval = cFit$neval,
    trsName = paste("multifit", cFit$algorithm, sep = "/"),
    stopCode = cFit$conv,
    stopMessage = cFit$status
  )

  nls.out <- list(m = m, data = c(list(y = y), list(...)), convInfo = convInfo, call = match.call())
  nls.out$call$algorithm <- algorithm
  nls.out$call$trace <- trace
  nls.out$call$formula <- nls.out$m$formula() ## needed for external generics
  if(trace) {
    nls.out$partrace <- cFit$partrace[seq_len(cFit$niter + 1L), , drop = FALSE]
    nls.out$devtrace <- cFit$ssrtrace[seq_len(cFit$niter + 1L)]
  }
  nls.out$control <- .ctrl
  if(!missing(weights))
    nls.out$weights <- weights
  if(!is.null(.lupars)) {
    nls.out$lower <- .lupars[1L, ]
    nls.out$upper <- .lupars[2L, ]
  }
  if(!is.null(cFit$irls)) {
    nls.out$irls <- cFit$irls
    nls.out$irls$irls_conv <- as.logical(!cFit$irls$irls_conv)
  }
  class(nls.out) <- "gsl_nls"

  return(nls.out)

}

#' Tunable Nonlinear Least Squares iteration parameters
#'
#' Allow the user to tune the characteristics of the \code{\link{gsl_nls}} and \code{\link{gsl_nls_large}}
#' nonlinear least squares algorithms.
#'
#' @param maxiter positive integer, termination occurs when the number of iterations reaches \code{maxiter}.
#' @param scale character, scaling method or damping strategy determining the diagonal scaling matrix D. The following options
#' are supported:
#' \itemize{
#' \item \code{"more"} Moré rescaling (default). This method makes the problem scale-invariant and has
#' been proven effective on a large class of problems.
#' \item \code{"levenberg"} Levenberg rescaling. This method has also proven effective on a large class of problems,
#' but is not scale-invariant. It may perform better for problems susceptible to \emph{parameter evaporation} (parameters going to infinity).
#' \item \code{"marquardt"} Marquardt rescaling. This method is scale-invariant, but it is generally
#' considered inferior to both the Levenberg and Moré strategies.
#' }
#' @param solver character, method used to solve the linear least squares system resulting as a subproblem in each iteration.
#' For large-scale problems fitted with \code{\link{gsl_nls_large}}, the Cholesky solver (\code{"cholesky"}) is always selected
#' and this parameter is not used. For least squares problems fitted with \code{\link{gsl_nls}} the following choices are supported:
#' \itemize{
#' \item \code{"qr"} QR decomposition of the Jacobian (default). This method will produce reliable solutions in cases
#' where the Jacobian is rank deficient or near-singular but does require more operations than the Cholesky method.
#' \item \code{"cholesky"} Cholesky decomposition of the Jacobian. This method is faster than the QR approach, however
#' it is susceptible to numerical instabilities if the Jacobian matrix is rank deficient or near-singular.
#' \item \code{"svd"} SVD decomposition of the Jacobian. This method will produce the most reliable solutions for
#' ill-conditioned Jacobians but is also the slowest.
#' }
#' @param fdtype character, method used to numerically approximate the Jacobian and/or second-order derivatives
#' when geodesic acceleration is used. Either \code{"forward"} for forward finite differencing or \code{"center"}
#' for centered finite differencing. For least squares problems solved with \code{\link{gsl_nls_large}}, numerical
#' approximation of the Jacobian matrix is not available and this parameter is only used to numerically approximate
#' the second-order derivatives (if geodesic acceleration is used).
#' @param factor_up numeric factor by which to increase the trust region radius when a search step is accepted.
#' Too large values may destabilize the search, too small values slow down the search, defaults to 2.
#' @param factor_down numeric factor by which to decrease the trust region radius when a search step is rejected.
#' Too large values may destabilize the search, too small values slow down the search, defaults to 3.
#' @param avmax numeric value, the ratio of the acceleration term to the velocity term when using geodesic acceleration to
#' solve the nonlinear least squares problem. Any steps with a ratio larger than \code{avmax} are rejected, defaults to 0.75.
#' For problems which experience difficulty converging, this threshold could be lowered.
#' @param h_df numeric value, the step size for approximating the Jacobian matrix with finite differences, defaults to \code{sqrt(.Machine$double.eps)}.
#' @param h_fvv numeric value, the step size for approximating the second directional derivative when geodesic acceleration
#' is used to solve the nonlinear least squares problem, defaults to 0.02. This is only used if no analytic second
#' directional derivative (\code{fvv}) is specified in \code{\link{gsl_nls}} or \code{\link{gsl_nls_large}}.
#' @param xtol numeric value, termination occurs when the relative change in parameters between iterations is \code{<= xtol}.
#' A general guideline for selecting the step tolerance is to choose \code{xtol = 10^(-d)} where \code{d} is the number of accurate
#' decimal digits desired in the parameters, defaults to \code{sqrt(.Machine$double.eps)}.
#' @param ftol numeric value, termination occurs when the relative change in sum of squared residuals between iterations is \code{<= ftol},
#' defaults to \code{sqrt(.Machine$double.eps)}.
#' @param gtol numeric value, termination occurs when the relative size of the gradient of the sum of squared residuals is \code{<= gtol},
#' indicating a local minimum, defaults to \code{sqrt(.Machine$double.eps)}
#' @param mstart_n positive integer, number of quasi-random points drawn in each major iteration, parameter \code{N} in Hickernell and Yuan (1997). Default is 30.
#' @param mstart_p positive integer, number of iterations of inexpensive local search to concentrate the sample, parameter \code{p} in Hickernell and Yuan (1997). Default is 5.
#' @param mstart_q positive integer, number of points retained in the concentrated sample, parameter \code{q} in Hickernell and Yuan (1997). Default is \code{mstart_n \%/\% 10}..
#' @param mstart_r positive integer, scaling factor of number of stationary points determining when the multi-start algorithm terminates, parameter \code{r} in Hickernell and Yuan (1997). Default is 4.
#' If the starting ranges for one or more parameters are unbounded and updated dynamically, \code{mstart_r} is multiplied by a factor 10 to avoid early termination.
#' @param mstart_s positive integer, minimum number of iterations a point needs to be retained before starting an efficient local search, parameter \code{s} in Hickernell and Yuan (1997). Default is 2.
#' @param mstart_tol numeric value, multiplicative tolerance \code{(1 + mstart_tol)} used as criterion to start an efficient local search (epsilon in Algorithm 2.1, Hickernell and Yuan (1997)).
#' @param mstart_maxiter positive integer, maximum number of iterations in the efficient local search algorithm (Algorithm B, Hickernell and Yuan (1997)), defaults to 10.
#' @param mstart_maxstart positive integer, minimum number of major iterations (Algorithm 2.1, Hickernell and Yuan (1997)) before the multi-start algorithm terminates, defaults to 250.
#' @param mstart_minsp positive integer, minimum number of detected stationary points before the multi-start algorithm terminates, defaults to 1.
#' @param irls_maxiter positive integer, maximum number of IRLS iterations, defaults to 50. Only used in case of a non-default loss function (\code{loss != "default"}) optimized by IRLS.
#' @param irls_xtol numeric value, termination of the IRLS procedure occurs when the relative change in parameters between IRLS iterations is \code{<= irls_xtol}, defaults to \code{.Machine$double.eps^(1/4)}.
#' Only used in case of a non-default loss function (\code{loss != "default"}) optimized by IRLS.
#' @param ... any additional arguments (currently not used).
#' @importFrom stats nls.control
#' @seealso \code{\link[stats]{nls.control}}
#' @seealso \url{https://www.gnu.org/software/gsl/doc/html/nls.html#tunable-parameters}
#' @note \code{ftol} is disabled in some versions of the GSL library.
#' @examples
#' ## default tuning parameters
#' gsl_nls_control()
#' @return A \code{list} with exactly twenty-three components:
#' \itemize{
#' \item maxiter
#' \item scale
#' \item solver
#' \item fdtype
#' \item factor_up
#' \item factor_down
#' \item avmax
#' \item h_df
#' \item h_fvv
#' \item xtol
#' \item ftol
#' \item gtol
#' \item mstart_n
#' \item mstart_p
#' \item mstart_q
#' \item mstart_r
#' \item mstart_s
#' \item mstart_tol
#' \item mstart_maxiter
#' \item mstart_maxstart
#' \item mstart_minsp
#' \item irls_maxiter
#' \item irls_xtol
#' }
#' with meanings as explained under 'Arguments'.
#' @references M. Galassi et al., \emph{GNU Scientific Library Reference Manual (3rd Ed.)}, ISBN 0954612078.
#' @references Hickernell, F.J. and Yuan, Y. (1997) \emph{“A simple multistart algorithm for global optimization”}, OR Transactions, Vol. 1 (2).
#' @export
gsl_nls_control <- function(maxiter = 100, scale = "more", solver = "qr",
                            fdtype = "forward", factor_up = 2, factor_down = 3, avmax = 0.75,
                            h_df = sqrt(.Machine$double.eps), h_fvv = 0.02, xtol = sqrt(.Machine$double.eps),
                            ftol = sqrt(.Machine$double.eps), gtol = sqrt(.Machine$double.eps),
                            mstart_n = 30, mstart_p = 5, mstart_q = mstart_n %/% 10, mstart_r = 4, mstart_s = 2,
                            mstart_tol = 0.25, mstart_maxiter = 10, mstart_maxstart = 250, mstart_minsp = 1,
                            irls_maxiter = 50, irls_xtol = .Machine$double.eps^.25, ...) {

  scale <- match.arg(scale, c("more", "levenberg", "marquardt"))
  solver <- match.arg(solver, c("qr", "cholesky", "svd"))
  fdtype <- match.arg(fdtype, c("forward", "center"))

  stopifnot(
    is.numeric(maxiter), length(maxiter) == 1, maxiter >= 1,
    is.numeric(factor_up), length(factor_up) == 1, factor_up > 0,
    is.numeric(factor_down), length(factor_down) == 1, factor_down > 0,
    is.numeric(avmax), length(avmax) == 1, avmax > 0,
    is.numeric(h_df), length(h_df) == 1, h_df > 0,
    is.numeric(h_fvv), length(h_fvv) == 1, h_fvv > 0,
    is.numeric(xtol), length(xtol) == 1, xtol > 0,
    is.numeric(ftol), length(ftol) == 1, ftol > 0,
    is.numeric(gtol), length(gtol) == 1, gtol > 0,
    is.numeric(mstart_n), length(mstart_n) == 1, mstart_n >= 1,
    is.numeric(mstart_p), length(mstart_p) == 1, mstart_p >= 1,
    is.numeric(mstart_q), length(mstart_q) == 1, mstart_q >= 1,
    is.numeric(mstart_r), length(mstart_r) == 1, mstart_r > 1,
    is.numeric(mstart_s), length(mstart_s) == 1, mstart_s >= 1,
    is.numeric(mstart_tol), length(mstart_tol) == 1, mstart_tol > 0,
    is.numeric(mstart_maxiter), length(mstart_maxiter) == 1, mstart_maxiter >= 1,
    is.numeric(mstart_maxstart), length(mstart_maxstart) == 1, mstart_maxstart >= 1,
    is.numeric(mstart_minsp), length(mstart_minsp) == 1, mstart_minsp >= 1,
    is.numeric(irls_maxiter), length(irls_maxiter) == 1, irls_maxiter >= 1,
    is.numeric(irls_xtol), length(irls_xtol) == 1, irls_xtol > 0
  )

  list(maxiter = as.integer(maxiter), scale = scale, solver = solver, fdtype = fdtype,
       factor_up = factor_up, factor_down = factor_down, avmax = avmax,
       h_df = h_df, h_fvv = h_fvv, xtol = xtol, ftol = ftol, gtol = gtol,
       mstart_n = as.integer(mstart_n), mstart_p = as.integer(mstart_p), mstart_q = as.integer(mstart_q),
       mstart_r = mstart_r, mstart_s = as.integer(mstart_s), mstart_tol = mstart_tol,
       mstart_maxiter = as.integer(mstart_maxiter), mstart_maxstart = as.integer(mstart_maxstart),
       mstart_minsp = as.integer(mstart_minsp), irls_maxiter = as.integer(irls_maxiter), irls_xtol = irls_xtol)

}

nlsModel <- function(form, data, cFit, swts, jac, upper=NULL) {
  ## thisEnv <- environment() # shared by all functions in the 'm' list; variable no longer needed
  env <- new.env(hash = TRUE, parent = environment(form))
  start <- cFit$par
  for(i in names(data)) env[[i]] <- data[[i]]
  ind <- as.list(start)
  parLength <- 0L
  for(i in names(ind) ) {
    temp <- start[[i]]
    storage.mode(temp) <- "double"
    env[[i]] <- temp
    ind[[i]] <- parLength + seq_along(temp)
    parLength <- parLength + length(temp)
  }
  getPars.noVarying <- function() unlist(mget(names(ind), env))
  getPars <- getPars.noVarying
  internalPars <- getPars()

  if(!is.null(upper)) upper <- rep_len(upper, parLength)  ## currently not used
  useParams <- rep_len(TRUE, parLength)
  lhs <- eval(form[[2L]], envir = env)
  rhs <- eval(form[[3L]], envir = env)
  .swts <- if(!missing(swts) && length(swts)) swts else rep_len(1, length(rhs))
  env$.swts <- .swts
  resid <- -cFit$resid
  dev <- cFit$ssr
  if(is.null(attr(rhs, "gradient"))) {
    getRHS.noVarying <- function() numericDeriv(form[[3L]], names(ind), env)
    getRHS <- getRHS.noVarying
    rhs <- getRHS()
  } else {
    getRHS.noVarying <- function() eval(form[[3L]], envir = env)
    getRHS <- getRHS.noVarying
  }
  dimGrad <- dim(attr(rhs, "gradient"))
  marg <- length(dimGrad)
  if(marg > 0L) {
    gradSetArgs <- vector("list", marg + 1L)
    for(i in 2L:marg)
      gradSetArgs[[i]] <- rep_len(TRUE, dimGrad[i-1L])
    useParams <- rep_len(TRUE, dimGrad[marg])
  } else {
    gradSetArgs <- vector("list", 2L)
    useParams <- rep_len(TRUE, length(attr(rhs, "gradient")))
  }
  npar <- length(useParams)
  gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
  gradCall <- switch(length(gradSetArgs) - 1L,
                     call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], drop = FALSE),
                     call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                          drop = FALSE),
                     call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                          gradSetArgs[[3L]], drop = FALSE),
                     call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                          gradSetArgs[[3L]], gradSetArgs[[4L]], drop = FALSE))
  getRHS.varying <- function() {
    ans <- getRHS.noVarying()
    attr(ans, "gradient") <- eval(gradCall)
    ans
  }
  if(length(gr <- attr(rhs, "gradient")) == 1L)
    attr(rhs, "gradient") <- gr <- as.vector(gr)

  if(!any(is.na(gr) | is.infinite(gr))) {
    QR <- if(is.matrix(.swts)) qr(.swts %*% gr) else qr(.swts * gr)
    qrDim <- min(dim(QR$qr))
    if(QR$rank < qrDim)
      warning("singular gradient matrix at parameter estimates")
    if(!is.null(cFit$irls)) {
      ## MM correction [Yohai, 1987. Thrm. 4.1]
      tau <- mean(cFit$irls$irls_psi^2) / mean(cFit$irls$irls_dpsi)^2
      QR$qr <- QR$qr / sqrt(tau)
    }
  } else {
    qrDim <- 0L
    stop("NA/NaN/Inf in gradient matrix at parameter estimates")
  }
  getPars.varying <- function() unlist(mget(names(ind), env))[useParams]
  setPars.noVarying <- function(newPars){
    internalPars <<- newPars # envir = thisEnv
    for(i in names(ind))
      env[[i]] <- unname(newPars[ ind[[i]] ])
  }
  setPars.varying <- function(newPars) {
    internalPars[useParams] <<- newPars
    for(i in names(ind))
      env[[i]] <- unname(internalPars[ ind[[i]] ])
  }
  setPars <- setPars.noVarying

  convCrit <- function() {
    if(npar == 0) return(0)
    rr <- qr.qty(QR, c(resid)) # rotated residual vector
    sqrt(sum(rr[1L:npar]^2) / sum(rr[-(1L:npar)]^2))
  }
  ## newdata gradient
  if(is.function(jac)) {
    gcall <- do.call(call, args = c(".jac", sapply(c(names(formals(jac)[1]),
                                                     intersect(c(".data", names(data)), names(formals(jac)))), as.name)), quote = TRUE)
    gcall[[2]] <- do.call(call, args = c(ifelse(is.list(start), "list", "c"),
                                         sapply(names(start), as.name)), quote = TRUE)
    env$.jac <- jac
  } else {
    gcall <- form[[3L]]
  }
  on.exit(remove(i, data, parLength, start, temp, m, gr,
                 marg, dimGrad, qrDim, gradSetArgs))
  m <- list(resid = function() resid,
            fitted = function() rhs,
            formula = function() form,
            deviance = function() dev,
            lhs = function() lhs,
            gradient = function() {
              if(is.matrix(.swts)) {
                .swts %*% attr(rhs, "gradient")
              } else {
                .swts * attr(rhs, "gradient")
              }
            },
            gradient1 = function(newdata = list()) {
              rho <- new.env(hash = TRUE, parent = env)
              for(i in names(newdata)) rho[[i]] <- newdata[[i]]
              if(is.function(env$.jac)) {
                rho$.data <- as.list(newdata)
                grad <- eval(gcall, envir = rho)
              } else {
                drv <- stats::numericDeriv(form[[3L]], names(ind), rho = rho)
                grad <- attr(drv, "gradient")
                colnames(grad) <- names(ind)
              }
              grad
            },
            conv = function() convCrit(),
            incr = function() qr.coef(QR, resid),
            setVarying = function(vary = rep_len(TRUE, np)) {
              np <- length(useParams)
              useParams <<- useP <-
                if(is.character(vary)) {
                  temp <- logical(np)
                  temp[unlist(ind[vary])] <- TRUE
                  temp
                } else if(is.logical(vary) && length(vary) != np)
                  stop("setVarying : 'vary' length must match length of parameters")
              else
                vary # envir = thisEnv
              gradCall[[length(gradCall) - 1L]] <<- useP
              if(all(useP)) {
                setPars <<- setPars.noVarying
                getPars <<- getPars.noVarying
                getRHS  <<-  getRHS.noVarying
                npar    <<- length(useP)
              } else {
                setPars <<- setPars.varying
                getPars <<- getPars.varying
                getRHS  <<-  getRHS.varying
                npar    <<- sum(useP)
              }
            },
            setPars = function(newPars) {
              setPars(newPars)
              resid0 <- (lhs - (rhs <<- getRHS()))
              resid <<- if(is.matrix(.swts)) c(.swts %*% resid0) else .swts * resid0 # envir = thisEnv {2 x}
              dev   <<- sum(resid^2) # envir = thisEnv
              if(length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
              QR <<- if(is.matrix(.swts)) qr(.swts %*% gr) else qr(.swts * gr) # envir = thisEnv
              (QR$rank < min(dim(QR$qr))) # to catch the singular gradient matrix
            },
            getPars = function() getPars(),
            getAllPars = function() getPars(),
            getEnv = function() env,
            trace = function() {
              cat(format(dev), ": ", format(getPars()))
              cat("\n")
            },
            Rmat = function() qr.R(QR),
            predict = function(newdata = list(), qr = FALSE)
              eval(form[[3L]], envir = as.list(newdata), enclos = env)
  )
  class(m) <- "nlsModel"
  m
}

gslModel <- function(fn, lhs, cFit, start, swts, jac, ...) {
  env <- new.env(hash = TRUE, parent = environment(fn))
  env$fn <- fn
  data <- list(...)
  for(i in names(data)) env[[i]] <- data[[i]]
  if(is.list(start)) {
    pars <- as.list(cFit$par)
  } else {
    pars <- cFit$par
  }
  parName <- names(formals(fn)[1])
  env[[parName]] <- pars
  fcall <- do.call(call, args = c("fn", sapply(c(parName, names(data)), as.name)), quote = TRUE)
  if(is.function(jac)) {
    gcall <- do.call(call, args = c(".gradient", sapply(c(names(formals(jac)[1]), names(data)), as.name)), quote = TRUE)
    env$.gradient <- jac
  } else {
    gcall <- fcall
  }
  gcall[[2]] <- do.call(call, args = c(ifelse(is.list(start), "list", "c"), sapply(names(start), as.name)), quote = TRUE)
  if(is.null(cFit$irls)) {
    gr <- cFit$grad
  } else {
    .swts <- if(!missing(swts) && length(swts)) swts else rep_len(1, dim(cFit$grad)[1])
    gr <- if(is.matrix(.swts)) .swts %*% cFit$grad else .swts * cFit$grad  ## weighted gradient
    gr <- gr / sqrt(cFit$irls$irls_weights)  ## gradient w/o irls weights
  }
  if(!any(is.na(gr) | is.infinite(gr))) {
    QR <- qr(gr)
    qrDim <- min(dim(QR$qr))
    if(QR$rank < qrDim)
      warning("singular gradient matrix at parameter estimates")
    if(!is.null(cFit$irls)) {
      ## MM correction [Yohai, 1987. Thrm. 4.1]
      tau <- mean(cFit$irls$irls_psi^2) / mean(cFit$irls$irls_dpsi)^2
      QR$qr <- QR$qr / sqrt(tau)
    }
  } else {
    QR <- structure(list(qr = gr), class = "qr")
    QR$qr[] <- NA_real_
    qrDim <- 0L
    warning("NA/NaN/Inf in gradient matrix at parameter estimates")
  }
  on.exit(remove(i, data, m, qrDim))
  m <- list(resid = function() -cFit$resid,
            fitted = function() eval(fcall, envir = env),
            formula = function() fn,
            deviance = function() cFit$ssr,
            lhs = function() lhs,
            gradient = function() gr,
            gradient1 = function(newdata = list()) {
              rho <- new.env(hash = TRUE, parent = env)
              for(i in names(newdata)) rho[[i]] <- newdata[[i]]
              for(i in names(pars)) rho[[i]] <- pars[[i]]
              if(is.function(env$.gradient)) {
                grad <- eval(gcall, envir = rho)
              } else {
                drv <- stats::numericDeriv(gcall, names(pars), rho = rho)
                grad <- attr(drv, "gradient")
                colnames(grad) <- names(pars)
              }
              grad
            },
            getPars = function() unlist(pars),
            getAllPars = function() unlist(pars),
            getEnv = function() env,
            Rmat = function() qr.R(QR),
            predict = function(newdata = list())
              eval(fcall, envir = as.list(newdata), enclos = env)
  )
  m
}
