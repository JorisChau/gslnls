#' GSL Large-scale Nonlinear Least Squares fitting
#'
#' Determine the nonlinear least-squares estimates of the parameters of a large
#' nonlinear model system using the \code{gsl_multilarge_nlinear} module present in
#' the GNU Scientific Library (GSL).
#'
#' @inheritParams gsl_nls
#' @param algorithm character string specifying the algorithm to use. The following choices are supported:
#' \itemize{
#' \item \code{"cgst"} Steihaug-Toint Conjugate Gradient algorithm (default), a generalization of the dogleg algorithm
#' that avoids solving for the Gauss-Newton step directly, instead using an iterative conjugate gradient algorithm.
#' The method performs well at points where the Jacobian is singular, and is also suitable for large-scale problems
#' where factoring the Jacobian matrix is prohibitively expensive.
#' \item \code{"lm"} Levenberg-Marquadt algorithm
#' \item \code{"lmaccel"} Levenberg-Marquadt algorithm with geodesic acceleration.
#' Can be faster than \code{"lm"} but less stable. Stability is controlled by the
#' \code{avmax} parameter in \code{control}, setting \code{avmax} to zero is analogous
#' to not using geodesic acceleration.
#' \item \code{"dogleg"} Powell's dogleg algorithm
#' \item \code{"ddogleg"} Double dogleg algorithm, an improvement over \code{"dogleg"}
#' by including information about the Gauss-Newton step while the iteration is still
#' far from the minimum.
#' \item \code{"subspace2D"} 2D generalization of the dogleg algorithm. This method
#' searches a larger subspace for a solution, it can converge more quickly than \code{"dogleg"}
#' on some problems.
#' }
#' @param jac a \link{function} returning the \code{n} by \code{p} dimensional Jacobian matrix of
#' the nonlinear model \code{fn}, where \code{n} is the number of observations and \code{p} the
#' number of parameters. The first argument must be the vector of parameters of length \code{p}.
#' Can also be \code{TRUE}, in which case \code{jac} is derived symbolically with \code{\link[stats]{deriv}},
#' this only works if \code{fn} is defined as a (non-selfstarting) formula. If \code{fn} is a \code{\link{selfStart}} model,
#' the Jacobian specified in the \code{"gradient"} attribute of the self-start model is used instead.
#' @param fvv a \link{function} returning an \code{n} dimensional vector containing
#' the second directional derivatives of the nonlinear model \code{fn}, with \code{n} the number of observations.
#' This argument is only used if geodesic acceleration is enabled (\code{algorithm = "lmaccel"}).
#' The first argument must be the vector of parameters of length \code{p} and the second argument must be the velocity vector
#' also of length \code{p}. Can also be \code{TRUE}, in which case \code{fvv} is derived
#' symbolically with \code{\link[stats]{deriv}}, this only works if \code{fn} is defined as a (non-selfstarting) formula.
#' If the model \link{function} in \code{fn} also returns a \code{"hessian"} attribute (similar to the \code{"gradient"} attribute
#' in a \code{selfStart} model), this Hessian matrix is used to evaluate the second directional derivatives instead.
#' @param trace logical value indicating if a trace of the iteration progress should be printed.
#' Default is \code{FALSE}. If \code{TRUE}, the residual (weighted) sum-of-squares,
#' the squared (Euclidean) norm of the current parameter estimates and the condition number of the Jacobian
#' are printed after each iteration.
#' @return
#' If \code{fn} is a \code{formula} returns a list object of class \code{nls}.
#' If \code{fn} is a \code{function} returns a list object of class \code{gsl_nls}.
#' See the individual method descriptions for the structures of the returned lists and the generic functions
#' applicable to objects of both classes.
#' @importFrom Matrix nnzero
#' @seealso \code{\link{gsl_nls}}
#' @seealso \url{https://www.gnu.org/software/gsl/doc/html/nls.html}
#' @references M. Galassi et al., \emph{GNU Scientific Library Reference Manual (3rd Ed.)}, ISBN 0954612078.
#' @examples
#' # Large NLS example
#' # (https://www.gnu.org/software/gsl/doc/html/nls.html#large-nonlinear-least-squares-example)
#'
#' ## number of parameters
#' p <- 250
#'
#' ## model function
#' f <- function(theta) {
#'   c(sqrt(1e-5) * (theta - 1), sum(theta^2) - 0.25)
#' }
#'
#' ## jacobian function
#' jac <- function(theta) {
#'   rbind(diag(sqrt(1e-5), nrow = length(theta)), 2 * t(theta))
#' }
#'
#' ## dense Levenberg-Marquadt
#' gsl_nls_large(
#'   fn = f,                       ## model
#'   y = rep(0, p + 1),            ## (dummy) responses
#'   start = 1:p,                  ## start values
#'   algorithm = "lm",             ## algorithm
#'   jac = jac,                    ## jacobian
#'   control = list(maxiter = 250)
#' )
#'
#' ## dense Steihaug-Toint conjugate gradient
#' gsl_nls_large(
#'   fn = f,                       ## model
#'   y = rep(0, p + 1),            ## (dummy) responses
#'   start = 1:p,                  ## start values
#'   jac = jac,                    ## jacobian
#'   algorithm = "cgst"            ## algorithm
#' )
#'
#' ## sparse Jacobian function
#' jacsp <- function(theta) {
#'   rbind(Matrix::Diagonal(x = sqrt(1e-5), n = length(theta)), 2 * t(theta))
#' }
#'
#' ## sparse Levenberg-Marquadt
#' gsl_nls_large(
#'   fn = f,                       ## model
#'   y = rep(0, p + 1),            ## (dummy) responses
#'   start = 1:p,                  ## start values
#'   algorithm = "lm",             ## algorithm
#'   jac = jacsp,                  ## sparse jacobian
#'   control = list(maxiter = 250)
#' )
#'
#' ## sparse Steihaug-Toint conjugate gradient
#' gsl_nls_large(
#'   fn = f,                       ## model
#'   y = rep(0, p + 1),            ## (dummy) responses
#'   start = 1:p,                  ## start values
#'   jac = jacsp,                  ## sparse jacobian
#'   algorithm = "cgst"            ## algorithm
#' )
#'
#' @export
gsl_nls_large <- function (fn, ...) {
  UseMethod("gsl_nls_large")
}

#' @describeIn gsl_nls_large
#' If \code{fn} is a \code{formula}, the returned list object is of classes \code{gsl_nls} and \code{nls}.
#' Therefore, all generic functions applicable to objects of class \code{nls}, such as \code{anova}, \code{coef}, \code{confint},
#' \code{deviance}, \code{df.residual}, \code{fitted}, \code{formula}, \code{logLik}, \code{nobs}, \code{predict}, \code{print}, \code{profile},
#' \code{residuals}, \code{summary}, \code{vcov} and \code{weights} are also applicable to the returned list object.
#' In addition, a method \code{confintd} is available for inference of derived parameters.
#' @export
gsl_nls_large.formula <- function(fn, data = parent.frame(), start,
                                  algorithm = c("cgst", "lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"),
                                  control = gsl_nls_control(), jac, fvv, trace = FALSE,
                                  subset, weights, na.action, model = FALSE, ...) {

  ## adapted from src/library/stats/nls.R
  formula <- as.formula(fn)
  algorithm <- match.arg(algorithm, c("cgst", "lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"))

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
    } else
      names(start)

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
        warning("No starting values specified for some parameters.\n",
                "Initializing ", paste(sQuote(nnn), collapse=", "),
                " to '1.'.\n",
                "Consider specifying 'start' or using a selfStart model", domain = NA)
        start <- setNames(rep_len(1., length(nnn)), nnn)
        varNames <- varNames[i <- is.na(match(varNames, nnn))]
        n <- n[i]
      }
      else                        # has 'start' but forgot some
        stop(gettextf("parameters without starting value in 'data': %s",
                      paste(nnn, collapse=", ")), domain = NA)
    }
  }
  else { ## length(varNames) == 0
    stop("no parameters to fit and/or no data variables present")
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
      for (i in names(start))
        startEnv[[i]] <- start[[i]]
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
        mf$fn <- mf$jac <- mf$fvv <- NULL
      ## need stats:: for non-standard evaluation
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights) model.weights(mf) else NULL
    }
    if (!is.null(wts) && any(wts < 0 | is.na(wts)))
      stop("missing or negative weights not allowed")
  }
  else {
    stop("no data variables present")
  }

  ## set up iteration
  if(missing(start))
    start <- getInitial(formula, data=mf, control=control, trace=trace)
  for(var in varNames[!varIndex])
    mf[[var]] <- eval(as.name(var), data, env)
  varNamesRHS <- varNamesRHS[ varNamesRHS %in% varNames[varIndex] ]

  ## function call
  .fn <- function(par, .data = mf) eval(formula[[3L]], envir = c(as.list(par), .data))
  .fcall <- tryCatch(.fn(start), error = function(err) err)

  if(inherits(.fcall, "error"))
    stop(sprintf("failed to evaluate 'fn' at starting values: %s", .fcall$message))
  if(!is.numeric(.fcall))
    stop("'fn' failed to return a numeric vector at starting values")
  if(any(is.na(.fcall)))
    stop("missing values returned by 'fn' at starting values")

  .lhs <- eval(formula[[2L]], envir = mf)

  ## n > p
  if(length(.lhs) < length(start)) {
    stop("zero or less residual degrees of freedom, cannot fit a model with less observations than parameters")
  }

  ## jac call
  if(missing(jac)) {
    jac <- NULL
  }
  if(!is.function(jac) && !is.null(attr(.fcall, "gradient"))) {
    jac <- function(par, .data = mf) attr(.fn(par, .data), "gradient")
  } else if(isTRUE(jac)) {
    .exprjac <- tryCatch(stats::deriv(formula[[3L]], namevec = names(start)), error = function(err) err)
    if(inherits(.exprjac, "error")) {
      stop(sprintf("failed to symbolically derive 'jac': %s", .exprjac$message))
    } else if(is.expression(.exprjac)){
      jac <- function(par, .data = mf) {
        grad <- eval(.exprjac, envir = c(as.list(par), .data))
        attr(grad, "gradient")
      }
    }
  }

  if(is.function(jac)) {
    .jac <- function(par) jac(par, ...)
    .dfcall <- tryCatch(.jac(start), error = function(err) err)
    if(inherits(.dfcall, "error"))
      stop(sprintf("failed to evaluate 'jac' at starting values: %s", .dfcall$message))
    if(!identical(dim(.dfcall), c(length(.lhs), length(start))))
      stop("'jac' failed to return a matrix of expected dimensions at starting values")
    if(!((is.numeric(.dfcall) && is.matrix(.dfcall)) || (inherits(.dfcall, "dMatrix") && inherits(.dfcall, "generalMatrix"))))
      stop("only \"matrix\", \"dgCMatrix\", \"dgRMatrix\", \"dgTMatrix\", and \"dgeMatrix\" classes are supported for 'jac'")
    if(any(is.na(.dfcall)))
      stop("missing values returned by 'jac' at starting values")
  } else {
    stop("analytic Jacobian function 'jac' is required, but none is available")
  }

  ## fvv call
  if(identical(algorithm, "lmaccel")) {
    if(missing(fvv)) {
      fvv <- NULL
    }
    if(!is.function(fvv) && !is.null(attr(.fcall, "hessian"))) {
      fvv <- function(par, v) {
        hess <- attr(.fn(par), "hessian")
        c(matrix(hess, nrow = nrow(hess), ncol = ncol(hess) * ncol(hess)) %*% c(outer(v, v)))
      }
    } else if(isTRUE(fvv)) {
      .exprfvv <- tryCatch(stats::deriv(formula[[3L]], namevec = names(start), hessian = TRUE), error = function(err) err)
      if(inherits(.exprfvv, "error")) {
        stop(sprintf("failed to symbolically derive 'fvv': %s", .exprfvv$message))
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
      .fvvcall <- tryCatch(.fvv(start, structure(rep(1, length(start)), names = names(start))), error = function(err) err)
      if(inherits(.fvvcall, "error"))
        stop(sprintf("failed to evaluate 'fvv' at starting values: %s", .fvvcall$message))
      if(!is.numeric(.fvvcall) || !identical(length(.fvvcall), length(.lhs)))
        stop("'fvv' failed to return a numeric vector equal in length to 'y' at starting values")
      if(any(is.na(.fvvcall)))
        stop("missing values returned by 'fvv' at starting values")
    } else {
      stop("analytic second derivative function 'fvv' is required, but none is available")
    }
  } else {
    .fvv <- NULL
  }

  ## control arguments
  trace <- isTRUE(trace)
  .ctrl <- gsl_nls_control()
  if(!missing(control)) {
    control <- as.list(control)
    .ctrl[names(control)] <- control
  }
  .ctrl$scale <- match.arg(.ctrl$scale, c("more", "levenberg", "marquadt"))
  .ctrl$solver <- "cholesky"  ## fixed
  .ctrl$fdtype <- match.arg(.ctrl$fdtype, c("forward", "center"))
  stopifnot(
    is.numeric(.ctrl$maxiter), length(.ctrl$maxiter) == 1, .ctrl$maxiter > 0,
    is.numeric(.ctrl$factor_up), length(.ctrl$factor_up) == 1, .ctrl$factor_up > 0,
    is.numeric(.ctrl$factor_down), length(.ctrl$factor_down) == 1, .ctrl$factor_down > 0,
    is.numeric(.ctrl$avmax), length(.ctrl$avmax) == 1, .ctrl$avmax > 0,
    is.numeric(.ctrl$h_df), length(.ctrl$h_df) == 1, .ctrl$h_df > 0,
    is.numeric(.ctrl$h_fvv), length(.ctrl$h_fvv) == 1, .ctrl$h_fvv > 0,
    is.numeric(.ctrl$xtol), length(.ctrl$xtol) == 1, .ctrl$xtol > 0,
    is.numeric(.ctrl$ftol), length(.ctrl$ftol) == 1, .ctrl$ftol > 0,
    is.numeric(.ctrl$gtol), length(.ctrl$gtol) == 1, .ctrl$gtol > 0
  )
  .ctrl_int <- c(
    maxiter = as.integer(.ctrl$maxiter),
    trace = isTRUE(trace),
    algorithm = match(algorithm, c("cgst", "lm", "lmaccel", "dogleg", "ddogleg", "subspace2D")) - 1L,
    scale = match(.ctrl$scale, c("more", "levenberg", "marquadt")) - 1L,
    fdtype = match(.ctrl$fdtype, c("forward", "center")) - 1L,
    jacclass = -2L,
    jacnz = 0L
  )

  ## sparse jacobian
  if(inherits(.dfcall, "sparseMatrix")) {
    .ctrl_int[7] <- as.integer(nnzero(.dfcall, na.counted = FALSE))
  }
  if(inherits(.dfcall, "dgTMatrix")) {
    .ctrl_int[6] <- 0L
  } else if(inherits(.dfcall, "dgCMatrix")) {
    .ctrl_int[6] <- 1L
  } else if(inherits(.dfcall, "dgRMatrix")) {
    .ctrl_int[6] <- 2L
  } else if(inherits(.dfcall, "dgeMatrix")) {
    .ctrl_int[6] <- -1L
  }

  .ctrl_dbl <- unlist(.ctrl[c("factor_up", "factor_down", "avmax", "h_df", "h_fvv", "xtol", "ftol", "gtol")])
  .ctrl_tol <- unlist(.ctrl[c("xtol", "ftol", "gtol")])

  ## optimize
  cFit <- .Call(C_nls_large, .fn, .lhs, .jac, .fvv, environment(), start, wts, .ctrl_int, .ctrl_dbl, PACKAGE = "gslnls")

  ## convert to nls object
  m <- nlsModel(formula, mf, cFit$par, wts, jac)

  convInfo <- list(
    isConv = as.logical(!cFit$conv),
    finIter = cFit$niter,
    finTol = cFit$ssrtol,
    nEval = cFit$neval,
    trsName = paste("multilarge", cFit$algorithm, sep = "/"),
    stopCode = cFit$conv,
    stopMessage = cFit$status
  )

  nls.out <- list(m = m, data = substitute(data), convInfo = convInfo, call = match.call())

  nls.out$call$algorithm <- algorithm
  nls.out$call$control <- nls.control() ## needed for profiler
  nls.out$call$trace <- trace
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
  class(nls.out) <- c("gsl_nls", "nls")

  return(nls.out)

}

#' @describeIn gsl_nls_large
#' If \code{fn} is a \code{function}, the first argument must be the vector of parameters and
#' the function should return a numeric vector containing the nonlinear model evaluations at
#' the provided parameter and predictor or covariate vectors. In addition, the argument \code{y}
#' needs to contain the numeric vector of observed responses, equal in length to the numeric
#' vector returned by \code{fn}. The returned list object is (only) of class \code{gsl_nls}.
#' Although the returned object is not of class \code{nls}, the following generic functions remain
#' applicable for an object of class \code{gsl_nls}: \code{anova}, \code{coef}, \code{confint}, \code{deviance},
#' \code{df.residual}, \code{fitted}, \code{formula}, \code{logLik}, \code{nobs}, \code{predict}, \code{print},
#' \code{residuals}, \code{summary}, \code{vcov} and \code{weights}. In addition, a method \code{confintd}
#' is available for inference of derived parameters.
#' @export
gsl_nls_large.function <- function(fn, y, start,
                                   algorithm = c("cgst", "lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"),
                                   control = gsl_nls_control(), jac, fvv, trace = FALSE,
                                   weights, ...) {

  ## algorithm
  algorithm <- match.arg(algorithm, c("cgst", "lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"))

  ## starting values
  if(missing(start)) {
    stop("starting values need to be provided if 'fn' is defined as a function")
  }

  ## function call
  if(!is.numeric(y))
    stop("'y' should be a numeric response vector")

  ## n > p
  if(length(y) < length(start)) {
    stop("zero or less residual degrees of freedom, cannot fit a model with less observations than parameters")
  }

  .fn <- function(par) fn(par, ...)
  .fcall <- tryCatch(.fn(start), error = function(err) err)
  if(inherits(.fcall, "error"))
    stop(sprintf("failed to evaluate 'fn' at starting values: %s", .fcall$message))
  if(!is.numeric(.fcall) || !identical(length(.fcall), length(y)))
    stop("'fn' failed to return a numeric vector equal in length to 'y' at starting values")
  if(any(is.na(.fcall)))
    stop("missing values returned by 'fn' at starting values")

  ## jac call
  if(missing(jac)) {
    jac <- NULL
  }
  if(!is.function(jac) && !is.null(attr(.fcall, "gradient"))) {
    jac <- function(par, ...) attr(fn(par, ...), "gradient")
  }
  if(is.function(jac)) {
    .jac <- function(par) jac(par, ...)
    .dfcall <- tryCatch(.jac(start), error = function(err) err)
    if(inherits(.dfcall, "error"))
      stop(sprintf("failed to evaluate 'jac' at starting values: %s", .dfcall$message))
    if(!identical(dim(.dfcall), c(length(y), length(start))))
      stop("'jac' failed to return a matrix of expected dimensions at starting values")
    if(!((is.numeric(.dfcall) && is.matrix(.dfcall)) || (inherits(.dfcall, "dMatrix") && inherits(.dfcall, "generalMatrix"))))
      stop("only \"matrix\", \"dgCMatrix\", \"dgRMatrix\", \"dgTMatrix\", and \"dgeMatrix\" classes are supported for 'jac'")
    if(any(is.na(.dfcall)))
      stop("missing values returned by 'jac' at starting values")
  } else {
    stop("analytic Jacobian function 'jac' is required, but none is available")
  }

  ## fvv call
  if(identical(algorithm, "lmaccel")) {
    if(missing(fvv)) {
      fvv <- NULL
    }
    if(!is.function(fvv) && !is.null(attr(.fcall, "hessian"))) {
      fvv <- function(par, v, ...) {
        hess <- attr(fn(par, ...), "hessian")
        c(matrix(hess, nrow = nrow(hess), ncol = ncol(hess) * ncol(hess)) %*% c(outer(v, v)))
      }
    }
    if(is.function(fvv)) {
      .fvv <- function(par, v) fvv(par, v, ...)
      .fvvcall <- tryCatch(.fvv(start, structure(rep(1, length(start)), names = names(start))), error = function(err) err)
      if(inherits(.fvvcall, "error"))
        stop(sprintf("failed to evaluate 'fvv' at starting values: %s", .fvvcall$message))
      if(!is.numeric(.fvvcall) || !identical(length(.fvvcall), length(y)))
        stop("'fvv' failed to return a numeric vector equal in length to 'y' at starting values")
      if(any(is.na(.fvvcall)))
        stop("missing values returned by 'fvv' at starting values")
    } else {
      stop("analytic second derivative function 'fvv' is required, but none is available")
    }
  } else {
    .fvv <- NULL
  }

  ## control arguments
  trace <- isTRUE(trace)
  .ctrl <- gsl_nls_control()
  if(!missing(control)) {
    control <- as.list(control)
    .ctrl[names(control)] <- control
  }
  .ctrl$scale <- match.arg(.ctrl$scale, c("more", "levenberg", "marquadt"))
  .ctrl$fdtype <- match.arg(.ctrl$fdtype, c("forward", "center"))
  stopifnot(
    is.numeric(.ctrl$maxiter), length(.ctrl$maxiter) == 1, .ctrl$maxiter > 0,
    is.numeric(.ctrl$factor_up), length(.ctrl$factor_up) == 1, .ctrl$factor_up > 0,
    is.numeric(.ctrl$factor_down), length(.ctrl$factor_down) == 1, .ctrl$factor_down > 0,
    is.numeric(.ctrl$avmax), length(.ctrl$avmax) == 1, .ctrl$avmax > 0,
    is.numeric(.ctrl$h_df), length(.ctrl$h_df) == 1, .ctrl$h_df > 0,
    is.numeric(.ctrl$h_fvv), length(.ctrl$h_fvv) == 1, .ctrl$h_fvv > 0,
    is.numeric(.ctrl$xtol), length(.ctrl$xtol) == 1, .ctrl$xtol > 0,
    is.numeric(.ctrl$ftol), length(.ctrl$ftol) == 1, .ctrl$ftol > 0,
    is.numeric(.ctrl$gtol), length(.ctrl$gtol) == 1, .ctrl$gtol > 0
  )
  .ctrl_int <- c(
    maxiter = as.integer(.ctrl$maxiter),
    trace = isTRUE(trace),
    algorithm = match(algorithm, c("cgst", "lm", "lmaccel", "dogleg", "ddogleg", "subspace2D")) - 1L,
    scale = match(.ctrl$scale, c("more", "levenberg", "marquadt")) - 1L,
    fdtype = match(.ctrl$fdtype, c("forward", "center")) - 1L,
    jacclass = -2L,
    jacnz = 0L
  )

  ## sparse jacobian
  if(inherits(.dfcall, "sparseMatrix")) {
    .ctrl_int[7] <- as.integer(nnzero(.dfcall, na.counted = FALSE))
  }
  if(inherits(.dfcall, "dgTMatrix")) {
    .ctrl_int[6] <- 0L
  } else if(inherits(.dfcall, "dgCMatrix")) {
    .ctrl_int[6] <- 1L
  } else if(inherits(.dfcall, "dgRMatrix")) {
    .ctrl_int[6] <- 2L
  } else if(inherits(.dfcall, "dgeMatrix")) {
    .ctrl_int[6] <- -1L
  }

  .ctrl_dbl <- unlist(.ctrl[c("factor_up", "factor_down", "avmax", "h_df", "h_fvv", "xtol", "ftol", "gtol")])
  .ctrl_tol <- unlist(.ctrl[c("xtol", "ftol", "gtol")])

  ## weights
  if(!missing(weights)) {
    if(!is.numeric(weights) || !identical(length(weights), length(y)))
      stop("'weights' should be numeric equal in length to 'y'")
    if (any(weights < 0 | is.na(weights)))
      stop("missing or negative weights not allowed")
  } else {
    weights <- NULL
  }

  ## optimize
  cFit <- .Call(C_nls_large, .fn, y, .jac, .fvv, environment(), start, weights, .ctrl_int, .ctrl_dbl, PACKAGE = "gslnls")

  m <- gslModel(fn, y, cFit, start, weights, jac, ...)

  ## mimick nls object
  convInfo <- list(
    isConv = as.logical(!cFit$conv),
    finIter = cFit$niter,
    finTol = cFit$ssrtol,
    nEval = cFit$neval,
    trsName = paste("multilarge", cFit$algorithm, sep = "/"),
    stopCode = cFit$conv,
    stopMessage = cFit$status
  )

  nls.out <- list(m = m, data = c(list(y = y), list(...)), convInfo = convInfo, call = match.call())
  nls.out$call$algorithm <- algorithm
  nls.out$call$trace <- trace
  if(trace) {
    nls.out$partrace <- cFit$partrace[seq_len(cFit$niter + 1L), , drop = FALSE]
    nls.out$devtrace <- cFit$ssrtrace[seq_len(cFit$niter + 1L)]
  }
  nls.out$control <- .ctrl
  if(!is.null(weights))
    nls.out$weights <- weights
  class(nls.out) <- "gsl_nls"

  return(nls.out)

}
