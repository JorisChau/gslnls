#' Extract model coefficients
#' @description Returns the fitted model coefficients from a \code{"gsl_nls"} object.
#' \code{coefficients} can also be used as an alias.
#' @param object An object inheriting from class \code{"gsl_nls"}.
#' @param ... At present no optional arguments are used.
#' @return Named numeric vector of fitted coefficients similar to \code{\link[stats]{coef}}
#' @seealso \code{\link[stats]{coef}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' coef(obj)
#' @export
coef.gsl_nls <- function(object, ...) {
  object$m$getAllPars()
}

#' Extract model fitted values
#' @description Returns the fitted responses from a \code{"gsl_nls"} object. \code{fitted.values}
#' can also be used as an alias.
#' @inheritParams coef.gsl_nls
#' @seealso \code{\link[stats]{fitted}}
#' @return Numeric vector of fitted responses similar to \code{\link[stats]{fitted}}.
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' fitted(obj)
#' @export
fitted.gsl_nls <- function(object, ...) {
  if(inherits(object, "nls")) {
    NextMethod()
  } else {
    val <- object$m$fitted()
    attr(val, "label") <- "Fitted values"
    val
  }
}

#' Extract the number of observations
#' @description Returns the number of \emph{observations} from a \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @return Integer number of observations similar to \code{\link[stats]{nobs}}
#' @seealso \code{\link[stats]{nobs}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' nobs(obj)
#' @export
nobs.gsl_nls <- function(object, ...) {
  if (is.null(w <- object$weights))
    length(object$m$resid())
  else
    sum(w != 0)
}

#' Model deviance
#' @description Returns the deviance of a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @return Numeric deviance value similar to \code{\link[stats]{deviance}}
#' @seealso \code{\link[stats]{deviance}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' deviance(obj)
#' @export
deviance.gsl_nls <- function(object, ...) {
  object$m$deviance()
}

#' Residual standard deviation
#' @description Returns the estimated (unweighted) residual standard deviation of a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @return Numeric residual standard deviation value similar to \code{\link[stats]{sigma}}
#' @seealso \code{\link[stats]{sigma}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' sigma(obj)
#' @export
sigma.gsl_nls <- function(object, ...) {
  if(is.null(object$irls)) {
    NextMethod()
  } else {
    object$irls$irls_sigma
  }
}

#' Extract model formula
#' @description Returns the model formula from a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @param x An object inheriting from class \code{"gsl_nls"}.
#' @return If the object inherits from class \code{"nls"} returns the fitted model as a \link{formula} similar
#' to \code{\link[stats]{formula}}. Otherwise returns the fitted model as a \link{function}.
#' @seealso \code{\link[stats]{formula}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' formula(obj)
#' @export
formula.gsl_nls <- function(x, ...) {
  x$m$formula()
}

#' Model summary
#' @description Returns the model summary for a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @param correlation logical; if \code{TRUE}, the correlation matrix of the estimated
#' parameters is returned and printed.
#' @param symbolic.cor logical; if \code{TRUE}, print the correlations in a symbolic form
#' (see \code{\link[stats]{symnum}}) rather than as numbers.
#' @return List object of class \code{"summary.gsl_nls"} similar to \code{\link[stats]{summary.nls}}
#' @seealso \code{\link[stats]{summary.nls}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' summary(obj)
#' @export
summary.gsl_nls <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  r <- as.vector(object$m$resid())  ## these are weighted residuals
  w <- object$weights
  n <- if (!is.null(w)) sum(w > 0) else length(r)
  param <- coef(object)
  pnames <- names(param)
  p <- length(param)
  rdf <- n - p
  if(is.null(object$irls)) {
    sigma <- if (rdf <= 0) NaN else sqrt(deviance(object)/rdf)
  } else {
    sigma <- object$irls$irls_sigma
  }
  XtXinv <- chol2inv(object$m$Rmat())
  dimnames(XtXinv) <- list(pnames, pnames)
  se <- sigma * sqrt(diag(XtXinv))
  tval <- param/se
  param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error",
                                    "t value", "Pr(>|t|)"))
  if(inherits(object, "nls")) {
    frm <- formula(object)
  } else {
    frm <- as.formula(sprintf("y ~ fn(%s)", paste(names(formals(formula(object))), collapse = ", ")))
  }
  ans <- list(formula = frm, residuals = r, sigma = sigma,
              df = c(p, rdf), cov.unscaled = XtXinv, call = object$call,
              convInfo = object$convInfo, control = object$control,
              na.action = object$na.action, coefficients = param, parameters = param,
              irls = object$irls)
  if (correlation && rdf > 0) {
    ans$correlation <- (XtXinv * sigma^2) / outer(se, se)
    ans$symbolic.cor <- symbolic.cor
  }
  class(ans) <- "summary.gsl_nls"
  ans
}

#' Print model object
#' @description Print method for a fitted \code{"gsl_nls"} object
#' @param x An object inheriting from class \code{"gsl_nls"}.
#' @param digits Minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... At present no optional arguments are used.
#' @return Returns the object \code{x} \emph{invisibly} (via \code{\link{invisible}}).
#' @noRd
#' @export
print.gsl_nls <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Nonlinear regression model\n")
  if(inherits(formula(x), "formula")) {
    cat("  model: ", deparse(formula(x)), "\n", sep = "")
    if(is.symbol(x$data)) {
      cat("   data: ", deparse(x$data), "\n", sep = "")
    }
  } else {
    cat("  model: ", sprintf("y ~ fn(%s)", paste(names(formals(formula(x))), collapse = ", ")), "\n", sep = "")
  }
  print(x$m$getAllPars(), digits = digits, ...)
  cat(" ", if(!is.null(x$weights) || !is.null(x$irls)) "weighted ",
      "residual sum-of-squares: ", format(x$m$deviance(), digits = digits),
      "\n", sep = "")
  convInfo(x, digits = digits)
  invisible(x)
}

#' Print model summary
#' @description Print method for a summary \code{"summary.gsl_nls"} object
#' @param x An object inheriting from class \code{"summary.gsl_nls"}
#' @param digits Minimal number of significant digits, see \code{\link{print.default}}
#' @inheritParams summary.nls
#' @return Returns the object \code{x} \emph(invisibly) (via \code{\link{invisible}}).
#' @noRd
#' @export
print.summary.gsl_nls <- function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
                                   signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nFormula: ", paste(deparse(x$formula), sep = "\n",
                           collapse = "\n"), "\n", sep = "")
  df <- x$df
  rdf <- df[2L]
  cat("\nParameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               ...)
  cat("\nResidual standard error:", format(signif(x$sigma,
                                                  digits)), "on", rdf, "degrees of freedom")
  cat("\n")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Parameter Estimates:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2L,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  convInfo(x, digits = digits)
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Calculate model predicted values
#' @description Returns predicted values for the expected response from a fitted \code{"gsl_nls"} object.
#' Asymptotic confidence or prediction (tolerance) intervals at a given \code{level} can be evaluated
#' by specifying the appropriate \code{interval} argument.
#' @inheritParams coef.gsl_nls
#' @param newdata A named list or data.frame in which to look for variables with which to predict. If
#' \code{newdata} is missing, the predicted values at the original data points are returned.
#' @param scale A numeric scalar or vector. If it is set, it is used as the residual standard deviation
#' (or vector of residual standard deviations) in the computation of the standard errors, otherwise
#' this information is extracted from the model fit.
#' @param interval A character string indicating if confidence or prediction (tolerance) intervals
#' at the specified level should be returned.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for the intervals (if any)
#' to be calculated.
#' @return If \code{interval = "none"} (default), a vector of predictions for the mean response. Otherwise,
#' a matrix with columns \code{fit}, \code{lwr} and \code{upr}. The first column (\code{fit}) contains
#' predictions for the mean response. The other two columns contain lower (\code{lwr}) and upper (\code{upr})
#' confidence or prediction bounds at the specified \code{level}.
#' @seealso \code{\link[stats]{predict.nls}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' predict(obj)
#' predict(obj, newdata = data.frame(x = 1:(2 * n) / n))
#' predict(obj, interval = "confidence")
#' predict(obj, interval = "prediction", level = 0.99)
#' @export
predict.gsl_nls <- function(object, newdata, scale = NULL, interval = c("none", "confidence", "prediction"), level = 0.95, ...) {
  interval <- match.arg(interval, c("none", "confidence", "prediction"))

  if (missing(newdata)) {
    fit <- as.vector(fitted(object))
    if(interval != "none") {
      Fdot <- object$m$gradient()  ## weighted gradient
      if(!is.null(object$weights))
        Fdot <- Fdot / sqrt(object$weights)
    }
  } else {
    fit <- object$m$predict(newdata)
    if(interval != "none") {
      Fdot <- object$m$gradient1(newdata)  ## unweighted gradient
    }
  }
  if(interval != "none") {
    if(is.null(scale))
      scale <- sigma(object)
    a <- c((1 - level) / 2, (1 + level) / 2)
    ses <- scale * sqrt(1 * (interval == "prediction") + rowSums(Fdot %*% chol2inv(object$m$Rmat()) * Fdot))
    ci <- fit + ses %o% qt(a, df.residual(object))
    cimat <- cbind(fit = fit, lwr = ci[, 1], upr = ci[, 2])
    return(cimat)
  } else {
    return(fit)
  }
}

#' Extract model residuals
#' @description Returns the model residuals from a fitted \code{"gsl_nls"} object.
#' \code{resid} can also be used as an alias.
#' @inheritParams coef.gsl_nls
#' @param type character; if \code{"response"} the raw residuals are returned, if \code{"pearson"}
#' the Pearson are returned, i.e. the raw residuals divided by their standard error.
#' @return Numeric vector of model residuals similar to \code{\link[stats]{residuals}}.
#' @seealso \code{\link[stats]{residuals}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' residuals(obj)
#' @export
residuals.gsl_nls <- function(object, type = c("response", "pearson"), ...) {
  if(inherits(object, "nls")) {
    NextMethod()
  } else {
    type <- match.arg(type)
    if (type == "pearson") {
      val <- as.vector(object$m$resid())
      std <- sqrt(sum(val^2)/(length(val) - length(coef(object))))
      val <- val/std
      attr(val, "label") <- "Standardized residuals"
    } else {
      val <- as.vector(object$m$lhs() - object$m$fitted())
      attr(val, "label") <- "Residuals"
    }
    val
  }
}

#' Extract model log-likelihood
#' @description Returns the model log-likelihood of a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @param REML logical value; included for compatibility reasons only, should not be used.
#' @return Numeric object of class \code{"logLik"} identical to \code{\link[stats]{logLik}}.
#' @seealso \code{\link[stats]{logLik}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' logLik(obj)
#' @export
logLik.gsl_nls <- function(object, REML = FALSE, ...) {
  if (REML)
    stop("cannot calculate REML log-likelihood for \"nls\" objects")
  res <- object$m$resid() # These are weighted residuals.
  N <- length(res)
  w <- if(!is.null(object$weights)) object$weights else rep_len(1, N)
  ## Note the trick for zero weights
  zw <- w == 0
  N <- sum(!zw)
  val <-  -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(res^2)))/2
  ## the formula here corresponds to estimating sigma^2.
  attr(val, "df") <- 1L + length(coef(object))
  attr(val, "nobs") <- attr(val, "nall") <- N
  class(val) <- "logLik"
  val
}

#' Residual degrees-of-freedom
#' @description Returns the residual degrees-of-freedom from a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @return Integer residual degrees-of-freedom similar to \code{\link[stats]{df.residual}}.
#' @seealso \code{\link[stats]{df.residual}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' df.residual(obj)
#' @export
df.residual.gsl_nls <- function(object, ...) {
  w <- object$weights
  n <- if(!is.null(w)) sum(w != 0) else length(object$m$resid())
  n - length(coef(object))
}

#' Calculate variance-covariance matrix
#' @description Returns the estimated variance-covariance matrix of the model parameters
#' from a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @return A matrix containing the estimated covariances between the parameter estimates similar
#' to \code{\link[stats]{vcov}} with row and column names corresponding to the parameter names given by \code{\link{coef.gsl_nls}}.
#' @seealso \code{\link[stats]{vcov}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' vcov(obj)
#' @export
vcov.gsl_nls <- function(object, ...) {
  pnames <- names(coef(object))
  XtXinv <- chol2inv(object$m$Rmat())
  dimnames(XtXinv) <- list(pnames, pnames)
  sigma(object)^2 * XtXinv
}

#' Calculate leverage values
#' @description Returns leverage values (hat values) from a fitted \code{"gsl_nls"} object based on the estimated
#' variance-covariance matrix of the model parameters.
#' @inheritParams coef.gsl_nls
#' @param model An object inheriting from class \code{"gsl_nls"}.
#' @return Numeric vector of leverage values similar to \code{\link[stats]{hatvalues}}.
#' @seealso \code{\link[stats]{hatvalues}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'  x = (1:n) / n,
#'  y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' hatvalues(obj)
#' @export
hatvalues.gsl_nls <- function(model, ...) {
  Fdot <- model$m$gradient()
  rowSums(Fdot %*% chol2inv(model$m$Rmat()) * Fdot)
}

#' Calculate Cook's distance
#' @description Returns Cook's distance values from a fitted \code{"gsl_nls"} object based on the estimated
#' variance-covariance matrix of the model parameters.
#' @inheritParams coef.gsl_nls
#' @param model An object inheriting from class \code{"gsl_nls"}.
#' @return Numeric vector of Cook's distance values similar to \code{\link[stats]{cooks.distance}}.
#' @seealso \code{\link[stats]{cooks.distance}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'  x = (1:n) / n,
#'  y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' cooks.distance(obj)
#' @export
cooks.distance.gsl_nls <- function(model, ...) {
  df <- df.residual(model)
  res <- model$m$resid()        # weighted residuals
  ssr <- deviance(model)        # weighted ssr
  p <- length(res) - df
  h <- hatvalues(model)
  res^2 / (p * (ssr / df)) * (h / (1 - h)^2)
}

#' Anova tables
#' @description Returns the analysis of variance (or deviance) tables for two or
#' more fitted \code{"gsl_nls"} objects.
#' @inheritParams coef.gsl_nls
#' @param ... Additional objects inheriting from class \code{"gsl_nls"}.
#' @return A data.frame object of class \code{"anova"} similar to \code{\link[stats]{anova}} representing the
#' analysis-of-variance table of the fitted model objects when printed.
#' @seealso \code{\link[stats]{anova}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + 1 + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj1 <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#' obj2 <- gsl_nls(fn = y ~ A * exp(-lam * x) + b, data = xy,
#'     start = c(A = 1, lam = 1, b = 0))
#'
#' anova(obj1, obj2)
#' @export
anova.gsl_nls <- function(object, ...) {
  if(length(list(object, ...)) > 1L) {
    objects <- list(object, ...)
    responses <- vapply(objects, function(x)
      if(inherits(formula(x), "formula")) as.character(formula(x)[[2L]]) else "y",
      character(1))
    models <- lapply(objects, function(x)
      if(inherits(formula(x), "formula")) {
        formula(x)
      } else {
        sprintf("y ~ fn(%s)", paste(names(formals(formula(x))), collapse = ", "))
      })
    models <- as.character(models)
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
      objects <- objects[sameresp]
      warning(gettextf("models with response %s removed because response differs from model 1",
                       sQuote(deparse(responses[!sameresp]))),
              domain = NA)
    }
    ## calculate the number of models
    nmodels <- length(objects)
    if (nmodels == 1L)
      stop("'anova' is only defined for sequences of \"nls\" objects")
    ## extract statistics
    df.r <- unlist(lapply(objects, df.residual))
    ss.r <- unlist(lapply(objects, deviance))
    df <- c(NA, -diff(df.r))
    ss <- c(NA, -diff(ss.r))
    ms <- ss/df
    f <- p <- rep_len(NA_real_, nmodels)
    for(i in 2:nmodels) {
      if(df[i] > 0) {
        f[i] <- ms[i]/(ss.r[i]/df.r[i])
        p[i] <- pf(f[i], df[i], df.r[i], lower.tail = FALSE)
      }
      else if(df[i] < 0) {
        f[i] <- ms[i]/(ss.r[i-1]/df.r[i-1])
        p[i] <- pf(f[i], -df[i], df.r[i-1], lower.tail = FALSE)
      }
      else {                          # df[i] == 0
        ss[i] <- 0
      }
    }
    table <- data.frame(df.r,ss.r,df,ss,f,p)
    dimnames(table) <- list(1L:nmodels, c("Res.Df", "Res.Sum Sq", "Df",
                                          "Sum Sq", "F value", "Pr(>F)"))
    ## construct table and title
    title <- "Analysis of Variance Table\n"
    topnote <- paste0("Model ", format(1L:nmodels), ": ", models,
                      collapse = "\n")
    ## calculate test statistic if needed
    structure(table, heading = c(title, topnote),
              class = c("anova", "data.frame")) # was "tabular"
  } else {
    stop("anova is only defined for sequences of \"gsl_nls\" objects")
  }
}

#' Confidence interval for model parameters
#' @description Returns asymptotic or profile likelihood confidence intervals for the parameters in a
#' fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @param parm A character vector of parameter names for which to evaluate confidence intervals, defaults
#' to all parameters.
#' @param level A numeric scalar between 0 and 1 giving the level of the parameter confidence intervals.
#' @param method Method to be used, either \code{"asymptotic"} for asymptotic confidence intervals or
#' \code{"profile"} for profile likelihood confidence intervals. The latter is only available for
#' \code{"gsl_nls"} objects that are also of class \code{"nls"}.
#' @return A matrix with columns giving the lower and upper confidence limits for each parameter.
#' @details
#' Method \code{"asymptotic"} assumes (approximate) normality of the errors in the model and calculates
#' standard asymptotic confidence intervals based on the quantiles of a t-distribution. Method \code{"profile"}
#' calculates profile likelihood confidence intervals using the \code{\link[MASS:confint]{confint.nls}} method
#' in the \CRANpkg{MASS} package and for this reason is only available for \code{"gsl_nls"} objects that
#' are \emph{also} of class \code{"nls"}.
#' @seealso \code{\link[stats]{confint}}, \code{\link[MASS:confint]{confint.nls}} in package \CRANpkg{MASS}.
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' ## asymptotic ci's
#' confint(obj)
#' \dontrun{
#' ## profile ci's (requires MASS)
#' confint(obj, method = "profile")
#' }
#' @export
confint.gsl_nls <- function(object, parm, level = 0.95, method = c("asymptotic", "profile"), ...) {
  method <- match.arg(method)
  if(identical(method, "profile")) {
    if(inherits(object, "nls")) {
      NextMethod()
    } else {
      stop("method 'profile' can only be used for \"nls\" objects")
    }
  } else {
    ## from confint.default
    cf <- coef(object)
    pnames <- names(cf)
    if(missing(parm))
      parm <- pnames
    else if(is.numeric(parm))
      parm <- pnames[parm]
    a <- c((1 - level) / 2, (1 + level) / 2)
    ses <- sqrt(diag(vcov(object)))[parm]
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3L), "%")
    ci <- array(NA_real_, dim = c(length(parm), 2), dimnames = list(parm, pct))
    ci[] <- cf[parm] + ses %o% qt(a, df.residual(object))
    ## artifically limit intervals
    if(!is.null(object$lower) && any(ci[parm, 1L] < object$lower)) {
      ci[parm, 1L] <- pmax(ci[parm, 1L], object$lower)
    }
    if(!is.null(object$upper) && any(ci[parm, 2L] > object$upper)) {
      ci[parm, 2L] <- pmin(ci[parm, 2L], object$upper)
    }
    return(ci)
  }
}

#' Confidence intervals for derived parameters
#' @description \code{confintd} is a generic function to compute confidence intervals for continuous functions
#' of the parameters in a fitted model. The function invokes particular \emph{methods} which depend on the
#' \code{\link{class}} of the first argument.
#' @param object A fitted model object.
#' @param expr An expression or character vector that can be transformed to an \code{\link{expression}}
#' giving the function(s) of the parameters to be evaluated. Each expression should evaluate to a numeric scalar.
#' @param level A numeric scalar between 0 and 1 giving the level of the derived parameter confidence intervals.
#' @param ... Additional argument(s) for methods
#' @seealso \code{\link[stats]{confint}}
#' @return A matrix with columns giving the fitted values and lower and upper confidence limits for
#' each derived parameter. The row names list the individual derived parameter expressions.
#' @export
confintd <- function(object, expr, level = 0.95, ...) {
  UseMethod("confintd")
}

#' Confidence intervals for derived parameters
#' @description Returns fitted values and confidence intervals for continuous functions of parameters
#' from a fitted \code{"gsl_nls"} object.
#' @inheritParams confintd
#' @param dtype A character string equal to \code{"symbolic"} for \emph{symbolic} differentiation of
#' \code{expr} with \code{\link[stats]{deriv}}, or \code{"numeric"} for \emph{numeric} differentiation
#' of \code{expr} with \code{\link[stats]{numericDeriv}} using forward finite differencing.
#' @return A matrix with columns giving the fitted values and lower and upper confidence limits for
#' each derived parameter. The row names list the individual derived parameter expressions.
#' @details
#' This method assumes (approximate) normality of the errors in the model and confidence intervals are
#' calculated using the \emph{delta method}, i.e. a first-order Taylor approximation of the (continuous)
#' function of the parameters. If \code{dtype = "symbolic"} (the default), \code{expr} is differentiated
#' with respect to the parameters using symbolic differentiation with \code{\link[stats]{deriv}}. As such,
#' each expression in \code{expr} must contain only operators that are known to \code{\link[stats]{deriv}}.
#' If \code{dtype = "numeric"}, \code{expr} is differentiated using numeric differentiation with
#' \code{\link[stats]{numericDeriv}}, which should be used if \code{expr} cannot be derived symbolically
#' with \code{\link[stats]{deriv}}.
#' @seealso \code{\link[stats]{confint}}
#' @examples
#' ## data
#' set.seed(1)
#' n <- 25
#' xy <- data.frame(
#'   x = (1:n) / n,
#'   y = 2.5 * exp(-1.5 * (1:n) / n) + rnorm(n, sd = 0.1)
#' )
#' ## model
#' obj <- gsl_nls(fn = y ~ A * exp(-lam * x), data = xy, start = c(A = 1, lam = 1))
#'
#' ## delta method ci's
#' confintd(obj, expr = c("log(lam)", "A / lam"))
#' @export
confintd.gsl_nls <- function(object, expr, level = 0.95, dtype = "symbolic", ...) {
  ## prepare expression
  if(is.character(expr)) {
    exprstr <- expr
    expr <- str2expression(expr)
  } else {
    expr <- as.expression(expr)
    exprstr <- as.character(expr)
  }
  ## standard errors
  p <- coef(object)
  plist <- as.list(p)
  a <- c((1 - level) / 2, (1 + level) / 2)
  fit <- vapply(expr, eval, numeric(1), envir = plist)
  ## derive expression
  dtype <- match.arg(dtype, choices = c("symbolic", "numeric"))
  if(identical(dtype, "symbolic")) {
    dexpr <- lapply(expr, stats::deriv, namevec = names(p))
    dfit <- lapply(dexpr, eval, envir = plist)
  } else {
    dfit <- lapply(expr, stats::numericDeriv, theta = names(p), rho = list2env(plist))
  }
  grad <- t(vapply(dfit, attr, p, which = "gradient"))
  ses <- sigma(object) * sqrt(rowSums(grad %*% chol2inv(object$m$Rmat()) * grad))
  ci <- fit + ses %o% qt(a, df.residual(object))
  cimat <- cbind(fit = fit, lwr = ci[, 1], upr = ci[, 2])
  rownames(cimat) <- exprstr
  return(cimat)
}

convInfo <- function(x, digits, show. = getOption("show.nls.convergence", TRUE)) {
  if(!inherits(x, "summary.gsl_nls")) {
    cat(sprintf("\nAlgorithm: %s, (scaling: %s, solver: %s)\n", x$convInfo$trsName, x$control$scale, x$control$solver))
  }
  if(!is.null(x$irls)) {
    with(x$irls, {
      if (!irls_conv || show.) {
        cat("\nNumber of IRLS iterations", if (irls_conv)
          "to convergence:"
          else "till stop:", irls_niter, "\nAchieved IRLS tolerance:",
          format(irls_tol, digits = digits))
        cat("\n")
      }
      if (!irls_conv) {
        cat("Reason stopped:", irls_status)
        cat("\n")
      }
    })
  }
  if(is.null(x$irls) || x$irls$irls_conv) {
    with(x$convInfo, {
      if (!isConv || show.) {
        cat("\nNumber of", if(!is.null(x$irls)) "NLS iterations" else "iterations",
            if (isConv) "to convergence:" else "till stop:", finIter,
            "\nAchieved", if(!is.null(x$irls)) "NLS tolerance:" else "convergence tolerance:",
            format(finTol, digits = digits))
        cat("\n")
      }
      if (!isConv) {
        cat("Reason stopped:", stopMessage)
        cat("\n")
      }
    })
  }
  invisible()
}
