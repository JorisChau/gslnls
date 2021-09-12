#' Extract model coefficients
#' @description Returns the fitted model coefficients from a \code{"gsl_nls"} object.
#' \code{coefficients} can also be used as an alias.
#' @param object An object inheriting from class \code{"gsl_nls"}.
#' @param ... At present no optional arguments are used.
#' @seealso \code{\link[stats]{coef}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
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
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
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
#' @seealso \code{\link[stats]{nobs}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
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
#' @seealso \code{\link[stats]{deviance}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
#' deviance(obj)
#' @export
deviance.gsl_nls <- function(object, ...) {
  object$m$deviance()
}

#' Extract model formula
#' @description Returns the model formula from a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @param x An object inheriting from class \code{"gsl_nls"}.
#' @seealso \code{\link[stats]{formula}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
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
#' @seealso \code{\link{summary}}
#' @examples 
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
#' summary(obj)
#' @export
summary.gsl_nls <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  if(inherits(object, "nls")) {
    NextMethod()
  } else {
    r <- as.vector(object$m$resid()) # These are weighted residuals.
    w <- object$weights
    n <- if (!is.null(w)) sum(w > 0) else length(r)
    param <- coef(object)
    pnames <- names(param)
    p <- length(param)
    rdf <- n - p
    resvar <- if(rdf <= 0) NaN else deviance(object)/rdf
    XtXinv <- chol2inv(object$m$Rmat())
    dimnames(XtXinv) <- list(pnames, pnames)
    se <- sqrt(diag(XtXinv) * resvar)
    tval <- param/se
    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <-
        list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans <- list(formula = as.formula(sprintf("y ~ fn(%s)", 
                paste(names(formals(formula(object))), collapse = ", "))), 
        residuals = r, sigma = sqrt(resvar),
        df = c(p, rdf), cov.unscaled = XtXinv,
        call = object$call,
        convInfo = object$convInfo,
        control = object$control,
        na.action = object$na.action,
        coefficients = param,
        parameters = param)# never documented, for back-compatibility
    if(correlation && rdf > 0) {
      ans$correlation <- (XtXinv * resvar)/outer(se, se)
      ans$symbolic.cor <- symbolic.cor
    }
    ## if(identical(object$call$algorithm, "port"))
    ##     ans$message <- object$message
    class(ans) <- "summary.nls"
    ans
  }
}

#' Print model object
#' @noRd
#' @export
print.gsl_nls <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Nonlinear regression model\n")
  if(inherits(formula(x), "formula")) {
    cat("  model: ", deparse(formula(x)), "\n", sep = "")
    cat("   data: ", deparse(x$data), "\n", sep = "")
  } else {
    cat("  model: ", sprintf("y ~ fn(%s)", paste(names(formals(formula(x))), collapse = ", ")), "\n", sep = "")
  }  
  print(x$m$getAllPars(), digits = digits, ...)
  cat(" ", if(!is.null(x$weights) && diff(range(x$weights))) "weighted ",
      "residual sum-of-squares: ", format(x$m$deviance(), digits = digits),
      "\n", sep = "")
  convInfo(x, digits = digits)
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
#' @seealso \code{\link[stats]{predict}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
#' predict(obj)
#' predict(obj, newdata = data.frame(x = 1:20))
#' predict(obj, interval = "confidence")
#' predict(obj, interval = "prediction", level = 0.99)
#' @export
predict.gsl_nls <- function(object, newdata, scale = NULL, interval = c("none", "confidence", "prediction"), level = 0.95, ...) {
  interval <- match.arg(interval, c("none", "confidence", "prediction"))
  if (missing(newdata)) {
    fit <- as.vector(fitted(object))
    if(interval != "none")
      Fdot <- object$m$gradient()
  } else {
    fit <- object$m$predict(newdata)
    if(interval != "none") {
      Fdot <- object$m$gradient1(newdata)
      if(!is.null(object$weights))
        warning("unweighted Jacobian matrix used to calculate standard errors, evaluate predictions without 'newdata' argument to use weighted Jacobian.")
    }
  }
  if(interval != "none") {
    if(is.null(scale))
      scale <- sigma(object)
    a <- c((1 - level) / 2, (1 + level) / 2)
    ses <- scale * sqrt(1 * (interval == "prediction") + rowSums(Fdot %*% solve(t(Fdot) %*% Fdot) * Fdot)) 
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
#' @seealso \code{\link[stats]{residuals}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
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
      lab <- "Residuals"
      attr(val, "label") <- lab
    }
    val
  }
}

#' Extract model log-likelihood
#' @description Returns the model log-likelihood of a fitted \code{"gsl_nls"} object.
#' @inheritParams coef.gsl_nls
#' @param REML logical value; included for compatibility reasons only, should not be used.
#' @seealso \code{\link[stats]{logLik}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
#' residuals(obj)
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
#' @seealso \code{\link[stats]{df.residual}}
#' @examples 
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
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
#' @seealso \code{\link[stats]{vcov}}
#' @examples
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## model
#' obj <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
#' vcov(obj)
#' @export
vcov.gsl_nls <- function(object, ...) {
  sm <- summary(object)
  sm$cov.unscaled * sm$sigma^2
}

#' Anova tables
#' @description Returns the analysis of variance (or deviance) tables for two or 
#' more fitted \code{"gsl_nls"} objects.
#' @inheritParams coef.gsl_nls
#' @param ... Additional objects inheriting from class \code{"gsl_nls"}.
#' @seealso \code{\link[stats]{anova}}
#' @examples 
#' ## data
#' set.seed(1)
#' x <- 1:10
#' y <- 1 + 0.5 * (1:10) + rnorm(10, sd = 0.1)
#' ## models
#' obj1 <- gsl_nls(fn = y ~ b * x, data = data.frame(x = x, y = y), start = c(b = 1))
#' obj2 <- gsl_nls(fn = y ~ a + b * x, data = data.frame(x = x, y = y), start = c(a = 0, b = 1))
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

convInfo <- function(x, digits, show. = getOption("show.nls.convergence", TRUE)) {
  with(x$convInfo, {
        cat(sprintf("\nAlgorithm: %s, (scaling: %s, solver: %s)\n", trsName, x$control$scale, x$control$solver))
        if(!isConv || show.) {
          cat("\nNumber of iterations",
              if(isConv) "to convergence:" else "till stop:", finIter,
              "\nAchieved convergence tolerance:",
              format(finTol, digits = digits))
          cat("\n")
        }
        if(!isConv) {
          cat("Reason stopped:", stopMessage)
          cat("\n")
        }
      })
  
  invisible()
}
