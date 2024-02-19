require(gslnls)
require(Matrix)

cat("Running gsnls unit tests...\n")

## test helper function
dotest <- function(itest, observed, expected) {
  if(!identical(observed, expected)) stop(sprintf("Test %s failed", itest), call. = FALSE)
}
dotest_tol <- function(itest, observed, expected, tol = (.Machine$double.eps)^0.25) {
  if(any(abs(observed - expected) > tol)) stop(sprintf("Test %s failed", itest), call. = FALSE)
}

## check example problems
nls_examples <- nls_test_list()
dotest("1.1", dim(nls_examples), c(59L, 5L))
dotest("1.2", nls_examples[c(1L, 37L), ], structure(list(name = c("Misra1a", "Rosenbrock"), class = c("formula", "function"),
                                                         p = c(2L, 2L), n = c(14L, 2L), check = c("p, n fixed", "p, n fixed")),
                                                    row.names = c(1L, 37L), class = "data.frame"))
## loop through all test formulas
for(i in 1:33) {
  nm <- nls_examples[i, "name"]
  n <- nls_examples[i, "n"]
  p <- nls_examples[i, "p"]
  nls_problem <- nls_test_problem(nm)
  dotest(sprintf("1.3.1: %s", nm), inherits(nls_problem, "nls_test_formula") && is.data.frame(nls_problem[["data"]]) &&
           inherits(nls_problem[["fn"]], "formula"), TRUE)
  dotest(sprintf("1.3.2: %s", nm), c(nrow(nls_problem[["data"]]), length(nls_problem[["start"]]), length(nls_problem[["target"]])), c(n, p, p))
}

## loop through all test functions
for(i in 34:59) {
  nm <- nls_examples[i, "name"]
  n <- nls_examples[i, "n"]
  p <- nls_examples[i, "p"]
  nls_problem <- nls_test_problem(nm)
  dotest(sprintf("1.4.1: %s", nm), inherits(nls_problem, "nls_test_function") && is.function(nls_problem[["fn"]]) && is.function(nls_problem[["jac"]]), TRUE)
  dotest(sprintf("1.4.2: %s", nm), c(length(nls_problem[["start"]]), length(nls_problem[["target"]])), c(p, p))
  fval <- nls_problem[["fn"]](x = nls_problem[["start"]])
  jval <- nls_problem[["jac"]](x = nls_problem[["start"]])
  dotest(sprintf("1.4.3: %s", nm), identical(length(fval), n) && all(is.finite(fval)), TRUE)
  dotest(sprintf("1.4.4: %s", nm), identical(dim(jval), c(n, p)) && all(is.finite(jval)), TRUE)
}

## gsl_nls

## formula input
misra1a <- nls_test_problem("Misra1a")

misra1a_fit1 <- with(misra1a, gsl_nls(fn = fn, data = data, start = as.list(start), trace = TRUE, model = TRUE))
dotest_tol("1.5.1", coef(misra1a_fit1), misra1a[["target"]])
misra1a_fit2 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "dogleg", jac = TRUE, control = list(scale = "levenberg")))
dotest_tol("1.5.2", coef(misra1a_fit2), misra1a[["target"]])
misra1a_fit3 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "lmaccel", fvv = TRUE, control = list(scale = "marquardt")))
dotest_tol("1.5.3", coef(misra1a_fit3), misra1a[["target"]])
misra1a_fit4 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "lm", weights = rep(100.0, 14L), control = list(solver = "cholesky")))
dotest_tol("1.5.4", coef(misra1a_fit4), misra1a[["target"]])
misra1a_fit5 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "ddogleg", lower = c(b1 = 100, b2 = 0)), control = list(solver = "svd"))
dotest_tol("1.5.5", coef(misra1a_fit5), misra1a[["target"]])
misra1a_fit6 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "subspace2D", upper = list(b1 = 500, b2 = 1)), control = list(fdtype = "center"))
dotest_tol("1.5.6", coef(misra1a_fit6), misra1a[["target"]])
misra1a_fit7 <- with(misra1a, gsl_nls(fn = fn, data = data, start = c(b1 = 300, b2 = 0), algorithm = "lm", jac = TRUE, lower = list(b1 = 250), upper = c(b2 = 1)))
dotest_tol("1.5.7", coef(misra1a_fit7), c(b1 = 250, misra1a[["target"]]["b2"]))

tools::assertError(gsl_nls(fn = ~ b1 * (1 - exp(-b2 * x)), data = misra1a$data))
tools::assertError(gsl_nls(fn = misra1a$fn, data = misra1a$data, start = c(b1 = 0)))
tools::assertError(gsl_nls(fn = misra1a$fn, data = list(x = misra1a$data$x, y = misra1a$data$y[1]), start = misra1a$start, lower = c(b1 = -Inf), upper = c(b1 = Inf)))

## function input
linear1 <- nls_test_problem("Linear, full rank", n = 5, p = 5)
linear1[["start"]] <- c(x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0)
linear1[["target"]] <- c(x1 = -1, x2 = -1, x3 = -1, x4 = -1, x5 = -1)

linear1_fit1 <- with(linear1, gsl_nls(fn = fn, y = y, start = start, algorithm = "lmaccel",  jac = jac, fvv = TRUE, trace = TRUE))
dotest_tol("1.5.8", coef(linear1_fit1), linear1[["target"]])
linear1_fit2 <- with(linear1, gsl_nls(fn = fn, y = y, start = start, algorithm = "dogleg", weights = rep(100.0, 5L)))
dotest_tol("1.5.9", coef(linear1_fit2), linear1[["target"]])
tools::assertWarning(linear1_fit3 <- with(linear1, gsl_nls(fn = fn, y = y, start = start, algorithm = "ddogleg", lower = start, upper = as.list(start))))
dotest_tol("1.5.10", coef(linear1_fit3), linear1[["start"]])

linear1_fn <- function(pars) {  ## access parameters as list
  fn <- numeric(5L)
  nms <- paste0("x", 1:5)
  for(i in 1:5) {
    fn[i] <- pars[[i]] - 2 * (pars$x1 + pars$x2 + pars$x3 + pars$x4 + pars$x5) / 5 - 1
  }
  attr(fn, "gradient") <- matrix(-2 / 5, nrow = 5L, ncol = 5L, dimnames = list(NULL, nms)) + diag(rep(1, 5L))
  attr(fn, "hessian") <- array(0, dim = c(5L, 5L, 5L), dimnames = list(NULL, nms, nms))
  return(fn)
}

linear1_fit4 <- with(linear1, gsl_nls(fn = linear1_fn, y = y, start = as.list(start), algorithm = "lmaccel"))
dotest_tol("1.5.11", coef(linear1_fit4), linear1[["target"]])
linear1_fit5 <- with(linear1, gsl_nls(fn = linear1_fn, y = y, start = as.list(start), algorithm = "subspace2D", lower = target / 2, upper = as.list(start)))
dotest_tol("1.5.12", coef(linear1_fit5), linear1[["target"]] / 2)

## gsl_nls_large

## formula input
misra1a_fit6 <- with(misra1a, gsl_nls_large(fn = fn, data = data, start = as.list(start), jac = TRUE, trace = TRUE))
dotest_tol("1.6.1", coef(misra1a_fit6), misra1a[["target"]])
misra1a_fit7 <- with(misra1a, gsl_nls_large(fn = fn, data = data, start = start, algorithm = "dogleg", jac = TRUE, control = list(scale = "levenberg")))
dotest_tol("1.6.2", coef(misra1a_fit7), misra1a[["target"]])
misra1a_fit8 <- with(misra1a, gsl_nls_large(fn = fn, data = data, start = start, algorithm = "lmaccel", jac = TRUE, fvv = TRUE, control = list(scale = "marquardt")))
dotest_tol("1.6.3", coef(misra1a_fit8), misra1a[["target"]])
misra1a_fit9 <- with(misra1a, gsl_nls_large(fn = fn, data = data, start = start, algorithm = "lm", weights = rep(1.0, 14L), jac = TRUE))
dotest_tol("1.6.4", coef(misra1a_fit9), misra1a[["target"]])

tools::assertError(gsl_nls_large(fn =  ~ b1 * (1 - exp(-b2 * x)), data = misra1a$data))
tools::assertError(gsl_nls_large(fn = misra1a$fn, data = misra1a$data, start = c(b1 = 0)))
tools::assertError(gsl_nls_large(fn = misra1a$fn, data = list(x = misra1a$data$x, y = misra1a$data$y[1]), start = misra1a$start))

## function input
linear1_fit6 <- with(linear1, gsl_nls_large(fn = fn, y = y, start = start, jac = jac, trace = TRUE, control = list(solver = "cholesky")))
dotest_tol("1.6.5", coef(linear1_fit6), linear1[["target"]])
linear1_fit7 <- with(linear1, gsl_nls_large(fn = fn, y = y, start = start, jac = jac, algorithm = "subspace2D", control = list(solver = "svd")))
dotest_tol("1.6.6", coef(linear1_fit7), linear1[["target"]])
linear1_fit8 <- with(linear1, gsl_nls_large(fn = fn, y = y, start = start, jac = jac, algorithm = "cgst"))
dotest_tol("1.6.7", coef(linear1_fit8), linear1[["target"]])
linear1_fit9 <- with(linear1, gsl_nls_large(fn = linear1_fn, y = y, start = as.list(start), algorithm = "lmaccel"))
dotest_tol("1.6.8", coef(linear1_fit9), linear1[["target"]])
linear1_fit10 <- with(linear1, gsl_nls_large(fn = fn, y = y, start = start, jac = jac, weights = rep(1.0, 5L)))
dotest_tol("1.6.9", coef(linear1_fit10), linear1[["target"]])

## multistart
boxbod <- nls_test_problem("BoxBOD")
madsen <- nls_test_problem("Madsen example")

boxbod_fit1 <- with(boxbod, gsl_nls(fn = fn, data = data, start = list(b1 = c(200, 250), b2 = c(0, 1)), trace = TRUE, control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.1", coef(boxbod_fit1), boxbod[["target"]])
boxbod_fit2 <- with(boxbod, gsl_nls(fn = fn, data = data, start = cbind(b1 = c(200, 250), b2 = c(0, 1)), algorithm = "lmaccel", fvv = TRUE, control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.2", coef(boxbod_fit2), boxbod[["target"]])
boxbod_fit3 <- with(boxbod, gsl_nls(fn = fn, data = data, start = list(b1 = c(200, 250), b2 = 1), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.3", coef(boxbod_fit3), boxbod[["target"]])
boxbod_fit4 <- with(boxbod, gsl_nls(fn = fn, data = data, start = c(b1 = 200, b2 = NA), lower = c(b2 = 0), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.4", coef(boxbod_fit4), boxbod[["target"]])
boxbod_fit5 <- with(boxbod, gsl_nls(fn = fn, data = data, start = list(b1 = NA, b2 = 0.5), jac = TRUE, upper = c(b2 = 1), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.5", coef(boxbod_fit5), boxbod[["target"]])
boxbod_fit6 <- with(boxbod, gsl_nls(fn = fn, data = data, start = cbind(b1 = NA, b2 = c(0, 1)), lower = list(b1 = 0), upper = list(b1 = 250), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.6", coef(boxbod_fit6), boxbod[["target"]])
boxbod_fit7 <- with(boxbod, gsl_nls(fn = fn, data = data, start = list(b1 = c(200, 250), b2 = NA), jac = TRUE, lower = list(b2 = 0), upper = c(b2 = 1), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.7", coef(boxbod_fit7), boxbod[["target"]])

madsen_fit1 <- with(madsen, gsl_nls(fn = fn, y = y, start = cbind(x1 = c(-1, 1), x2 = c(0, 1)), trace = TRUE, control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.8", coef(madsen_fit1), madsen[["target"]])
madsen_fit2 <- with(madsen, gsl_nls(fn = fn, y = y, start = c(x1 = 0, x2 = NA), jac = jac, lower = c(x2 = 0), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.9", coef(madsen_fit2), madsen[["target"]])
madsen_fit3 <- with(madsen, gsl_nls(fn = fn, y = y, start = c(x1 = NA, x2 = 0), upper = c(x2 = 1), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.10", coef(madsen_fit3), madsen[["target"]])
madsen_fit4 <- with(madsen, gsl_nls(fn = fn, y = y, start = c(x1 = NA, x2 = NA), jac = jac, control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.11", coef(madsen_fit4), madsen[["target"]])
madsen_fit5 <- with(madsen, gsl_nls(fn = fn, y = y, start = cbind(x1 = c(-0.5, -0.5), x2 = c(1, 1)), jac = jac, lower = c(x1 = -Inf)))
dotest_tol("1.7.12", coef(madsen_fit4), madsen[["target"]])

linear1_fit10 <- with(linear1, gsl_nls(fn = linear1_fn, y = y, start = list(x1 = NA, x2 = 0, x3 = 0, x4 = 0, x5 = 0), algorithm = "lmaccel", lower = list(x1 = -5), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.13", coef(linear1_fit10), linear1[["target"]])
linear1_fit11 <- with(linear1, gsl_nls(fn = linear1_fn, y = y, start = list(x1 = c(-5, 0), x2 = 0, x3 = 0, x4 = 0, x5 = NA), lower = list(x1 = -5, x5 = -5), control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.14", coef(linear1_fit11), linear1[["target"]])

## miscellaneous

x <- c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25,
       1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625,
       2.75, 2.875, 3)
y <- c(5.84338654731442, 5.19105642195752, 4.22753924085235, 4.24773432418906,
       3.44420970665891, 2.75291103735449, 2.74511959989887, 2.53031291992822,
       2.25959613865552, 1.84855990274743, 2.14472012633735, 1.73313947376304,
       1.37168597767387, 0.883220159182727, 1.64343151470703, 1.28903993722273,
       1.24488777606458, 1.44233369960031, 1.37639589033285, 1.29031441255294,
       1.3473330721821, 1.29301855572576, 1.09945871878213, 0.569662114359861,
       1.21050141966489)

## self-start
ss_fit <- gsl_nls(
  fn =  y ~ SSasymp(x, Asym, R0, lrc),                   ## model formula
  data = data.frame(y = y, x = x),                       ## model fit data
)
ss_fit_large <- gsl_nls_large(
  fn =  y ~ SSasymp(x, Asym, R0, lrc),                   ## model formula
  data = data.frame(x = x, y = y)                        ## model fit data
)
## non self-start
SSasymp2 <- deriv(
  ~Asym+(R0-Asym)*exp(-exp(lrc)*input),
  namevec = c("Asym", "R0", "lrc"),
  function.arg = c("input", "Asym", "R0", "lrc"),
  hessian = TRUE
)
ss_fit2 <- gsl_nls(
  fn = y ~ SSasymp2(x, Asym, R0, lrc),
  data = data.frame(x = x, y = y),
  start = c(Asym = 1, R0 = 6, lrc = 0.25),
  algorithm = "lmaccel"
)
ss_fit2_large <- gsl_nls_large(
  fn = y ~ SSasymp2(x, Asym, R0, lrc),
  data = data.frame(x = x, y = y),
  start = c(Asym = 1, R0 = 6, lrc = 0.25),
  algorithm = "lmaccel"
)

dotest("1.8.1", inherits(ss_fit, "gsl_nls"), TRUE)
dotest("1.8.2", inherits(ss_fit_large, "gsl_nls"), TRUE)
dotest("1.8.3", inherits(ss_fit2, "gsl_nls"), TRUE)
dotest("1.8.4", inherits(ss_fit2_large, "gsl_nls"), TRUE)

## sparse jacobian
misra1a_fit_sp <- gsl_nls_large(
  fn = misra1a$fn,
  data = misra1a$data,
  start = misra1a$start,
  jac = function(par, x) {
    .expr3 <- exp(-par[["b2"]] * x)
    .expr4 <- 1 - .expr3
    .grad <- Matrix::Matrix(c(.expr4, par[["b1"]] * (.expr3 * x)), nrow = length(x), ncol = 2L,
                            dimnames = list(NULL, c("b1", "b2")), sparse = TRUE)
    return(.grad)
  },
  x = misra1a$data$x
)

dotest_tol("1.8.5", coef(misra1a_fit_sp), misra1a[["target"]])

penalty_fit_dgC <- gsl_nls_large(
  fn = function(theta) c(sqrt(1e-5) * (theta - 1), sum(theta^2) - 0.25),
  y = rep(0, 11L),
  start = rep(0.15, 10L),
  jac = function(theta) rbind(Matrix::Diagonal(x = sqrt(1e-5), n = length(theta)), 2 * t(theta))
)
penalty_fit_dgR <- gsl_nls_large(
  fn = function(theta) c(sqrt(1e-5) * (theta - 1), sum(theta^2) - 0.25),
  y = rep(0, 11L),
  start = rep(0.15, 10L),
  jac = function(theta) as(rbind(Matrix::Diagonal(x = sqrt(1e-5), n = length(theta)), 2 * t(theta)), Class = "RsparseMatrix")
)
penalty_fit_dgT <- gsl_nls_large(
  fn = function(theta) c(sqrt(1e-5) * (theta - 1), sum(theta^2) - 0.25),
  y = rep(0, 11L),
  start = rep(0.15, 10L),
  jac = function(theta) as(rbind(diag(rep(sqrt(1e-5), length(theta))), 2 * t(theta)), Class = "dgTMatrix")
)
penalty_fit_dge <- gsl_nls_large(
  fn = function(theta) c(sqrt(1e-5) * (theta - 1), sum(theta^2) - 0.25),
  y = rep(0, 11L),
  start = rep(0.15, 10L),
  jac = function(theta) as(rbind(diag(rep(sqrt(1e-5), length(theta))), 2 * t(theta)), Class = "dgeMatrix")
)

dotest_tol("1.8.6", inherits(penalty_fit_dgC, "gsl_nls"), TRUE)
dotest_tol("1.8.7", inherits(penalty_fit_dgR, "gsl_nls"), TRUE)
dotest_tol("1.8.8", inherits(penalty_fit_dgT, "gsl_nls"), TRUE)
dotest_tol("1.8.9", inherits(penalty_fit_dge, "gsl_nls"), TRUE)

## gsl_nls methods

dotest("1.9.1", length(fitted(misra1a_fit1)), 14L)
dotest("1.9.2", nobs(misra1a_fit1), 14L)
dotest_tol("1.9.3", deviance(misra1a_fit1), .1245514)
dotest_tol("1.9.4a", sigma(misra1a_fit1), .1018788)
dotest_tol("1.9.4b", sigma(misra1a_fit4), .1018788)
dotest("1.9.5", names(summary(misra1a_fit1)), c("formula", "residuals", "sigma", "df", "cov.unscaled", "call",
                                                "convInfo", "control", "na.action", "coefficients", "parameters"))
dotest("1.9.6", dim(predict(misra1a_fit1, interval = "prediction")), c(14L, 3L))
dotest("1.9.7", dim(predict(misra1a_fit1, newdata = misra1a$data, interval = "confidence")), c(14L, 3L))
dotest("1.9.8", length(residuals(misra1a_fit1)), 14L)
dotest("1.9.9", inherits(logLik(misra1a_fit1), "logLik"), TRUE)
dotest("1.9.10", df.residual(misra1a_fit1), 12L)
dotest("1.9.11", dim(vcov(misra1a_fit1)), c(2L, 2L))
dotest("1.9.12", names(anova(misra1a_fit1, misra1a_fit2)), c("Res.Df", "Res.Sum Sq", "Df", "Sum Sq", "F value", "Pr(>F)"))
if(requireNamespace("MASS")) {
  dotest("1.9.13", dim(confint(misra1a_fit1, method = "profile")), c(2L, 2L))
}
dotest("1.9.14", dim(confintd(misra1a_fit1, expr = "b1 + b2")), c(1L, 3L))
dotest("1.9.15", capture.output(misra1a_fit1)[c(1, 2)], c("Nonlinear regression model", "  model: y ~ b1 * (1 - exp(-b2 * x))"))

dotest("1.10.1", length(fitted(madsen_fit1)), 3L)
dotest("1.10.2", nobs(madsen_fit1), 3L)
dotest_tol("1.10.3", deviance(madsen_fit1), 0.7731991)
dotest_tol("1.10.4", sigma(madsen_fit1), 0.8793174)
dotest("1.10.5", names(summary(madsen_fit1)), c("formula", "residuals", "sigma", "df", "cov.unscaled", "call",
                                                "convInfo", "control", "na.action", "coefficients", "parameters"))
dotest("1.10.6", dim(predict(madsen_fit1, interval = "prediction")), c(3L, 3L))
dotest("1.10.7", dim(predict(madsen_fit1, newdata = data.frame(), interval = "confidence")), c(3L, 3L))
dotest("1.10.8a", length(residuals(madsen_fit1, type = "pearson")), 3L)
dotest("1.10.8b", length(residuals(madsen_fit1, type = "response")), 3L)
dotest("1.10.9", inherits(logLik(madsen_fit1), "logLik"), TRUE)
dotest("1.10.10", df.residual(madsen_fit1), 1L)
dotest("1.10.11", dim(vcov(madsen_fit1)), c(2L, 2L))
dotest("1.10.12", names(anova(madsen_fit1, madsen_fit2)), c("Res.Df", "Res.Sum Sq", "Df", "Sum Sq", "F value", "Pr(>F)"))
dotest("1.10.13", dim(confint(madsen_fit1, parm = c(1L, 2L))), c(2L, 2L))
dotest("1.10.14", dim(confintd(madsen_fit1, expr = quote(x1 - x2), dtype = "numeric")), c(1L, 3L))

cat("Completed gslnls unit tests\n")
