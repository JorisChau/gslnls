require(gslnls)

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
dotest("1.1", dim(nls_examples), c(53L, 5L))
dotest("1.2", nls_examples[c(1L, 31L), ], structure(list(name = c("Misra1a", "Rosenbrock"), class = c("formula", "function"),
                                                         p = c(2L, 2L), n = c(14L, 2L), check = c("p, n fixed", "p, n fixed")),
                                                    row.names = c(1L, 31L), class = "data.frame"))
## loop through all test formulas
for(i in 1:27) {
  nm <- nls_examples[i, "name"]
  n <- nls_examples[i, "n"]
  p <- nls_examples[i, "p"]
  nls_problem <- nls_test_problem(nm)
  dotest(sprintf("1.3.1: %s", nm), inherits(nls_problem, "nls_test_formula") && is.data.frame(nls_problem[["data"]]) &&
           inherits(nls_problem[["fn"]], "formula"), TRUE)
  dotest(sprintf("1.3.2: %s", nm), c(nrow(nls_problem[["data"]]), length(nls_problem[["start"]]), length(nls_problem[["target"]])), c(n, p, p))
}

## loop through all test functions
for(i in 28:53) {
  nm <- nls_examples[i, "name"]
  n <- nls_examples[i, "n"]
  p <- nls_examples[i, "p"]
  nls_problem <- nls_test_problem(nm)
  dotest(sprintf("1.4.1: %s", nm), inherits(nls_problem, "nls_test_function") && is.function(nls_problem[["fn"]]) && is.function(nls_problem[["jac"]]), TRUE)
  dotest(sprintf("1.4.2: %s", nm), c(length(nls_problem[["start"]]), length(nls_problem[["target"]])), c(p, p))
  fval <- nls_problem[["fn"]](x = nls_problem[["start"]])
  jval <- nls_problem[["jac"]](x = nls_problem[["start"]])
  dotest(sprintf("1.4.3: %s", nm), is.numeric(fval) && !any(is.na(fval)) && length(fval) == n, TRUE)
  dotest(sprintf("1.4.4: %s", nm), is.numeric(jval) && !any(is.na(jval)), TRUE)
  dotest(sprintf("1.4.5: %s", nm), dim(jval), c(n, p))
}

## gsl_nls

control <- matrix(
  c("lm", "more", "qr", "forward",
    "dogleg", "levenberg", "cholesky", "forward",
    "ddogleg", "marquardt", "svd", "center",
    "subspace2D", "more", "qr", "center"),
  byrow = TRUE,
  ncol = 4
)

## formula input
misra1a <- nls_test_problem("Misra1a")

misra1a_fit1 <- with(misra1a, gsl_nls(fn = fn, data = data, start = as.list(start), trace = TRUE))
dotest_tol("1.5.2", coef(misra1a_fit1), misra1a[["target"]])

misra1a_fit2 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "dogleg", jac = TRUE, control = list(scale = "levenberg")))
dotest_tol("1.5.3", coef(misra1a_fit2), misra1a[["target"]])

misra1a_fit3 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "lmaccel", fvv = TRUE, control = list(scale = "marquardt")))
dotest_tol("1.5.4", coef(misra1a_fit3), misra1a[["target"]])

misra1a_fit4 <- with(misra1a, gsl_nls(fn = fn, data = data, start = start, algorithm = "ddogleg", weights = rep(100.0, 14L), control = list(solver = "cholesky")))
dotest_tol("1.5.5", coef(misra1a_fit4), misra1a[["target"]])

## function input
linear1 <- nls_test_problem("Linear, full rank")

linear1_fit1 <- with(linear1, gsl_nls(fn = fn, y = y, start = start, jac = jac, algorithm = "subspace2D", trace = TRUE, control = list(solver = "svd")))
dotest_tol("1.5.6", coef(linear1_fit1), linear1[["target"]])

linear1_fit2 <- with(linear1, gsl_nls(fn = fn, y = y, start = start, weights = rep(100.0, 10L), control = list(fdtype = "center")))
dotest_tol("1.5.7", coef(linear1_fit2), linear1[["target"]])

## warning
# tools::assertWarning(with(linear1, gsl_nls(fn = fn, y = y, start = target, jac = jac, weights = rep(100.0, 10L))))

## gsl_nls_large

## formula input
misra1a_fit6 <- with(misra1a, gsl_nls_large(fn = fn, data = data, start = as.list(start), jac = TRUE, trace = TRUE))
dotest_tol("1.6.1", coef(misra1a_fit6), misra1a[["target"]])

misra1a_fit7 <- with(misra1a, gsl_nls_large(fn = fn, data = data, start = start, algorithm = "dogleg", jac = TRUE, control = list(scale = "levenberg")))
dotest_tol("1.6.2", coef(misra1a_fit7), misra1a[["target"]])

misra1a_fit8 <- with(misra1a, gsl_nls_large(fn = fn, data = data, start = start, algorithm = "lmaccel", jac = TRUE, fvv = TRUE, control = list(scale = "marquardt")))
dotest_tol("1.6.3", coef(misra1a_fit8), misra1a[["target"]])

## function input
linear1_fit3 <- with(linear1, gsl_nls_large(fn = fn, y = y, start = start, jac = jac, algorithm = "ddogleg", trace = TRUE, control = list(solver = "cholesky")))
dotest_tol("1.6.4", coef(linear1_fit3), linear1[["target"]])

linear1_fit4 <- with(linear1, gsl_nls_large(fn = fn, y = y, start = start, jac = jac, algorithm = "subspace2D", control = list(solver = "svd")))
dotest_tol("1.6.5", coef(linear1_fit4), linear1[["target"]])

linear1_fit5 <- with(linear1, gsl_nls_large(fn = fn, y = y, start = start, jac = jac, algorithm = "cgst"))
dotest_tol("1.6.6", coef(linear1_fit5), linear1[["target"]])

## multistart
boxbod <- nls_test_problem("BoxBOD")
boxbod3 <- nls_test_problem("Box 3-dimensional")

mstart <- list(b1 = c(200, 250), b2 = c(0, 1))
boxbod_fit <- with(boxbod, gsl_nls(fn = fn, data = data, start = mstart, trace = TRUE, control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.1", coef(boxbod_fit), boxbod[["target"]])

mstart <- list(x1 = c(1, 1), x2 = c(0, 20), x3 = c(1, 1))
boxbod3_fit <-  with(boxbod3, gsl_nls(fn = fn, y = y, start = mstart, trace = TRUE, jac = jac, control = list(mstart_n = 5, mstart_r = 1.1)))
dotest_tol("1.7.2", coef(boxbod3_fit), boxbod3[["target"]])

## miscellaneous

## self-start
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

ss_fit <- gsl_nls(
  fn =  y ~ SSasymp(x, Asym, R0, lrc),                   ## model formula
  data = data.frame(x = x, y = y)                        ## model fit data
)
dotest("1.8.1", inherits(ss_fit, "gsl_nls"), TRUE)

## sparse jacobian
p <- 10
sp_fit <- gsl_nls_large(
  fn = function(theta) c(sqrt(1e-5) * (theta - 1), sum(theta^2) - 0.25),
  y = rep(0, p + 1),
  start = rep(0.15, p),
  jac = function(theta) rbind(Matrix::Diagonal(x = sqrt(1e-5), n = length(theta)), 2 * t(theta))
)
dotest("1.8.2", inherits(sp_fit, "gsl_nls"), TRUE)

## gsl_nls methods

dotest("1.9.1", length(fitted(misra1a_fit1)), 14L)
dotest("1.9.2", nobs(misra1a_fit1), 14L)
dotest_tol("1.9.3", deviance(misra1a_fit1), .1245514)
dotest_tol("1.9.4", sigma(misra1a_fit1), .1018788)
dotest("1.9.5", names(summary(misra1a_fit1)), c("formula", "residuals", "sigma", "df", "cov.unscaled", "call",
                                               "convInfo", "control", "na.action", "coefficients", "parameters"))
dotest("1.9.6", dim(predict(misra1a_fit1, interval = "prediction")), c(14L, 3L))
dotest("1.9.7", length(predict(misra1a_fit1, newdata = misra1a$data)), 14L)
dotest("1.9.8", length(residuals(misra1a_fit1)), 14L)
dotest("1.9.9", inherits(logLik(misra1a_fit1), "logLik"), TRUE)
dotest("1.9.10", df.residual(misra1a_fit1), 12L)
dotest("1.9.11", dim(vcov(misra1a_fit1)), c(2L, 2L))
dotest("1.9.12", names(anova(misra1a_fit1, misra1a_fit2)), c("Res.Df", "Res.Sum Sq", "Df", "Sum Sq", "F value", "Pr(>F)"))
dotest("1.9.13", dim(confint(misra1a_fit1)), c(2L, 2L))
dotest("1.9.14", dim(confintd(misra1a_fit1, expr = "b1 + b2")), c(1L, 3L))
dotest("1.9.15", capture.output(misra1a_fit1)[c(1, 2)], c("Nonlinear regression model", "  model: y ~ b1 * (1 - exp(-b2 * x))"))

dotest("1.10.1", length(fitted(linear1_fit1)), 10L)
dotest("1.10.2", nobs(linear1_fit1), 10L)
dotest_tol("1.10.3", deviance(linear1_fit1), 5.0)
dotest_tol("1.10.4", sigma(linear1_fit1), 1.0)
dotest("1.10.5", names(summary(linear1_fit1)), c("formula", "residuals", "sigma", "df", "cov.unscaled", "call",
                                               "convInfo", "control", "na.action", "coefficients", "parameters"))
dotest("1.10.6", dim(predict(linear1_fit1, interval = "prediction")), c(10L, 3L))
dotest("1.10.7", length(residuals(linear1_fit1, type = "pearson")), 10L)
dotest("1.10.8", inherits(logLik(linear1_fit1), "logLik"), TRUE)
dotest("1.10.9", df.residual(linear1_fit1), 5L)
dotest("1.10.10", dim(vcov(linear1_fit1)), c(5L, 5L))
dotest("1.10.11", names(anova(linear1_fit1, linear1_fit2)), c("Res.Df", "Res.Sum Sq", "Df", "Sum Sq", "F value", "Pr(>F)"))
dotest("1.10.12", dim(confint(linear1_fit1)), c(5L, 2L))
dotest("1.10.13", dim(confintd(linear1_fit1, expr = "x1 - x2")), c(1L, 3L))

cat("Completed gslnls unit tests\n")
