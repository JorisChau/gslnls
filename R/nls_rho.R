#' Loss functions with tunable parameters
#'
#' Allow the user to tune the coefficient(s) of the loss functions
#' available in \code{\link{gsl_nls}} and \code{\link{gsl_nls_large}}.
#'
#' @param rho character loss function name, one of \code{"huber", "bisquare", "hampel", "ggw", "lqq", "optimal", "welsh"}.
#' @param cc numeric vector of loss function tuning parameters, with length depends on the chosen loss function.
#' If \code{NULL}, the default tuning parameters are returned.
#' @seealso \url{https://CRAN.R-project.org/package=robustbase/vignettes/psi_functions.pdf}
#' @examples
#' ## huber loss function with default tuning parameter
#' huber_loss()
#' @return A \code{list} with two components:
#' \itemize{
#' \item rho
#' \item cc
#' }
#' with meanings as explained under 'Arguments'.
#' @export
gsl_nls_loss <- function(rho = c("default", "huber", "barron", "bisquare", "welsh", "optimal", "hampel", "ggw", "lqq"), cc = NULL) {

  rho <- match.arg(rho, c("default", "huber", "barron", "bisquare", "welsh", "optimal", "hampel", "ggw", "lqq"))

  ## default tuning parameters
  cc_default <- switch(rho,
       default = numeric(0),   ## not used
       huber = 1.345,
       barron = c(0.0, 1.345),  ## cauchy loss
       bisquare = 4.685061,
       welsh = 2.11,
       optimal = 1.060158,
       hampel = c(1.352413, 3.155630, 7.212868),
       ggw = c(0.648, 1.0, 1.694),
       lqq = c(1.473, 0.982, 1.5)
  )

  if(is.null(cc) || identical(rho, "default")) {
    cc <- cc_default
  } else if(length(cc) != length(cc_default)) {
    stop(sprintf("'cc' must be of length %d for function '%s'", length(cc_default), rho))
  } else {
    cc <- as.numeric(cc)
    if(identical(rho, "barron") && cc[1] > 2) {
      warning("Robustness parameter (alpha) in Barron loss function cannot be larger than 2")
      cc[1] <- 2
    }
  }

  list(rho = rho, cc = cc)

}
