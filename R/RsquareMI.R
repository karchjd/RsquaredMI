#' @export
print.RsquaredPooled <- function(x, ...) {
  cat("R-squared SP:", "\n")
  print(x$rtotal)
  cat("\n")
  cat("Beta Coefficients SP:", "\n")
  print(x$total)
  if (!is.null(x$alternative_adj_R2)) {
    cat("\n")
    cat("Alternative adjusted R^2 estimates:", "\n")
    print(x$alternative_adj_R2)
  }
}
#' Calculate R-squared with Standardized Predictors
#'
#' This function calculates the R-squared value for a linear model applied to a
#' multiply imputed dataset, along with standardized regression coefficients.
#' Optionally, it can also return the confidence intervals of the standardized
#' regression coefficients and the zero-order correlations.
#'
#' @param object The results of a regression on a multiply imputed dataset of
#' class `mira` from the `mice` package.
#' @param conf Logical. If `TRUE`, the function returns the confidence intervals
#' of the standardized regression coefficients.
#' @param cor Logical. If `TRUE`, the function returns the zero-order
#' correlations between the outcome and each predictor.
#' @param conf.level A real number between 0 and 1 specifying the confidence level
#' of the confidence intervals.
#' @param alternative_adj_R2 Logical. If `TRUE`, the function returns alternative
#' estimates of adjusted R^2, as described in the references
#' @return A list of class `RsquaredMI` containing the following elements:
#' \item{r_squared}{The R-squared value calculated using standardized predictors.}
#' \item{r}{The square root of the R-squared value, or the multiple correlation R.}
#' \item{rtotal}{A vector containing both the R-squared and R.}
#' \item{beta}{The standardized regression coefficients.}
#' \item{lower}{The lowerbound of the condidence intervals of the standardized
#' regression coefficients  (if `conf = TRUE`).}
#' \item{upper}{The upperbound of the condidence intervals of the standardized
#' regression coefficients  (if `conf = TRUE`).}
#' \item{dfe}{The error degrees of freedom of the condidence intervals of
#' the standardized regression coefficients  (if `conf = TRUE`).}
#' \item{zero}{The zero-order correlations between the outcome and each predictor}
#' \item{total}{A matrix containing the betas and optionally (if `cor = TRUE`), the error degrees of
#' freedom, confidence intervals, and zero-order correlations.}
#'
#'
#' @details The function first completes the imputed datasets using 'mice::complete'.
#' It then calculates the linear model on each imputed dataset and averages
#' the standardized coefficients and correlations across imputations.
#' The final R-squared value is computed as the sum of the products of the
#' averaged standardized coefficients and averaged correlations.
#' The confidence intervals of the standardized regression coefficients are
#' calculated under the assumption that the variables are multivarate normally
#' distributed
#'
#'
#' @examples
#' library(mice)
#' imp <- mice(nhanes, print = FALSE, seed = 16117)
#' fit <- with(imp, lm(chl ~ age + hyp + bmi))
#' RsquareSP(fit)
#' @references
#' Van Ginkel, J.R., & Karch, J.D. (2024). A comparison of different measures of
#' the proportion of explained variance in multiply imputed data sets.
#' British Journal of Mathematical and Statistical Psychology. \doi{10.1111/bmsp.12344}
#'
#' Karch, J.D. (2024). Improving on Adjusted R-squared. Collabra: Psychology. \doi{10.1525/collabra.343}
#'
#' Van Ginkel, J.R. (2020). Standardized regression coefficients and newly proposed
#' estimators for R^2 in multiply imputed data. Psychometrika. \doi{10.1007/s11336-020-09696-4}

#' @export
RsquareSP <- function(object,
                      cor = FALSE,
                      conf = FALSE,
                      conf.level = 0.95,
                      alternative_adj_R2 = FALSE) {
  ## input checks
  if (!mice::is.mira(object)) {
    stop("The object must have class 'mira'")
  }
  if ((m <- length(object$analyses)) < 2) {
    stop("At least two imputations are needed for pooling.\n")
  }
  if (class((object$analyses[[1]]))[1] != "lm") {
    stop("r^2 can only be calculated for results of the 'lm' modeling function")
  }

  ## init
  alpha <- 1 - conf.level
  NumberOfImp <- length(object$analyses)
  DFE <- object$analyses[[1]]$df.residual
  vars <- colnames(object$analyses[[1]]$model)
  outcome <- vars[1]
  predictors <- vars[2:length(vars)]
  meanbeta <- meancor <- Umeanbeta <- cjmean <- Sxjsquaremean <-
    Sxjysquaremean <- rep(0, times = length(predictors))
  Sesquaremean <- bjsquaremean <- bSxbmean <- Sysquaremean <- 0
  meanbetam <- matrix(0, NumberOfImp, length(predictors))

  ## main loop
  for (m in 1:NumberOfImp) {
    datasetm <- object$analyses[[m]]$model
    model <- object$analyses[[1]]$call
    modelb <- object$analyses[[m]]
    modelbeta <- lm.beta::lm.beta(modelb)
    meanbetam[m, ] <- modelbeta$standardized.coefficients[
      2:length(modelbeta$standardized.coefficients)
    ]
    meanbeta <- meanbeta + meanbetam[m, ]
    meancor <- meancor + cor(datasetm)[1, 2:ncol(datasetm)]
    # CALCULATING BETA SES
    if (conf) {
      covxy <- stats::cov(datasetm)
      cj <- diag(solve(covxy[predictors, predictors]))
      cjmean <- cjmean + cj
      Sxjsquare <- diag(covxy)[predictors]
      Sxjsquaremean <- Sxjsquaremean + Sxjsquare
      Sxjysquare <- covxy[1, 2:ncol(covxy)]
      Sxjysquaremean <- Sxjysquaremean + Sxjysquare
      Sesquare <- as.vector((t(modelb$residuals) %*% modelb$residuals) / DFE)
      Sesquaremean <- Sesquaremean + Sesquare
      Sysquare <- var(datasetm[, outcome])
      Sysquaremean <- Sysquaremean + Sysquare
      bjsquare <- (modelb$coefficients)^2
      bjsquare <- bjsquare[2:length(bjsquare)]
      bjsquaremean <- bjsquaremean + bjsquare
      bSxb <- t(modelb$coefficients)[-1] %*% covxy[predictors, predictors] %*%
        (modelb$coefficients)[-1]
      bSxbmean <- bSxbmean + bSxb
      SEbeta <- sqrt((Sxjsquare * cj * Sesquare) / ((nrow(datasetm) - 3) * Sysquare) +
        (bjsquare * (Sxjsquare * as.vector(bSxb) - Sxjsquare * Sesquare - Sxjysquare)) /
        ((nrow(datasetm) - 3) * (sqrt(Sysquare))^4))
      SEb <- sqrt(diag(stats::vcov(modelb)))
      SEb <- SEb[2:length(SEb)]
      Umeanbeta <- Umeanbeta + SEbeta^2
    }
  }

  ## pool results
  vars <- c(
    "meanbeta", "meancor"
  )
  if (conf) {
    vars <- c(vars, "Sxjsquaremean", "Sxjysquaremean",
      "Sysquaremean", "bjsquaremean", "bSxbmean", "Umeanbeta", "cjmean",
      "Sesquaremean"
    )
  }
  for (var in vars) {
    assign(var, get(var) / NumberOfImp)
  }
  if (conf) {
    Bmeanbeta <- matrixStats::colVars(meanbetam)
    Tmeanbeta <- Umeanbeta + (1 + 1 / NumberOfImp) * Bmeanbeta
    ## CALCULATING DEGREES OF FREEDOM (REITER, 2007)
    vcom <- ((DFE + 1) / (DFE + 3)) * DFE
    gammameanbeta <- (1 + 1 / NumberOfImp) * Bmeanbeta / Tmeanbeta
    vmmeanbeta <- (NumberOfImp - 1) * gammameanbeta^(-2)
    vobsmeanbeta <- (1 - gammameanbeta) * vcom
    DFEmeanbeta <- (1 / vmmeanbeta + 1 / vobsmeanbeta)^(-1)
    lowermeanbeta <- meanbeta + stats::qt(alpha / 2, DFEmeanbeta) * sqrt(Tmeanbeta)
    uppermeanbeta <- meanbeta + stats::qt(1 - alpha / 2, DFEmeanbeta) * sqrt(Tmeanbeta)
  }

  ### adjusted stuff
  r_squared <- sum(meanbeta * meancor)
  N <- nrow(datasetm)
  p <- ncol(datasetm) - 1
  if (alternative_adj_R2) {
    alt_adjusted_rs <- altR2::estimate_adj_R2(r_squared, N = nrow(datasetm), p = ncol(datasetm) - 1)
    r_adj <- alt_adjusted_rs["Ezekiel"]
    alt_adjusted_rs <- alt_adjusted_rs[c(
      "Olkin_Pratt_Exact", "Pratt", "Claudy",
      "Wherry", "Smith", "Maximum_Likelihood",
      "Olkin_Pratt_K_1", "Olkin_Pratt_K_2",
      "Olkin_Pratt_K_5"
    )]
  } else {
    r_adj <- 1 - (N - 1) / (N - p - 1) * (1 - r_squared)
  }

  ## propagate results
  ### rsquared and relatives
  results <- list()
  class(results) <- "RsquaredPooled"
  results$r_squared <- r_squared
  results$r <- sqrt(results$r_squared)
  results$r_adj <- r_adj
  results$rtotal <- c(results$r_squared, results$r, results$r_adj)
  names(results$rtotal) <- c("R^2", "R", "adj. R^2")

  ## beta coefs
  results$beta <- meanbeta
  names(results$beta) <- predictors
  results$total <- as.matrix(results$beta)
  Names <- "Beta"

  if (conf) {
    results$lower <- lowermeanbeta
    results$upper <- uppermeanbeta
    results$dfe <- DFEmeanbeta
    results$total <- cbind(results$total, results$dfe, results$lower, results$upper)
    Names <- c(Names, "df", "Lowerbound", "Upperbound")
  }
  if (cor) {
    results$zero <- meancor
    names(results$zero) <- predictors
    results$total <- cbind(results$total, results$zero)
    Names <- c(Names, "Zero Order")
  }
  colnames(results$total) <- Names

  if (alternative_adj_R2) {
    results$alternative_adj_R2 <- alt_adjusted_rs
  }
  return(results)
}
