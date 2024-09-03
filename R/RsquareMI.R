print.Output <- function(x) {
  cat("R-squared SP:", "\n")
  print(x$Rtotal)
  cat("\n")
  cat("Beta Coefficients SP:", "\n")
  print(x$total)
  cat("\n")
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
#' of the standardized regression coefficients. Default is `FALSE`.
#' @param cor Logical. If `TRUE`, the function returns the zero-order
#' correlations between the outcome and each predictor. Default is `FALSE`.
#' @param alpha A real number between 0 and 1 specifying the significance level
#' of the confidence intervals. Default is `0.05`.
#'
#' @return A list of class `RsquaredMI` containing the following elements:
#' \item{R_squared}{The R-squared value calculated using standardized predictors.}
#' \item{R}{The square root of the R-squared value, or the multiple correlation R.}
#' \item{Rtotal}{A vector containing both the R-squared and R.}
#' \item{beta}{The standardized regression coefficients.}
#' \item{lower}{The lowerbound of the condidence intervals of the standardized
#' regression coefficients  (if `conf = TRUE`).}
#' \item{upper}{The upperbound of the condidence intervals of the standardized
#' regression coefficients  (if `conf = TRUE`).}
#' \item{dfe}{The error degrees of freedom of the condidence intervals of
#' the standardized regression coefficients  (if `conf = TRUE`).}
#' \item{zero}{The zero-order correlations between the outcome and each
#' predictor (if `cor = TRUE`).}
#'
#' @details The function first completes the imputed datasets using `mice::complete`.
#' It then calculates the linear model on each imputed dataset and averages
#' the standardized coefficients and correlations across imputations.
#' The final R-squared value is computed as the sum of the products of the
#' averaged standardized coefficients and averaged correlations.
#' The confidence intervals of the standardized regression coefficients are
#' calculated under the assumption that the variables are multivarate normally
#' distributed
#'
#' @import mice
#' @importFrom lm.beta lm.beta
#'
#' @examples
#' \dontrun{
#' # Assuming `imp` is a multiply imputed dataset from `mice`
#' model <- as.formula(y ~ x1 + x2)
#' result <- RsquareSP(dataset = imp, model = model, beta = TRUE, cor = TRUE)
#' print(result$R_squared)
#' print(result$beta)
#' print(result$zero)
#' }
#'
#' @export
RsquareSP <- function(object,
                      cor = FALSE,
                      conf = FALSE,
                      alpha = 0.05) {

  ## input checks
  if (!is.mira(object)) {
    stop("The object must have class 'mira'")
  }
  if (is.mira(object)) {
    if ((m <- length(object$analyses)) < 2) {
      stop("At least two imputations are needed for pooling.\n")
    }
    if (class((object$analyses[[1]]))[1] != "lm") {
      stop("r^2 can only be calculated for results of the 'lm' modeling function")
    }
  }

  ## init variables
  results <- list(R_squared = NULL, R = NULL, Rtotal = NULL, Beta = NULL,
                  Lower = NULL, Upper = NULL, Dfe = NULL, Zero = NULL,
                  Total = NULL)
  class(results) <- "Output"
  NumberOfImp <- length(object$analyses)
  datasetm <- object$analyses[[1]]$model
  DFE <- object$analyses[[1]]$df.residual
  vars <- colnames(datasetm)
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
      2:length(modelbeta$standardized.coefficients)]
    meanbeta <- meanbeta + meanbetam[m, ]
    meancor <- meancor + cor(datasetm)[1, 2:ncol(datasetm)]

    # CALCULATING BETA SES
    covxy <- cov(datasetm)
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
    SEb <- sqrt(diag(vcov(modelb)))
    SEb <- SEb[2:length(SEb)]
    Umeanbeta <- Umeanbeta + SEbeta^2
  }

  ## pool results
  meanbeta <- meanbeta / NumberOfImp
  meancor <- meancor / NumberOfImp
  Sxjsquaremean <- Sxjsquaremean / NumberOfImp
  Sxjysquaremean <- Sxjysquaremean / NumberOfImp
  Sysquaremean <- Sysquaremean / NumberOfImp
  bjsquaremean <- bjsquaremean / NumberOfImp
  bSxbmean <- bSxbmean / NumberOfImp
  Umeanbeta <- Umeanbeta / NumberOfImp
  cjmean <- cjmean / NumberOfImp
  Sesquaremean <- Sesquaremean / NumberOfImp
  Bmeanbeta <- matrixStats::colVars(meanbetam)
  Tmeanbeta <- Umeanbeta + (1 + 1 / NumberOfImp) * Bmeanbeta

  ## CALCULATING DEGREES OF FREEDOM (REITER, 2007)
  vcom <- ((DFE + 1) / (DFE + 3)) * DFE
  gammameanbeta <- (1 + 1 / NumberOfImp) * Bmeanbeta / Tmeanbeta
  vmmeanbeta <- (NumberOfImp - 1) * gammameanbeta^(-2)
  vobsmeanbeta <- (1 - gammameanbeta) * vcom
  DFEmeanbeta <- (1 / vmmeanbeta + 1 / vobsmeanbeta)^(-1)

  lowermeanbeta <- meanbeta + qt(alpha / 2, DFEmeanbeta) * sqrt(Tmeanbeta)
  uppermeanbeta <- meanbeta + qt(1 - alpha / 2, DFEmeanbeta) * sqrt(Tmeanbeta)

  ## propagate results
  results$R_squared <- sum(meanbeta * meancor)
  results$R <- sqrt(results$R_squared)
  results$Rtotal <- c(results$R_squared, results$R)
  names(results$Rtotal) <- c("R^2", "R")
  results$beta <- meanbeta
  results$lower <- lowermeanbeta
  results$upper <- uppermeanbeta
  results$dfe <- DFEmeanbeta
  names(results$beta) <- predictors
  results$total <- as.matrix(results$beta)
  Names <- "Beta"
  colnames(results$total) <- Names
  if (conf == TRUE) {
    results$lower <- lowermeanbeta
    results$upper <- uppermeanbeta
    results$dfe <- DFEmeanbeta
    results$total <- cbind(results$total, results$dfe, results$lower, results$upper)
    Names <- c(Names, "df", "Lowerbound", "Upperbound")
    colnames(results$total) <- Names
  }
  if (cor == TRUE) {
    results$zero <- meancor
    names(results$zero) <- predictors
    results$total <- cbind(results$total, results$zero)
    colnames(results$total) <- c(Names, "Zero Order")
  }
  return(results)
}
