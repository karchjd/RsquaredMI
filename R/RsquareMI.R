print.RsquaredMI <- function(x) {
  cat("R-squared SP:", "\n", x$R_squared, "\n", "\n")
  if (!is.null(x$beta)) {
    cat("Beta Coefficients SP:", "\n")
    print(x$beta)
    cat("\n")
  }
  if (!is.null(x$zero)) {
    cat("Zero Order SP:", "\n")
    print(x$zero)
    cat("\n")
  }
}

#' Calculate R-squared with Standardized Predictors
#'
#' This function calculates the R-squared value for a linear model applied to a multiply imputed dataset.
#' Optionally, it can also return the standardized regression coefficients and zero-order correlations.
#'
#' @param dataset A multiply imputed dataset of class `mids` from the `mice` package.
#' @param model A linear model formula specifying the model to be fitted.
#' @param beta Logical. If `TRUE`, the function returns the standardized regression coefficients. Default is `FALSE`.
#' @param cor Logical. If `TRUE`, the function returns the zero-order correlations between the outcome and each predictor. Default is `FALSE`.
#'
#' @return A list of class `RsquaredMI` containing the following elements:
#' \item{R_squared}{The R-squared value calculated using standardized predictors.}
#' \item{beta}{The standardized regression coefficients (if `beta = TRUE`).}
#' \item{zero}{The zero-order correlations between the outcome and each predictor (if `cor = TRUE`).}
#'
#' @details The function first completes the imputed datasets using `mice::complete`.
#' It then calculates the linear model on each imputed dataset and averages the standardized coefficients and correlations across imputations.
#' The final R-squared value is computed as the sum of the products of the averaged standardized coefficients and averaged correlations.
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
RsquareSP <- function(dataset,
                      model,
                      beta = FALSE,
                      cor = FALSE) {
  results <- list(R_squared = NULL, Beta = NULL, Zero = NULL)
  class(results) <- "RsquaredMI"

  NumberOfImp <- dataset$m
  dataset <- mice::complete(dataset, "long", inc = TRUE)
  modeltemp <- lm(model, dataset)
  vars <- names(modeltemp$model)
  outcome <- vars[1]
  predictors <- vars[2:length(vars)]

  meanbeta <- meancor <- rep(0, times = length(predictors))
  for (m in 1:NumberOfImp) {
    datasetm <- dataset[dataset[, 1] == m, vars]
    modelb <- lm(model, datasetm)
    modelbeta <- lm.beta::lm.beta(modelb)
    meanbeta <- meanbeta + modelbeta$standardized.coefficients[2:length(modelbeta$standardized.coefficients)]
    meancor <- meancor + cor(datasetm)[1, 2:ncol(datasetm)]
  }
  meanbeta <- meanbeta / NumberOfImp
  meancor <- meancor / NumberOfImp
  results$R_squared <- sum(meanbeta * meancor)

  if (beta) {
    results$beta <- meanbeta
    names(results$beta) <- predictors
  }
  if (cor) {
    results$zero <- meancor
    names(results$zero) <- predictors
  }

  return(results)
}

RsquareSP <- function(dataset,
                      model,
                      beta = FALSE,
                      cor = FALSE) {
  results <- list(R_squared = NULL, Beta = NULL, Zero = NULL)
  class(results) <- "RsquaredMI"

  NumberOfImp <- dataset$m
  dataset <- mice::complete(dataset, "long", inc = TRUE)
  modeltemp <- lm(model, dataset)
  vars <- names(modeltemp$model)
  outcome <- vars[1]
  predictors <- vars[2:length(vars)]

  meanbeta <- meancor <- rep(0, times = length(predictors))
  for (m in 1:NumberOfImp) {
    datasetm <- dataset[dataset[, 1] == m, vars]
    modelb <- lm(model, datasetm)
    modelbeta <- lm.beta::lm.beta(modelb)
    meanbeta <- meanbeta + modelbeta$standardized.coefficients[2:length(modelbeta$standardized.coefficients)]
    meancor <- meancor + cor(datasetm)[1, 2:ncol(datasetm)]
  }
  meanbeta <- meanbeta / NumberOfImp
  meancor <- meancor / NumberOfImp
  results$R_squared <- sum(meanbeta * meancor)

  if (beta) {
    results$beta <- meanbeta
    names(results$beta) <- predictors
  }
  if (cor) {
    results$zero <- meancor
    names(results$zero) <- predictors
  }

  return(results)
}
