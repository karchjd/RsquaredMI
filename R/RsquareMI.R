print.Output <- function(x) {
  cat('R-squared SP:', '\n', x$R_squared,'\n','\n')
  if (!is.null(x$beta)) {
    cat('Beta Coefficients SP:', '\n')
    print(x$beta)
    cat('\n')
  }
  if (!is.null(x$zero)) {
    cat('Zero Order SP:', '\n')
    print(x$zero)
    cat('\n')
  }
}
RsquareSP = function(dataset,
                     model,
                     beta = FALSE,
                     cor = FALSE){
  library(lm.beta)
  results <- list(R_squared=NULL,Beta=NULL,Zero=NULL)
  class(results) <- "Output"
  NumberOfImp = dataset$m
  dataset <- complete(dataset, "long", inc = TRUE)
  modeltemp <- lm(model,dataset)
  vars <- names(modeltemp$model)
  outcome <- vars[1]
  predictors <- vars[2:length(vars)]
  meanbeta <- meancor <- rep(0, times=length(predictors))
  for (m in 1:NumberOfImp){    
    datasetm <- dataset[dataset[,1]==m,vars]
    modelb <- lm(model,datasetm)
    modelbeta <- lm.beta(modelb)
    meanbeta <- meanbeta + modelbeta$standardized.coefficients[2:length(modelbeta$standardized.coefficients)]
    meancor <- meancor + cor(datasetm)[1,2:ncol(datasetm)]
  }
  meanbeta <- meanbeta/NumberOfImp
  meancor <- meancor/NumberOfImp
  results$R_squared <- sum(meanbeta*meancor)

  if (beta == TRUE) {
     results$beta <- meanbeta
     names(results$beta) <- predictors
  }
  if (cor == TRUE) {
    results$zero <- meancor
    names(results$zero) <- predictors
  } 
  return(results)
}