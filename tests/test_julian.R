library(mice)
imp <- mice(nhanes, print = FALSE, seed = 16117)
fit <- with(imp, lm(chl ~ age + hyp + bmi))

# input: mira object
tmp <- RsquareSP(fit, conf = TRUE, cor = TRUE)
print(tmp)
