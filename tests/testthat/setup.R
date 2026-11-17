library(mice)

# NOTE: The fixture and all expected values are based on the OLD mice 3.19.0
# (CRAN), i.e. before the pmm tie-breaking fix in mice 3.20.0
# (https://github.com/amices/mice/issues/757). Newer mice versions produce
# different imputations for the same seed. The fixture was created with mice 3.19.0 via
#   imp <- mice(nhanes, print = FALSE, seed = 16117)
#   saveRDS(complete(imp, "long", include = TRUE), "tests/testthat/fixtures/nhanes_imp.rds")
imp <- as.mids(readRDS(test_path("fixtures", "nhanes_imp.rds")))
fit <- with(imp, lm(chl ~ age + hyp + bmi))
