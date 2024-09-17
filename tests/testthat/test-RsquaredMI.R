check_numbers <- function(res, exp) {
  expect_equal(res, exp, ignore_attr = TRUE)
}

check_normal_res <- function(res) {
  check_numbers(res$rtotal, c(0.42250072, 0.65000055, 0.34000082))
  check_numbers(res$total[, "Beta"], c(0.58470893, 0.02098065, 0.56689497))
}

check_cor <- function(res) {
  check_numbers(res$total[, "Zero Order"], c(0.38055967, 0.32691257, 0.34067199))
}

check_df <- function(res) {
  check_numbers(res$total[, "df"], c(10.0132074, 17.8762882, 10.7091772))
}

check_conf <- function(res) {
  check_numbers(res$total[, "Lowerbound"], c(0.01044635, -0.42610359, 0.07355896))
  check_numbers(res$total[, "Upperbound"], c(1.1589715, 0.4680649, 1.0602310))
}

check_conf_90 <- function(res) {
  check_numbers(res$total[, "Lowerbound"], c(0.1175587887822, -0.3479899966129, 0.1646934226470))
  check_numbers(res$total[, "Upperbound"], c(1.0518590698629, 0.3899513041215, 0.9690965223859))
}

check_alternatives <- function(res) {
  check_numbers(
    res$alternative_adj_R2,
    c(
      0.3613807710744, 0.3595300248822, 0.3886422964785, 0.3700007801723,
      0.3437508126795, 0.3586601016493, 0.3646193027210, 0.3618120605415,
      0.3613837129718
    )
  )
}

expect_null <- function(res, fields) {
  present_fields <- fields[fields %in% names(res)]
  if (length(present_fields) > 0) {
    stop(paste("Fields", paste(present_fields, collapse = ", "), "are present in 'res'."))
  }
  invisible(TRUE)
}

test_that("normal run", {
  res <- RsquareSP(fit)
  check_normal_res(res)
  expect_null(res, c("lower", "upper", "dfe", "alternative_adj_R2", "zero"))
  print_out <- capture.output(print(res))
  exp_out <- c(
    "R-squared SP: ",
    "      R^2         R  adj. R^2 ",
    "0.4225007 0.6500006 0.3400008 ",
    "",
    "Beta Coefficients SP: ",
    "          Beta",
    "age 0.58470893",
    "hyp 0.02098065",
    "bmi 0.56689497"
  )
  expect_identical(print_out, exp_out)
})

test_that("cor true", {
  res <- RsquareSP(fit, cor = TRUE)
  check_normal_res(res)
  check_cor(res)
  expect_null(res, c("lower", "upper", "dfe", "alternative_adj_R2"))
})

test_that("conf true", {
  res <- RsquareSP(fit, conf = TRUE)
  check_normal_res(res)
  check_conf(res)
  # expect_null(res, c("zero", "alternative_adj_R2"))
})


test_that("conf changed to .9", {
  res <- RsquareSP(fit, conf = TRUE, conf.level = 0.9)
  check_normal_res(res)
  check_df(res)
  check_conf_90(res)
})

test_that("alternatives adjusted RSquared", {
  res <- RsquareSP(fit, alternative_adj_R2 = TRUE)
  check_normal_res(res)
  check_alternatives(res)

  print_out <- capture.output(print(res))
  exp_out <- c(
    "R-squared SP: ",
    "      R^2         R  adj. R^2 ",
    "0.4225007 0.6500006 0.3400008 ",
    "",
    "Beta Coefficients SP: ",
    "          Beta",
    "age 0.58470893",
    "hyp 0.02098065",
    "bmi 0.56689497",
    "",
    "Alternative adjusted R^2 estimates: ",
    " Olkin_Pratt_Exact              Pratt             Claudy             Wherry ",
    "         0.3613808          0.3595300          0.3886423          0.3700008 ",
    "             Smith Maximum_Likelihood    Olkin_Pratt_K_1    Olkin_Pratt_K_2 ",
    "         0.3437508          0.3586601          0.3646193          0.3618121 ",
    "   Olkin_Pratt_K_5 ",
    "         0.3613837 "
  )
  expect_identical(exp_out, print_out)
  expect_null(res, c("lower", "upper", "dfe", "zero"))
})

test_that("all on", {
  res <- RsquareSP(fit, cor = TRUE, conf = TRUE, alternative_adj_R2 = TRUE)
  check_normal_res(res)
  check_conf(res)
  check_cor(res)
  check_df(res)
  check_alternatives(res)
})


test_that("sanity unidimensional", {
  fit <- with(imp, lm(chl ~ age))
  res <- RsquareSP(fit, cor = TRUE, conf = TRUE, alternative_adj_R2 = TRUE)
  check_numbers(res$rtotal["R"], res$beta[1])
})

test_that("incorrect input", {
  expect_error(RsquareSP(5), "The object must have class 'mira'")
})

test_that("only one imputation", {
  imp <- mice(nhanes, print = FALSE, seed = 16117, m = 1)
  fit <- with(imp, lm(chl ~ age + hyp + bmi))
  expect_error(RsquareSP(fit), "At least two imputations are needed for pooling")
})

test_that("non lm input", {
  imp <- mice(nhanes, print = FALSE, seed = 16117)
  fit <- with(imp, t.test(chl ~ hyp))
  expect_error(RsquareSP(fit), "r\\^2 can only be calculated for results of the 'lm' modeling function")
})


