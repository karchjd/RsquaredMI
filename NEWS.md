# RSquaredMI 0.2.1

* Tests now use stored imputed datasets instead of calling `mice()`, so they
  no longer depend on the installed 'mice' version. This fixes test failures
  with 'mice' 3.20.0, whose fix to predictive mean matching changes the
  imputations for a given seed. No user-facing changes.

# RSquaredMI 0.2.0

# RSquaredMI 0.1.2

* Added support for categorical predictors
