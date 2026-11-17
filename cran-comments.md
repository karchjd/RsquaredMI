## Submission

This is a maintenance release requested by the 'mice' maintainer ahead of the
'mice' 3.20.0 release. 'mice' 3.20.0 fixes predictive mean matching, which
changes the imputations for a given seed and made our tests fail in the
reverse dependency check. The tests now use stored imputed datasets, so they
no longer depend on the 'mice' version. There are no changes to the package
code.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no reverse dependencies.
