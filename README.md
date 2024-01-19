# fastglmpca

[![R-CMD-check](https://github.com/stephenslab/fastglmpca/workflows/R-CMD-check/badge.svg)](https://github.com/stephenslab/fastglmpca/actions)
[![Build status](https://ci.appveyor.com/api/projects/status/o8xpogal3pfca7ub?svg=true)](https://ci.appveyor.com/project/pcarbo/fastglmpca)

Fast algorithms for estimating GLM-PCA models to count data.

To install the package, clone or download the git repository, then run

```R
remotes::install_local("fastglmpca")
```

assuming your R working directory contains a copy of the `fastglmpca`
repository.

If you have issues downloading or running the package, please post a
github issue or email ericweine15@gmail.com.

## Developer notes

This is the command used to check the package before submitting to
CRAN:

```r
library(rhub)
check_for_cran(".",show_status = TRUE,
  env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false",
               `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"))
```
