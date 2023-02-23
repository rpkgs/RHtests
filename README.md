
# RHtests

<!-- badges: start -->
[![R-CMD-check](https://github.com/rpkgs/RHtests/workflows/R-CMD-check/badge.svg)](https://github.com/rpkgs/RHtests/actions)
[![codecov](https://codecov.io/gh/rpkgs/RHtests/branch/master/graph/badge.svg)](https://codecov.io/gh/rpkgs/RHtests)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/RHtests)](https://cran.r-project.org/package=RHtests)
<!-- badges: end -->


## Installation

``` r
devtools::install_github("rpkgs/RHtest")
```

## 使用步骤

数据均一化一站式处理：

```r
FUN_FindU  <- if (has_ref) FindU.wRef else FindU
FUN_FindUD <- if (has_ref) FindUD.wRef else FindUD

RHtests_read(data, data.ref)
U  <- FUN_FindU(output = prefix, is_plot = is_plot)
if (is_empty(U$TP)) return(NULL)

UD <- FUN_FindUD(InCs = U$TP, output = prefix, is_plot = is_plot)
if (is_empty(UD$TP)) return(NULL)

TP  <- UD$TP
TP2 <- TP_adjustByMeta(TP, metadata, maxgap = maxgap)

r <- RHtests_stepsize(data = NULL, data.ref = NULL, TP2, has_ref,
    prefix, is_plot, verbose)
```

see details in `RHtests_process`.

## 版权

Please do not fork！

若愿共同维护，可私信获取授权。
