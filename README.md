
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

数据均一化，一站式处理：

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


## References

1. Xu, W., Li, Q., Wang, X. L., Yang, S., Cao, L., & Feng, Y. (2013). Homogenization of Chinese daily surface air temperatures and analysis of trends in the extreme temperature indices: HOMOGENIZE CHINESE DAILY TEMPERATURES. Journal of Geophysical Research: Atmospheres, 118(17), 9708–9720. <https://doi.org/10.1002/jgrd.50791>

2. 朱亚妮, 曹丽娟, 唐国利, & 周自江. (2015). 中国地面相对湿度非均一性检验及订正. **气候变化研究进展**, 11(6), 379–386.

3. Li, Z., Yan, Z., Zhu, Y., Freychet, N., & Tett, S. (2020). Homogenized Daily Relative Humidity Series in China during 1960–2017. Advances in Atmospheric Sciences, 37(4), 318–327. <https://doi.org/10.1007/s00376-020-9180-0>
