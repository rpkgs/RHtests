#' 示例数据与临界值表
#'
#' 供 [FindU()]、[StepSize()] 等函数的示例与测试使用。
#'
#' @format
#' - `dat_example01`: 636 × 4，`year`、`month`、`day`（月序列 `day = 0`）、`data`
#' - `dat_example02`: 612 × 5，同上，另有参考序列列 `data.ref`
#' - `dat_prcp`: 9131 × 4，日降水，缺测为 `NA`
#' - `dat_temp`: 22645 × 3，`date`、`Tavg`、`ref`
#' - `PFtable`: 792 × 8 矩阵，PMF 临界值表，首列为序列长度
#' - `PTtable`: 1232 × 8 矩阵，PMT 临界值表，结构同上
#'
#' @name RHtests-data
#' @aliases dat_example01 dat_example02 dat_prcp dat_temp PFtable PTtable
NULL
