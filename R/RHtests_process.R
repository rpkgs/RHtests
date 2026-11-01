RHtests_read <- function(data, data.ref = NULL, plev = 0.95) {
  if (is.null(data.ref)) Read(data, plev = plev) else Read.wRef(data, data.ref, plev = plev)
}


#' RHtests_stepsize
#' 
#' stepsize will not alter TP2, only delte non-significant record.
#'
#' @param data 基准序列；仅当 `.state` 为 `NULL` 时使用。
#' @param data.ref 参考序列，可为 `NULL`。
#' @param TP2 待检验的突变点表。
#' @param has_ref 是否使用参考序列流程；默认依据 `data.ref` 是否为空。
#' @param prefix 输出文件名前缀。
#' @param is_plot 是否绘图。
#' @param verbose 是否输出迭代信息。
#' @param ... 保留参数。
#' @param debug 是否写出检验过程统计量。
#' @param .state 上游传入的 `RHState` 对象；提供后不再读取输入数据。
#' @export
RHtests_stepsize <- function(
    data = NULL, data.ref = NULL, TP2,
    has_ref = !is.null(data.ref),
    prefix = "./OUTPUT/example03", 
    is_plot = FALSE, verbose = TRUE, ..., 
    debug = TRUE, .state = NULL)
{
  mkdir(dirname(prefix))

  # 每次迭代共享同一份状态，避免重复读取和污染全局环境。
  state <- if (is.null(.state)) RHtests_read(data, data.ref) else .state
  FUN_step <- if (has_ref) StepSize.wRef else StepSize

  if (nrow(TP2) == 0) {
    return(NULL)
  }
  r <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot, plev = state@plev, .state = state)

  times <- 1
  while (times < nrow(TP2)) {
    if (!is.list(r) || length(r$TP) == 0) {
      return(NULL)
    }
    TP2 <- TP_remove_least(r$TP)

    if (nrow(TP2) < nrow(r$TP)) {
      if (verbose) {
        ok(sprintf("[stepwize] %dth iteration:", times))
        print(TP2)
      }

      times <- times + 1
      r <- FUN_step(InCs = TP2, output = prefix, is_plot = is_plot, debug = debug,
                    plev = state@plev, .state = state)
    } else {
      break
    }
  }
  r %>% set_class("stepsize")
}

#' process RHtests without reference series
#' 
#' @param data A data.frame or data.table, with the columns of c('year',
#' 'month', 'day', 'data')
#' 
#' @param data.ref A data.frame or data.table, with the columns of c('year',
#' 'month', 'day', 'data'). The reference series of `data`, should have the same
#' length as `data`.
#' 
#' @param metadata A data.frame or data.table, with the columns of
#' `TurningPoint` date.
#' @param prefix 输出文件名前缀；为 `NULL` 时用 `outdir`。
#' @param outdir 输出目录。
#' @param maxgap 突变点与元数据日期的最大容差（天）。
#' @param ... 保留参数。
#' @param debug 是否写出检验过程统计量。
#' @param is_plot 是否绘图。
#' @param verbose 是否输出进度信息。
#'
#' @example R/examples/run_ex01.R
#' @example R/examples/run_ex02.R
#' 
#' @export
RHtests_process <- function(data, data.ref = NULL, metadata, 
    prefix = NULL, outdir = "OUTPUT/example01", 
    maxgap = 366, 
    ..., 
    debug = FALSE, 
    is_plot = TRUE, verbose = TRUE)
{
  mkdir(outdir)
  
  if(is.null(prefix)) prefix <- outdir
  prefix <- paste0(outdir, "/", basename(prefix))
  
  has_ref = !is.null(data.ref)
  FUN_FindU  <- if (has_ref) FindU.wRef else FindU
  FUN_FindUD <- if (has_ref) FindUD.wRef else FindUD
  
  state <- RHtests_read(data, data.ref)
  U <- FUN_FindU(output = prefix, is_plot = is_plot, debug = debug, .state = state)
  if (is_empty(U$TP)) return(NULL)

  UD <- FUN_FindUD(InCs = U$TP, output = prefix, is_plot = is_plot, .state = state)
  if (is_empty(UD$TP)) return(NULL)

  TP  <- UD$TP
  TP2 <- TP_adjustByMeta(TP, metadata, maxgap = maxgap)
  
  r <- RHtests_stepsize(data = NULL, data.ref = NULL, TP2, has_ref,
      prefix, is_plot, verbose, debug = debug, .state = state)
  # if null returned,
  if (is.null(r$TP) || nrow(r$TP) == 0) return(NULL)
  r$TP %<>% merge_metainfo(metadata)
  r
}
