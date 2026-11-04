#' RHtests 计算状态
#'
#' 保存观测序列、季节周期及检验临界值。由 [Read] 或 [Read.wRef] 构造，
#' 再传给 [FindU]、[FindUD]、[FindU.wRef]、[FindUD.wRef]、[StepSize]、
#' [StepSize.wRef] 复用，避免各步骤重复读数据或依赖全局变量。
#'
#' 各属性的来源：
#'
#' - `series`：`Read()`/`Read.wRef()` 内部 `ReadDLY.g()` 解码后整理的序列字段，
#'   无参考序列时含 `Y0`（观测值）、`IY0`（观测日期 YYYYMMDD）、`Ti`（时间序号，
#'   从首年起按 `Nt` 计）、`Imd`（月日 mmdd）、`IY0flg`（相邻观测是否连续，0/1）、
#'   `itable`（有效观测表）、`ori.itable`（原始表）、`ooflg`（观测在原始表中的行号）；
#'   有参考序列时改含 `Y0`（基准减参考的差值）、`TiB`、`IYBflg`、`bdata`、
#'   `ori.bdata`、`owflg`。
#' - `season`：同上的 `Icy`（季节周期内的月日序列）与 `Nt`（周期长度，
#'   年/月/日序列分别为 1、12、365/366）。
#' - `critical_values`：`Read()` 调 `readPFtable(N, plev)`、`Read.wRef()` 调
#'   `readPTtable(N, plev)` 生成，含 `phi`（自相关取值网格）与
#'   `PFmax`/`PTmax`（行数为序列长度 `N` 的临界值表）。
#' - `reference`：`Read.wRef()` 置 `TRUE`，`Read()` 置 `FALSE`。
#' - `plev`：构造时传入的置信水平，供 `check_rh_state()` 校验下游参数一致。
#'
#' @keywords internal
RHState <- S7::new_class(
  "RHState",
  properties = list(
    # Read()/Read.wRef() 整理后的序列字段；两种流程字段不同，校验时分别检查
    series = S7::class_list,
    # Read()/Read.wRef() 的 Icy、Nt
    season = S7::class_list,
    # readPFtable()/readPTtable() 的返回值：phi 与 PFmax 或 PTmax
    critical_values = S7::class_list,
    # Read.wRef() 为 TRUE，Read() 为 FALSE
    reference = S7::class_logical,
    # Read()/Read.wRef() 的 plev 实参
    plev = S7::class_numeric
  ),
  validator = function(self) {
    if (length(self@reference) != 1L || is.na(self@reference)) {
      return("reference must be a single non-missing logical value")
    }
    if (length(self@plev) != 1L || is.na(self@plev) ||
      !self@plev %in% c(0.75, 0.8, 0.9, 0.95, 0.99, 0.9999)) {
      return("plev must be a supported confidence level")
    }

    series <- self@series
    season <- self@season
    values <- self@critical_values
    required <- c("Y0", "IY0", "IY0flg", "Ti")
    required <- c(required, if (self@reference) {
      c("TiB", "IYBflg", "bdata", "ori.bdata", "owflg")
    } else {
      c("Imd", "itable", "ori.itable", "ooflg")
    })
    missing <- setdiff(required, names(series))
    if (length(missing)) {
      return(paste("Missing series fields: ", paste(missing, collapse = ", ")))
    }

    if (!is.numeric(series$Y0) || !is.numeric(series$IY0) || !is.numeric(series$Ti) ||
      !(is.numeric(series$IY0flg) || is.logical(series$IY0flg)) ||
      anyNA(series$IY0flg) || any(!series$IY0flg %in% c(0, 1))) {
      return("Values, dates and time indices must be numeric; continuity flags must be 0 or 1")
    }
    n <- length(series$Y0)
    if (n == 0L || any(lengths(series[c("IY0", "IY0flg", "Ti")]) != n)) {
      return("Values, dates, time indices and continuity flags must have the same nonzero length")
    }
    if (!is.numeric(season$Nt) || length(season$Nt) != 1L || is.na(season$Nt) ||
      season$Nt < 1L || length(season$Icy) != season$Nt) {
      return("Season length must match the number of cycle indices")
    }
    if (self@reference) {
      if (!is.data.frame(series$bdata) || !is.data.frame(series$ori.bdata) ||
        length(series$TiB) != nrow(series$bdata) ||
        length(series$IYBflg) != nrow(series$bdata) ||
        !is.logical(series$owflg) || length(series$owflg) != nrow(series$ori.bdata) ||
        anyNA(series$owflg) || sum(series$owflg) != nrow(series$bdata)) {
        return("Base series, time indices and observation flags have inconsistent lengths")
      }
    } else if (!is.data.frame(series$itable) || !is.data.frame(series$ori.itable) ||
      nrow(series$itable) != n || length(series$Imd) != n ||
      !is.numeric(series$ooflg) || length(series$ooflg) != n || anyNA(series$ooflg) ||
      any(series$ooflg < 1L | series$ooflg > nrow(series$ori.itable))) {
      return("Base series, cycle indices and observation positions have inconsistent lengths")
    }

    table <- if (self@reference) values$PTmax else values$PFmax
    if (!is.numeric(values$phi) || !is.matrix(table) ||
      ncol(table) != length(values$phi) || nrow(table) < n) {
      return("Critical value table does not match autocorrelation values or series length")
    }
    NULL
  }
)

# 由 Read()/Read.wRef() 调用：把读取结果整理为带类型与一致性校验的状态对象。
# series 中的 Icy、Nt 被拆到 season，其余字段原样进入 series。
new_rh_state <- function(series, critical_values, reference = FALSE, plev = 0.95) {
  season <- series[c("Icy", "Nt")]
  series[c("Icy", "Nt")] <- NULL
  RHState(
    series = series, season = season, critical_values = critical_values,
    reference = reference, plev = plev
  )
}

# 由 FindU/FindUD/StepSize（及 wRef 版）在计算前调用：
# 校验状态类型、是否带参考序列，以及 plev 是否与 Read()/Read.wRef() 建表时一致。
check_rh_state <- function(state, reference, plev) {
  if (!S7::S7_inherits(state, RHState) || state@reference != reference) {
    kind <- if (reference) "with" else "without"
    stop("A state ", kind, " a reference series is required", call. = FALSE)
  }
  if (!identical(as.numeric(state@plev), as.numeric(plev))) {
    stop("The confidence level must match the input state", call. = FALSE)
  }
  invisible(state)
}
