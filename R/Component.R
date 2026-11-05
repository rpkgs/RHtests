# y 为观测值，t_mmdd 为对应的月日（MMDD）。
# 返回按月日升序排列的季节均值 mu，以及去季节序列 anomaly。
RemoveSeasonal <- function(y, t_mmdd) {
  inds <- sort(unique(t_mmdd)) # 月尺度最多 12 个，日尺度最多 366 个
  mu <- numeric(length(inds))
  anomaly <- numeric(length(y))
  for (i in seq_along(inds)) {
    ind <- which(t_mmdd == inds[i])
    mu[i] <- mean(y[ind], na.rm = TRUE)
    anomaly[ind] <- y[ind] - mu[i]
  }
  list(mu = mu, anomaly = anomaly)
}

#' Multiple breakpoints piecewise Regression (same slope, different interception)
#'
#' @param Y 响应值（观测序列）。
#' @param X 自变量（通常为时间序列）。
#' @param Ips 各段末尾的位置，末位为序列长度。
#' @param ... 传给 `lm_solve()`。
#'
#' @return
#' - `sig`   : coefficients
#' - `fitted`: fitted value
#' - `resi`  : residual of fitted value
LS_MultiBreaks <- function(Y, X, Ips, ...) {
  n_breaks <- length(Ips) - 1
  design <- cbind(1, X - mean(X))
  for (i in seq_len(n_breaks)) {
    segment <- numeric(length(Y))
    segment[(Ips[i] + 1):Ips[i + 1]] <- 1
    design <- cbind(design, segment)
  }
  lm_solve(design, Y, ...)
}

#' single-breakpoint piecewise regression (same slope, different interception)
#'
#' @param Y The response vector
#' @param T_anorm The predictor vector
#'
#' @param Ic the position of change point, `[1, Ic]` and `[Ic + 1, end]` are
#' corresponding to the perid
#' @param ... 传给 `lm_solve()`。
#'
#' @return
#' - `sig`   : coefficients
#' - `fitted`: fitted value
#' - `resi`  : residual of fitted value
#'
#' - `SSE`   : sum of square error
#' @export
LS_OneBreak <- function(Y, X, Ic, ...) {
  n <- length(Y)
  # T_anorm = T - mean(T)
  design <- cbind(1, X)
  if (!is.na(Ic)) {
    design <- cbind(design, c(rep(0, Ic), rep(1, n - Ic)))
  }
  lm_solve(design, Y, ...)
}

lm_solve <- function(design, response, only.SSE = FALSE, ...) {
  l <- .lm.fit(design, as.matrix(response))
  # sig  <- as.matrix(l$coefficients)
  fitted <- design %*% l$coefficients
  # resi <- response-fitted
  SSE <- sum(l$residuals^2)
  if (!only.SSE) {
    list(
      sig = l$coefficients,
      fitted = as.vector(fitted),
      resi = l$residuals,
      SSE = SSE
    )
  } else {
    list(SSE = SSE)
  }
}

# microbenchmark::microbenchmark(
#   sig  = solve(t(design) %*% design) %*% t(design) %*% response,
#   sig3 = as.matrix(.lm.fit(design, response)$coefficients),
#   sig2 = chol2inv(chol(t(design) %*% design)) %*% t(design) %*% response,
#   times = 100
# )
