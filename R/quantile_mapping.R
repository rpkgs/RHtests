#' 分位数匹配订正（简化版）
#'
#' 将各段按段内排序分为 `Mq` 个分位组，把各组均值订正到基准段，
#' 再按分位位置做样条插值得到逐点订正量。
#'
#' @param y 待订正序列。
#' @param breaks 突变点位置，即各段（除末段外）的末个观测；无突变点时为空。
#' @param ref 基准段序号，默认最后一段。
#' @param Mq 分位组数上限；依据窗口样本数自动减少，最少为 1 组。
#' @param min_n 每个分位组的最少样本数；窗口不足时仍使用 1 组。
#' @param window 突变点两侧用于估计订正量的最大观测数；0 表示整段。
#'
#' @return 订正后的序列。
#' @export
quantile_mapping <- function(y, breaks, ref = length(breaks) + 1, Mq = 10, min_n = 5, window = 0) {
  .quantile_mapping(y, breaks, ref, Mq, min_n, window)$y
}

# Dif 为 NULL 时使用目标值，否则使用目标与参考序列的差值
.quantile_mapping <- function(y, breaks, ref, Mq, min_n, window, Dif = NULL) {
  n <- diff(c(0, breaks, length(y)))
  seg <- rep(seq_along(n), n)
  with_ref <- !is.null(Dif)
  w <- if (window > 0) window else Inf
  if (!with_ref) Mq <- max(1, min(Mq, 100, min(n, w) %/% min_n))
  p <- ave(y, seg, FUN = \(x) rank(x, ties.method = if (with_ref) "last" else "first") / length(x))
  AdjM <- rep(0, length(n))
  if (!length(breaks)) {
    return(list(y = y, W = rep(NA_real_, length(y)), p = p,
      adj = matrix(0, 1, Mq + 2), Mq = Mq, AdjM = AdjM))
  }

  # 突变点两侧各分位组的均值
  if (with_ref) {
    groups <- .qm_ref_groups(y, Dif, c(breaks, length(y)), Mq, window, min_n, window > 0)
    Mq <- groups$Mq
    before <- groups$before
    after <- groups$after
    if (window == 0 && length(breaks) > 1) {
      after[seq_len(length(breaks) - 1), ] <- before[2:length(breaks), ]
    }
    shift <- c(0, cumsum(groups$mean_shift))
    AdjM <- shift[ref] - shift
    knots <- groups$breaks - groups$offset
  } else {
    starts <- c(0, breaks)
    stops <- c(breaks, length(y))
    means <- \(x) as.vector(tapply(x, ceiling(rank(x, ties.method = "first") * Mq / length(x)), mean))
    before <- rbind(do.call(rbind, lapply(seq_along(breaks), \(j)
      means(y[(max(starts[j], breaks[j] - w) + 1):breaks[j]]))), 0)
    after <- rbind(do.call(rbind, lapply(seq_along(breaks), \(j)
      means(y[(breaks[j] + 1):min(stops[j + 1], breaks[j] + w)]))), 0)
    knots <- (0:(Mq + 1) - .5) / Mq
  }

  # 累加跳变并插值到各观测；有参考序列时增加均值约束
  adj <- .qm_adjustments(before, after, ref)
  if (!with_ref) adj <- adj[, c(1, 1:Mq, Mq), drop = FALSE]
  W <- rep(NA_real_, length(y))
  adjusted <- y
  for (i in setdiff(seq_along(n), ref)) {
    j <- seg == i
    W[j] <- .qm_interpolate(knots, adj[i, ], p[j], Mq)
    if (with_ref) {
      remain <- AdjM[i] - sum(W[j]) / n[i]
      selected <- j & if (remain < 0) W < 0 else W > 0
      W[selected] <- W[selected] + remain * n[i] * W[selected] / sum(W[selected])
    }
    adjusted[j] <- y[j] + W[j]
  }
  list(y = adjusted, W = W, p = p, adj = adj, Mq = Mq, AdjM = AdjM)
}
