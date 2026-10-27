#' 分位数匹配订正（简化版）
#'
#' 将各段按段内排序分为 `Mq` 个分位组，把各组均值订正到基准段，
#' 再按分位位置做样条插值得到逐点订正量。
#'
#' @param y 待订正序列。
#' @param breaks 突变点位置，即各段（除末段外）的末个观测。
#' @param ref 基准段序号，默认最后一段。
#' @param Mq 分位组数上限；每组至少 `min_n` 个样本。
#' @param min_n 每个分位组的最少样本数。
#' @param window 突变点两侧用于估计订正量的最大观测数；0 表示整段。
#'
#' @return 订正后的序列。
#' @export
quantile_mapping <- function(y, breaks, ref = length(breaks) + 1, Mq = 10, min_n = 5, window = 0) {
  .quantile_mapping(y, breaks, ref, Mq, min_n, window)$y
}

# 返回订正后的序列 y、逐点订正量 W、段内分位位置 p、
# 各段分位组订正量 adj [Nseg, Mq + 2] 及实际分组数 Mq
.quantile_mapping <- function(y, breaks, ref, Mq, min_n, window) {
  n <- diff(c(0, breaks, length(y))) # 每段的长度
  seg <- rep(seq_along(n), n)        # 1:nseg
  Mq <- max(1, min(Mq, 100, min(n) %/% min_n))
  p <- ave(y, seg, FUN = \(x) rank(x, ties.method = "first") / length(x)) # 段内分位位置

  # 突变点前、后窗口内各分位组的均值 [Nbreak, Mq]
  w <- if (window > 0) window else Inf
  starts <- c(0, breaks)
  stops <- c(breaks, length(y))
  means <- \(x) as.vector(tapply(x, ceiling(rank(x, ties.method = "first") * Mq / length(x)), mean))
  before <- do.call(rbind, lapply(seq_along(breaks), \(j)
    means(y[(max(starts[j], breaks[j] - w) + 1):breaks[j]])))
  after <- do.call(rbind, lapply(seq_along(breaks), \(j)
    means(y[(breaks[j] + 1):min(stops[j + 1], breaks[j] + w)])))

  # 累加跳变，使各段订正到基准段 [Nseg, Mq]
  jump <- apply(rbind(0, after - before), 2, cumsum)
  adj <- -sweep(jump, 2, jump[ref, ])
  adj <- adj[, c(1, 1:Mq, Mq), drop = FALSE] # 两端复制相邻组

  knots <- (0:(Mq + 1) - .5) / Mq            # 组中心，两端各延拓一个节点
  W <- rep(NA_real_, length(y))
  for (i in setdiff(seq_along(n), ref)) {
    j <- seg == i
    W[j] <- splinefun(knots, adj[i, ], "natural")(p[j])
  }
  list(y = y + replace(W, is.na(W), 0), W = W, p = p, adj = adj, Mq = Mq)
}
