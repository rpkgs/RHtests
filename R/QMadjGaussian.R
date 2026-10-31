#' 分位数匹配订正
#'
#' @param y 待订正序列，按时间排序。
#' @param Dif 目标与参考序列的差值，与 `y` 对齐；缺失值可为 `NA`。
#' @param Ips 突变点位置及末尾观测的位置；突变点为前一段的末个观测。
#' @param Mq 分位数组数上限；1 表示均值订正，样本不足时自动减少。
#'   无参考序列时 `Mq <= 0`、有参考序列时 `Mq < 0` 自动确定组数。
#'   有参考序列时，逐次减少组数，直至各窗口的每组有效样本数满足要求；
#'   减至 1 组仍不足时给出警告，窗口无有效数据时报错。
#' @param Iseg.adj 基准段序号，该段不作订正。
#' @param Nadj 突变点两侧用于估计订正量的最大观测数；0 表示整段。
#'   有参考序列时，由 `Ny4a * Nt` 覆盖。
#' @param Nt 每年观测次数。
#' @param Ny4a 有参考序列时，突变点两侧用于估计订正量的年数；0 表示整段。
#'
#' @return 列表，`Nseg` 为段数：
#'   - `PA`：订正后序列。
#'   - `Mq`：实际分位组数。
#'   - `Nseg.mn`：最短段样本数；有参考序列时为各段有效差值的最小数量。
#'   - `osmean`：绘图矩阵 `[Mq + 2, 1 + 2 * Nseg]`，首列为行号，其后每段依次为分位位置和订正量。
#'   - `osp`：绘图矩阵 `[length(y), 3]`，依次为原始位置、排序后的分位位置和对应订正量；首列不随排序变化。
#'   - `W`：逐点订正量（无参考序列）。
#'   - `DB`：最后一段的订正量（有参考序列）。
#'   - `AdjM`：各段均值订正量（有参考序列）。
#'
#' @example R/examples/ex_QMadjGaussian.R
#' @export
QMadjGaussian <- function(y, Ips, Mq, Iseg.adj, Nadj, Nt) {
  min_n <- if (Nt <= 12) 5 else 20
  if (Mq <= 0) Mq <- 100
  r <- .quantile_mapping(y, breaks = head(Ips, -1), ref = Iseg.adj,
    Mq = Mq, min_n = min_n, window = Nadj)

  # 保留绘图数据的原有格式
  probs <- seq(0, by = 1 / r$Mq, length.out = r$Mq + 2)
  group_plot <- .qm_plot_means(probs, r$adj, .5 / r$Mq)
  segs <- .qm_segments(y, Ips)
  # point_plot[length(y), 3]：原始位置、段内排序后的分位位置及对应订正量；首列不随排序变化
  point_plot <- do.call(rbind, lapply(segs, \(idx) {
    pos <- r$quantile_pos[idx]
    adjustment <- r$adjustment[idx]
    order_p <- order(pos)
    cbind(idx, pos[order_p], adjustment[order_p])
  }))
  dimnames(point_plot) <- NULL
  listk(PA = r$y, W = r$adjustment, Mq = r$Mq, Nseg.mn = min(diff(c(0, Ips))),
    osmean = group_plot, osp = point_plot)
}

#' @rdname QMadjGaussian
#' @export
QMadjGaussian.wRef <- function(y, Dif, Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a) {
  segs <- .qm_segments(y, Ips)
  min_n <- if (Nt <= 12) 5 else 20
  Nseg.mn <- min(length(y), map_int(segs, \(idx) sum(!is.na(Dif[idx]))))
  Mq <- max(1, min(Nseg.mn %/% min_n, if (Mq < 0) 100 else Mq, 100))

  r <- .quantile_mapping(y, breaks = head(Ips, -1), ref = Iseg.adj, Mq = Mq,
    min_n = min_n, window = Ny4a * Nt, Dif = Dif)
  Mq <- r$Mq
  probs <- seq(0, by = 1 / Mq, length.out = Mq + 2)
  group_plot <- .qm_plot_means(probs, r$adj, .5 / Mq)

  # point_plot[length(y), 3] 同上；基准段的分位位置和订正量按旧版格式置零
  point_plot <- NULL
  for (i in seq_along(segs)) {
    idx <- segs[[i]]
    if (i == Iseg.adj) {
      point_plot <- rbind(point_plot, cbind(idx, 0, 0))
      next
    }
    quantile_pos <- r$quantile_pos[idx]
    adjustment <- r$adjustment[idx]
    order_p <- order(quantile_pos)
    point_plot <- rbind(point_plot, cbind(idx, quantile_pos[order_p], adjustment[order_p]))
    # 沿用两端闭区间，边界上的观测计入相邻两组
    bin_means <- map_dbl(seq_len(Mq), \(k) {
      mean(adjustment[quantile_pos >= (k - 1) / Mq & quantile_pos <= k / Mq])
    })
    # 每段占两列：分位位置、订正量
    group_plot[, i * 2 + 1] <- bin_means[c(1, seq_len(Mq), Mq)]
  }
  dimnames(point_plot) <- NULL
  listk(PA = r$y, DB = r$adjustment[segs[[length(segs)]]], Mq, Nseg.mn,
    osmean = group_plot, osp = point_plot, AdjM = r$AdjM)
}

# 按突变点将观测位置分段
.qm_segments <- function(y, Ips) {
  split(seq_along(y), rep(seq_along(Ips), diff(c(0, Ips))))
}
