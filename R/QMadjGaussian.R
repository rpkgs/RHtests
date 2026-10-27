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
#' @return 列表：`PA` 为订正序列，`Mq` 为实际分组数，`Nseg.mn` 为最短段
#'   样本数，`osmean`、`osp` 用于绘图。无参考版本的 `W` 为逐点订正量；
#'   有参考版本另含 `DB`（最后一段的订正量）和 `AdjM`（各段均值订正量）。
#'
#' @example R/examples/ex_QMadjGaussian.R
#' @export
QMadjGaussian <- function(y, Ips, Mq, Iseg.adj, Nadj, Nt) {
  min_n <- if (Nt <= 12) 5 else 20
  if (Mq <= 0) Mq <- 100
  r <- .quantile_mapping(y, head(Ips, -1), Iseg.adj, Mq, min_n, Nadj)

  # 保留绘图数据的原有格式
  breaks <- seq(0, by = 1 / r$Mq, length.out = r$Mq + 2)
  osmean <- .qm_plot_means(breaks, r$adj, .5 / r$Mq)
  segs <- split(seq_along(y), rep(seq_along(Ips), diff(c(0, Ips))))
  osp <- do.call(rbind, lapply(segs, \(idx) {
    o <- order(r$p[idx])
    cbind(idx, r$p[idx][o], r$W[idx][o])
  }))
  dimnames(osp) <- NULL
  list(PA = r$y, W = r$W, Mq = r$Mq, Nseg.mn = min(diff(c(0, Ips))),
    osmean = osmean, osp = osp)
}

#' @rdname QMadjGaussian
#' @export
QMadjGaussian.wRef <- function(y, Dif, Ips, Mq, Iseg.adj, Nadj, Nt, Ny4a) {
  segs <- split(seq_along(y), rep(seq_along(Ips), diff(c(0, Ips))))
  Nmin.Mq <- if (Nt <= 12) 5 else 20
  Nseg.mn <- min(length(y), vapply(segs, function(idx) sum(!is.na(Dif[idx])), integer(1)))
  Mq <- max(1, min(Nseg.mn %/% Nmin.Mq, if (Mq < 0) 100 else Mq, 100))

  r <- .quantile_mapping(y, head(Ips, -1), Iseg.adj, Mq, Nmin.Mq, Ny4a * Nt, Dif)
  Mq <- r$Mq
  breaks <- seq(0, by = 1 / Mq, length.out = Mq + 2)
  osmean <- .qm_plot_means(breaks, r$adj, .5 / Mq)

  # 兼容旧版绘图数据
  osp <- NULL
  for (i in seq_along(segs)) {
    idx <- segs[[i]]
    if (i == Iseg.adj) {
      osp <- rbind(osp, cbind(idx, 0, 0))
      next
    }
    p <- r$p[idx]
    DB <- r$W[idx]
    osp <- rbind(osp, cbind(idx, p[order(p)], DB[order(p)]))
    bin_means <- vapply(seq_len(Mq), function(k) {
      mean(DB[p >= (k - 1) / Mq & p <= k / Mq])
    }, numeric(1))
    osmean[, i * 2 + 1] <- bin_means[c(1, seq_len(Mq), Mq)]
  }
  dimnames(osp) <- NULL
  list(PA = r$y, DB = r$W[segs[[length(segs)]]], Mq = Mq, Nseg.mn = Nseg.mn,
    osmean = osmean, osp = osp, AdjM = r$AdjM)
}
