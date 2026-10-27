#' 分位数匹配订正
#'
#' @param y 待订正序列，按时间排序。
#' @param Dif 目标与参考序列的差值，与 `y` 对齐；缺失值可为 `NA`。
#' @param Ips 突变点位置及末尾观测的位置；突变点为前一段的末个观测。
#' @param Mq 分位数组数上限；1 表示均值订正，样本不足时自动减少。
#'   无参考序列时 `Mq <= 0`、有参考序列时 `Mq < 0` 自动确定组数。
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
  Ns <- length(Ips) - 1
  N <- length(y)
  starts <- c(0, head(Ips, -1))
  Nmin.Mq <- if (Nt <= 12) 5 else 20
  Nseg.mn <- min(N, vapply(seq_along(Ips), function(i) {
    sum(!is.na(Dif[(starts[i] + 1):Ips[i]]))
  }, integer(1)))
  Mq <- min(floor(Nseg.mn / Nmin.Mq), if (Mq < 0) 100 else Mq, 100)
  if (Mq == 0) Mq <- 1
  Nadj <- Ny4a * Nt

  # 估计突变前后的分位组均值，样本不足时减少组数
  groups <- .qm_ref_groups(y, Dif, Ips, Mq, Nadj, Nmin.Mq, Ny4a > 0)
  Mq <- groups$Mq
  breaks <- groups$breaks
  Fd <- groups$offset
  EBb <- groups$before
  EBa <- groups$after
  EYd <- groups$mean_shift

  AdjM <- rep(0, Ns + 1)
  if (Iseg.adj > 1) AdjM[seq_len(Iseg.adj - 1)] <- rev(cumsum(EYd[(Iseg.adj - 1):1]))
  if (Iseg.adj <= Ns) AdjM[(Iseg.adj + 1):(Ns + 1)] <- -cumsum(EYd[Iseg.adj:Ns])
  if (Ny4a == 0 && Ns > 1) EBa[seq_len(Ns - 1), ] <- EBb[2:Ns, ]
  Adj <- .qm_adjustments(EBb, EBa, Iseg.adj)

  osmean <- .qm_plot_means(breaks, Adj, Fd)

  BA <- y
  osp <- NULL
  for (i in seq_len(Ns + 1)) {
    idx <- (starts[i] + 1):Ips[i]
    Nseg <- length(idx)
    DB <- rep(NA, Nseg)
    if (i == Iseg.adj) {
      osp <- rbind(osp, cbind(idx, 0, 0))
      next
    }

    # 同值排序沿用原规则，避免改变分位位置
    iindx <- rev(sort(y[idx], index.return = TRUE, decreasing = TRUE)$ix)
    irank <- sort(iindx, index.return = TRUE)$ix
    Rs <- irank / Nseg
    DB <- .qm_interpolate(breaks - Fd, Adj[i, ], Rs, Mq)

    # 使分位数订正量的均值与均值跳变订正一致
    DB <- .qm_match_mean(DB, AdjM[i])
    BA[idx] <- y[idx] + DB

    ors <- sort(Rs, index.return = TRUE)$ix
    osp <- rbind(osp, cbind(idx, Rs[ors], DB[ors]))
    osmean <- .qm_plot_ref(osmean, i, Rs[ors], DB[ors], Mq)
  }
  dimnames(osp) <- NULL
  list(PA = BA, DB = DB, Mq = Mq, Nseg.mn = Nseg.mn,
    osmean = osmean, osp = osp, AdjM = AdjM)
}
