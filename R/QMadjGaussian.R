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
#' @param Nt 每年观测次数；无参考版本从全局环境读取。
#' @param Ny4a 有参考序列时，突变点两侧用于估计订正量的年数；0 表示整段。
#'
#' @return 列表：`PA` 为订正序列，`Mq` 为实际分组数，`Nseg.mn` 为最短段
#'   样本数，`osmean`、`osp` 用于绘图。无参考版本的 `W` 为逐点订正量；
#'   有参考版本另含 `DB`（最后一段的订正量）和 `AdjM`（各段均值订正量）。
#'
#' @example R/examples/ex_QMadjGaussian.R
#' @export
QMadjGaussian <- function(y, Ips, Mq, Iseg.adj, Nadj) {
  Ns <- length(Ips) - 1
  N <- length(y)
  starts <- c(0, head(Ips, -1))
  Nseg.mn <- min(N, Ips - starts)
  Ln <- if (Nt <= 12) 5 else 20
  Mq <- min(floor(Nseg.mn / Ln), if (Mq <= 0) 100 else Mq, 100)
  Mq <- max(Mq, 1)
  breaks <- seq(0, by = 1 / Mq, length.out = Mq + 2)
  EPb <- EPa <- matrix(0, Ns + 1, Mq + 2)

  # 估计突变前后各分位组的均值
  for (i in seq_len(Ns)) {
    left <- if (Nadj > 0) max(starts[i], Ips[i] - Nadj) else starts[i]
    EPb[i, ] <- .qm_means(y[(left + 1):Ips[i]], Mq, breaks)
    right <- if (Nadj > 0) min(Ips[i + 1], Ips[i] + Nadj) else Ips[i + 1]
    if (Nadj > 0 || i == Ns) {
      EPa[i, ] <- .qm_means(y[(Ips[i] + 1):right], Mq, breaks)
    }
  }
  if (Nadj == 0 && Ns > 1) EPa[seq_len(Ns - 1), ] <- EPb[2:Ns, ]
  Adj <- .qm_adjustments(EPb, EPa, Iseg.adj)

  # 保留绘图数据的原有列顺序
  osmean <- seq_len(Mq + 2)
  for (i in seq_len(Ns + 1)) {
    osmean <- cbind(osmean, breaks - .5 / Mq, Adj[i, ])
  }

  PA <- y
  W <- rep(NA, N)
  osp <- NULL
  for (i in seq_len(Ns + 1)) {
    idx <- (starts[i] + 1):Ips[i]
    ranks <- sort(sort(y[idx], index.return = TRUE)$ix, index.return = TRUE)$ix
    Rs <- ranks / length(idx)
    if (i != Iseg.adj) {
      W[idx] <- .qm_interpolate(breaks - .5 / Mq, Adj[i, ], Rs, Mq)
      PA[idx] <- y[idx] + W[idx]
    }
    ors <- sort(Rs, index.return = TRUE)$ix
    osp <- rbind(osp, cbind(idx, Rs[ors], W[idx][ors]))
  }
  dimnames(osp) <- NULL
  list(PA = PA, W = W, Mq = Mq, Nseg.mn = Nseg.mn, osmean = osmean, osp = osp)
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

  # 有效样本不足时，按原顺序减少组数并重算
  repeat {
    retry <- FALSE
    Fd <- .5 / floor(Mq)
    breaks <- seq(0, by = 1 / Mq, length.out = Mq + 2)
    EBa <- EBb <- matrix(0, Ns + 1, Mq + 2)
    EYd <- numeric(Ns)
    for (i in seq_len(Ns)) {
      left <- if (Ny4a > 0) max(starts[i], Ips[i] - Nadj) else starts[i]
      idx <- (left + 1):Ips[i]
      before <- .qm_ref_means(y[idx], Dif[idx], Mq, breaks, Nmin.Mq, keep_short = FALSE)
      EBb[i, ] <- before$bins
      Mq <- before$Mq
      retry <- retry || before$retry

      right <- if (Ny4a > 0) min(Ips[i + 1], Ips[i] + Nadj) else Ips[i + 1]
      idx <- (Ips[i] + 1):right
      after <- .qm_ref_means(y[idx], Dif[idx], Mq, breaks, Nmin.Mq, keep_short = TRUE)
      EBa[i, ] <- after$bins
      Mq <- after$Mq
      retry <- retry || after$retry
      EYd[i] <- after$mean - before$mean
    }
    if (!retry) break
  }

  AdjM <- rep(0, Ns + 1)
  if (Iseg.adj > 1) AdjM[seq_len(Iseg.adj - 1)] <- rev(cumsum(EYd[(Iseg.adj - 1):1]))
  if (Iseg.adj <= Ns) AdjM[(Iseg.adj + 1):(Ns + 1)] <- -cumsum(EYd[Iseg.adj:Ns])
  if (Ny4a == 0 && Ns > 1) EBa[seq_len(Ns - 1), ] <- EBb[2:Ns, ]
  Adj <- .qm_adjustments(EBb, EBa, Iseg.adj)

  osmean <- seq_len(Mq + 2)
  for (i in seq_len(Ns + 1)) {
    osmean <- cbind(osmean, breaks - Fd, Adj[i, ])
  }

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
    remain <- AdjM[i] - sum(DB) / Nseg
    selected <- if (remain < 0) DB < 0 else DB > 0
    DB[selected] <- DB[selected] + remain * Nseg * DB[selected] / sum(DB[selected])
    BA[idx] <- y[idx] + DB

    ors <- sort(Rs, index.return = TRUE)$ix
    osp <- rbind(osp, cbind(idx, Rs[ors], DB[ors]))
    ttmp <- cbind(Rs[ors], DB[ors])
    for (k in seq_len(Mq)) {
      osmean[k + 1, i * 2 + 1] <- mean(
        ttmp[ttmp[, 1] >= (k - 1) / Mq & ttmp[, 1] <= k / Mq, 2]
      )
    }
    osmean[1, i * 2 + 1] <- osmean[2, i * 2 + 1]
    osmean[Mq + 2, i * 2 + 1] <- osmean[Mq + 1, i * 2 + 1]
  }
  dimnames(osp) <- NULL
  list(PA = BA, DB = DB, Mq = Mq, Nseg.mn = Nseg.mn,
    osmean = osmean, osp = osp, AdjM = AdjM)
}
