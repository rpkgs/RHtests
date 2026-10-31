#' 绘制各段 QM 订正量的分布图
#'
#' @param QMout `QMadjGaussian` 的结果，含 `osp` 与 `osmean`。
#' @param Ips 各段末个观测的位置。
#' @param Iseg.adj 基准段序号，该段不绘制。
#' @param Ns 突变点个数。
#'
#' @return 无返回值，直接绘图。
plot_qm_adjustments <- function(QMout, Ips, Iseg.adj, Ns) {
  Mq <- QMout$Mq
  if (Ns < 1) return(invisible(NULL))

  osp <- QMout$osp
  osmean <- QMout$osmean
  par(mar = c(4, 5, 3, 2) + .1, cex = .8, mfrow = c(2, 2), mgp = c(1.2, .5, 0))
  Fd <- .5 / Mq
  ymax <- max(osp[, 2:3], na.rm = TRUE)
  ymin <- min(osp[, 2:3], na.rm = TRUE)

  col <- 0
  np <- 0
  for (i in seq_len(Ns + 1)) {
    I1 <- if (i == 1) 1 else Ips[i - 1] + 1
    I2 <- Ips[i]
    if (i == Iseg.adj) next
    np <- np + 1
    # 各段的分位位置与订正量
    pos <- osmean[[paste0("quantile_pos_", i)]]
    adj <- osmean[[paste0("adjustment_", i)]]
    if (col == 0) {
      col <- 2
      plot(osp[I1:I2, 2], osp[I1:I2, 3],
        xlim = c(0, 1), ylim = c(ymin, ymax),
        type = "l", lwd = 1, col = col,
        xlab = "Cumulative Frequency", ylab = "QM Adjustment"
      )
      title(cex.main = .9, main = paste("distribution of QM adjustments with Mq=", Mq), line = .5)
    } else {
      col <- col + 1
      lines(osp[I1:I2, 2], osp[I1:I2, 3], lwd = 1, col = col)
    }
    # 分位组箱体
    for (j in seq_len(Mq)) {
      lines(c(pos[j + 1] - Fd, pos[j + 1] + Fd), rep(adj[j + 1], 2),
        col = col, lty = 2, lwd = .5
      )
      if (j < Mq) {
        lines(rep(pos[j + 1] + Fd, 2), c(adj[j + 1], adj[j + 2]),
          col = col, lty = 2, lwd = .5
        )
      }
    }
    text(.15, ymax - np * (ymax - ymin) / (Ns * 3), paste("Seg.", i))
    lines(c(.25, .30), rep(ymax - np * (ymax - ymin) / (Ns * 3), 2), lwd = 2, col = col)
  }
}
