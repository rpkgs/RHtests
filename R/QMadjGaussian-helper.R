# 按目标值排序，计算各分位组的差值均值和有效样本数
.qm_ref_means <- function(y, y_diff, Mq) {
  sorted <- rev(sort(y, index.return = TRUE, decreasing = TRUE)$ix)
  group <- ceiling(seq_along(y) * Mq / length(y))
  values <- split(y_diff[sorted], factor(group, levels = seq_len(Mq)))
  list(means = vapply(values, mean, numeric(1), na.rm = TRUE),
    counts = vapply(values, function(x) sum(!is.na(x)), integer(1)))
}

# 逐次减少组数，直到各窗口、各分位组均满足最小有效样本数
.qm_ref_groups <- function(y, y_diff, Ips, Mq, Nadj, min_count, local_window) {
  starts <- c(0, head(Ips, -1))
  Ns <- length(Ips) - 1
  windows <- unlist(lapply(seq_len(Ns), function(i) {
    left <- if (local_window) max(starts[i], Ips[i] - Nadj) else starts[i]
    right <- if (local_window) min(Ips[i + 1], Ips[i] + Nadj) else Ips[i + 1]
    list((left + 1):Ips[i], (Ips[i] + 1):right)
  }), recursive = FALSE)
  repeat {
    groups <- lapply(windows, function(idx) .qm_ref_means(y[idx], y_diff[idx], Mq))
    enough <- all(vapply(groups, function(x) all(x$counts >= min_count), logical(1)))
    if (enough || Mq == 1) break
    Mq <- Mq - 1
  }
  if (!enough) {
    counts <- unlist(lapply(groups, `[[`, "counts"))
    if (any(counts == 0)) stop("参考序列在部分窗口中没有有效数据。", call. = FALSE)
    warning("参考序列有效数据不足，无法继续减少分位数组数。", call. = FALSE)
  }

  before <- after <- matrix(0, Ns + 1, Mq + 2)
  mean_shift <- numeric(Ns)
  for (i in seq_len(Ns)) {
    b <- 2 * i - 1
    a <- 2 * i
    before[i, ] <- groups[[b]]$means[c(1, seq_len(Mq), Mq)]
    after[i, ] <- groups[[a]]$means[c(1, seq_len(Mq), Mq)]
    mean_shift[i] <- mean(y_diff[windows[[a]]], na.rm = TRUE) -
      mean(y_diff[windows[[b]]], na.rm = TRUE)
  }
  list(before = before, after = after, mean_shift = mean_shift, Mq = Mq,
    breaks = seq(0, by = 1 / Mq, length.out = Mq + 2), offset = .5 / Mq)
}

#' 累加相邻突变的订正量
#' @param before,after 突变前、后的分位组均值矩阵。
#' @param ref_segment 基准段序号。
#' @return 各段的分位组订正量矩阵。
#' @example R/examples/ex_qm_adjustments.R
.qm_adjustments <- function(before, after, ref_segment) {
  idx <- seq_len(nrow(before) - 1)
  jumps <- rbind(0, after[idx, , drop = FALSE] - before[idx, , drop = FALSE])
  cumulative <- apply(jumps, 2, cumsum)
  -sweep(cumulative, 2, cumulative[ref_segment, ])
}

# 调整逐点订正量，使其均值与均值跳变订正一致
.qm_match_mean <- function(adjustment, target) {
  n <- length(adjustment)
  remain <- target - sum(adjustment) / n
  selected <- if (remain < 0) adjustment < 0 else adjustment > 0
  adjustment[selected] <- adjustment[selected] +
    remain * n * adjustment[selected] / sum(adjustment[selected])
  adjustment
}

# 整理分位组订正量，保留原有绘图数据列顺序
.qm_plot_means <- function(breaks, Adj, offset) {
  osmean <- seq_len(length(breaks))
  for (i in seq_len(nrow(Adj))) {
    osmean <- cbind(osmean, breaks - offset, Adj[i, ])
  }
  osmean
}

# 有参考序列时，用逐点订正量更新绘图用的分位组均值
.qm_plot_ref <- function(osmean, segment, prob, adjustment, Mq) {
  ttmp <- cbind(prob, adjustment)
  col <- segment * 2 + 1
  for (k in seq_len(Mq)) {
    osmean[k + 1, col] <- mean(
      ttmp[ttmp[, 1] >= (k - 1) / Mq & ttmp[, 1] <= k / Mq, 2]
    )
  }
  osmean[1, col] <- osmean[2, col]
  osmean[Mq + 2, col] <- osmean[Mq + 1, col]
  osmean
}

#' 将分位组订正量插值到各观测值
#'
#' @param knots 分位组中心及两端延拓的插值节点。
#' @param values 节点对应的订正量。
#' @param prob 待插值的概率位置。
#' @param n_quantiles 分位数组数。
#'
#' @return 各概率位置对应的订正量。
.qm_interpolate <- function(knots, values, prob, n_quantiles = length(knots) - 2) {
  if (n_quantiles == 1) values[1] <- values[2]
  adjustment <- splinefun(knots, values, "natural")(prob)
  # 样条节点之外保留旧版三次外推规则
  outside <- which(prob < min(knots) | prob > max(knots))
  if (length(outside)) {
    derivatives <- splineN(knots, values, 2E30, 2E30)
    adjustment[outside] <- vapply(prob[outside], function(p) {
      splintN(knots, values, derivatives, p)
    }, numeric(1))
  }
  adjustment
}
