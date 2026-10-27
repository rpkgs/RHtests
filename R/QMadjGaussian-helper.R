# 逐次减少组数，直到各窗口、各分位组均满足最小有效样本数
.qm_ref_groups <- function(y, y_diff, Ips, Mq, Nadj, min_count, local_window) {
  starts <- c(0, head(Ips, -1))
  Ns <- length(Ips) - 1
  # 先取各突变点前、后的窗口；Nadj 为窗口长度
  before_windows <- after_windows <- vector("list", Ns)
  for (i in seq_len(Ns)) {
    left <- if (local_window) max(starts[i], Ips[i] - Nadj) else starts[i]
    right <- if (local_window) min(Ips[i + 1], Ips[i] + Nadj) else Ips[i + 1]
    before_windows[[i]] <- (left + 1):Ips[i]
    after_windows[[i]] <- (Ips[i] + 1):right
  }
  windows <- c(before_windows, after_windows)

  # 按目标值排序差值，逐次减少组数直至有效样本足够
  sorted_diffs <- lapply(windows, function(idx) {
    order_y <- rev(sort(y[idx], index.return = TRUE, decreasing = TRUE)$ix)
    y_diff[idx][order_y]
  })
  repeat {
    groups <- lapply(sorted_diffs, function(dif) {
      group_id <- ceiling(seq_along(dif) * Mq / length(dif))
      split(dif, factor(group_id, levels = seq_len(Mq)))
    })
    counts <- unlist(lapply(groups, function(window) {
      vapply(window, function(group) sum(!is.na(group)), integer(1))
    }))
    if (all(counts >= min_count) || Mq == 1) break
    Mq <- Mq - 1
  }
  if (any(counts < min_count)) {
    if (any(counts == 0)) stop("参考序列在部分窗口中没有有效数据。", call. = FALSE)
    warning("参考序列有效数据不足，无法继续减少分位数组数。", call. = FALSE)
  }

  # 计算各窗口的分位组均值；两端延拓，末行留给最后一段
  group_means <- do.call(rbind, lapply(groups, function(window) {
    vapply(window, mean, numeric(1), na.rm = TRUE)
  }))
  group_means <- unname(group_means[, c(1, seq_len(Mq), Mq), drop = FALSE])
  before <- rbind(group_means[seq_len(Ns), , drop = FALSE], 0)
  after <- rbind(group_means[Ns + seq_len(Ns), , drop = FALSE], 0)
  window_means <- vapply(windows, function(idx) mean(y_diff[idx], na.rm = TRUE), numeric(1))
  mean_shift <- window_means[Ns + seq_len(Ns)] - window_means[seq_len(Ns)]
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

# 整理分位组订正量，保留原有绘图数据列顺序
.qm_plot_means <- function(breaks, Adj, offset) {
  osmean <- seq_len(length(breaks))
  for (i in seq_len(nrow(Adj))) {
    osmean <- cbind(osmean, breaks - offset, Adj[i, ])
  }
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
