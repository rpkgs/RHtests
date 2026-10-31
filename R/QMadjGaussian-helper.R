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
  sorted_diffs <- lapply(windows, \(idx) {
    order_y <- rev(sort(y[idx], index.return = TRUE, decreasing = TRUE)$ix)
    y_diff[idx][order_y]
  })
  repeat {
    groups <- lapply(sorted_diffs, \(dif) {
      group_id <- ceiling(seq_along(dif) * Mq / length(dif))
      split(dif, factor(group_id, levels = seq_len(Mq)))
    })
    counts <- unlist(lapply(groups, \(window) {
      map_int(window, \(group) sum(!is.na(group)))
    }))
    if (all(counts >= min_count) || Mq == 1) break
    Mq <- Mq - 1
  }
  if (any(counts < min_count)) {
    if (any(counts == 0)) stop("参考序列在部分窗口中没有有效数据。", call. = FALSE)
    warning("参考序列有效数据不足，无法继续减少分位数组数。", call. = FALSE)
  }

  # 计算各窗口的分位组均值；两端延拓，末行留给最后一段
  group_means <- lapply(groups, \(dif) map_dbl(dif, mean, na.rm = TRUE)) %>%
    do.call(rbind, .) %>% # [Nseg, Nq + 2]
    {
      unname(.[, c(1, seq_len(Mq), Mq)])
    }

  inds_before <- 1:Ns
  inds_after <- Ns + (1:Ns)
  before <- rbind(group_means[inds_before, , drop = FALSE], 0)
  after <- rbind(group_means[inds_after, , drop = FALSE], 0)

  window_means <- map_dbl(windows, \(idx) mean(y_diff[idx], na.rm = TRUE))
  mean_shift <- window_means[inds_after] - window_means[inds_before]
  list(
    before = before, after = after, mean_shift = mean_shift, Mq = Mq,
    breaks = seq(0, by = 1 / Mq, length.out = Mq + 2), offset = .5 / Mq
  )
}

#' 累加相邻突变的订正量
#' @param before,after 突变前、后的分位组均值矩阵。
#' @param ref_segment 基准段序号。
#' @return 各段的分位组订正量矩阵，维度与 `before` 相同：
#'   `[Nseg, Mq]`；含两端延拓值时为 `[Nseg, Mq + 2]`。
#' @example R/examples/ex_qm_adjustments.R
.qm_adjustments <- function(before, after, ref_segment) {
  idx <- seq_len(nrow(before) - 1)
  jumps <- rbind(0, after[idx, , drop = FALSE] - before[idx, , drop = FALSE])
  cumulative <- apply(jumps, 2, cumsum)
  -sweep(cumulative, 2, cumulative[ref_segment, ])
}

#' 整理分位组订正量，保留原有绘图数据列顺序
#'
#' @param breaks 分位组边界，长度为 `Mq + 2`。
#' @param Adj 订正量矩阵 `[Nseg, Mq + 2]`，行对应各段，列对应分位组及两端延拓值。
#' @param offset 分位组中心相对于边界的偏移量。
#' @return 绘图数据矩阵，维度为 `[Mq + 2, 1 + 2 * Nseg]`。
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
#' @return 各概率位置对应的订正量向量，长度为 `length(prob)`。
.qm_interpolate <- function(knots, values, prob, n_quantiles = length(knots) - 2) {
  if (n_quantiles == 1) values[1] <- values[2]
  adjustment <- splinefun(knots, values, "natural")(prob)
  # 样条节点之外保留旧版三次外推规则
  outside <- which(prob < min(knots) | prob > max(knots))
  if (length(outside)) {
    derivatives <- splineN(knots, values, 2E30, 2E30)
    adjustment[outside] <- map_dbl(prob[outside], \(p) {
      splintN(knots, values, derivatives, p)
    })
  }
  adjustment
}
