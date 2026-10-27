

# 按目标值分组，计算差值均值；保留原来的减组次序和边界容差
# y_diff = y - y_ref
#' @importFrom magrittr %>% %<>%
.qm_ref_means <- function(y, y_diff, n_quantiles, breaks, min_count = 5, keep_short = FALSE) {
  # 保留原来的同值排序规则
  sorted <- rev(sort(y, index.return = TRUE, decreasing = TRUE)$ix)
  group_means <- numeric(length(breaks))
  retry <- FALSE

  for (k in 2:(n_quantiles + 1)) {
    # 分位边界转换为排序位置，修正浮点数误差
    bounds <- floor(length(y)) * breaks[c(k - 1, k)]
    # 接近整数时先取整，避免浮点误差导致 floor 少取一个
    bounds <- ifelse(abs(bounds - round(bounds)) < 1E-9, round(bounds), floor(bounds))
    inds <- (bounds[1] + 1):bounds[2]
    yi_diff <- y_diff[sorted[inds]]
    enough <- sum(!is.na(yi_diff)) >= min_count
    if (!enough) {
      if (n_quantiles > 1) {
        retry <- TRUE
        n_quantiles <- n_quantiles - 1
      } else {
        warning("参考序列有效数据不足，无法继续减少分位数组数。", call. = FALSE)
      }
    }
    # 突变后侧即使样本不足也计算均值，沿用原逻辑
    if (enough || keep_short) group_means[k] <- mean(yi_diff, na.rm = TRUE)
  }

  # 两端复制相邻组，用于插值延拓
  group_means[1] <- group_means[2]
  group_means[n_quantiles + 2] <- group_means[n_quantiles + 1]
  list(bins = group_means, mean = mean(y_diff, na.rm = TRUE),
    Mq = n_quantiles, retry = retry)
}

# 估计各突变点两侧的分位组均值；有效样本不足时重新分组
.qm_ref_groups <- function(y, y_diff, Ips, Mq, Nadj, min_count, local_window) {
  Ns <- length(Ips) - 1
  starts <- c(0, head(Ips, -1))
  repeat {
    retry <- FALSE
    offset <- .5 / floor(Mq)
    breaks <- seq(0, by = 1 / Mq, length.out = Mq + 2)
    after <- before <- matrix(0, Ns + 1, Mq + 2)
    mean_shift <- numeric(Ns)
    for (i in seq_len(Ns)) {
      left <- if (local_window) max(starts[i], Ips[i] - Nadj) else starts[i]
      idx <- (left + 1):Ips[i]
      b <- .qm_ref_means(y[idx], y_diff[idx], Mq, breaks, min_count, keep_short = FALSE)
      before[i, ] <- b$bins
      Mq <- b$Mq
      retry <- retry || b$retry

      right <- if (local_window) min(Ips[i + 1], Ips[i] + Nadj) else Ips[i + 1]
      idx <- (Ips[i] + 1):right
      a <- .qm_ref_means(y[idx], y_diff[idx], Mq, breaks, min_count, keep_short = TRUE)
      after[i, ] <- a$bins
      Mq <- a$Mq
      retry <- retry || a$retry
      mean_shift[i] <- a$mean - b$mean
    }
    if (!retry) break
  }
  list(before = before, after = after, mean_shift = mean_shift,
    Mq = Mq, breaks = breaks, offset = offset)
}

#' 累加相邻突变的订正量
#' @param before,after 突变前、后的分位组均值矩阵。
#' @param ref_segment 基准段序号。
#' @return 各段的分位组订正量矩阵。
#' @example R/examples/ex_qm_adjustments.R
.qm_adjustments <- function(before, after, ref_segment) {
  n_changes <- nrow(before) - 1
  adjustments <- matrix(0, nrow(before), ncol(before))

  for (k in seq_len(ncol(before))) {
    # 基准段之前：累加从当前段到基准段的跳变
    if (ref_segment > 1) {
      for (i in seq_len(ref_segment - 1)) {
        changes <- i:(ref_segment - 1)
        adjustments[i, k] <- sum(after[changes, k]) - sum(before[changes, k])
      }
    }
    # 基准段之后：反向累加跳变
    if (ref_segment <= n_changes) {
      for (i in (ref_segment + 1):(n_changes + 1)) {
        changes <- ref_segment:(i - 1)
        adjustments[i, k] <- sum(before[changes, k]) - sum(after[changes, k])
      }
    }
  }
  adjustments
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
  derivatives <- splineN(knots, values, 2E30, 2E30)
  adjustment <- numeric(length(prob))
  for (i in seq_along(prob)) {
    adjustment[i] <- splintN(knots, values, derivatives, prob[i])
  }
  adjustment
}
