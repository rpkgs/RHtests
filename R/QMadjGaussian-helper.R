
# 计算各分位组均值，并延拓两端
.qm_means <- function(y, breaks) {
  n_quantiles <- length(breaks) - 2
  sorted <- sort(y, index.return = TRUE)$ix
  bounds <- floor(length(y) * breaks)
  values <- numeric(n_quantiles + 2)

  for (k in 2:(n_quantiles + 1)) {
    inds <- (bounds[k - 1] + 1):bounds[k]
    values[k] <- mean(y[sorted[inds]])
  }
  values[1] <- values[2]
  values[n_quantiles + 2] <- values[n_quantiles + 1]
  values
}

# 按目标值分组，计算差值均值；保留原来的减组次序和边界容差
.qm_ref_means <- function(y, difference, n_quantiles, breaks, min_count, keep_short) {
  # 保留原来的同值排序规则
  sorted <- rev(sort(y, index.return = TRUE, decreasing = TRUE)$ix)
  group_means <- numeric(length(breaks))
  retry <- FALSE

  for (k in 2:(n_quantiles + 1)) {
    # 分位边界转换为排序位置，修正浮点数误差
    bounds <- floor(length(y)) * breaks[c(k - 1, k)]
    near_integer <- abs(bounds - round(bounds)) < 1E-9
    bounds[near_integer] <- round(bounds[near_integer])
    bounds <- floor(bounds)
    inds <- (bounds[1] + 1):bounds[2]
    values <- difference[sorted[inds]]
    enough <- sum(!is.na(values)) >= min_count
    if (!enough) {
      if (n_quantiles > 1) {
        retry <- TRUE
        n_quantiles <- n_quantiles - 1
      } else {
        warning("参考序列有效数据不足，无法继续减少分位数组数。", call. = FALSE)
      }
    }
    # 突变后侧即使样本不足也计算均值，沿用原逻辑
    if (enough || keep_short) group_means[k] <- mean(values, na.rm = TRUE)
  }

  # 两端复制相邻组，用于插值延拓
  group_means[1] <- group_means[2]
  group_means[n_quantiles + 2] <- group_means[n_quantiles + 1]
  list(bins = group_means, mean = mean(difference, na.rm = TRUE),
    Mq = n_quantiles, retry = retry)
}

# 累加相邻突变的订正量，保持原有求和顺序
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

# 将分位组订正量插值到各观测值
.qm_interpolate <- function(knots, values, prob, n_quantiles) {
  if (n_quantiles == 1) values[1] <- values[2]
  derivatives <- splineN(knots, values, 2E30, 2E30)
  adjustment <- numeric(length(prob))
  for (i in seq_along(prob)) {
    adjustment[i] <- splintN(knots, values, derivatives, prob[i])
  }
  adjustment
}
