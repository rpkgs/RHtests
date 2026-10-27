
# 计算各分位组均值，并延拓两端
.qm_means <- function(y, Mq, breaks) {
  sorted <- sort(y, index.return = TRUE)$ix
  bounds <- floor(length(y) * breaks)
  values <- numeric(Mq + 2)
  for (k in 2:(Mq + 1)) {
    values[k] <- mean(y[sorted[(bounds[k - 1] + 1):bounds[k]]])
  }
  values[1] <- values[2]
  values[Mq + 2] <- values[Mq + 1]
  values
}

# 按目标值分组，计算差值均值；保留原来的减组次序和边界容差
.qm_ref_means <- function(y, Dif, Mq, breaks, min_count, keep_short) {
  iindx <- rev(sort(y, index.return = TRUE, decreasing = TRUE)$ix)
  bins <- numeric(length(breaks))
  retry <- FALSE
  for (k in 2:(Mq + 1)) {
    bounds <- floor(length(y)) * breaks[c(k - 1, k)]
    near <- abs(bounds - round(bounds)) < 1E-9
    bounds[near] <- round(bounds[near])
    bounds <- floor(bounds)
    values <- Dif[iindx[(bounds[1] + 1):bounds[2]]]
    enough <- sum(!is.na(values)) >= min_count
    if (!enough) {
      if (Mq > 1) {
        retry <- TRUE
        Mq <- Mq - 1
      } else {
        ErrorMSG <<- paste(
          "!!!! Warning - The Ref series has too many missing values to do QM-adjustments.\n",
          "Only mean-adjustments were done. You may want to use a better Ref series if possible!\n",
          get("ErrorMSG", envir = .GlobalEnv), "\n"
        )
      }
    }
    # 突变后侧即使样本不足也计算均值，沿用原逻辑
    if (enough || keep_short) bins[k] <- mean(values, na.rm = TRUE)
  }
  bins[1] <- bins[2]
  bins[Mq + 2] <- bins[Mq + 1]
  list(bins = bins, mean = mean(Dif, na.rm = TRUE), Mq = Mq, retry = retry)
}

# 累加相邻突变的订正量，保持原有求和顺序
.qm_adjustments <- function(before, after, ref) {
  Ns <- nrow(before) - 1
  Adj <- matrix(0, nrow(before), ncol(before))
  for (k in seq_len(ncol(before))) {
    if (ref > 1) {
      for (i in seq_len(ref - 1)) {
        Adj[i, k] <- sum(after[i:(ref - 1), k]) - sum(before[i:(ref - 1), k])
      }
    }
    if (ref <= Ns) {
      for (i in (ref + 1):(Ns + 1)) {
        Adj[i, k] <- sum(before[ref:(i - 1), k]) - sum(after[ref:(i - 1), k])
      }
    }
  }
  Adj
}

# 将分位组订正量插值到各观测值
.qm_interpolate <- function(dx, fdx, prob, Mq) {
  if (Mq == 1) fdx[1] <- fdx[2]
  fdx2 <- splineN(dx, fdx, 2E30, 2E30)
  adjustment <- numeric(length(prob))
  for (i in seq_along(prob)) {
    adjustment[i] <- splintN(dx, fdx, fdx2, prob[i])
  }
  adjustment
}
