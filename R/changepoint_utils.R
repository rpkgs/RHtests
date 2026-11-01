# 选择订正基准段：0 表示最长段，超出范围时选择最后一段。
adjustment_segment <- function(Ips, Ns, Iadj) {
  if (Ns > 0) {
    Nsegs <- Ips - c(0, Ips[1:Ns])
    longest <- sort(Nsegs, index.return = TRUE, decreasing = TRUE)$ix[1]
  } else {
    longest <- 0
  }

  if (Iadj > (Ns + 1) | longest == 0) Ns + 1
  else if (Iadj == 0) longest
  else Iadj
}

# 优先匹配日期，否则取不晚于该日期的最后一个位置。
match_changepoint_dates <- function(dates, IY0) {
  N <- length(IY0)
  Ips <- c(rep(0, length(dates)), N)
  for (i in seq_along(dates)) {
    it <- match(dates[i], IY0)
    Ips[i] <- if (!is.na(it)) it else max(c(1:N)[IY0 <= dates[i]])
  }
  Ips
}
