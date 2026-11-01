# 将突变日期映射到当天或之前最近的观测位置。
locate_changepoints <- function(InCs, dates) {
  TP <- if (is.character(InCs)) fread(InCs) else InCs
  if (is.Date(TP$date)) TP$date <- date2num(TP$date)
  Ns <- nrow(TP)
  if (is.null(Ns)) Ns <- 0

  N <- length(dates)
  Ips <- c(rep(0, Ns), N)
  for (i in seq_len(Ns)) {
    it <- match(TP$date[i], dates)
    Ips[i] <- if (!is.na(it)) it else max(seq_len(N)[dates <= TP$date[i]])
  }
  list(Ips = Ips, Ids = TP$kind, Ns = Ns)
}

# 将各段订正量加到对应的观测值。
add_segment_adjustment <- function(x, segments, adjustment) {
  x + rep(adjustment, times = segments$end - segments$start + 1L)
}

# 选择订正所用的基准段。
select_base_segment <- function(breaks, Iadj) {
  Ns <- length(breaks) - 1
  longest <- if (Ns > 0) {
    sort(breaks - c(0, breaks[seq_len(Ns)]), index.return = TRUE, decreasing = TRUE)$ix[1]
  } else 0

  if (Iadj > Ns + 1 || longest == 0) Ns + 1
  else if (Iadj == 0) longest
  else Iadj
}

# 根据置信界限判定 0 类和 1 类突变点。
classify_changepoint <- function(Id, probL, probU, plev, statistic, lower, upper) {
  if (Id == 0) {
    if (statistic >= upper) return("Yes ")
    if (probU < plev) return("No  ")
    if (probL < plev) return("?   ")
    return("YifD")
  }
  if (Id == 1) {
    if (statistic < lower) return("No  ")
    if (statistic < upper) return("?   ")
    return("Yes ")
  }
  stop("Unknown changepoint type: ", Id)
}

# 恢复拟合时剔除的闰日数据。
restore_leap_days <- function(odata, original, break_dates, Rb, Ti, trend, QMout,
                          Mq, Iseg.adj, Nadj, Nt, fitted_cols, adjusted_cols) {
  md <- original$month * 100 + original$day
  observed <- !is.na(original$data)
  leap <- which(observed & md == 229)
  if (!length(leap)) return(odata)

  Ns <- length(break_dates) - 1
  if (Ns > 0) {
    tdata <- original[observed, ]
    dates <- ymd2num(tdata)
    Ips <- c(match(break_dates[seq_len(Ns)], dates), length(dates))
    leap_fit <- which(tdata$month * 100 + tdata$day == 229)
    Rb1 <- tdata$data
    Rb1[-leap_fit] <- Rb
    Ti1 <- rep(NA, length(dates))
    Ti1[-leap_fit] <- Ti
    for (i in leap_fit) {
      Rb1[i] <- tdata$data[i] + Rb1[i - 1] - tdata$data[i - 1]
      Ti1[i] <- Ti1[i - 1]
    }
    if (QMout$Mq > 1) {
      B <- QMadjGaussian(Rb1, Ips, Mq, Iseg.adj, Nadj, Nt)$PA + trend * Ti1
      odata$QM_adjusted[leap] <- round(B[leap_fit], 4)
    }
  } else {
    odata$QM_adjusted[md == 229] <- odata$base[md == 229]
  }

  for (i in leap) {
    neighbor <- if (Ns > 0 && odata$date[i - 1] %in% break_dates) i + 1 else i - 1
    for (column in fitted_cols) odata[[column]][i] <- odata[[column]][neighbor]
    for (column in adjusted_cols) {
      odata[[column]][i] <- odata$base[i] - (odata$base[neighbor] - odata[[column]][neighbor])
    }
  }
  odata
}
