# Convert year, month and day columns to YYYYMMDD.
.make_date <- function(data) {
  data$year * 10000 + data$month * 100 + data$day
}

# Map changepoint dates to the last observation on or before each date.
stepsize_points <- function(InCs, dates) {
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

# Get inclusive boundaries for each segment.
stepsize_bounds <- function(Ips) {
  list(start = c(1, Ips[-length(Ips)] + 1), end = Ips)
}

# Choose the reference segment for adjustment.
stepsize_segment <- function(Ips, Iadj) {
  Ns <- length(Ips) - 1
  longest <- if (Ns > 0) {
    sort(Ips - c(0, Ips[seq_len(Ns)]), index.return = TRUE, decreasing = TRUE)$ix[1]
  } else 0

  if (Iadj > Ns + 1 || longest == 0) Ns + 1
  else if (Iadj == 0) longest
  else Iadj
}

# Classify type-0 and type-1 changepoints using their confidence bounds.
stepsize_classify <- function(Id, probL, probU, plev, statistic, lower, upper) {
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

# Place fitted values back into the original observation positions.
stepsize_expand <- function(x, rows, n) {
  result <- rep(NA_real_, n)
  result[rows] <- round(x, 4)
  result
}

# Restore leap-day values omitted from the fitting series.
stepsize_leap <- function(odata, original, break_dates, Rb, Ti, trend, QMout,
                          Mq, Iseg.adj, Nadj, Nt, fitted_cols, adjusted_cols) {
  md <- original$month * 100 + original$day
  observed <- !is.na(original$data)
  leap <- which(observed & md == 229)
  if (!length(leap)) return(odata)

  Ns <- length(break_dates) - 1
  if (Ns > 0) {
    tdata <- original[observed, ]
    dates <- .make_date(tdata)
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

# Write quantile-mapping diagnostics for the selected reference segment.
stepsize_qm_log <- function(QMout, segments, Iseg.adj, Ny4a, adjustment, file) {
  fprintf("Nseg_shortest = %d ; Mq =  %d ; Ny4a =  %g \n",
    QMout$Nseg.mn, QMout$Mq, Ny4a, file = file, append = TRUE
  )
  fprintf("\n Adjust to segment %d : from %d to %d \n",
    Iseg.adj, segments$start[Iseg.adj], segments$end[Iseg.adj], file = file, append = TRUE
  )
  if (QMout$Mq <= 1) return(invisible(NULL))

  Ns <- length(segments$end) - 1
  oline <- paste("#Fcat: frequency category boundaries\n",
    "#DP: Difference in the category means\n#", sep = ""
  )
  for (i in seq_len(Ns)) oline <- paste(oline, paste0("Fcat.Seg", i), paste0("DP.Seg", i))
  cat(paste(oline, "\n"), file = file, append = TRUE)
  write.table(round(QMout$osmean, 4), file = file, append = TRUE, row.names = FALSE, col.names = FALSE)
  for (i in seq_along(segments$end)) {
    if (i != Iseg.adj) {
      fprintf("Seg. %d: mean of QM-adjustments =%.4f\n", i, adjustment(i), file = file, append = TRUE)
    }
  }
}
